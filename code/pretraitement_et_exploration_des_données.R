###############################################################################
#  CHAPITRE 2 — COMPRÉHENSION, PRÉTRAITEMENT ET EXPLORATION DES DONNÉES
#  (analyse sur la vitesse brute, pas de transformation avancée ici)
###############################################################################

## ───────────────────────── 0 • Librairies & thème ──────────────────────────
library(dplyr)          # manipulation
library(tibble)         # tibble
library(ggplot2)        # graphiques
library(momentuHMM)     # classe de données (pas de modélisation ici)
library(circular)       # statistiques directionnelles
library(fitdistrplus)   # ajustement Gamma
library(data.table)     # data.table
library(hexbin)         # nuages hexbin
theme_set(theme_minimal())

## ───────────────────────── 0bis • Dossier « figures » ──────────────────────
figdir <- file.path(getwd(), "figures")
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

# ferme tout PDF resté ouvert (au cas où)
while (!is.null(dev.list())) dev.off()

## ───────────────────────── 1 • Description des données ─────────────────────
data <- read.table(file.choose(), header = TRUE, sep = ",", stringsAsFactors = FALSE)
str(data)
summary(data)

## ───────────────────────── 2 • Préparation & nettoyage ─────────────────────
wrap_pi <- function(th) ((th + pi) %% (2 * pi)) - pi       # ramène l’angle dans (–π, π]

compute_step_angle <- function(df) {
  df %>% 
    arrange(t) %>% 
    mutate(
      dx    = x - lag(x),
      dy    = y - lag(y),
      dt    = t - lag(t),
      dir   = atan2(dy, dx),
      step  = sqrt(dx^2 + dy^2) / dt,
      angle = wrap_pi(dir - lag(dir))
    )
}

df_HMM <- data %>% 
  mutate(ID = id) %>% 
  dplyr::select(ID, x, y, t) %>%            # on force dplyr::select
  group_by(ID) %>% 
  compute_step_angle() %>% 
  ungroup() %>% 
  filter(!is.na(step), !is.na(angle))

class(df_HMM) <- c("momentuHMMData", "data.frame")

cat("Observations finales :", nrow(df_HMM), "\n")
cat("Trajectoires (ID)    :", dplyr::n_distinct(df_HMM$ID), "\n")

## ───────────────────────── 3 • Trajectoires TOP-10 ─────────────────────────
top10_ids <- df_HMM %>% 
  count(ID, name = "n_pas") %>% 
  arrange(desc(n_pas)) %>% 
  slice_head(n = 10) %>% 
  pull(ID)

# points début / fin
start_end <- df_HMM %>% 
  filter(ID %in% top10_ids) %>% 
  group_by(ID) %>% 
  slice(c(1, n())) %>%            # première et dernière ligne
  mutate(point_type = c("Début", "Fin")) %>% 
  ungroup()

g_traj <- df_HMM %>% 
  filter(ID %in% top10_ids) %>% 
  mutate(ID = factor(ID, levels = top10_ids)) %>% 
  ggplot(aes(x, y, group = ID)) +
  geom_path(colour = "steelblue", linewidth = .4, alpha = .8) +
  geom_point(data = start_end,
             aes(x = x, y = y, colour = point_type, shape = point_type),
             size = 1.8) +
  scale_colour_manual(values = c("Début" = "red", "Fin" = "darkgreen")) +
  scale_shape_manual(values = c("Début" = 16, "Fin" = 17)) +
  coord_equal() +
  facet_wrap(~ ID, ncol = 5) +
  labs(title = "Trajectoires des 10 fourmis les plus suivies",
       x = "X (mm)", y = "Y (mm)", colour = "", shape = "")

print(g_traj)
ggsave(file.path(figdir, "02_traj_top10.pdf"), g_traj, width = 8, height = 4.5)



## ───────────────────────── 4 • Parcours le plus long ───────────────────────
long_ID <- names(which.max(table(df_HMM$ID)))
df_long <- filter(df_HMM, ID == long_ID)

g_pluslong <- ggplot(df_long, aes(x = x, y = y)) +
  geom_point(size = .3) +
  coord_equal() +
  labs(title = paste0("Trajectoire complète — ID ", long_ID),
       x = "X (mm)", y = "Y (mm)")

print(g_pluslong)
ggsave(file.path(figdir, "02_traj_plus_long.pdf"),
       g_pluslong, width = 8, height = 5.5)

## ───────────────────────── 5 • Distribution de la vitesse ──────────────────
fit_gamma <- fitdist(df_long$step, dist = "gamma")

g_step <- ggplot(df_long, aes(step)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue4", alpha = .6) +
  geom_density(colour = "red", linewidth = 1.1) +
  stat_function(fun  = dgamma,
                args = list(shape = fit_gamma$estimate["shape"],
                            rate  = fit_gamma$estimate["rate"]),
                colour = "black", linewidth = 1.1) +
  labs(title    = paste("Distribution de la vitesse — ID", long_ID),
       subtitle = sprintf("Gamma(shape = %.2f, rate = %.2f) · AIC = %.1f",
                          fit_gamma$estimate["shape"],
                          fit_gamma$estimate["rate"],
                          fit_gamma$aic),
       x = "step (mm·s⁻¹)", y = "Densité")

print(g_step)
ggsave(file.path(figdir, "02_hist_gamma_step.pdf"),
       g_step, width = 8, height = 5.5)

pdf(file.path(figdir, "02_diag_gamma_step.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2))
denscomp(fit_gamma, legendtext = "Gamma ajustée")
qqcomp(fit_gamma, main = "Q–Q plot")
ppcomp(fit_gamma, main = "P–P plot")
cdfcomp(fit_gamma, legendtext = "Gamma ajustée")
dev.off()


#___________________________distrubution de step^2________________________
## ───────── 5bis • Distribution de la vitesse transformée (step^2) ─────────

# Ajout de la transformation
df_long <- df_long %>%
  mutate(step2 = step^2)

# Ajustement d'une Gamma sur step^2
fit_gamma2 <- fitdist(df_long$step2, dist = "gamma")

# Histogramme + densité noyau + densité Gamma ajustée
g_step2 <- ggplot(df_long, aes(step2)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue4", alpha = .6) +
  geom_density(colour = "red", linewidth = 1.1) +
  stat_function(fun  = dgamma,
                args = list(shape = fit_gamma2$estimate["shape"],
                            rate  = fit_gamma2$estimate["rate"]),
                colour = "black", linewidth = 1.1) +
  labs(
    title    = paste("Distribution de la vitesse transformée — ID", long_ID),
    subtitle = sprintf("Gamma(shape = %.2f, rate = %.2f) · AIC = %.1f",
                       fit_gamma2$estimate["shape"],
                       fit_gamma2$estimate["rate"],
                       fit_gamma2$aic),
    x = "step² (mm²·s⁻²)",
    y = "Densité"
  )

print(g_step2)
ggsave(file.path(figdir, "02_hist_gamma_step2.pdf"),
       g_step2, width = 8, height = 5.5)

# Planche de diagnostics pour step^2
pdf(file.path(figdir, "02_diag_gamma_step2.pdf"), width = 8, height = 8)
par(mfrow = c(2, 2))
denscomp(fit_gamma2, legendtext = "Gamma ajustée")
qqcomp(fit_gamma2, main = "Q–Q plot (step²)")
ppcomp(fit_gamma2, main = "P–P plot (step²)")
cdfcomp(fit_gamma2, legendtext = "Gamma ajustée")
dev.off()



## ─────  6 • Turning angle — stats & rose plot pour la fourmi id 243 ─────
library(circular)

## 1) Angles de la trajectoire la plus longue (déjà en radians)
ang <- na.omit(df_long$angle)

## 2) Objet circulaire : turning angles dans (-π, π]
y <- as.circular(
  ang,
  units  = "radians",
  modulo = "pi"   # pour des turning angles (−π, π]
)


## 3) Statistiques directionnelles
mu   <- mean.circular(y) |> as.numeric()   # moyenne directionnelle
Rbar <- rho.circular(y)  |> as.numeric()   # longueur du vecteur résultant

## 4) Figure
pdf(file.path(figdir, "02_rose_angle_id243.pdf"), width = 5, height = 5)
par(mar = c(1, 1, 1, 1))

## a) Tiges / points empilés sur le cercle
plot(
  y,
  stack    = TRUE,
  cex      = 0.7,
  tcl.text = 0.2
)

## b) Histogramme circulaire (rose) par-dessus
rose.diag(
  y,
  bins     = 30,
  col      = "grey",
  prop     = 1.8,
  add      = TRUE,
  cex      = 0.7,
  tcl.text = 0.2
)

## c) Flèche de la moyenne directionnelle
arrows(
  x0 = 0, y0 = 0,
  x1 = Rbar * cos(mu),
  y1 = Rbar * sin(mu),
  col    = "red",
  lwd    = 2,
  length = 0.12
)

dev.off()

## ───────────────────────── 7 • Vitesse vs turning angle ────────────────────
library(data.table)
library(ggplot2)
library(hexbin)

DT <- as.data.table(df_HMM)[!is.na(step) & !is.na(angle)]

# turning angle signé en degrés (-180, 180]
DT[, theta_deg := angle * 180 / pi]

# Corrélation de Spearman (approximation à cause des ex aequo)
spearman_turning <- cor.test(DT$theta_deg, DT$step,
                             method = "spearman",
                             exact  = FALSE)
rho_hat <- unname(spearman_turning$estimate)  # ≈ -0.005

p_turning <- ggplot(DT, aes(x = theta_deg, y = step)) +
  stat_binhex(bins = 120) +
  scale_fill_gradient(trans = "log10", name = "N") +
  geom_smooth(method = "gam",
              formula = y ~ s(x),
              se = FALSE,
              colour = "white",
              linewidth = 0.7) +
  labs(
    x = "turning angle (°)",
    y = "vitesse (mm·s\u207b\u00b9)",      # mm·s^{-1}
    title = "Vitesse vs turning angle",
    subtitle = bquote(rho[Spearman] == .(round(rho_hat, 3)) ~
                        "(p < 10"^{-5}*")")
  ) +
  theme_bw()

print(p_turning)

ggsave(file.path(figdir, "02_step_vs_turning_angle.pdf"),
       p_turning, width = 7, height = 4.5)

###############################################################################
#  FIN CHAPITRE 2
###############################################################################




