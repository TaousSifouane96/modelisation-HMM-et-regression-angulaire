###############################################################################
#  CHAPITRE 0 — LIBRAIRIES
###############################################################################
library(dplyr)
library(tibble)
library(purrr)
library(momentuHMM)
library(fitdistrplus)
library(circular)
library(CircStats)
library(ggplot2)
library(moments)         # skewness()
library(progressr)       # progress bars
library(data.table)      # étude inter-parcours
handlers(global = TRUE)

###############################################################################
#  CHAPITRE 1 — CHARGEMENT & PRÉ-TRAITEMENT
###############################################################################

## 1) S'assurer que df_HMM contient bien la colonne step_original ----------
##    Si elle n'existe pas, on la reconstruit à partir de step

if (!"step_original" %in% names(df_HMM)) {
  # on suppose que la colonne 'step' contient les pas bruts
  df_HMM$step_original <- df_HMM$step
}

## 2) Choisir automatiquement la trajectoire la plus longue -----------------

library(data.table)
df_HMM_dt <- as.data.table(df_HMM)

long_ID <- df_HMM_dt[, .N, by = ID][order(-N)][1, ID]
cat("Trajectoire la plus longue : ID =", long_ID, "\n")

## 3) Extraire cette trajectoire et nettoyer ---------------------------------

df_HMM_1 <- df_HMM_dt[ID == long_ID &
                        is.finite(step_original) &
                        step_original > 0]

cat("Nombre de pas dans df_HMM_1 :", nrow(df_HMM_1), "\n")

if (nrow(df_HMM_1) <= 1) {
  stop("df_HMM_1 ne contient pas assez de pas (>1) pour ajuster une Gamma.")
}

## 4) Construire les quatre vecteurs de vitesse -----------------------------

step      <- df_HMM_1$step_original                 # pas brut   d
step2     <- step^2                                 # d²
step_p15  <- (1 + step)^1.5 - 1                     # (1+d)^{1.5}-1
step_p2   <- (1 + step)^2  - 1                      # (1+d)^{2}-1

## Vérifications rapides
cat("Longueur step     :", length(step),    "\n")
cat("Longueur step2    :", length(step2),   "\n")
cat("Longueur step_p15 :", length(step_p15),"\n")
cat("Longueur step_p2  :", length(step_p2), "\n")
############################################################################
# 2) Ajustement Gamma par ML ----------------------------------------------
############################################################################
library(fitdistrplus)

fit_step     <- fitdist(step,      "gamma", method = "mle")
fit_step2    <- fitdist(step2,     "gamma", method = "mle")
fit_step_p15 <- fitdist(step_p15,  "gamma", method = "mle")
fit_step_p2  <- fitdist(step_p2,   "gamma", method = "mle")

summary(fit_step)
summary(fit_step2)
summary(fit_step_p15)
summary(fit_step_p2)
############################################################################
############################################################################
# 0) S’assurer qu’aucun device n’est encore actif --------------------------
############################################################################
while (!is.null(dev.list())) dev.off()   # ferme tout

############################################################################
# 1) AFFICHAGE INTERACTIF (Plots RStudio) ---------------------------------
############################################################################
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(fit_step,      histo = TRUE, dens = TRUE,
     qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = "Gamma — pas brut (d)",
      xlab = "d (mm/s)", ylab = "Densité")

plot(fit_step2,     histo = TRUE, dens = TRUE,
     qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = "Gamma — d\u00B2",
      xlab = expression(d^2), ylab = "Densité")

plot(fit_step_p15,  histo = TRUE, dens = TRUE,
     qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = expression(Gamma~~(1+d)^{1.5}-1),
      xlab = "(1+d)^{1.5}-1", ylab = "Densité")

plot(fit_step_p2,   histo = TRUE, dens = TRUE,
     qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = expression(Gamma~~(1+d)^{2}-1),
      xlab = "(1+d)^{2}-1", ylab = "Densité")

par(mfrow = c(1, 1))   # reset

############################################################################
# 2) (OPTION) ENREGISTRER LA MÊME FIGURE -----------------------------------
############################################################################
png("figures/compare_gamma_base.png", width = 1000, height = 800)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(fit_step,      histo = TRUE, dens = TRUE, qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = "Gamma — pas brut (d)",      xlab = "d (mm/s)",      ylab = "Densité")

plot(fit_step2,     histo = TRUE, dens = TRUE, qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = "Gamma — d\u00B2",             xlab = expression(d^2), ylab = "Densité")

plot(fit_step_p15,  histo = TRUE, dens = TRUE, qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = expression(Gamma~~(1+d)^{1.5}-1),
      xlab = "(1+d)^{1.5}-1", ylab = "Densité")

plot(fit_step_p2,   histo = TRUE, dens = TRUE, qq = FALSE, pp = FALSE, cdf = FALSE)
title(main = expression(Gamma~~(1+d)^{2}-1),
      xlab = "(1+d)^{2}-1",  ylab = "Densité")

dev.off()                        # IMPORTANT : termine le fichier


####

###############################################################################
#  CHAPITRE 3 — FONCTION HMM (SRLR)
###############################################################################
HMM_SRLR <- function(
    data, nbStates = 2,
    dist = list(step="gamma", angle="vm"),
    n_short = 50, maxiter_short = 20, maxiter_long = 1000,
    eps_diag = 1e-3, jitter_sd = 0.4, nll_tol = 2, seed = 123, ...
){
  if(!is.null(seed)) set.seed(seed)
  gfit <- MASS::fitdistr(data$step, densfun="gamma")$estimate
  vfit <- circular::mle.vonmises(circular(data$angle))
  base <- list(
    step  = c(gfit["shape"], 1/gfit["rate"]),
    angle = vfit$kappa
  )
  scalers <- seq(0.5, 1.5, length.out = nbStates)
  base$step  <- c(base$step[1]*scalers, base$step[2]*scalers)
  base$angle <- base$angle*scalers
  jitter <- function(x) x*rlnorm(length(x),-0.5*jitter_sd^2,jitter_sd)
  gen0   <- function() list(step=jitter(base$step),angle=jitter(base$angle))
  short_tbl <- progressr::with_progress({
    p <- progressr::progressor(n_short)
    map_dfr(1:n_short, function(i){
      p()
      fi <- try(fitHMM(data, nbStates, dist, Par0=gen0(),
                       control=list(maxit=maxiter_short), angleMean=FALSE, ...),
                silent=TRUE)
      ok <- !inherits(fi,"try-error") &&
        all(diag(fi$mle$gamma)>eps_diag & diag(fi$mle$gamma)<1-eps_diag)
      tibble(run=i,valid=ok,nll=if(ok) fi$mod$minimum else Inf,fit=list(if(ok)fi))
    })
  })
  best <- short_tbl %>% filter(valid,nll<=min(nll)+nll_tol) %>% slice_min(nll, n=1)
  Par0 <- list(
    step  = c(best$fit[[1]]$mle$step["mean",],best$fit[[1]]$mle$step["sd",]),
    angle = best$fit[[1]]$mle$angle["concentration",]
  )
  long <- fitHMM(data, nbStates, dist, Par0=Par0,
                 control=list(maxit=maxiter_long), angleMean=FALSE, ...)
  list(long_run=long)
}

###############################################################################
#  CHAPITRE 4 — TRANSFORMATIONS & MODÈLES
###############################################################################
crit_table <- tibble(transfo=character(),n_states=integer(),
                     logLik=numeric(),AIC=numeric(),BIC=numeric(),skew_res=numeric())

run_block <- function(df, label, n_states){
  key <- paste0(label,"_",n_states)
  m   <- HMM_SRLR(df, nbStates=n_states)
  png(paste0("pseudo_residuals_",key,".png"),800,600); plotPR(m$long_run); dev.off()
  ll  <- -m$long_run$mod$minimum
  k   <- length(m$long_run$mod$estimate)
  n   <- nrow(df)
  crit_table <<- add_row(crit_table,
                         transfo=label, n_states=n_states,
                         logLik=ll, AIC=-2*ll+2*k, BIC=-2*ll+log(n)*k,
                         skew_res=skewness(pseudoRes(m$long_run)$stepRes,na.rm=TRUE))
  model_store[[key]] <<- m
}

model_store <- list()

## 4.1 step brut
df_raw <- mutate(df_HMM_1, step=step_original)
run_block(df_raw,"step",2); run_block(df_raw,"step",3)

## 4.2 step^2
df_sq  <- mutate(df_HMM_1, step=step_original^2)
run_block(df_sq,"step2",2); run_block(df_sq,"step2",3)

## 4.3 (1+step)^1.5 – 1
df_p15 <- mutate(df_HMM_1, step=(1+step_original)^1.5-1)
run_block(df_p15,"step_p15",2); run_block(df_p15,"step_p15",3)

## 4.4 (1+step)^2 – 1
df_p2  <- mutate(df_HMM_1, step=(1+step_original)^2-1)
run_block(df_p2,"step_p2",2);  run_block(df_p2,"step_p2",3)

print(crit_table)
write.csv(crit_table,"table_criteres_HMM.csv",row.names=FALSE)
################
#########
###############################################################################
# 2) Résumé + figures pour chaque modèle -------------------------------------
###############################################################################
dir.create("figures", showWarnings = FALSE)

imap(model_store, ~{
  fit <- .x$long_run           # objet fitHMM
  nom <- .y                    # ex. "step_2", "step2_3", …
  
  ## 2-a. résumé paramétrique
  capture.output(
    summary(fit),
    file = file.path("figures", sprintf("summary_%s.txt", nom))
  )
  
  ## 2-b. plot() général
  png(file.path("figures", sprintf("plot_%s.png", nom)),
      width = 1800, height = 1200, res = 200)
  plot(fit, ask = FALSE)
  dev.off()
  
  ## 2-c. diagnostics pseudo-résidus
  png(file.path("figures", sprintf("plotPR_%s.png", nom)),
      width = 1800, height = 1200, res = 200)
  plotPR(fit)                  # plus d’arguments ask/plot
  dev.off()
})


message("\n✅  Tableau 'crit_table_complet.csv' + figures dans le dossier 'figures/'")
#  RÉSUMÉ « console » complet pour chaque modèle
###############################################################################

imap(model_store, ~{
  cat("\n============================================\n")
  cat(" Résumé général du modèle :", .y, "\n")
  cat("============================================\n\n")
  
  print(.x$long_run)   
})
################
###############
###############################################################################
#  CHAPITRE 5 — PARAMÈTRES DU MODÈLE FINAL (exemple : (1+step)^2-1, 3 états)
###############################################################################
fit_final <- model_store[["step_p2_3"]]$long_run
plot(fit_final)


gamma_from_mean_sd <- function(mu,sd){data.frame(mean=mu,sd=sd,
                                                 shape=mu^2/sd^2,scale=sd^2/mu)}
Rbar_from_kappa <- function(k) besselI(k,1)/besselI(k,0)
stationary_dist <- function(G){v<-Re(eigen(t(G))$vectors[,1]);v/sum(v)}
mean_run_length <- function(G) 1/(1-diag(G))

pars_step  <- fit_final$mle$step
pars_angle <- fit_final$mle$angle
Gamma      <- fit_final$mle$gamma

tab_step  <- gamma_from_mean_sd(pars_step["mean",],pars_step["sd",]) %>%
  mutate(state=factor(seq_along(mean))) %>% relocate(state)
tab_angle <- data.frame(state=factor(seq_along(pars_angle["concentration",])),
                        mean=pars_angle["mean",],
                        kappa=pars_angle["concentration",],
                        Rbar=Rbar_from_kappa(pars_angle["concentration",]))
tab_trans <- as.data.frame(Gamma) %>% mutate(from=factor(1:nrow(Gamma))) %>% relocate(from)
tab_pi    <- data.frame(state=factor(1:nrow(Gamma)),
                        pi_stationnaire=stationary_dist(Gamma),
                        duree_moyenne=mean_run_length(Gamma))
print(tab_angle)
print(tab_pi)
print(tab_step)
print(tab_trans)

write.csv(tab_step ,"params_step_gamma.csv",row.names=FALSE)
write.csv(tab_angle,"params_angle_vm.csv",row.names=FALSE)
write.csv(tab_trans,"transition_matrix.csv",row.names=FALSE)
write.csv(tab_pi   ,"stationary_durations.csv",row.names=FALSE)

###############################################################################
#  CHAPITRE 6 — ÉTUDE INTER-PARCOURS (version propre avec test LRT sur κ)
###############################################################################

library(data.table)
library(dplyr)
library(purrr)
library(circular)
library(readr)
library(ggplot2)

## 6.1 Sélection des 10 trajectoires les plus longues -------------------------

setDT(df_HMM)

top_ids <- df_HMM[, .N, by = ID][order(-N)][1:10, ID]

df_HMM_top <- df_HMM[ID %in% top_ids]

## Transformation de la vitesse selon le modèle retenu : (1+step)^2 - 1
df_HMM_top[, step := (1 + step)^2 - 1]
df_HMM_top <- df_HMM_top[is.finite(step) & is.finite(angle) & step > 0]

## 6.2 Ajustement d’un HMM (3 états) pour chaque fourmi -----------------------

hmm_models_par_id <- list()

for(id in unique(df_HMM_top$ID)){
  cat("Ajustement ID :", id, "\n")
  
  df_id <- df_HMM_top[ID == id]
  class(df_id) <- c("momentuHMMData", "data.frame")
  
  if (nrow(df_id) < 50) next
  
  m <- try(HMM_SRLR(df_id, nbStates = 3), silent = TRUE)
  
  if (!inherits(m, "try-error")) {
    hmm_models_par_id[[as.character(id)]] <- m
  }
}

## Vérifier que des modèles existent :
print(hmm_models_par_id |> names())
cat("Nombre de modèles ajustés :", length(hmm_models_par_id), "\n")

## 6.3 Fonctions utilitaires --------------------------------------------------

# π stationnaire d’une matrice de transition G
pi_stationnaire <- function(G){
  v <- Re(eigen(t(G))$vectors[, 1])
  v / sum(v)
}

# ordonner les états par vitesse moyenne croissante
ordre_etats <- function(long_run_fit){
  order(long_run_fit$mle$step["mean", ])
}

## 6.4 Tables descriptives : kappa, Gamma, π, durées --------------------------

tab_kappa <- list()
tab_step  <- list()
tab_pi    <- list()
tab_G     <- list()

for(id in names(hmm_models_par_id)){
  
  fit_long <- hmm_models_par_id[[id]]$long_run
  ord      <- ordre_etats(fit_long)
  
  ## paramètres angle (von Mises)
  kap <- fit_long$mle$angle["concentration", ord]
  Rb  <- Rbar_from_kappa(kap)
  
  tab_kappa[[id]] <- data.frame(
    ID    = as.numeric(id),
    state = 1:3,
    Kappa = kap,
    Rbar  = Rb
  )
  
  ## paramètres vitesse (Gamma) — à partir des moyennes et écarts-types
  m  <- fit_long$mle$step["mean", ord]
  sd <- fit_long$mle$step["sd",   ord]
  
  tab_step[[id]] <- gamma_from_mean_sd(m, sd) %>%
    mutate(ID = as.numeric(id),
           state = 1:3) %>%
    dplyr::relocate(ID, state)
  
  ## matrice de transition et π stationnaire
  G      <- fit_long$mle$gamma[ord, ord]
  pi_hat <- pi_stationnaire(G)
  dur    <- mean_run_length(G)
  
  tab_pi[[id]] <- data.frame(
    ID              = as.numeric(id),
    state           = 1:3,
    pi_stationnaire = pi_hat,
    duree_moy       = dur
  )
  
  tab_G[[id]] <- cbind(
    data.frame(ID = as.numeric(id), from = 1:3),
    setNames(as.data.frame(G), paste0("to_", 1:3))
  )
}

tab_kappa <- bind_rows(tab_kappa)
tab_step  <- bind_rows(tab_step)
tab_pi    <- bind_rows(tab_pi)
tab_G     <- bind_rows(tab_G)

## 6.5 TEST SUR LES KAPPA : TEST LRT (H0 : κ commun entre fourmis) -----------

# log-vraisemblance von Mises (mu = 0 fixé)
vm_loglik_mu0 <- function(alpha, kappa) {
  sum(dvonmises(
    circular(alpha),
    mu    = circular(0),
    kappa = kappa,
    log   = TRUE
  ))
}

# Test LRT d'homogénéité des κ pour un état donné
# H0 : même κ pour toutes les fourmis
# H1 : κ_i spécifique à chaque fourmi
test_kappa_state <- function(df_state) {
  # df_state doit contenir : ID, angle (radians), state
  
  df_state <- df_state %>%
    dplyr::filter(!is.na(angle), !is.na(ID))
  
  by_id <- split(df_state$angle, df_state$ID)
  M     <- length(by_id)        # nombre de fourmis
  
  if (M < 2) {
    stop("Il faut au moins deux fourmis pour tester l'homogénéité.")
  }
  
  ## H1 : un κ_i par fourmi (mu = 0 fixé)
  fit_H1 <- lapply(by_id, function(a)
    mle.vonmises(circular(a), mu = circular(0)))
  
  kappa_H1 <- sapply(fit_H1, function(f) f$kappa)
  names(kappa_H1) <- names(by_id)
  
  ll_H1 <- sum(mapply(function(a, k) vm_loglik_mu0(a, k),
                      by_id, kappa_H1))
  
  ## H0 : un κ commun (mu = 0 fixé)
  fit_H0   <- mle.vonmises(circular(df_state$angle),
                           mu = circular(0))
  kappa_H0 <- fit_H0$kappa
  ll_H0    <- vm_loglik_mu0(df_state$angle, kappa_H0)
  
  ## Statistique LRT
  LRT <- 2 * (ll_H1 - ll_H0)
  df  <- M - 1
  p   <- 1 - pchisq(LRT, df)
  
  list(
    LRT       = LRT,
    df        = df,
    p.value   = p,
    kappa_H0  = kappa_H0,
    kappa_H1  = kappa_H1
  )
}

## 6.6 Construction du tableau (ID, state, angle) pour les tests --------------

angles_states_list <- vector("list", length(hmm_models_par_id))
names(angles_states_list) <- names(hmm_models_par_id)

i <- 1
for(id in names(hmm_models_par_id)){
  fit_id <- hmm_models_par_id[[id]]$long_run
  
  # données brutes pour cet ID
  df_id <- df_HMM_top[ID == as.numeric(id)]
  
  # états décodés (1,2,3) pour cet ID
  z_id  <- viterbi(fit_id)
  
  # sécurité : on tronque à la même longueur si besoin
  n <- min(length(z_id), nrow(df_id))
  
  angles_states_list[[id]] <- data.frame(
    ID    = as.numeric(id),
    state = z_id[1:n],
    angle = df_id$angle[1:n]   # turning angle en radians
  )
}

angles_states <- dplyr::bind_rows(angles_states_list)
angles_states <- dplyr::as_tibble(angles_states)

## Vérification rapide
print(colnames(angles_states))
str(angles_states)
dplyr::count(angles_states, state)

## 6.7 Application du test LRT pour chaque état -------------------------------

res_tests_kappa <- angles_states %>%
  dplyr::group_by(state) %>%
  dplyr::group_modify(~{
    res <- test_kappa_state(.x)
    tibble(
      state    = unique(.x$state),
      LRT      = res$LRT,
      df       = res$df,
      p_value  = res$p.value,
      kappa_H0 = res$kappa_H0
    )
  }) %>%
  dplyr::ungroup()

print(res_tests_kappa)
write.csv(res_tests_kappa, "tests_kappa_LRT.csv", row.names = FALSE)

## 6.8 Sauvegarde des tableaux descriptifs ------------------------------------

fwrite(tab_kappa, "kappa_by_id.csv")
fwrite(tab_step,  "gamma_step_by_id.csv")
fwrite(tab_pi,    "stationary_pi_dur_by_id.csv")
fwrite(tab_G,     "gamma_matrix_by_id.csv")

## 6.9 Boxplots des kappa par état -------------------------------------------

kapp <- tab_kappa   # déjà en mémoire

p_kappa <- ggplot(kapp, aes(x = factor(state), y = Kappa)) +
  geom_boxplot(fill = "grey90") +
  geom_jitter(width = .2, height = 0,
              size = 1.4, alpha = .7) +
  labs(x = "État", y = "Concentration κ") +
  theme_minimal()

print(p_kappa)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/box_kappa_states.png",
       plot   = p_kappa,
       width  = 5,
       height = 4)
