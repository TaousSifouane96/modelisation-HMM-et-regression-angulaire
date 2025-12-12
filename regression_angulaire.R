############################################################
# 0) LIBRAIRIES
############################################################
if (!requireNamespace("CircularRegression", quietly = TRUE)) {
  stop("Installe d'abord le package CircularRegression depuis GitHub.")
}

library(CircularRegression)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(grid)
library(gridExtra)
############################################################
# 1) CHARGER / DÉFINIR TES DONNÉES DE FOURMIS
############################################################
# Remet un angle dans (-pi, pi]
wrap_angle <- function(a) atan2(sin(a), cos(a))

###############################################################################
## 1 • Lecture du jeu de données + calcul des angles de virage
##    → adapte juste le file.choose() si tu veux mettre un chemin fixe
###############################################################################

library(dplyr)
library(tidyr)

df_all <- read.table(file.choose(), header = TRUE, sep = ",",
                     stringsAsFactors = FALSE) |>
  dplyr::rename(ID = id) |>
  dplyr::arrange(ID, t) |>
  dplyr::group_by(ID) |>
  dplyr::mutate(
    dx    = x - dplyr::lag(x),
    dy    = y - dplyr::lag(y),
    dir   = atan2(dy, dx),
    angle = wrap_angle(dir - dplyr::lag(dir))  # <--- ici on utilise wrap_angle
  ) |>
  dplyr::ungroup() |>
  tidyr::drop_na(angle) |>
  dplyr::select(ID, t, angle)   # on garde t (pour l’ordre) et angle

message("• Angles restants : ", nrow(df_all),
        " | individus : ", dplyr::n_distinct(df_all$ID))

cat("Structure des données utilisées pour la régression angulaire :\n")
str(df_all)

############################################################
# 2) FONCTION consensus (modèle de Rivest / Nicosia / Baillargeon)
############################################################

consensus <- function(formula, data, model = "simplified", weights = NULL, 
                      initbeta = NULL, control = list()) {
  
  call  <- mfcall <- match.call()
  model <- model[1]
  
  # extraire les données selon la formule
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  
  nobs     <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  nterms   <- length(nomterms)
  
  p <- if (model == "simplified") nterms - 1 else nterms
  nparam <- if (model == "simplified") p + 1 else p + 2
  
  paramname <- nomterms
  if (model == "complete") {
    paramname <- c(paramname, paste0("beta", p + 1))
  }
  
  # y = angle réponse
  y <- as.vector(mf[, 1])
  
  # variables explicatives
  noms <- strsplit(nomterms, split = ":")
  noms <- do.call(rbind, noms)
  
  if (model == "simplified") {
    x0   <- mf[, noms[1, 1]]
    noms <- noms[-1, , drop = FALSE]
  }
  
  matx <- as.matrix(mf[, noms[, 1], drop = FALSE])
  
  if (ncol(noms) == 1) {
    matz <- matrix(1, ncol = ncol(matx), nrow = nrow(matx))
  } else {
    matz <- as.matrix(mf[, noms[, 2], drop = FALSE])
    matz[, noms[, 2] == noms[, 1]] <- 1
  }
  
  weight <- rep(1, nobs) * (is.null(weights)) + (!is.null(weights)) * weights
  
  # log-vraisemblance
  LL <- function(param) {
    angleref <- if (model == "simplified") x0 else rep(param[p + 2], nobs)
    
    sinmu <- param[1] * sin(angleref) +
      (matz * sin(matx)) %*% param[2:(p + 1)]
    cosmu <- param[1] * cos(angleref) +
      (matz * cos(matx)) %*% param[2:(p + 1)]
    
    long <- as.vector(sqrt(sinmu^2 + cosmu^2))
    mui  <- as.vector(atan2(sinmu, cosmu))
    
    term1 <- param[1] * cos(y - angleref) +
      (matz * cos(y - matx)) %*% param[2:(p + 1)]
    
    LL_val <- sum(term1) - sum(log(besselI(long, 0, expon.scaled = FALSE)))
    
    list(LL = LL_val, long = long, mui = mui)
  }
  
  # mise à jour des paramètres
  paramUpdate <- function(paramk, long, mui) {
    angleref <- if (model == "simplified") x0 else rep(paramk[p + 2], nobs)
    matx0    <- cbind(angleref, matx)
    matz0    <- cbind(rep(1, nobs), matz)
    
    Along <- as.vector(
      besselI(long, 1, expon.scaled = FALSE) /
        besselI(long, 0, expon.scaled = FALSE)
    )
    
    matu <- matz0 * (cos(y - matx0) - cos(matx0 - mui) * Along)
    
    if (model == "complete") {
      matu <- cbind(
        matu,
        paramk[1] * sin(y - angleref) - sin(mui - angleref) * Along
      )
    }
    
    vecs <- colSums(matu)
    names(vecs) <- paramname
    
    Xc <- matz0 * cos(matx0 - mui)
    Xs <- matz0 * sin(matx0 - mui)
    
    if (model == "complete") {
      Xc <- cbind(Xc, paramk[1] * sin(mui - paramk[p + 2]))
      Xs <- cbind(Xs, paramk[1] * cos(mui - paramk[p + 2]))
    }
    
    Dc <- diag(1 - Along/long - Along^2, nrow = nobs, ncol = nobs)
    Ds <- diag(Along/long,              nrow = nobs, ncol = nobs)
    
    matI <- t(Xc) %*% Dc %*% Xc + t(Xs) %*% Ds %*% Xs
    colnames(matI) <- rownames(matI) <- paramname
    
    dparam  <- as.vector(solve(matI, vecs))
    paramk1 <- paramk + dparam
    
    list(paramk1 = paramk1, dparam = dparam, matu = matu, matI = matI)
  }
  
  # valeurs initiales
  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 1000 else control$pginit
    pg     <- round(pginit^(1/nparam))
    
    possparam <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, pg + 2)]),
                     p + 1)
    
    if (model == "complete") {
      possparam[[nparam]] <- seq(0, 2 * pi, length.out = pg + 2)[-c(1, pg + 2)]
    }
    
    possVal <- cbind(expand.grid(possparam), NA)
    colnames(possVal) <- c(paramname, "LL")
    
    maxLL <- function(param) LL(param = param)$LL
    possVal[, nparam + 1] <- apply(possVal[, 1:nparam], 1, maxLL)
    
    paramk <- unlist(
      possVal[which.max(possVal[, nparam + 1]), 1:nparam]
    )
  } else {
    if (length(initbeta) != nparam) {
      stop("initbeta doit être de longueur ", nparam)
    }
    paramk <- initbeta
  }
  
  calcul <- LL(param = paramk)
  maxLLk <- calcul$LL
  long   <- calcul$long
  mui    <- calcul$mui
  
  iter    <- 0
  iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 1e-06 else control$mindiff
  conv    <- FALSE
  
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 3)
  colnames(iter.detail) <- c(paramname, "maxLL", "iter", "nitersh")
  iter.detail[1, ] <- c(paramk, maxLLk, iter, iter.sh)
  
  while (!conv && iter <= maxiter) {
    maj     <- paramUpdate(paramk = paramk, long = long, mui = mui)
    paramk1 <- maj$paramk1
    dparam  <- maj$dparam
    
    calcul  <- LL(param = paramk1)
    maxLLk1 <- calcul$LL
    long    <- calcul$long
    mui     <- calcul$mui
    
    iter.sh <- 0
    while (maxLLk1 < maxLLk) {
      iter.sh <- iter.sh + 1
      paramk1 <- paramk + dparam / (2^iter.sh)
      calcul  <- LL(param = paramk1)
      maxLLk1 <- calcul$LL
      long    <- calcul$long
      mui     <- calcul$mui
      if (iter.sh >= maxiter) break
    }
    
    if (maxLLk1 < maxLLk) {
      conv <- FALSE
      warning("l'algorithme n'a pas convergé (log-vraisemblance a diminué)")
      break
    } else {
      conv <- if (maxLLk1 - maxLLk > mindiff) FALSE else TRUE
      paramk <- paramk1
      maxLLk <- maxLLk1
      iter   <- iter + 1
      iter.detail[iter + 1, ] <- c(paramk, maxLLk, iter, iter.sh)
    }
  }
  
  if (iter > maxiter + 1) {
    warning("l'algorithme n'a pas convergé (maxiter atteint)")
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }
  
  if (maxLLk == maxLLk1) {
    maj <- paramUpdate(paramk = paramk, long = long, mui = mui)
  }
  matu <- maj$matu
  matI <- maj$matI
  
  v1 <- solve(matI)
  
  mid <- matrix(0, ncol = nparam, nrow = nparam)
  for (i in 1:nobs) {
    mid <- mid + t(matu[i, , drop = FALSE]) %*% matu[i, , drop = FALSE]
  }
  v2 <- v1 %*% mid %*% v1
  
  paramb <- paramk[2:(p + 1)] / paramk[1]
  
  matDeriv <- rbind(
    -paramk[2:(p + 1)] / paramk[1]^2,
    diag(1/paramk[1], nrow = p, ncol = p)
  )
  
  vb  <- t(matDeriv) %*% v1[1:(p + 1), 1:(p + 1)] %*% matDeriv
  vb2 <- t(matDeriv) %*% v2[1:(p + 1), 1:(p + 1)] %*% matDeriv
  
  names(paramb) <- colnames(vb) <- rownames(vb) <-
    colnames(vb2) <- rownames(vb2) <- paramname[-1]
  
  zvalue <- abs(paramk) / sqrt(diag(v2))
  pval   <- round(2 * pnorm(zvalue, lower.tail = FALSE), 5)
  
  parameters <- cbind(paramk, sqrt(diag(v2)), zvalue, pval)
  colnames(parameters) <- c("estimate", "robust std", "z value", "P(|z|>.)")
  rownames(parameters) <- paramname
  
  zvaluebeta <- abs(paramb) / sqrt(diag(vb2))
  pbeta      <- round(2 * pnorm(zvaluebeta, lower.tail = FALSE), 5)
  
  parambeta <- cbind(paramb, sqrt(diag(vb2)), zvaluebeta, pbeta)
  colnames(parambeta) <- c("estimate", "stderr", "z value", "P(|z|>.)")
  rownames(parambeta) <- names(paramb)
  
  res <- sin(y - mui)
  autocorr <- acf(res, plot = FALSE)
  
  out <- list(
    MaxLL        = maxLLk,
    parameters   = parameters,
    varcov1      = v1,
    varcov2      = v2,
    parambeta    = parambeta,
    varcovbeta1  = vb,
    varcovbeta2  = vb2,
    autocorr     = autocorr,
    matx         = matx,
    matz         = matz,
    y            = y,
    long         = long,
    mui          = mui,
    res          = res,
    iter.detail  = iter.detail,
    call         = call
  )
  class(out) <- "consensus"
  out
}

plot.consensus <- function(x, ...) {
  cat("\nGoodness of Fit:\n")
  cat("\nCorrelation between |e_i| and 1/l_i :\n")
  print.default(
    cor(x = abs(x$y - x$mui), y = 1/x$long, method = "spearman"),
    print.gap = 2, quote = FALSE, right = TRUE, ...
  )
  cat("\n")
  gg_qq(fitted = x$mui, resid = sqrt(x$long) * (x$y - x$mui))
  invisible(x)
}

print.consensus <- function(x, ...) {
  cat("\nMaximum log-likelihood :", x$MaxLL, "\n")
  cat("\nKappa Parameters:\n")
  print.default(x$parameters, print.gap = 2, quote = FALSE, right = TRUE, ...)
  cat("\n")
  invisible(x)
}

############################################################
# 3) FONCTIONS DE DIAGNOSTIC
############################################################

gg_qq <- function(fitted, resid, distribution = "norm", ...,
                  line.estimate = NULL, conf = 0.95,
                  labels = NULL) {
  
  # scatter résidus vs valeurs ajustées
  df <- data.frame(predicted.mean = fitted, residual = resid)
  p1 <- ggplot(df, aes(x = predicted.mean, y = residual)) +
    geom_point(colour = "blue") +
    geom_abline(intercept = 0, slope = 0) +
    labs(x = "Predicted Mean", y = "Residual")
  
  # histogramme des résidus
  residual <- data.frame(resid = resid)
  gg <- ggplot(residual, aes(x = resid))
  range.residual <- range(residual$resid)
  gg <- gg + geom_histogram(
    binwidth = (range.residual[2] - range.residual[1]) / 10,
    colour = "black", aes(y = ..density.., fill = ..count..)
  )
  p2 <- gg + stat_function(
    fun  = dnorm,
    color = "red",
    args = list(
      mean = mean(residual$resid),
      sd   = sd(residual$resid)
    )
  )
  
  # QQ-plot
  x <- resid
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x   <- na.omit(x)
  ord <- order(x)
  n   <- length(x)
  P   <- ppoints(n)
  df2 <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if (is.null(line.estimate)) {
    Q.x <- quantile(df2$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b   <- diff(Q.x) / diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf) / 2)
  SE <- (coef[2] / d.function(df2$z)) * sqrt(P * (1 - P) / n)
  fit.value <- coef[1] + coef[2] * df2$z
  df2$upper <- fit.value + zz * SE
  df2$lower <- fit.value - zz * SE
  
  p3 <- ggplot(df2, aes(x = z, y = ord.x)) +
    geom_point(colour = "blue") +
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Theoretical Quantiles", y = "Residuals")
  
  # boxplot des résidus
  residual$boxplot <- as.factor(rep(1, length(residual$resid)))
  p4 <- ggplot(residual, aes(x = boxplot, y = resid)) +
    geom_boxplot(colour = "blue")
  
  multiplot(p1, p2, p3, p4, cols = 2)
  grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
               top = "Diagnostic plots for scaled residuals")
}

hist.resid <- function(resid) {
  residual <- data.frame(resid = resid)
  gg <- ggplot(residual, aes(x = resid))
  gg <- gg + geom_histogram(
    bins = 12, colour = "black",
    aes(y = ..density.., fill = ..count..)
  )
  gg <- gg + scale_fill_gradient("Count")
  gg <- gg + stat_function(
    fun   = dnorm,
    color = "red",
    args  = list(
      mean = mean(residual$resid),
      sd   = sd(residual$resid)
    )
  )
  gg
}

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols,
                     nrow = ceiling(numPlots / cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

############################################################
# 4) FONCTIONS SPÉCIFIQUES AU MODÈLE AUTO-RÉGRESSIF
############################################################

# construit les lags angle_tm1, angle_tm2, ..., angle_tmp
make_AR_for_p <- function(df_id, response = "angle", p) {
  dat <- df_id
  for (k in 1:p) {
    lag_name <- paste0(response, "_tm", k)
    dat[[lag_name]] <- dplyr::lag(dat[[response]], k)
  }
  lag_cols <- paste0(response, "_tm", 1:p)
  dat %>% tidyr::drop_na(dplyr::all_of(lag_cols))
}

# calcule L_p = moyenne de cos(y - mu)
compute_Lp_from_fit <- function(fit) {
  mean(cos(fit$y - fit$mui))
}

# ajuste AR(1)...AR(p_max) pour UNE fourmi (ID)
fit_AR_consensus_one_id <- function(df_all, id, p_max = 10) {
  df_id <- df_all %>%
    dplyr::filter(ID == id) %>%
    dplyr::arrange(row_number()) %>%
    tidyr::drop_na(angle)
  
  n_total <- nrow(df_id)
  cat("\n===================== ID =", id, "=====================\n")
  cat("Nombre de points :", n_total, "\n")
  
  fits   <- vector("list", p_max)
  n_vec  <- integer(p_max)
  Lp_vec <- rep(NA_real_, p_max)
  
  for (p in 1:p_max) {
    dat_p <- make_AR_for_p(df_id, response = "angle", p = p)
    n_vec[p] <- nrow(dat_p)
    
    if (n_vec[p] <= p + 1) {
      cat("  p =", p, ": trop peu d'observations, on saute.\n")
      next
    }
    
    rhs_terms <- paste0("angle_tm", 1:p)
    form_p <- as.formula(
      paste("angle ~", paste(rhs_terms, collapse = " + "))
    )
    
    fit_p <- try(
      consensus(formula = form_p,
                data    = dat_p,
                model   = "simplified"),
      silent = TRUE
    )
    
    if (inherits(fit_p, "try-error")) {
      cat("  p =", p, "→ échec de convergence.\n")
      fits[[p]] <- NULL
      Lp_vec[p] <- NA_real_
    } else {
      fits[[p]] <- fit_p
      Lp_vec[p] <- compute_Lp_from_fit(fit_p)
      cat("  p =", p, "→ L_p ≈", round(Lp_vec[p], 4), "\n")
    }
  }
  
  score_tbl <- tibble::tibble(
    p   = 1:p_max,
    n   = n_vec,
    L_p = Lp_vec
  )
  cat("\nTableau L_p pour ID =", id, ":\n")
  print(score_tbl, n = p_max)
  
  ok_idx <- which(!is.na(score_tbl$L_p))
  if (length(ok_idx) == 0L) {
    warning("Aucun modèle valide pour ID = ", id)
    return(list(
      ID       = id,
      n        = n_total,
      score    = score_tbl,
      best_p   = NA_integer_,
      best_Lp  = NA_real_,
      best_fit = NULL,
      rho1     = NA_real_
    ))
  }
  
  L_max <- max(score_tbl$L_p[ok_idx], na.rm = TRUE)
  best_row <- score_tbl %>%
    dplyr::filter(!is.na(L_p), abs(L_p - L_max) < 1e-10) %>%
    dplyr::arrange(p) %>%
    dplyr::slice(1)
  
  best_p   <- best_row$p
  best_Lp  <- best_row$L_p
  best_fit <- fits[[best_p]]
  
  rho1 <- NA_real_
  if (!is.null(best_fit)) {
    res_circ <- best_fit$res
    if (length(res_circ) >= 2) {
      rho1 <- cor(res_circ[-1], res_circ[-length(res_circ)])
    }
  }
  
  cat("\n>>> ID =", id, "\n")
  cat(">>> p* =", best_p,
      " (maximise L_p = 1/n Σ cos(y - μ))\n")
  cat("    L_p(p*) ≈", round(best_Lp, 4), "\n")
  cat("    Corrélation lag-1 des résidus (rho1) ≈",
      round(rho1, 4), "\n")
  
  list(
    ID       = id,
    n        = n_total,
    score    = score_tbl,
    best_p   = best_p,
    best_Lp  = best_Lp,
    best_fit = best_fit,
    rho1     = rho1
  )
}

############################################################
# 5) ANALYSE DES 10 TRAJECTOIRES LES PLUS LONGUES
############################################################

# ATTENTION : s'assurer que df_all est défini plus haut
# et contient au moins : ID, angle

# sélectionner les 10 fourmis ayant le plus de points
top10_ids <- df_all %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(n_points = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_points)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(ID)

cat("\nIDs des 10 trajectoires les plus longues :\n")
print(top10_ids)

p_max <- 10

# ajuster AR(1)...AR(10) pour chaque fourmi
res_list <- lapply(top10_ids, function(id) {
  fit_AR_consensus_one_id(df_all = df_all, id = id, p_max = p_max)
})

# résumé global
summary_top10 <- purrr::map_dfr(res_list, function(r) {
  tibble::tibble(
    ID       = r$ID,
    n_points = r$n,
    best_p   = r$best_p,
    best_Lp  = r$best_Lp,
    rho1     = r$rho1
  )
}) %>%
  dplyr::arrange(dplyr::desc(n_points))

cat("\nRésumé global pour les 10 trajectoires les plus longues :\n")
print(summary_top10)

# afficher print(best_fit) pour chaque trajectoire
for (r in res_list) {
  cat("=====================================================\n")
  cat("Résumé du modèle consensus AR(p*) pour ID =", r$ID, "\n")
  cat("n_points =", r$n,
      " | p* =", r$best_p,
      " | L_p(p*) ≈", round(r$best_Lp, 4),
      " | rho1 ≈", round(r$rho1, 4), "\n\n")
  
  if (!is.null(r$best_fit)) {
    print(r$best_fit)  # utilise print.consensus
  } else {
    cat("Aucun modèle ajusté pour cet ID.\n\n")
  }
}



for (i in seq_along(res_list)) {
  id_i   <- res_list[[i]]$ID
  fit_i  <- res_list[[i]]$best_fit
  
  cat("\nHistogramme des résidus pour ID =", id_i, "\n")
  print(hist.resid(fit_i$res))
}


# 1. Préparer les lags pour une fourmi
df_243 <- df_all %>%
  filter(ID == 243) %>%
  arrange(t) %>%
  mutate(
    angle_tm1  = lag(angle, 1),
    angle_tm2  = lag(angle, 2),
    angle_tm3  = lag(angle, 3),
    angle_tm4  = lag(angle, 4),
    angle_tm5  = lag(angle, 5),
    angle_tm6  = lag(angle, 6),
    angle_tm7  = lag(angle, 7),
    angle_tm8  = lag(angle, 8),
    angle_tm9  = lag(angle, 9),
    angle_tm10 = lag(angle, 10)
  ) %>%
  drop_na()   # enlève les premières lignes sans lags

# 2. Ajuster le modèle consensus AR(10)
fit_243 <- consensus(
  angle ~ angle_tm1 + angle_tm2 + angle_tm3 + angle_tm4 + angle_tm5 +
    angle_tm6 + angle_tm7 + angle_tm8 + angle_tm9 + angle_tm10,
  data  = df_243,
  model = "simplified"
)

# 3. Résumé des paramètres
print(fit_243)

plot(fit_243)

dir.create("figures_diag", showWarnings = FALSE)

for (i in seq_along(res_list)) {
  id_i  <- res_list[[i]]$ID
  fit_i <- res_list[[i]]$best_fit
  
  png(file.path("figures_diag",
                paste0("diag_residus_ID", id_i, ".png")),
      width = 1000, height = 500, res = 120)
  plot(fit_i)      # les 4 graphes
  dev.off()
}
