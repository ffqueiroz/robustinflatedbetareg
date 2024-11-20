## Functions to fit robust inflated beta regression

# This function fit two robust estimators for inflated beta regression models (M-LSE and M-LME)
# The tuning constants are selected via the data-driven algorithm.

fit.betainflated <- function(y, S, X, Z, linkmu="logit", 
                             linkphi="log", linkvartheta = "logit"){
  
  cl <- match.call()
  
  if ((min(y) == 0 & max(y) == 1))
    stop("Invalid dependent variable, all observations must be in [0, 1) or (0,1].")
  
  index0 <- y != 0 & y != 1 # y \in (0,1)
  
  y.1 <- y[index0]
  X.1 <- as.matrix(add_intercept(X[index0, ]))
  Z.1 <- as.matrix(add_intercept(Z[index0, ]))
  S <- as.matrix(add_intercept(S))
  
  y.2 <- ifelse(index0, 0, 1) # I(y=0) or I(y=1)
  
  fit.LSMLE <- robustbetareg:::LSMLE.fit(y.1, X.1, Z.1, link = linkmu, link.phi = linkphi)
  fit.LMDPDE <- robustbetareg:::LMDPDE.fit(y.1, X.1, Z.1, link = linkmu, link.phi = linkphi)  
  fit.MDPDE.disc <- MDPDE_BIN(y.2, S, linkvartheta = linkvartheta)
  
  Upsilon.M_LSE <- c(fit.MDPDE.disc$kappa, fit.LSMLE$coefficients$mean, fit.LSMLE$coefficients$precision)
  Upsilon.M_LME <- c(fit.MDPDE.disc$kappa, fit.LMDPDE$coefficients$mean, fit.LMDPDE$coefficients$precision)
  
  
  X <- as.matrix(add_intercept(X))
  Z <- as.matrix(add_intercept(Z))
  
  se.M_LSE <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_LSE, "M_LSE", 
                                                         fit.LSMLE$Tuning, S = S, X = X, Z = Z, y = y))
  se.M_LME <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_LME, "M_LME", fit.LMDPDE$Tuning, S = S, X = X, Z = Z, y = y))
  
  # Coefficients M-LME
  cf.MLME <- cbind(Upsilon.M_LME, se.M_LME, Upsilon.M_LME/se.M_LME, 2 * pnorm(-abs(Upsilon.M_LME/se.M_LME)))
  cf.MLSE <- cbind(Upsilon.M_LSE, se.M_LSE, Upsilon.M_LSE/se.M_LSE, 2 * pnorm(-abs(Upsilon.M_LSE/se.M_LSE)))
  
  p0 <- ncol(S)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  
  colnames(cf.MLME) <- colnames(cf.MLSE) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  row.names(cf.MLME) <- row.names(cf.MLSE) <- c(colnames(S), row.names(cf.MLSE)[-(1:p0)])
  
  
  cf.MLME <- list(probability = cf.MLME[seq.int(length.out = p0), , drop = FALSE],
             mean = cf.MLME[seq.int(length.out = p1) + p0, , drop = FALSE],
             precision = cf.MLME[seq.int(length.out = p2) + p0 + p1,  , drop = FALSE])
  
  cf.MLSE <- list(probability = cf.MLSE[seq.int(length.out = p0), , drop = FALSE],
                  mean = cf.MLSE[seq.int(length.out = p1) + p0, , drop = FALSE],
                  precision = cf.MLSE[seq.int(length.out = p2) + p0 + p1, , drop = FALSE])
  
  object <- NULL
  
  c <- ifelse(any(y == 0), 0, ifelse(any(y == 1), 1, NA))
  
  object$coefficients$MLME <- cf.MLME
  object$coefficients$MLSE <- cf.MLSE
  
  object$tuning$discrete <- fit.MDPDE.disc$alpha_value
  object$tuning$LSMLE <- fit.LSMLE$Tuning
  object$tuning$LMDPDE <- fit.LMDPDE$Tuning
  
  cat("\nCall:", deparse(cl, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  digits <- 4
  
  cat("\n\n Output for the M-LSE \n")
  
  if(NROW(object$coefficients$MLSE$mean)) {
    cat(paste("\nCoefficients (mean model with ", linkmu, " link):\n", sep = ""))
    printCoefmat(object$coefficients$MLSE$mean, digits = digits, signif.legend = FALSE)
  } else cat("\nNo coefficients (in mean model)\n")
  
  if(NROW(object$coefficients$MLSE$precision)) {
    cat(paste("\nSigma coefficients (precision model with ", linkphi, " link):\n", sep = ""))
    printCoefmat(object$coefficients$MLSE$precision, digits = digits, signif.legend = FALSE)
  } else cat("\nNo coefficients (in precision model)\n")
  
  if(NROW(object$coefficients$MLSE$probability)) {
    cat(paste("\nVartheta coefficients (Probability of ",c ," model with ", linkvartheta, " link):\n", sep = ""))
    printCoefmat(object$coefficients$MLSE$probability, digits = digits, signif.legend = FALSE)
  } else cat("\nNo coefficients (in probability of ",c ," model)\n")
  
  
  aux <- object$coefficients$MLSE[c("mean", "precision", "probability")]
  
  if(getOption("show.signif.stars") & any(do.call("rbind", aux)[, 4L] < 0.1))
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  
  cat("\nTuninf constant:", object$tuning$discrete, "(discrete part) and ", object$tuning$LSMLE , "(continuous part).")
  
  invisible(object)
}



