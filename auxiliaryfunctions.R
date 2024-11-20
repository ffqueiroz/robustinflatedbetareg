add_intercept <- function(data) {
  if (!all(data[, 1] == 1)) {
    data <- cbind(Intercept = 1, data)
  }
  return(data)
}

make.dalpha.deta <- function(linkstr){
  switch(linkstr, "logit" = {
    logit_link <- make.link("logit")
    function(eta) logit_link$mu.eta(eta) * (1 - 2 * logit_link$linkinv(eta))
  },
  "probit" = function(eta) -eta * pmax(dnorm(eta), .Machine$double.eps),
  "cauchit" = function(eta) -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps),
  "cloglog" = function(eta) pmax((1 - exp(eta)) * exp(eta - exp(eta)), .Machine$double.eps),
  "loglog" = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
  "log" = function(eta) pmax(exp(eta), .Machine$double.eps),
  "sqrt" = function(eta) rep.int(2, length(eta)),
  "1/alpha^2" = function(eta) 3/(4 * eta^2.5),
  "inverse" = function(eta) 2/(eta^3)
  )
}

MDPDE_BIN <- function(y, S, alphaoptimal=TRUE, alpha0=0, m=3, L=0.02, alphamax=1, spac=0.05, method="BFGS", maxit=200, linkvartheta = "logit", weights=FALSE){
  #****Packages required****#
  if(!suppressWarnings(require(matlib))){suppressWarnings(install.packages("matlib"));suppressWarnings(library(matlib))}#used for Ginv function
  if(!suppressWarnings(require(BBmisc))){suppressWarnings(install.packages("BBmisc"));suppressWarnings(library(BBmisc))} #used for is.error function
  if(!suppressWarnings(require(robust))){suppressWarnings(install.packages("robust"));suppressWarnings(library(robust))} #used for is.error function
  
  #**********link functions**********#
  
  vartheta_link_functions <- function(kappa, S, linkvartheta="logit"){
    eta <- as.vector(S%*%kappa)	 
    linkstr <- linkvartheta
    if (linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      linkobj$dalpha.deta <- make.dalpha.deta(linkstr)
    } else {
      linkobj <- structure(list(linkfun = function(alpha) -log(-log(alpha)),
                                linkinv = function(eta){
                                  pmax(pmin(exp(-exp(-eta)),
                                            1 - .Machine$double.eps), .Machine$double.eps)
                                },
                                mu.eta = function(eta) {
                                  eta <- pmin(eta, 700)
                                  pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                                },
                                dalpha.deta = function(eta) pmax(exp(-exp(-eta) -
                                                                       eta) * expm1(-eta), .Machine$double.eps),
                                valideta = function(eta) TRUE, name = "loglog"),
                           class = "link-glm")
    }
    vartheta <- linkobj$linkinv(eta)
    T_1 <- diag(linkobj$mu.eta(eta))
    
    results <- list(vartheta= vartheta, T_1 = T_1)
    return(results)
  }#ends mean link function
  
  rho <- function(t, alpha_const){
    (1 + (1/alpha_const))*(1 - exp(-alpha_const*t))
  }
  Gt <- function(t, alpha_const){
    t^(alpha_const + 1)
  }
  
  #*************Function to be minimized**********************#
  log_liksurrogate <- function(kappa) 
  {
    kk1 <- ncol(S)
    vartheta <- vartheta_link_functions(kappa = kappa, S=S, linkvartheta=linkvartheta)$vartheta
    
    deviance <- suppressWarnings(-log(dbinom(y, 1, vartheta)))
    phi      <- suppressWarnings(rho(deviance, alpha_const) + Gt(vartheta, alpha_const) + Gt(1 - vartheta, alpha_const))
    ll       <- suppressWarnings(phi)
    if (any(!is.finite(ll)))
      NaN
    else sum(ll)
  }
  
  #********************Score-type function*******************#
  
  psi <- function(t, alpha_const){
    (alpha_const + 1)*exp(-alpha_const*t)
  }
  dpsi <- function(t, alpha_const){
    -alpha_const*(alpha_const + 1)*exp(-alpha_const*t)
  }
  
  Score <- function(kappa) { 
    avScore = numeric(0) 
    vartheta  <- as.vector(vartheta_link_functions(kappa = kappa,  S=S, linkvartheta=linkvartheta)$vartheta)
    T_1 <- diag(vartheta_link_functions(kappa = kappa, S=S, linkvartheta=linkvartheta)$T_1)
    
    
    Mt         <- (1 - vartheta)*psi(-log(vartheta), alpha_const) + vartheta*psi(-log(1 - vartheta), alpha_const)
    dMt.dalpha <- psi(-log(1 - vartheta), alpha_const) - psi(-log(vartheta), alpha_const) +
      (1/(vartheta*(1 - vartheta)))*((vartheta^2)*dpsi(-log(1 - vartheta), alpha_const) -
                                       ((1 - vartheta)^2)*dpsi(-log(vartheta), alpha_const))
    
    Psi <- as.vector((vartheta - y)*Mt*T_1/(vartheta*(1-vartheta)))
    
    #****vector type Score*****#
    avScore <- t(S)%*%Psi
    return(avScore)
  }#ends score-type function
  
  #********************Covariance matrix****************#
  V_matrix <- function(kappa) {  
    
    vartheta <- vartheta_link_functions(kappa = kappa,  S=S, linkvartheta=linkvartheta)$vartheta
    t_1 <- vartheta_link_functions(kappa = kappa,  S=S, linkvartheta=linkvartheta)$T_1
    
    M <- (alpha_const + 1)*diag((1-vartheta)*exp(alpha_const*log(vartheta)) + vartheta*exp(alpha_const*log(1 - vartheta)))
    vartheta_star <- diag(1/(vartheta*(1-vartheta)))
    
    
    A <- t(S)%*%M%*%vartheta_star%*%(t_1^2)%*%S
    B <- t(S)%*%(M^2)%*%vartheta_star%*%(t_1^2)%*%S
    
    
    Vq <- tryCatch(solve(A)%*%(B)%*%t(solve(A)), error=function(e) {e})  #asymptotic covariance matrix
    if(is.error(Vq)){
      Vq <- Ginv(A)%*%(B)%*%t(Ginv(A))
    }
    return(Vq)
  }#ends Covariance matrix function
  
  #**********************Starting Values function*****************#
  Starting_values <- function(y, S, linkvartheta){
    n <- nrow(S)
    fit_mle <- glm(y ~ S[,-1], family="binomial")
    initials_mle <- thetaStart <- as.numeric(coef(fit_mle))
    
    thetaStart <- suppressWarnings(tryCatch(as.numeric(robust::glmRob(y ~ S[,-1], family = "binomial",  method = "mallows")$coef)
                                            , error=function(e) {e}))#maximization 
    if(is.error(thetaStart)){
      thetaStart <- as.numeric(fit_mle$coefficients)
    }
    
    results <- list(thetaStart=thetaStart, MLEfit= fit_mle)
    return(results)
  }#ends function of starting values
  
  #----------------------------------------------------------------------------------------------------------------#
  #------------------------------------->Estimation procedure starts here<-----------------------------------------#
  #----------------------------------------------------------------------------------------------------------------#
  
  #***Evaluating MLE and starting values for SMLE***#
  n <- nrow(S)
  SV <- Starting_values(y=y, S=S, linkvartheta=linkvartheta)
  thetaStart <- SV$thetaStart 
  fit_mle <- SV$MLEfit
  theta_mle <- as.numeric(coef(fit_mle))
  se_mle <- as.numeric(sqrt(diag(vcov(fit_mle))))
  
  trace_mle <-  sum(se_mle^2)	
  z_mle <- as.matrix(theta_mle/se_mle)
  pvalue_mle <- apply(z_mle,1,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
  etahat_mle <- S%*%theta_mle
  varthetahat_mle <- vartheta_link_functions(kappa = theta_mle, S=S, linkvartheta=linkvartheta)$vartheta
  
  #********RESULTS for alpha0=1********#
  results_mle_w <- list(kappa = theta_mle, se_kappa = se_mle,
                        z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                        alpha_value = 0, weights= rep(1,n)) #results with weights
  results_mle <- list(kappa = theta_mle, se_kappa = se_mle,
                      z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                      alpha_value = 0) #results without weights
  
  results_NC_w <- list(kappa = theta_mle, se_kappa = se_mle,
                       z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                       alpha_value = 0, weights= rep(1,n), Warning = "Non-convergence with alpha=alpha0; alpha=0 is used.")
  results_NC <- list(kappa = theta_mle, se_kappa = se_mle,
                     z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                     alpha_value = 0, Warning = "Non-convergence with alpha=alpha0; alpha=0 is used.")                 
  
  #------------------------------->If q is fixed<----------------------------------#
  if(alphaoptimal==FALSE){#starts option "alphaoptimal==FALSE"
    if(alpha0==0){#******alpha0 equal to 0 (MLE)*******#
      if(weights==T){
        results <- results_mle_w
      }
      else{
        results <- results_mle
      }#closes else of weights
      return(results) #returns results based on MLE
    }#closes if for alpha0 = 0 
    else{ #opens if for alpha0 > 0 
      alpha_const <- alpha0
      #-------->maximization with q fixed<------------#
      res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1, maxit = maxit),
                                             gr = Score ,method=method), error=function(e) {e}))#maximization 
      if(is.error(res)){#starts else in case of error in the estimation procedure
        if(weights==T){
          results <- results_NC_w
        }
        else{
          results <- results_NC
        }
        return(results)
      }   
      else{
        if(res$convergence == 0||res$convergence == 1){ 
          estimates <- res$par
          log.lik=res$value
          Vq <- V_matrix(estimates)	
          se_smle <-  suppressWarnings(t(sqrt(diag(Vq))))	 
          trace_smle <-  sum(diag(Vq))	
          z_smle <- as.matrix(estimates/se_smle)
          pvalue_smle <- apply(z_smle,2,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
          etahat_smle <- S%*%estimates
          varthetahat_smle <- vartheta_link_functions(kappa = estimates, S=S, linkvartheta=linkvartheta)$vartheta
          if(weights==T){
            w <- (dbinom(y, 1, varthetahat_smle))^(1.0 - alpha0)
            results <- list(kappa = estimates, se_kappa = se_smle,
                            z_kappa = z_smle, 
                            pvalue_kappa = pvalue_smle, alpha_value = alpha0,
                            weights = as.numeric(t(w)), Initial_Values = thetaStart)
          }
          else{
            results <- list(kappa = estimates, se_kappa = se_smle,
                            z_kappa = z_smle, 
                            pvalue_kappa = pvalue_smle, alpha_value = alpha0,
                            Initial_Values = thetaStart)
          }
          return(results)
        }
        else {# in case of no convergence, take the MLE
          if(weights==T){
            results <- results_NC_w
          }
          else{
            results <- results_NC
          }#closes else of weights 
          return(results)
        }
      }
    }#closes else of alpha0 < 1 
    
    #------------------------------->If q is not fixed<----------------------------------#
    
  }else {#starts option "alphaoptimal==TRUE"
    
    #***************Searching for the optimal q starts here****************#
    kk1 <- ncol(S)
    alpha_gridI <- sort(seq(from=0,to=0.5,by=spac), decreasing = FALSE) #first grid
    Size_gridI <- length(alpha_gridI)
    thetahat_I <- matrix(numeric(0),nrow=kk1, ncol= Size_gridI) 
    se_smleI <- matrix(numeric(0),nrow=kk1, ncol= Size_gridI)
    trace_smleI <- matrix(numeric(0), nrow=1, ncol= Size_gridI)
    SQVI <- numeric(Size_gridI-1)
    thetahat_I[,1] <-  theta_mle
    se_smleI[,1] <-  se_mle
    
    for(j in 1:(Size_gridI-1)){
      alpha_const <- alpha_gridI[j+1]
      resT <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1, maxit=maxit), 
                                              gr = Score ,method=method), error=function(e) {e}))#maximization
      if(is.error(resT)){ #"non-convergence" 
        thetahat_I[,j + 1.0] <- NA   
        se_smleI[,j + 1.0] <- NA
        trace_smleI[,j+1.0] <- NA
        SQVI[j] <- NA
      }
      else{# if estimation procedure converged
        estimatesI <- resT$par
        thetahat_I[,j + 1.0] <- estimatesI
        VqI <- V_matrix(thetahat_I[,j + 1.0])	
        se_smleI[,j + 1.0] =  suppressWarnings(t(sqrt(diag(VqI))))
        trace_smleI[,j+1.0] <- sum(diag(VqI)) 	  
        
        #**** calculate SQV ****#
        SQVI[j] <- round((1.0/(kk1))*sqrt(sum( (thetahat_I[,j]/(sqrt(n)*se_smleI[,j]) - 
                                                  thetahat_I[,j + 1.0]/(sqrt(n)*se_smleI[,j + 1.0]))^2)),5); 
      }                                                     
    }  
    
    alphavaluesout <- alpha_gridI[-1][SQVI>=L] 
    alphastart <- ifelse(suppressWarnings(max(alphavaluesout)) == 0.5, 0.5,
                         suppressWarnings(min(alpha_gridI[alpha_gridI>max(alphavaluesout)])))
    seq_pQS <- 1:Size_gridI
    position_alphastart <- seq_pQS[alpha_gridI==alphastart]
    
    if(suppressWarnings(min(alphavaluesout, na.rm=T))==Inf){# stability of estimates for q in [0.8, 1]
      alpha_optimal <- 0 
      if(weights==T){ #returns MLE with weights
        results <- list(kappa = theta_mle, se_kappa = se_mle,
                        z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                        alpha_value = 0, weights = rep(1,n), Initial_Values = thetaStart)
      }
      else{ #returns MLE without weights
        results <-list(kappa = theta_mle, se_kappa = se_mle,
                       z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                       alpha_value = 0, Initial_Values = thetaStart)
      }
      
    }else{ #if there is SQV>=L in [0.8, 1]
      
      #***grid of values for q starting in alphastart***#
      alpha_values <- sort(seq(from=alphastart,to=alphamax,by=spac), decreasing = FALSE) 
      nalpha <- length(alpha_values)
      SQV <- numeric(nalpha-1)
      thetahat_n <- matrix(numeric(0),nrow=kk1,ncol= nalpha)
      se_smle <- matrix(numeric(0),nrow=kk1,ncol= nalpha)
      trace_smle <- matrix(numeric(0),nrow=1,ncol= nalpha)		  
      thetahat_n[,1] =  thetahat_I[,position_alphastart]
      se_smle[,1] =  se_smleI[,position_alphastart]
      trace_smle[,1] =  sum(se_smle[,1])	  	 
      counter <- 1
      grid <- m;  alpha_const <- alphastart
      f1 <- 0.0; fc <- 0; f1e <- 0.0; cfailure <- 0.0
      alpha_valuesF <- sort(seq(from=0,to=alphamax,by=spac), decreasing = FALSE)
      SQVF <- rep(NA, length(alpha_valuesF) - 1)
      SQVF[1:length(SQVI)] <- SQVI
      signal1 <- 0; signal2 <- 0
      
      thetahat_all <- matrix(numeric(0),nrow=kk1,ncol= length(alpha_valuesF))
      se_smle_all <- matrix(numeric(0),nrow=kk1,ncol= length(alpha_valuesF))
      trace_smle_all <- matrix(numeric(0),nrow=1,ncol= length(alpha_valuesF))	
      
      thetahat_all[,1:Size_gridI] <-  thetahat_I 
      se_smle_all[,1:Size_gridI] <-  se_smleI  
      trace_smle_all[,1:Size_gridI] <- trace_smleI   
      
      #Results MLE
      resultsMLE_OW <- list(kappa = theta_mle, se_kappa = se_mle,
                            z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                            alpha_value = 0, weights = rep(1,n), Initial_Values = thetaStart, Warning= "Lack of stability")
      
      resultsMLE_O <- list(kappa = theta_mle, se_kappa = se_mle,
                           z_kappa = z_mle, pvalue_kappa = pvalue_mle,
                           alpha_value = 0, Initial_Values = thetaStart, Warning= "Lack of stability")
      
      #***************Searching for the optimal q from alphastart****************#    
      while(grid > 0 && alpha_const <= alphamax){
        alpha_const <- alpha_values[1.0 + counter]
        res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1, maxit=maxit), 
                                               gr = Score ,method=method), error=function(e) {e}))#maximization
        
        if(is.error(res) && alpha_const < alphamax){
          counter <- counter + 1.0
          grid = m # if no convergence, take one more q for the grid
          f1= f1 + 1.0  # sum of failures
          if(f1 == nalpha - 1.0) {#if the number of failures is equal to the grid size, take SMLE = MLE  
            grid = 0;	  #grid = 0 means that the search for the optimal value of q is over
            alpha_optimal =	0	 #alpha_optimal = 1 means that SMLE = MLE 
          }
          
          if(f1>=3){
            signal1 <- 1
            grid <- 0
            qM <- as.numeric(na.omit(alpha_valuesF[-length(alpha_values)][SQVF<L]))
            if(length(qM)<=3){
              alpha_optimal <- 0 
              if(weights==T){ #returns MLE with weights
                results <- resultsMLE_OW
              }
              else{ #returns MLE without weights
                results <- resultsMLE_O
              }
            }else{
              gride <- m
              qtest <- 1.0
              counter_test <- 1
              while(gride>0 && qtest>min(qM)){
                if(is.na(SQVF[counter_test])){
                  gride = m
                }else{
                  if(SQVF[counter_test] < L){
                    gride = gride - 1.0;
                    if(gride==0) alpha_optimal <- alpha_valuesF[counter_test-m+1.0] 
                  }  else{  
                    gride = m
                  } 
                } 
                counter_test <- counter_test + 1.0 
                qtest <- alpha_valuesF[counter_test]                       
              }
              if(gride>0 && qtest==min(qM)){
                alpha_optimal <- 0     
                if(weights==T){ #returns MLE with weights
                  results <- resultsMLE_OW
                }
                else{ #returns MLE without weights
                  results <- resultsMLE_O
                }
              }   
            }
            
          }#closes if f1>=3
        }
        else{# if estimation procedure converged
          estimates <- res$par
          log.lik=res$value
          thetahat_n[,counter + 1.0] <- estimates
          Vq <- V_matrix(thetahat_n[,counter + 1.0])	
          se_smle[,counter + 1.0] =  suppressWarnings(t(sqrt(diag(Vq))))	  
          trace_smle[,counter + 1.0] =  sum(diag(Vq))	 
          
          #****checking stability of the estimates****#
          SQV[counter] <- round((1.0/(kk1))*sqrt(sum( (thetahat_n[,counter]/(sqrt(n)*se_smle[,counter]) - 
                                                         thetahat_n[,counter+1.0]/(sqrt(n)*se_smle[,counter + 1.0]))^2)),5) 
          
          SQVF[length(alpha_valuesF) - length(alpha_values) + 1.0 + counter] <- SQV[counter]
          thetahat_all[,length(alpha_valuesF)-length(alpha_values)+ 1.0 + counter] <- estimates
          se_smle_all[,length(alpha_valuesF)-length(alpha_values)+ 1.0 + counter] <- suppressWarnings(t(sqrt(diag(Vq))))
          trace_smle_all[,length(alpha_valuesF)-length(alpha_values)+ 1.0 + counter] <- sum(diag(Vq))
          
          #***In case of NaN in SQV***#
          if(is.nan(SQV[counter])){
            grid = m # if SQV[counter]=NaN, take one more q for the grid
            f1= f1 + 1.0  # sum of failures
          }	else{  	
            if(SQV[counter] < L){# if stability condition is satisfied
              grid = grid - 1.0 # subtract 1 in the grid	
              if(grid==0) alpha_optimal = alpha_values[counter-m+1.0] #if grid = 0 take the maxivarthetam m (i.e, take the q-max of the grid) 
            }  else{  # if  stability condition is not satisfied  
              grid = m
              f1e = f1e + 1.0       
            }
          }                                                             
          counter = counter + 1.0
          
        }#closes else of "the estimation procedure converged"
        
        if((grid>0)&&alpha_const == alphamax){#if alphamax is reached and there is no stability of the estimates
          grid <- 0; signal2 <- 1
          qME  <- as.numeric(na.omit(alpha_valuesF[-length(alpha_values)][SQVF<L]))
          alpha_optimal <- 0
          if(length(qME)<=3){
            alpha_optimal <- 0;  
            if(weights==T){ #returns MLE with weights
              results <- resultsMLE_OW
            }
            else{ #returns MLE without weights
              results <- resultsMLE_O
            }
          }else{
            
            grideE <- m
            qtestE <- 0
            counter_testE <- 1
            while(grideE>0&&qtestE>min(qME)){
              if(is.na(SQVF[counter_testE])){
                grideE = m
              }else{
                if(SQVF[counter_testE] < L){
                  grideE = grideE - 1.0;
                  if(grideE==0) alpha_optimal <- alpha_valuesF[counter_testE-m+1.0] 
                }  else{  
                  grideE = m
                } 
              } 
              counter_testE <- counter_testE + 1.0 
              qtestE <- alpha_valuesF[counter_testE]                      
            }
            if(grideE>0&&qtestE==min(qME)){
              alpha_optimal <- 0     
              if(weights==T){ #returns MLE with weights
                results <- resultsMLE_OW
              }
              else{ #returns MLE without weights
                results <-resultsMLE_O
              }
            }  
          }
        } ## closes if alphamax is reached and there is no stability of the estimates
      }
      
      
      #***************Searching for the optimal q ends here****************# 
      
      #****Selecting estimates corresponding to the optimal q****#
      
      if(alpha_optimal==0){
        if(weights==TRUE){
          results <- resultsMLE_OW
        }
        else{
          results <- resultsMLE_O
        }
      }else{
        if(signal1==1||signal2==1){
          seq <- 1:length(alpha_valuesF)
          index_op <- seq[alpha_valuesF==alpha_optimal]
          thehat_n_optimal <- thetahat_all[,index_op]
          se_smle_otimo <- se_smle_all[,index_op]
          trace_smle_optimal <- trace_smle_all[,index_op]
          z_smle <- as.matrix(thetahat_all[,index_op]/se_smle_all[,index_op])
        }else{
          seq <- 1:nalpha
          index_op <- seq[alpha_values==alpha_optimal]
          thehat_n_optimal <- thetahat_n[,index_op]
          se_smle_otimo <- se_smle[,index_op]
          trace_smle_optimal <- trace_smle[,index_op]
          z_smle <- as.matrix(thetahat_n[,index_op]/se_smle[,index_op])
        } 
        
        pvalue_smle <- apply(z_smle,1,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
        etahat_optimal <- S%*%thehat_n_optimal
        varthetahat_optimal <- vartheta_link_functions(kappa = thehat_n_optimal[1:kk1], S=S, linkvartheta=linkvartheta)$vartheta
        if(weights==T){
          w <- (dbinom(y, 1, varthetahat_optimal))^(1.0 - alpha_optimal)
          results <- list(kappa = thehat_n_optimal, se_kappa = se_smle_otimo,
                          z_kappa = z_smle, pvalue_kappa = pvalue_smle,
                          alpha_value = alpha_optimal, weights=as.numeric(t(w)), Initial_Values = thetaStart)
        }
        else{
          results <- list(kappa = thehat_n_optimal, se_kappa = se_smle_otimo,
                          z_kappa = z_smle, pvalue_kappa = pvalue_smle,
                          alpha_value = alpha_optimal, Initial_Values = thetaStart)
        }
        
        
      }
    } #closes else "#if there is SQV>=L in [0.8, 1]"
    return(results)
  }# ends option "qoptimal==TRUE"
}# ends function

vcov_RobInfBeta <- function(Upsilon, type, alpha_const, S, X, Z, y, linkvartheta = "logit", linkmu="logit", linkphi="log"){
  mean_link_functions <- function(beta, y, X, linkmu="logit"){
    kk1 <- ncol(X);  n <- nrow(X);  eta <- as.vector(X%*%beta)	 
    if(linkmu == "logit"){
      mu <- exp(eta)/(1.0+exp(eta))
      T_1 <- diag(mu*(1.0-mu))
      yg1 <- log(y/(1-y))
    }
    if(linkmu == "probit"){
      mu <- pnorm(eta) 
      T_1 <- diag(dnorm(qnorm(mu)))
      yg1 <- qnorm(y)
    }
    if(linkmu == "cloglog"){
      mu <- 1.0 - exp(-exp(eta)) 
      T_1 <- diag((mu - 1)*log(1 - mu))
      yg1 <- log(-log(1-y))
    }
    if(linkmu == "log"){
      mu <- exp(eta) 
      T_1 <- diag(mu)
      yg1 <- log(y)
    }
    if(linkmu == "loglog"){
      mu <- exp(-exp(-eta)) 
      T_1 <- diag(-mu*log(mu))
      yg1 <- -log(-log(y))
    }
    if(linkmu == "cauchit"){
      mu <- (pi^(-1))*atan(eta) + 0.5   
      T_1 <- diag((1/pi)*((cospi(mu-0.5))^2))     
      yg1 <- tan(pi*(y-0.5))                          
    }
    results <- list(mu= mu, T_1 = T_1, yg1=yg1)
    return(results)
  }#ends mean link function
  precision_link_functions <- function(gama, Z, linkphi="log"){
    kk2 <- ncol(Z);  n <- nrow(Z);  delta <- as.vector(Z%*%gama) 
    if(linkphi == "log"){
      phi <- exp(delta) 
      T_2 <- diag(phi)
    }
    if(linkphi == "identify"){
      phi <- delta 
      T_2 <- diag(rep(1,n))
    }
    if(linkphi == "sqrt"){
      phi <- delta^2 
      T_2 <- diag(2*sqrt(phi))
    }
    results <- list(phi=phi, T_2=T_2)
    return(results)
  }#ends precision link function
  vartheta_link_functions <- function(kappa, y, S, linkvartheta="logit"){
    eta <- as.vector(S%*%kappa)	 
    linkstr <- linkvartheta
    if (linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      linkobj$dalpha.deta <- make.dalpha.deta(linkstr)
    } else {
      linkobj <- structure(list(linkfun = function(alpha) -log(-log(alpha)),
                                linkinv = function(eta){
                                  pmax(pmin(exp(-exp(-eta)),
                                            1 - .Machine$double.eps), .Machine$double.eps)
                                },
                                mu.eta = function(eta) {
                                  eta <- pmin(eta, 700)
                                  pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                                },
                                dalpha.deta = function(eta) pmax(exp(-exp(-eta) -
                                                                       eta) * expm1(-eta), .Machine$double.eps),
                                valideta = function(eta) TRUE, name = "loglog"),
                           class = "link-glm")
    }
    vartheta <- linkobj$linkinv(eta)
    T_1 <- diag(linkobj$mu.eta(eta))
    
    results <- list(vartheta= vartheta, T_1 = T_1)
    return(results)
  }#ends mean link function
  
  kk0 <- ncol(S)
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  kappa_p <- Upsilon[1:kk0]
  beta_p <- Upsilon[(kk0 + 1):(kk0 + kk1)]
  gama_p <- Upsilon[(kk0 + kk1 + 1.0):(kk0 + kk1 + kk2)]
  
  theta <- c(beta_p, gama_p)
  vartheta <- vartheta_link_functions(kappa = kappa_p, S=S, linkvartheta=linkvartheta)$vartheta
  
  if(type == "M_LSE"){
    muhat <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$mu
    phihat_alpha <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$phi 
    
    ahat_alpha <- muhat*phihat_alpha
    bhat_alpha <- (1.0 - muhat)*phihat_alpha
    
    T_1       <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$T_1
    T_2_alpha <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$T_2
    
    phihat_n <- phihat_alpha/(1 - alpha_const)
    muhat_n <- 	muhat
    
    ahat_n <- muhat_n*phihat_n
    bhat_n <- (1.0 - muhat_n)*phihat_n
    
    phihat_1_alpha <- phihat_n*(1 + alpha_const)	#expression of phi_(1 + alpha)
    
    a1_alphahat <- muhat_n*phihat_1_alpha
    b1_alphahat <- (1.0 - muhat_n)*phihat_1_alpha
    
    m_phi <- diag(phihat_n)
    psi1_n <- psigamma(ahat_n, 1.0) # psi(mu*phi)
    psi2_n <- psigamma(bhat_n, 1.0) # psi((1-mu)*phi)
    psi3_n <- psigamma(phihat_n, 1.0) #psi(phi)
    
    A <- diag(1 - vartheta)
    V <- diag(psi1_n + psi2_n)	# ok
    B1 <- diag(exp((1 - alpha_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - (lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha)))) #ok
    B2 <- diag(exp(lgamma(a1_alphahat) + lgamma(b1_alphahat) - lgamma(phihat_1_alpha) - (2.0*(alpha_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n))  +
                                                                                           lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha))))	# ok
    
    C <- diag(phihat_n*(muhat_n*psi1_n - (1.0 - muhat_n)*psi2_n))	# ok
    D <- diag((muhat_n^2.0)*psi1_n + ((1.0 - muhat_n)^2.0)*psi2_n - psi3_n) # ok
    
    psi1_1_alpha <- psigamma(a1_alphahat, 1.0)	
    psi2_1_alpha <- psigamma(b1_alphahat, 1.0)
    psi3_1_alpha <- psigamma(phihat_1_alpha, 1.0)	
    
    V_1_alpha <- diag(psi1_1_alpha + psi2_1_alpha) # ok
    C_1_alpha <- diag(phihat_n*(muhat_n*psi1_1_alpha - (1.0 - muhat_n)*psi2_1_alpha)) # ok
    D_1_alpha <- diag((muhat_n^2.0)*psi1_1_alpha + ((1.0 - muhat_n)^2.0)*psi2_1_alpha - psi3_1_alpha) # ok	
    
    Jalpha_betabeta <- as.matrix((1 - alpha_const)*t(X)%*%B1%*%A%*%(T_1^2.0)%*%(m_phi^2.0)%*%V%*%X)
    Jalpha_betagamma <- as.matrix(t(X)%*%B1%*%T_1%*%T_2_alpha%*%A%*%C%*%Z)
    Jalpha_gammagamma <- as.matrix(((1 - alpha_const)^(-1))*t(Z)%*%B1%*%A%*%(T_2_alpha^2.0)%*%D%*%Z)  
    
    Jalpha <- matrix(numeric(0), kk1+kk2, kk1+kk2) 
    Jalpha[1:kk1,1:kk1] <- Jalpha_betabeta
    Jalpha[1:kk1,(kk1+1):(kk1+kk2)] <- Jalpha_betagamma
    Jalpha[(kk1+1):(kk1+kk2),1:kk1] <- t(Jalpha_betagamma)
    Jalpha[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Jalpha_gammagamma)
    Jalpha <- -Jalpha
    
    Kalpha_betabeta <- as.matrix(t(X)%*%B2%*%(T_1^2.0)%*%A%*%(m_phi^2.0)%*%V_1_alpha%*%X) 
    Kalpha_betagamma <- as.matrix(((1 - alpha_const)^(-1))*t(X)%*%B2%*%T_1%*%A%*%T_2_alpha%*%C_1_alpha%*%Z)	 
    Kalpha_gammagamma <- as.matrix(((1 - alpha_const)^(-2))*t(Z)%*%B2%*%A%*%(T_2_alpha^2.0)%*%D_1_alpha%*%Z)
    Kalpha <- matrix(numeric(0),kk1+kk2,kk1+kk2)
    Kalpha[1:kk1,1:kk1] <- Kalpha_betabeta
    Kalpha[1:kk1,(kk1+1):(kk1+kk2)] <- Kalpha_betagamma
    Kalpha[(kk1+1):(kk1+kk2),1:kk1] <- t(Kalpha_betagamma)
    Kalpha[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Kalpha_gammagamma)
    Vq <- tryCatch( solve(Jalpha)%*%Kalpha%*%t(solve(Jalpha)), error=function(e) {e})  #asymptotic covariance matrix
    if(is.error(Vq)){
      Vq <- Ginv(Jalpha)%*%Kalpha%*%t(Ginv(Jalpha)) 
    }
    return(sqrt(diag(Vq)))
  }
  if(type == "M_ME"){
    q_const <- 1 - alpha_const
    kk1 <- ncol(X)
    kk2 <- ncol(Z)
    n <- nrow(X)
    beta <- theta[1:kk1]
    gama <- theta[(kk1+1.0):(kk1+kk2)]    	   
    mu <- mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$mu
    phi <- precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$phi
    a <- mu*phi						 
    b <- (1.0 - mu)*phi
    t_1 <- diag(mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$T_1)
    t_2 <- diag(precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$T_2)
    mustar <- psigamma(a, 0) - psigamma(b, 0)	  
    mudagger <- psigamma(b, 0) - psigamma(a + b, 0) 	
    m_phi <- diag(phi)	
    psi1 <- psigamma(a, 1.0)
    psi2 <- psigamma(b, 1.0) 
    psi3 <- psigamma(a + b, 1.0) 
    a_q <- (2.0 - q_const)*(a - 1.0) + 1.0						 
    b_q <- (2.0 - q_const)*(b - 1.0) + 1.0
    psi1_q <- psigamma(a_q, 1.0)
    psi2_q <- psigamma(b_q, 1.0) 
    psi3_q <- psigamma(a_q + b_q, 1.0) 
    a2_q <- (3.0 - 2.0*q_const)*(a - 1.0) + 1.0 
    b2_q <- (3.0 - 2.0*q_const)*(b - 1.0) + 1.0  
    psi1_2q <- psigamma(a2_q, 1.0)
    psi2_2q <- psigamma(b2_q, 1.0) 
    psi3_2q <- psigamma(a2_q + b2_q, 1.0)  
    mustar_q <- psigamma(a_q, 0) - psigamma(b_q, 0)
    mustar_2q <- psigamma(a2_q, 0) - psigamma(b2_q, 0)
    mudagger_q <-  psigamma(b_q, 0) - psigamma(a_q + b_q, 0)
    mudagger_2q <-  psigamma(b2_q, 0) - psigamma(a2_q + b2_q, 0)
    K <-  exp(lgamma(a_q) + lgamma(b_q) - lgamma(a_q + b_q) - (2.0 - q_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))							  
    K2 <-  exp(lgamma(a2_q) + lgamma(b2_q) - lgamma(a2_q + b2_q) - (3.0 - 2.0*q_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))   
    gama11_q <- diag(K*(phi^2.0)*(t_1^2.0)*(psi1_q + psi2_q + (mustar_q - mustar)^2.0))
    gama12_q <- diag(K*phi*t_1*t_2*(mu*(psi1_q +  psi2_q + (mustar_q - mustar)^2.0) - psi2_q + (mustar_q - mustar)*(mudagger_q - mudagger)))
    gama22_q <- diag(K*(t_2^2.0)*((mu^2.0)*psi1_q + ((1.0 - mu)^2.0)*psi2_q - psi3_q + (mu*(mustar_q - mustar) + mudagger_q - mudagger)^2.0))
    gama11_2q <- diag(K2*(phi^2.0)*(t_1^2.0)*(psi1_2q + psi2_2q + (mustar_2q - mustar)^2.0))
    gama12_2q <- diag(K2*phi*t_1*t_2*(mu*(psi1_2q + psi2_2q + (mustar_2q - mustar)^2.0) - psi2_2q + (mustar_2q - mustar)*(mudagger_2q - mudagger)))
    gama22_2q <- diag(K2*(t_2^2.0)*((mu^2.0)*psi1_2q + ((1.0 - mu)^2.0)*psi2_2q - psi3_2q + (mu*(mustar_2q - mustar) + mudagger_2q - mudagger)^2.0))
    E1q <- diag(K*phi*t_1*(mustar_q - mustar))
    E2q <- diag(K*t_2*(mu*(mustar_q - mustar) + mudagger_q - mudagger))
    A <- diag(1 - vartheta)
    
    Psin_betabeta <- -(2.0 - q_const)*as.matrix(t(X)%*%A%*%gama11_q%*%X)
    Psin_betagamma <- -(2.0 - q_const)*as.matrix(t(X)%*%A%*%gama12_q%*%Z)
    Psin_gammagamma <- -(2.0 - q_const)*as.matrix(t(Z)%*%A%*%gama22_q%*%Z) 
    Psin <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Psin[1:kk1,1:kk1] <- Psin_betabeta
    Psin[1:kk1,(kk1+1):(kk1+kk2)] <- Psin_betagamma 
    Psin[(kk1+1):(kk1+kk2),1:kk1] <- t(Psin_betagamma)
    Psin[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Psin_gammagamma 	
    
    Omegan_betabeta <- ((2.0 - q_const)^2.0)*as.matrix(t(X)%*%A%*%(gama11_2q - E1q^2.0)%*%X)
    Omegan_betagamma <- ((2.0 - q_const)^2.0)*as.matrix(t(X)%*%A%*%(gama12_2q - E1q*E2q)%*%Z)
    Omegan_gammagamma <- ((2.0 - q_const)^2.0)*as.matrix(t(Z)%*%A%*%(gama22_2q - E2q^2.0)%*%Z)  
    Omegan <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Omegan[1:kk1,1:kk1] <- Omegan_betabeta
    Omegan[1:kk1,(kk1+1):(kk1+kk2)] <- Omegan_betagamma 
    Omegan[(kk1+1):(kk1+kk2),1:kk1] <- t(Omegan_betagamma)
    Omegan[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Omegan_gammagamma		
    Vq <- tryCatch(solve(Psin)%*%(Omegan)%*%t(solve(Psin)), error=function(e) {e})  #asymptotic covariance matrix
    if(is.error(Vq)){
      Vq <- Ginv(Psin)%*%(Omegan)%*%t(Ginv(Psin))
    }
    return(sqrt(diag(Vq)))
  }
  if(type == "M_SE"){
    q_const <- 1 - alpha_const
    kk1 <- ncol(X)
    kk2 <- ncol(Z)
    beta_p <- theta[1:kk1]
    gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
    muhat_q <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$mu
    phihat_q <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$phi  
    ahat_q <- muhat_q*phihat_q
    bhat_q <- (1.0 - muhat_q)*phihat_q
    T_1_q <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$T_1
    T_2_q <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$T_2
    phihat_n <-   (1.0/q_const)*(phihat_q - 2.0) + 2.0
    muhat_n <- 	((1.0/q_const)*(muhat_q*phihat_q  - 1.0) + 1.0)/phihat_n
    ahat_n <- muhat_n*phihat_n
    bhat_n <- (1.0 - muhat_n)*phihat_n
    phihat_2_q <- (2.0 - q_const)*(phihat_n - 2.0) + 2.0;	#expression of phi_(2-q)
    muhat_2_q <- ((2.0 - q_const)*(ahat_n - 1.0) + 1.0)/phihat_2_q; #expression of mu_(2-q)
    a2_qhat <- muhat_2_q*phihat_2_q
    b2_qhat <- (1.0 - muhat_2_q)*phihat_2_q
    mustarhat_n <- psigamma(ahat_n, 0) - psigamma(bhat_n, 0)
    mudaggerhat_n <-  psigamma(bhat_n, 0) - psigamma(phihat_n, 0)	
    mustarhat_2_q <- psigamma(a2_qhat, 0) - psigamma(b2_qhat, 0)
    mudaggerhat_2_q <-  psigamma(b2_qhat, 0) - psigamma(phihat_2_q, 0)	
    muhat_d_2_q <- muhat_q*(mustarhat_2_q - mustarhat_n) + mudaggerhat_2_q - mudaggerhat_n
    m_phiq <- diag(phihat_q)
    psi1_n <- psigamma(ahat_n, 1.0) 
    psi2_n <- psigamma(bhat_n, 1.0) 
    psi3_n <- psigamma(phihat_n, 1.0) 
    V_n <- diag(psi1_n + psi2_n)	
    B1 <- diag(exp(q_const*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - (lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
    B2 <- diag(exp(lgamma(a2_qhat) + lgamma(b2_qhat) - lgamma(phihat_2_q) - (2.0*(1.0-q_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n))  +
                                                                               lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))	
    C_q_0 <- diag(phihat_q*(muhat_q*psi1_n - (1.0 - muhat_q)*psi2_n))	
    D_q_0 <- diag((muhat_q^2.0)*psi1_n + ((1.0 - muhat_q)^2.0)*psi2_n - psi3_n) 
    psi1_2_q <- psigamma(a2_qhat, 1.0)	
    psi2_2_q <- psigamma(b2_qhat, 1.0)
    psi3_2_q <- psigamma(phihat_2_q, 1.0)		
    V_2_q <- diag(psi1_2_q + psi2_2_q)
    C_q_2_q <- diag(phihat_q*(muhat_q*psi1_2_q - (1.0 - muhat_q)*psi2_2_q))
    D_q_2_q <- diag((muhat_q^2.0)*psi1_2_q + ((1.0 - muhat_q)^2.0)*psi2_2_q - psi3_2_q)	
    M1 <- diag(mustarhat_2_q - mustarhat_n)
    M2 <- diag(muhat_d_2_q)
    M3 <- diag((mustarhat_2_q - mustarhat_n)*muhat_d_2_q)	
    A <- diag(1 - vartheta)
    Jq_betabeta <- as.matrix(t(X)%*%B1%*%A%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%V_n%*%X)
    Jq_betagamma <- as.matrix(t(X)%*%B1%*%A%*%T_1_q%*%T_2_q%*%C_q_0%*%Z)
    Jq_gammagamma <- as.matrix(t(Z)%*%B1%*%A%*%(T_2_q^2.0)%*%D_q_0%*%Z)  
    Jq <- matrix(numeric(0), kk1+kk2, kk1+kk2) 
    Jq[1:kk1,1:kk1] <- Jq_betabeta
    Jq[1:kk1,(kk1+1):(kk1+kk2)] <- Jq_betagamma
    Jq[(kk1+1):(kk1+kk2),1:kk1] <- t(Jq_betagamma)
    Jq[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Jq_gammagamma)
    Jq <- -(q_const^(-1.0))*Jq
    Kq_betabeta <- as.matrix(t(X)%*%B2%*%A%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%(V_2_q+M1^2.0)%*%X) 
    Kq_betagamma <- as.matrix(t(X)%*%B2%*%A%*%T_1_q%*%T_2_q%*%m_phiq%*%((C_q_2_q*(1/phihat_q)) + M3)%*%Z)	 
    Kq_gammagamma <- as.matrix(t(Z)%*%B2%*%A%*%(T_2_q^2.0)%*%(D_q_2_q+M2^2.0)%*%Z)
    Kq <- matrix(numeric(0),kk1+kk2,kk1+kk2)
    Kq[1:kk1,1:kk1] <- Kq_betabeta
    Kq[1:kk1,(kk1+1):(kk1+kk2)] <- Kq_betagamma
    Kq[(kk1+1):(kk1+kk2),1:kk1] <- t(Kq_betagamma)
    Kq[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Kq_gammagamma)
    Kq <- (q_const^(-2.0))*Kq	
    Vq <- tryCatch( solve(Jq)%*%Kq%*%t(solve(Jq)), error=function(e) {e})  #asymptotic covariance matrix
    if(is.error(Vq)){
      Vq <- Ginv(Jq)%*%Kq%*%t(Ginv(Jq)) 
    }
    return(sqrt(diag(Vq)))
  }
  if(type == "M_LME"){
    kk1 <- ncol(X)
    kk2 <- ncol(Z)
    beta <- theta[1:kk1]
    gama <- theta[(kk1+1.0):(kk1+kk2)]
    
    mu  <- mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$mu
    phi <- precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$phi 
    
    a <- mu*phi
    b <- (1 - mu)*phi
    
    t_1 <- diag(mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$T_1)
    t_2 <- diag(precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$T_2)
    
    mustar <- psigamma(a, 0) - psigamma(b, 0)	  
    mudagger <- psigamma(b, 0) - psigamma(a + b, 0) 	
    
    m_phi <- diag(phi)	
    
    psi1 <- psigamma(a, 1.0)
    psi2 <- psigamma(b, 1.0) 
    psi3 <- psigamma(a + b, 1.0) 
    
    a_alpha <- mu*phi*(1 + alpha_const)			 
    b_alpha <- (1 - mu)*phi*(1 + alpha_const)	
    
    psi1_alpha <- psigamma(a_alpha, 1.0)
    psi2_alpha <- psigamma(b_alpha, 1.0) 
    psi3_alpha <- psigamma(a_alpha + b_alpha, 1.0) 
    
    a2_alpha <- mu*phi*(1 + 2*alpha_const)		
    b2_alpha <- (1 - mu)*phi*(1 + 2*alpha_const)	 
    psi1_2alpha <- psigamma(a2_alpha, 1.0)
    psi2_2alpha <- psigamma(b2_alpha, 1.0) 
    psi3_2alpha <- psigamma(a2_alpha + b2_alpha, 1.0)  
    
    mustar_alpha <- psigamma(a_alpha, 0) - psigamma(b_alpha, 0)
    mustar_2alpha <- psigamma(a2_alpha, 0) - psigamma(b2_alpha, 0)
    
    mudagger_alpha <-  psigamma(b_alpha, 0) - psigamma(a_alpha + b_alpha, 0)
    mudagger_2alpha <-  psigamma(b2_alpha, 0) - psigamma(a2_alpha + b2_alpha, 0)
    
    K <-  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) - (1 + alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))							  
    K2 <-  exp(lgamma(a2_alpha) + lgamma(b2_alpha) - lgamma(a2_alpha + b2_alpha) - (1 + 2*alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))   
    
    gama11_alpha <- diag(K*(phi^2.0)*(t_1^2.0)*(psi1_alpha + psi2_alpha + (mustar_alpha - mustar)^2.0))
    gama12_alpha <- diag(K*phi*t_1*t_2*(mu*(psi1_alpha +  psi2_alpha + (mustar_alpha - mustar)^2.0) - psi2_alpha + (mustar_alpha - mustar)*(mudagger_alpha - mudagger)))
    gama22_alpha <- diag(K*(t_2^2.0)*((mu^2.0)*psi1_alpha + ((1.0 - mu)^2.0)*psi2_alpha - psi3_alpha + (mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)^2.0))
    
    gama11_2alpha <- diag(K2*(phi^2.0)*(t_1^2.0)*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0))
    gama12_2alpha <- diag(K2*phi*t_1*t_2*(mu*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0) - psi2_2alpha + (mustar_2alpha - mustar)*(mudagger_2alpha - mudagger)))
    gama22_2alpha <- diag(K2*(t_2^2.0)*((mu^2.0)*psi1_2alpha + ((1.0 - mu)^2.0)*psi2_2alpha - psi3_2alpha + (mu*(mustar_2alpha - mustar) + mudagger_2alpha - mudagger)^2.0))
    
    gama1_alpha <- diag(K*phi*t_1*(mustar_alpha - mustar))
    gama2_alpha <- diag(K*t_2*(mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger))
    A <- diag(1 - vartheta)
    
    Psin_betabeta <- (1 + alpha_const)*as.matrix(t(X)%*%A%*%gama11_alpha%*%X)
    Psin_betagamma <- (1 + alpha_const)*as.matrix(t(X)%*%A%*%gama12_alpha%*%Z)
    Psin_gammagamma <- (1 + alpha_const)*as.matrix(t(Z)%*%A%*%gama22_alpha%*%Z) 
    Psin <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Psin[1:kk1,1:kk1] <- Psin_betabeta
    Psin[1:kk1,(kk1+1):(kk1+kk2)] <- Psin_betagamma 
    Psin[(kk1+1):(kk1+kk2),1:kk1] <- t(Psin_betagamma)
    Psin[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Psin_gammagamma 	
    
    Omegan_betabeta <- ((1 + alpha_const)^2.0)*as.matrix(t(X)%*%A%*%(gama11_2alpha - gama1_alpha^2.0)%*%X)
    Omegan_betagamma <- ((1 + alpha_const)^2.0)*as.matrix(t(X)%*%A%*%(gama12_2alpha - gama1_alpha*gama2_alpha)%*%Z)
    Omegan_gammagamma <- ((1 + alpha_const)^2.0)*as.matrix(t(Z)%*%A%*%(gama22_2alpha - gama2_alpha^2.0)%*%Z)  
    Omegan <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Omegan[1:kk1,1:kk1] <- Omegan_betabeta
    Omegan[1:kk1,(kk1+1):(kk1+kk2)] <- Omegan_betagamma 
    Omegan[(kk1+1):(kk1+kk2),1:kk1] <- t(Omegan_betagamma)
    Omegan[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Omegan_gammagamma		
    Vq <- tryCatch(solve(Psin)%*%(Omegan)%*%t(solve(Psin)), error=function(e) {e})  #asymptotic covariance matrix
    if(is.error(Vq)){
      Vq <- Ginv(Psin)%*%(Omegan)%*%t(Ginv(Psin))
    }
    return(sqrt(diag(Vq)))
  }
  
}

# This function provides residuals for the fits of inflated beta regression models.
# The function provides randomized quantile residual (type = "overall") or residuals
# for the discrete and continuous part separately (type = "bypart"); 
# Espinheira et al. (2008) 'standardized weighted residual 2' for the continuous
# and deviance for the discrete part.
# The link function for mu and vartheta are logit, and for phi, log link
residuals_ZIBE <- function(y, S, X, Z, theta, type = "overall") { 
  
  S <- as.matrix(add_intercept(S))
  X <- as.matrix(add_intercept(X))
  Z <- as.matrix(add_intercept(Z))
  
  p0 <- ncol(S)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  n  <- nrow(X)
  
  kapa_p <- theta[1:p0]
  beta_p <- theta[(p0 + 1):(p0 + p1)]
  gama_p <- theta[(p0 + p1 + 1):(p0 + p1 + p2)] 
  eta0_hat <- as.vector(S%*%kapa_p)	
  eta1_hat <- as.vector(X%*%beta_p)	
  eta3_hat <- as.vector(Z%*%gama_p) 		
  
  # vartheta - logit link function
  varthetahat  <- exp(eta0_hat)/(1 + exp(eta0_hat))
  T_0hat <- diag((varthetahat*(1 - varthetahat)))
  
  # mu - logit link function
  muhat  <- exp(eta1_hat)/(1 + exp(eta1_hat))
  T_1hat <- diag((muhat*(1 - muhat)))
  
  # phi - log link function
  phihat <- exp(eta3_hat) 
  T_2hat <- diag(phihat)
  
  # quantile residual
  if(type == "overall"){ 
    cdf <- varthetahat + (1 - varthetahat)*pbeta(y, muhat*phihat, (1-muhat)*phihat)
    cdf <- ifelse(cdf == 1, 0.999, cdf)
    cdf <- ifelse(cdf == 0, 0.001, cdf)
    res <- ifelse(y == 0, qnorm(runif(length(y),0, varthetahat)), qnorm(cdf))
  }
  
  # Espinheira et al. (2008) 'standardized weighted residual 2' for the continuos
  # and deviance for the discrete component
  if(type == "bypart"){
    
    # Continuous part
    PhiM <- diag(phihat)
    mustarhat <- psigamma(muhat*phihat, 0) - psigamma((1.0-muhat)*phihat, 0)	
    psi1hat <- psigamma(muhat*phihat, 1.0) 
    psi2hat <- psigamma((1.0-muhat)*phihat, 1.0)  
    What <- diag(phihat*(psi1hat+psi2hat))%*%(T_1hat^2.0)
    tempinvhat <- solve(t(X)%*%PhiM%*%What%*%X)
    Hhat <- diag(((What%*%PhiM)^(0.5))%*%X%*%tempinvhat%*%t(X)*((PhiM%*%What)^(0.5)))
    vhat <- psi1hat + psi2hat		 
    r_p1 <- (log(y/(1-y)) - mustarhat)/sqrt(vhat) 	
    r_p2_aux <- 1.0/sqrt(1.0- Hhat)
    r_p2 <- r_p1*r_p2_aux
    r_continuous <- r_p2[y>0]
    
    # Discrete part
    index <- ifelse(y == 0, 1, 0)
    qi    <- varthetahat*(1 - varthetahat)
    Q     <- diag(qi)
    Hhat2 <- diag(((Q)^(0.5))%*%Z%*% solve(t(Z)%*%Q%*%Z)%*%t(Z)*((Q)^(0.5)))
    r_discrete <- (index - varthetahat)/sqrt(varthetahat*(1 - varthetahat)*(1 - Hhat2))
    
    res <- list(res.cont = r_continuous, res.disc = r_discrete)
  }
  return(res)
}#ends function

# This function provides envelope plots based on the MLE and M-LSE (specify estimator = "MLE" and 
#  estimator = "M-LSE", respectively). The randomized quantile residual is used.
# The link function for mu and vartheta are logit, and for phi, log link
envelope <- function(y, S, X, Z, theta, estimator,
                     rep = 40, conf = 0.95){
  
  S <- as.matrix(add_intercept(S))
  X <- as.matrix(add_intercept(X))
  Z <- as.matrix(add_intercept(Z))
  
  p0 <- ncol(S)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  n  <- nrow(X)
  kapa_p <- theta[1:p0]
  beta_p <- theta[(p0 + 1):(p0 + p1)]
  gama_p <- theta[(p0 + p1 + 1):(p0 + p1 + p2)] 
  eta0_hat <- as.vector(S%*%kapa_p)	
  eta1_hat <- as.vector(X%*%beta_p)	
  eta3_hat <- as.vector(Z%*%gama_p) 		
  
  # vartheta - logit link function
  varthetahat  <- exp(eta0_hat)/(1 + exp(eta0_hat))
  
  # mu - logit link function
  muhat  <- exp(eta1_hat)/(1 + exp(eta1_hat))
  
  # phi - log link function
  phihat <- exp(eta3_hat) 
  
  a = muhat*phihat
  b = (1 - muhat)*phihat
  
  resRid <- residuals_ZIBE(y, S, X, Z, theta, type = "overall") 
  
  resid_env <- matrix(0, n, rep)
  i <- 1
  bar <- txtProgressBar(min = 0, max = rep, initial = 0, width = 50, char = "+", style = 3)
  while(i <= rep){
    tryCatch({
      
      prob <- runif(n)
      y_env <-  suppressWarnings(ifelse((varthetahat >= prob), 
                                        0, qbeta((prob - varthetahat)/(1-varthetahat), a, b, lower.tail = TRUE, log.p = FALSE)))
      y_env <- as.vector(y_env)
      if(estimator == "MLE"){
        
        index0_env <- y_env != 0 & y_env != 1 # y \in (0,1)
        
        y.1_env <- y_env[index0_env]
        X.1_env <- as.matrix(add_intercept(X[index0_env, ]))
        Z.1_env <- as.matrix(add_intercept(Z[index0_env, ]))
        
        y.2_env <- ifelse(index0_env, 0, 1) # I(y=0) or I(y=1)
        
        #discrete part
        fit.ber.env <- glm.fit(S, y.2_env, family = binomial(link = "logit"))$coef
        
        #continuous part
        fit.beta.env <- betareg::betareg.fit(X.1_env, y.1_env, Z.1_env)$coef
        
        val <- c(fit.ber.env, fit.beta.env$mean, fit.beta.env$precision)
      }else{
        if(estimator == "MLSE"){
          val <- fit.betainflated(y_env, S, X, Z)$coefficients$MLSE
          val <- c(val$probability[,1], val$mean[,1], val$precision[,1])
        }
      }
      
      res_env <- residuals_ZIBE(y_env, S, X, Z, val, type = "overall") 
      
      resid_env[,i] <- sort(res_env)
      setTxtProgressBar(bar,i)
      i = i + 1
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  liml <- apply(resid_env, 1, quantile, prob=(1-conf)/2)
  limu <- apply(resid_env, 1, quantile, prob=(1-(1-conf)/2))
  mm   <- apply(resid_env, 1, median)
  
  faixaRid <- range(resRid , liml , limu)
  xlab <- "Quantile N(0,1)"
  ylab <- "Residuals"
  main <- ""
  envcol <- "black"
  ylim <- faixaRid
  
  qqnormInt <- function(y,IDENTIFY = TRUE){
    qqnorm(y, pch = "+", las = 1, ylim = ylim, xlab = xlab, ylab = ylab, main = main, cex = 0.8,
           lwd = 3, las = 1) -> X
    cat(paste("\nClick on points you would like to identify and press Esc."))
    if(IDENTIFY) return(identify(X, cex = 0.8))
    invisible(X)
  }
  qqnormInt(resRid)
  
  par(new = TRUE)
  qqnorm(liml, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
         col = envcol)
  par(new = TRUE)
  qqnorm(limu, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
         col = envcol)
  par(new = TRUE)
  qqnorm(mm, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 2, main = main)
}


