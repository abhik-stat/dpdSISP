## =============================================================================
## Standard maximum likelihood (ML), via lme4.
## -----------------------------------------------------------------------------
LMM_MLE <- function(Y, X, Z, ID, gam = 0.25, maxitr = 500) {
  fit <- suppressWarnings(lme4::lmer(Y ~ X[, -1] + (Z[, -1] | ID), REML = FALSE))

  RE  <- as.matrix(lme4::ranef(fit)$ID)
  Beta <- lme4::fixef(fit)
  Sig  <- stats::sigma(fit)
  R    <- as.matrix(lme4::VarCorr(fit)$ID)

  ## singular fits are flagged here so a caller can distinguish "no signal" 
  ## from "the fit didn't really work."
  singular <- lme4::isSingular(fit)

  list(Beta = Beta, Sig = Sig, R = R, RE = RE, Mu = NULL, weight = NULL,
      itr = 0L, singular = singular)
}


## -----------------------------------------------------------------------------
## Standard restricted maximum likelihood (REML), via lme4. 
## -----------------------------------------------------------------------------
LMM_REML <- function(Y, X, Z, ID, gam = 0.25, maxitr = 500) {
  fit <- suppressWarnings(lme4::lmer(Y ~ X[, -1] + (Z[, -1] | ID), REML = TRUE))

  RE  <- as.matrix(lme4::ranef(fit)$ID)
  Beta <- lme4::fixef(fit)
  Sig  <- stats::sigma(fit)
  R    <- as.matrix(lme4::VarCorr(fit)$ID)
  singular <- lme4::isSingular(fit)

  list(Beta = Beta, Sig = Sig, R = R, RE = RE, Mu = NULL, weight = NULL,
      itr = 0L, singular = singular)
}


## -----------------------------------------------------------------------------
## Minimum hierarchical gamma-divergence estimator (HGD), fixed tuning
## parameter, via an MM (minorize-maximize) algorithm. Needs lme4 (for the
## initial fit), mvtnorm (for dmvnorm), Matrix (for the sparse ZZ).
##
## -----------------------------------------------------------------------------
LMM_HGD <- function(Y, X, Z, ID, gam=0.5, bw=NULL, maxitr=100, init.fit=NULL){
  ## preparation 
  qq <- ncol(Z)    # dimension of random effects
  m <- length(unique(ID))

  ID_int <- match(ID, unique(ID))
  ni <- table(ID_int)
  N  <- sum(ni)
  tr <- function(mat) sum(diag(mat))
  
  ## initial fit (standard maximum likelihood)
  if(is.null(init.fit)){
    fit0 <- try(suppressWarnings(lme4::lmer(Y ~ X[, -1] + (Z[, -1] | ID_int), 
                                            REML = FALSE)), silent = TRUE)
    
    if (!inherits(fit0, "try-error")) {
      Beta <- lme4::fixef(fit0)
      V <- as.matrix(coef(fit0)$ID_int[,1:qq])
      V.mat <- as.matrix(cbind(V[,qq], V[,-qq]) )
      V <- as.vector(V.mat)
      
      Sig <- summary(fit0)$sigma
      R <- as.matrix(VarCorr(fit0)$ID_int)[1:qq, 1:qq] 
      R <- R + 0.01 * diag(qq)     # regularization to avoid numerical error
      } else {
        fit1 <- try(suppressWarnings(lme4::lmer(Y ~ X[, -1] + (1 | ID_int), 
                                                REML = FALSE)), silent = TRUE)
        
        if (!inherits(fit1, "try-error")) {
          Beta <- lme4::fixef(fit1)
          Sig  <- stats::sigma(fit1)
          icept_var <- as.numeric(lme4::VarCorr(fit1)$ID_int[1, 1])
          icept_re  <- as.numeric(lme4::ranef(fit1)$ID_int[[1]])
          V.mat <- matrix(0, m, qq); V.mat[, 1] <- icept_re
          V <- as.vector(V.mat)
          R <- diag(qq) * max(icept_var, 1e-3)   
          R[1, 1] <- max(icept_var, 1e-3)
          } else {
            ## final fallback: no lme4 call at all, pure data-driven guess
            Beta <- as.vector(stats::coef(stats::lm.fit(X, Y)))
            Sig  <- sqrt(max(stats::var(Y) / 2, 1e-3))
            V <- rep(0, m * qq)
            R <- diag(qq) * max(stats::var(Y) / 4, 1e-3)
          }
        }
    init.fit <- list(Beta = Beta, RE = V, V = R, Sig = Sig) 
    } else {
      Beta <- init.fit$Beta
      V <- init.fit$RE
      R <- init.fit$V
      Sig <- init.fit$Sig
    }
  
  
  ## design matrix for random effects (block per random-effects dimension,
  ## each block m columns wide, one column per cluster)
  ZZ <- matrix(0, N, m * qq)
  for (i in seq_len(m)) {
    for (k in seq_len(qq)) {
      ZZ[ID_int == i, m * (k - 1) + i] <- Z[ID_int == i, k]
    }
  }
  ZZ <- methods::as(ZZ, "sparseMatrix")
  
  ## weight (used for cluster bootstrap; 1 for every cluster otherwise)
  if (is.null(bw)) bw <- rep(1, m)
  bbw <- bw[ID_int]
  
  
  ##   MM algorithm   ##
  converged <- FALSE
  for (k in seq_len(maxitr)) {
    Beta0 <- Beta
    V0 <- V
    
    # weight
    mu <- as.vector(X%*%Beta+ZZ%*%V)
    val <- dnorm(Y, mu, Sig)^(gam)
    ww <- val/mean(val)*bbw
    
    # update beta
    resid <- as.vector(Y-ZZ%*%V)
    Beta <- as.vector( solve(t(ww*X)%*%X)%*%t(X)%*%(ww*resid) )
    
    # update sigma 
    resid <- as.vector(Y-X%*%Beta-ZZ%*%V)
    IS <- list()
    for(i in 1:m){
      sZ <- Z[ID_int==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      IS[[i]] <- solve( Sig^2*diag(ni[i]) + sZ%*%R%*%t(sZ) )
    }
    N2 <- Sig^2*sum(unlist(lapply(IS, tr))) - N*gam/(1+gam)
    Sig <- sqrt( sum(ww*resid^2)/N2 )
    
    # update random effects
    V.mat <- matrix(V, m, qq)
    log.uu <- gam*dmvnorm(V.mat, rep(0,qq), R, log=T)
    log.uu <- log.uu - max(log.uu)
    val <- exp(log.uu)
    uu <- val/mean(val)*bw
    resid <- as.vector(Y-X%*%Beta)
    V.mat <- matrix(NA, m, qq)
    for(i in 1:m){
      sZ <- Z[ID_int==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      V.mat[i,] <- as.vector(
        solve(t(ww[ID_int==i]*sZ)%*%sZ + uu[i]*Sig^2*solve(R) 
              + 10^(-5)*diag(qq))%*%t(sZ)%*%(ww[ID_int==i]*resid[ID_int==i]) )
    }
    V <- as.vector(V.mat)
    
    # update R
    mat <- 0
    for(i in 1:m){
      sZ <- Z[ID_int==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      IS <- solve( Sig^2*diag(ni[i]) + sZ%*%R%*%t(sZ) )
      mat <- mat + t(sZ)%*%IS%*%sZ
    }
    R <- (1+gam)*(t(uu*V.mat)%*%V.mat - R%*%mat%*%R + m*R)/m
    R <- R + 0.01*diag(qq)     # regularization to avoid numerical error
    
    # checking convergence
    dd <- sum(abs(Beta - Beta0)) / sum(abs(Beta0)+0.0001)
    if( dd<10^(-6) ){converged <- TRUE; break }
  }
  Mu <- as.vector(X%*%Beta+ZZ%*%V)
  
  ## Result
  Res <- list(Beta=Beta, Sig=Sig, R=R, RE=V, Mu=Mu, itr=k, converged = converged)
  return(Res)
}

