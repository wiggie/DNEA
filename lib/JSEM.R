#To compute the trace of a matrix
matTr <- function(z) sum(diag(z))


##--------------------------------------------\
## function: CGM_AHP_train.R
##--------------------------------------------\
## This is the code from Guo et al. (2011)
## trainX: data
## trainY: labels for categories (1, 2, 3,...)
## lambda_value: tuning parameter
## Make sure loading "glasso" package before using the code.
##--------------------------------------------\
CGM_AHP_train <- function(
  trainX, 
  trainY, 
  lambda_value,
  adaptive_weight = array(1, c(length(unique(trainY)), ncol(trainX), ncol(trainX))),
  eta = 0.01,
  limkappa = 1e+6  #limit for condition number of the sample cov 
){
  ## Set the general paramters
  K <- length(unique(trainY))
  p <- ncol(trainX)
  diff_value <- 1e+10
  count <- 0
  tol_value <- 0.01
  max_iter <- 30
  
  ## Set the optimizaiton parameters
  OMEGA <- array(0, c(K, p, p))
  S <- array(0, c(K, p, p))
  OMEGA_new <- array(0, c(K, p, p))
  nk <- rep(0, K)
  
  ## Initialize Omega
  for (k in seq(1, K)) {
    idx <- which(trainY == k)
    S[k, , ] <- cov(trainX[idx, ])
    if (kappa(S[k, , ]) > limkappa) {
      S[k, , ] <- S[k, , ] + eta * diag(p)
    }
    tmp <- solve(S[k, , ])
    OMEGA[k, , ] <- tmp
    nk[k] <- length(idx)
  }
  
  ## Start loop
  while ((count < max_iter) & (diff_value > tol_value)) {
    tmp <- apply(abs(OMEGA), c(2, 3), sum)
    tmp[abs(tmp) < 1e-10] <- 1e-10
    V <- 1/sqrt(tmp)
    
    for (k in seq(1, K)) {
      penalty_matrix <- lambda_value * adaptive_weight[k, , ] * V
      obj_glasso <- glasso(S[k, , ], penalty_matrix, maxit = 100)
      OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi))/2
      #OMEGA_new[k, , ] <- obj_glasso$wi
    }
    
    ## Check the convergence
    diff_value <- sum(abs(OMEGA_new - OMEGA))/sum(abs(OMEGA))
    count <- count + 1
    OMEGA <- OMEGA_new
    #cat(count, ', diff_value=', diff_value, '\n')
  }
  
  ## Filter the noise
  list.Omega <- NULL
  Theta <- vector("list", K)
  for (k in seq(1, K)) {
    ome <- OMEGA[k, , ]
    ww <- diag(ome)
    ww[abs(ww) < 1e-10] <- 1e-10
    ww <- diag(1/sqrt(ww))
    tmp <- ww %*% ome %*% ww
    ome[abs(tmp) < 1e-05] <- 0
    OMEGA[k, , ] <- ome
    list.Omega[[k]] <- ome
    Theta[[k]] <- diag(diag(list.Omega[[k]])^(-0.5)) %*% list.Omega[[k]] %*% diag(diag(list.Omega[[k]])^(-0.5))
  }
  
  output <- list()
  output$OMEGA <- list.Omega
  output$S <- S
  output$Theta <- Theta
  output$lambda <- lambda_value
  
  return(output)
}



#to select the optimal tuning parameter for Guo et al.
##--------------------------------------------\
CGM_AHP_tune <- function(
  trainX,   # training data
  testX,    # test data
  model,    # labels for models. 
  lambda,   # a vector of supplied lambda
  BIC=FALSE, # whether to compute the bic.score
  eps=1e-06,  
  eta=0.01,
  limkappa = 1e+6
){
  p = dim(trainX)[2]
  K = length(unique(model))
  N <- length(lambda)
  n <- rep(0, K)
  for (k in 1:K){ 
    n[k] <- length(which(model == k))
  }
  
  bic.score <- rep(0, N) 
  likelihood <- rep(0, N)
  for (j in 1:N){
    cat("The ", j, "-th step in tuning... \n")
    
    Omega.hat <- CGM_AHP_train(trainX = trainX, trainY = model, lambda_value = lambda[j], eta = eta)$OMEGA
    
    for (k in 1:K){
      data <- testX[which(model == k), ]      
      empcov <- cov(data) 
      if (kappa(empcov) > limkappa){
        empcov = empcov + eta * diag(p)      
      }   
      likelihood[j] = likelihood[j] + matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]]))
    }    
    
    if (BIC){
      for (k in 1:K){
        no.edge = sum(abs(Omega.hat[[k]])>eps) - p
        empcov <- cov(trainX[which(model==k),])
        bic.score[j] = bic.score[j] +  matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]])) + log(n[k]) * no.edge/(2*n[k])
      }
    }      
  }
  
  out <- list(BIC = bic.score, likelihood = likelihood)
  
  return(out)
}



##Stability selection for Guo's method conditional on lastar
##selected from cross validation
CGM_AHP_stabsel <- function(X, cnt, lastar, seed.base = 100, eta=0.01, limkappa=1e+6) {
  K = 1
  p = ncol(X)

  if (is.null(dim(X))) {
    K = length(X)
    p = ncol(X[[1]])
  }
  
  n = lapply(X, nrow)  
  
  X1 = vector("list", K)
  X2 = vector("list", K)
  sel.mat = vector("list", K)
  for (k in 1:K){
    sel.mat[[k]] = matrix(0, p, p)
  }
  count = 0
  
  for (i in 1:cnt) {
    set.seed(i+seed.base)
    cat("count is:", i, "\n")
    model.1 = NULL 
    model.2 = NULL 
    for (k in 1:K){
      ind.1 = sample(seq(1, n[[k]]), n[[k]]/2, F)
      ind.2 = seq(1, n[[k]])[match(seq(1, n[[k]]), ind.1, 0) == 0]
      X1[[k]] = X[[k]][ind.1, ]
      X2[[k]] = X[[k]][ind.2, ]
      model.1 = c(model.1, rep(k, length(ind.1)))
      model.2 = c(model.2, rep(k, length(ind.2)))
    }
    
    tmp.1 = try(CGM_AHP_train(trainX=do.call(rbind, X1), trainY=model.1, lambda_value=lastar, limkappa = limkappa, eta=eta))
    tmp.2 = try(CGM_AHP_train(trainX=do.call(rbind, X2), trainY=model.2, lambda_value=lastar, limkappa = limkappa, eta=eta))
    
    if (inherits(tmp.1, "try-error") || inherits(tmp.2, "try-error")){
      warning("There might be some error!")
      next;
    }
    
    for (k in 1:K){
      mat1 = tmp.1$OMEGA[[k]]
      mat1[which(abs(mat1)>1e-5)] = 1
      diag(mat1) = 0
      mat2 = tmp.2$OMEGA[[k]]
      mat2[which(abs(mat2)>1e-5)] = 1
      diag(mat2) = 0
      sel.mat[[k]] = sel.mat[[k]] + mat1 + mat2
    }
    
    count = count + 1
  }
  
  return(list(mat = sel.mat, count = count))
}

##--------------------------------------------\
# Tuning parameter selection:
# Given a grid of lambda, compute the bic.score and
# likelihood based on a test dataset
##--------------------------------------------\
glasso_tune <- function(
  trainX, # training data n by p
  testX,  # test data n by p
  rho,    # a vector of tuning parameters
  weights, #the weight matrix for penalization
  eps=1e-06,
  eta=0.01, 
  limkappa=1e+6
){
  p = ncol(trainX)
  n = nrow(trainX)
  bic.score <- rep(0, length(rho))
  for (j in 1:length(rho)){
    cat("The ", j, "-th step in glasso tuning... \n")
    
    fit <- glasso(s = cov(trainX), rho = rho[j]*weights, penalize.diagonal = FALSE)
    sigInv = (fit$wi + t(fit$wi))/2
    sigInv[abs(sigInv)<eps] = 0
    
    empcov <- cov(testX) #empirical cov
    if (kappa(empcov) > limkappa){
      empcov = empcov + eta * diag(p)
    }
    bic.score[j] = matTr(empcov %*% fit$wi) - log(det(fit$wi)) +  + log(n) * (sum((abs(sigInv) > eps))-p)/(2*n)
  }
  
  out <- list(BIC = bic.score)
  return(out)
}

require(glasso)
adjDGlasso <- function(
  X, #the n by p data matrix
  weights=1, #the weight for the penalty
  theta_star=NULL, #true precision matrix
  lambda = NULL, 
  FDR.type='BH', #FDR control procedure
  alpha = 0.1, #the significance level for FDR control
  quiet=TRUE
){
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center = T, scale = F)
  empcov <- (1/n) * (t(X) %*% X) #empirical cov
  if (is.null(lambda)){
    lambda <- sqrt(log(p)/n)
  }
  
  if (!quiet){print('fit glasso')}
  Theta.hat.from.Glasso <- glasso(s=empcov, rho=lambda*weights, penalize.diagonal=FALSE)$wi
  Theta.hat.from.Glasso <- (Theta.hat.from.Glasso + t(Theta.hat.from.Glasso))/2
  
  if (!quiet){print('done')}
  
  if (!quiet){print('de-biasing glasso')}
  ## T.hat = Theta.hat - Theta.hat * (Sigma.hat - Theta.hat^{-1}) * Theta.hat
  ## Theta.hat and Sigma.hat are both symmetric.
  temp.mat <- empcov - chol2inv(chol(Theta.hat.from.Glasso))
  temp.vec <- as.vector(Theta.hat.from.Glasso %*% temp.mat %*% t(Theta.hat.from.Glasso))
  T.hat <- as.vector(Theta.hat.from.Glasso) - temp.vec
  T.hat <- matrix(T.hat,nrow = p)
  
  if (!quiet){print('done')}
  
  sigma.hat2 <- array(0,c(p,p))
  for (i in 1:p){
    for (j in 1:p){
      sigma.hat2[i,j] <- Theta.hat.from.Glasso[i,j]^2+Theta.hat.from.Glasso[i,i]*Theta.hat.from.Glasso[j,j]
    }
  }
  
  test.stat <- sqrt(n)*T.hat/sqrt(sigma.hat2)
  std.statistic <- NULL
  if (!is.null(theta_star)){
    std.statistic <- sqrt(n)*(T.hat - theta_star)/sqrt(sigma.hat2)
  }
  
  pvals <- 2*(pnorm(abs(test.stat), lower.tail=FALSE))
  pvals.vec <- lowerTriangle(pvals,diag=FALSE)
  adjpvals.vec <- p.adjust(pvals.vec, FDR.type)
  coeff <- diag(1,p) - cov2cor(T.hat)
  # coeff <- - pmax(pmin(T.hat, 1), -1)
  Qmat <- matrix(0, p, p)
  lowerTriangle(Qmat, diag=FALSE) <- adjpvals.vec
  Qmat <- Qmat + t(Qmat)
  diag(Qmat) <- rep(1, p)
  Qmat.fdr <- (Qmat <= alpha)
  
  # Qmat is the p by p matrix of qvalues;
  # Qmat.fdr is the thresholded matrix of qvalues based on alpha
  return(list(Theta=coeff, pvalue=pvals, qvalue = Qmat, qvalue.fdr=Qmat.fdr, statistic=std.statistic))
}


##-----------------------------------\
##    StructDiff.full.R
##-----------------------------------\
### Purpose: to compare the estimated and the truth matrices. 
### Inputs: 
##    Ahat: a list of the estimated matrices. 
##    Amat: a list of the corresponding true matrices.
##  Outputs:
##   FPrate, FNrate, SHD, Floss, KL loss
##   The above measures are as defined in the paper, which are the average deviance
##   from all pairs of matrices.
##--------------------------------------------\
## The following is specifically for calculating the ROC curve
StructDiff <- function(Ahat, Amat, eps = 1e-08){
  K <- length(Amat)
  
  if (is.null(dim(Amat)) == F){
    # indicating there's only one matrix compared.
    K <- 1
    Ahat <- list(Ahat)
    Amat <- list(Amat)
  }
  
  p = dim(Amat[[1]])[1]
  TP <- rep(0, K)
  FP <- rep(0, K)
  TN <- rep(0, K)
  FN <- rep(0, K)
  SHD <- rep(0, K)
  FPrate <- rep(0, K)
  TPrate <- rep(0, K)
  FNrate <- rep(0, K)
  Pr <- rep(0, K)
  Re <- rep(0, K)
  F1 <- rep(0, K)
  Floss <- rep(0, K)
  for (k in 1:K){
    Floss[k] <- norm(Ahat[[k]] - Amat[[k]], "F")/norm(Amat[[k]], "F") 
    diag(Ahat[[k]]) = 0
    diag(Amat[[k]]) = 0
    
    TP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) > eps))
    TN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) <= eps))
    FP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) <= eps))
    FN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) > eps))  
    SHD[k] <- FP[k] + FN[k]
    
    P <- TP[k] + FN[k]
    N <- TN[k] + FP[k]
    TPrate[k] <- TP[k]/(P + eps)
    FPrate[k] <- FP[k]/(N + eps)
    FNrate[k] <- FN[k]/(P + eps)
    
    Re[k] <- TP[k]/(P + eps) ## Recall
    Pr[k] <- TP[k]/(TP[k] + FP[k] + eps) ## Precision
    
    F1[k] <- (2 * Pr[k] * Re[k])/(Pr[k] + Re[k] + eps)
    
  }
  
  dev = data.frame(TP = TP, FP = FP, TN = TN, FN = FN, 
                   TPrate = TPrate, FPrate = FPrate, FNrate = FNrate, 
                   SHD = SHD, Precision = Pr, Recall = Re, F1 = F1, FL = Floss)
  return(dev)
} 

