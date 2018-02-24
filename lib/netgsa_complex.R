
## NetGSA for complex experiments including
## (1) Estimation for directed and undirected networks
## (2) Testing for 2 or more conditions

NetGSA <-
  function(
    A,  #Adjacency matrix in a list	    	
    x, 	#The p x n data matrix		      	
    y,    #vector of class indicators of length n
    B, 	#indicator matrix for pathways (npath x p)	    	  	
    lklMethod = c("REML","ML"),
    directed = FALSE,          
    eta = 1e-1,           
    lim4kappa = 500       
  ){
    this.call <- match.call()
    lklMethod <- match.arg(lklMethod)
    # s2profile <- match.arg(s2profile)
    
    p = dim(x)[1] #No. of genes
    n = length(y) #No. of samples in total
    
    if (dim(x)[2] != n) {
      stop("The dimensions of the data matrix and class vector don't match.")
    }
    
    if (dim(B)[2] != p) {
      stop("The dimensions of the data matrix and indicator matrix don't match.")
    }
    
    if (length(unique(y)) < 2) {
      stop("There should be at least 2 unique classes in the class indicator!")
    }
    
    if (min(sapply(lapply(A, abs), sum))==0) {
      stop("No network interactions were found!")
    }
    
    ##-----------------
    ##Determine whether the network is DAG
    ##Assume A1 and A2 are of the same type (directed or undirected)
    isNetDAG = FALSE
    if (directed) {    
      gA <- graph.adjacency(A[[1]], mode="directed")
      isNetDAG = is.dag(gA)
    }  
    
    if (p > 5000) {
      warning("netGSA may be slow for datasets with large number of genes.")
    }
    
    ##-----------------
    ##setting up control parameters for the var estimation procedures
    varEstCntrl = list(lklMethod = lklMethod,                    
                       s2profile = "se",   
                       lb = 0.5,           
                       ub = 100,           
                       tol = 0.01)         
    
    ##-----------------
    ##Find influence matrices based on adjacency matrices A1 and A2
    ##Check if the influence matrices are well conditioned. Otherwise update eta.
    if (directed){
      D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
      
      tmp = min(sapply(D, kappa)) 
      while ((tmp> lim4kappa) && !isNetDAG) {
        eta = eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))
        D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        tmp = min(sapply(D, kappa)) 
      }
      
      DD = lapply(D, function(m) m %*% t(m))
      
      tmp = min(sapply(DD, kappa))    
      while ((tmp > lim4kappa) && !isNetDAG) {
        eta = eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
        D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        DD = lapply(D, function(m) m %*% t(m))
        tmp = min(sapply(DD, kappa))
      }
      
      # #     if (((tmp > lim4kappa) && !isNetDAG)) {
      # warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
      # }
    } else {
      #Undirected gaussian graphical model
      #Need to normalize the matrix differently
      Ip = diag( rep(1,p) )
      print('Making call in netgsa')
      D = lapply(A, function(m) t(chol(solve(Ip - m))) ) 
    }
    
    output = call.netgsa(D, x, y, B, varEstCntrl)
    
    return(output)
  }

call.netgsa <-
  function(
    D,    	            #Influence matrices in a list
    x, 			      	#the p x n data matrix
    y, 			      	#vector of class indicators of length n
    B, 		    	  	#indicator matix for pathways (npath x p)
    varEstCntrl           #parameters to pass for variance estimation
  ){
    p = nrow(x) #No. of genes
    n = length(y) #No. of samples in total
    npath = nrow(B) #No. of gene sets to consider
    
    y = as.integer(as.factor(y))
    ncond = length(unique(y))
    n_vec = as.numeric(table(y))
    
    ##The identity matrix 
    Ip = matrix(0, p, p)
    diag(Ip) = 1
    
    xx = vector("list", ncond)  
    for (i in 1:ncond){
      xx[[i]] = x[, (y==i)]
    }
    
    DInv = lapply(D, solve)
    DtD = lapply(1:ncond, function(k) D[[k]] %*% t(D[[k]]))
    tDD = lapply(1:ncond, function(k) t(D[[k]]) %*% D[[k]]) 
    DDInv = lapply(tDD, solve)
    
    #Initialzing pvalues of each test
    pvals = matrix(0, npath, 1)
    
    ##------------------
    ##ESTIMATION OF BETA
    ##------------------
    beta = vector("list", ncond)
    for (i in 1:ncond){
      beta[[i]] = (1/n_vec[i]) * DInv[[i]] %*% apply(xx[[i]], 1, sum)
    }
    
    ##-----------------
    ##ESTIMATION OF VAR COMPONENTS
    ##-----------------
    #First calculate residuals: 
    resid = vector("list", ncond)
    for (i in 1:ncond){
      resid[[i]] = xx[[i]] - replicate(n_vec[i], c(D[[i]] %*% beta[[i]]))
    }
    temp = do.call(cbind, resid)
    
    #initial estimates
    sg0 = (mean(apply(temp, 2, sd)))^2
    se0 = (sd(diag(temp)))^2

    ##Estimation of variance components is based on profile likelihood with Newton's method.
    ##It can be done using either ML or REML.
    if (p<5000){
      # approxEst = approxVarEst(se0, sg0, D, resid, n_vec, control = varEstCntrl) 
      # se0 = approxEst$s2e
      # sg0 = approxEst$s2g
      # if (min(c(se0,sg0)<0)){
      #   sg0 = (mean(apply(temp, 2, sd)))^2
      #   se0 = (sd(diag(temp)))^2
      # }
      # print(se0); print(sg0);
      S = profile.newton.se(sg0/se0, D, resid, control = varEstCntrl)
      s2e = S$s2e
      s2g = S$s2g  	
    } else {
      s2e = se0
      s2g = sg0  	
    } 
    
    if (ncond==2){
      ##Use the T test
      res = netgsa.ttest(s2g, s2e, D, DtD, DDInv,n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g)      
    } else {
      ##Use the F test
      res = netgsa.ftest(s2g, s2e, D, DtD, DDInv, n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g) 
    } 
    
    return(output)
  }

adj2inf <-
  function(AA, isDAG = FALSE, eta=0.1){
    p = dim(AA)[1]
    Ip = diag( rep(1,p) )
    
    if(isDAG){
      
      DD = solve(Ip - AA)
      return(DD)
      
    }else{
      tmp = eigen(AA, symmetric=FALSE, only.values=TRUE)$values
      
      if( max(abs(tmp)) < 1-eta ){
        DD = solve(Ip - AA)
      }else{
        # This is consistent with the ref #2
        tmp = apply( abs(AA), 1, sum )
        tmp = matrix(tmp+eta, ncol=p, nrow=p, byrow=FALSE) 
        tmp = AA / tmp
        DD = solve(Ip - tmp)
      }  
      
      return(DD)
    }
  }

edgelist2adj <-
  function(file, vertex.names, mode=c("directed", "undirected")) {
    this.call <- match.call()
    mode <- match.arg(mode)
    is.directed <- ifelse(mode=="directed", TRUE, FALSE)
    
    dat = read.table(file, header=TRUE) 
    el = as.matrix(dat) 
    el[,1] = as.character(el[,1])
    el[,2] = as.character(el[,2])
    g = graph.edgelist(el, directed=is.directed)
    
    extra.vertices <- setdiff(vertex.names, V(g)$name)
    if (length(extra.vertices)>0){
      g <- add.vertices(g, nv = length(extra.vertices), name = extra.vertices)
    }
    Adj = as.matrix(get.adjacency(g, type="both"))
    
    reOrder <- match(vertex.names,rownames(Adj))
    Adj <- Adj[reOrder, reOrder]
    
    return(Adj)
  }

graphlaplacian <-
  function(A, zeta=0.01){
    if (sum(A != t(A)) > 0){
      stop("This method only works for symmetric A!")
    }
    Adeg = apply(abs(A), 1, sum)   # a p-dim vector
    Adeg[Adeg==0] = 1
    AdegInv = (Adeg + zeta)^(-0.5)
    Lunnorm = diag(Adeg) - A
    Lnorm = diag(AdegInv) %*% A %*% diag(AdegInv)
    
    return(list(Lunnorm = Lunnorm, Lnorm = Lnorm))
  }


glmnet.soft <- function(x, y, lambda){
  require(glmnet)
  
  if(ncol(x) > 1){		## use glmnet
    fit = glmnet(x, y, family="gaussian", alpha=1, lambda=lambda)
    beta = as.matrix(fit$beta)
  }else{				    ## use soft thresholding
    #beta = cor(y,x)	## assuming data is centered and scaled
    beta = lm(y~x)$coef[2]
    beta = sign(beta) * max((abs(beta) - lambda/2),0)
  }
  return(beta)
} 


## Purpose: to select the covariance matrix based on external information about 0's and 1's.  
## Input: X, the data matrix of n by p 
##        zero: the given external information about 0's 
##        one:  the given external information about 1's 
##        lambda: the tuning parameter in lasso regression, only accept one at a time
##        eps:  the minimum threshold value for determining edges. By default, it's 1e-8. 
##        verbose: whether to display intermediate outputs. 
## Output:  
##    Adj, the adjacency matrix selected 
##	  infmat, the influence matrix
## Notes:
##       1) In this version, it is assumed that the columns of X are ordered accoring to 
##			a correct (Wald) causal order, such that no X_j is a parent of X_k k \le j.
##		 2) To estimate the network, first each node is regressed on the known edges (one). 
##			The reisdual obtained from this regression is then used to find the 
##			additional edges, among the nodes that could potentially interact with the given
##			node (those not in zero).
##		 3) Given the causal ordering of nodes, the resulting adjacency matrix is lower triangular
##			(see Shojaie & Michailidis, 2010). Thus, only lower triangular parts of zero and one 
##			are used in this function. For this reason, it is important that both of these matrices
##			are also ordered according to the causal order of the nodes in X.
##------------------------------------------------- /
netEst.dir <- function(X, zero=NULL, one=NULL, lambda, verbose=FALSE, eps=1e-08) {
  n = dim(X)[1]
  p = dim(X)[2]
  Adj = matrix(0, p, p)
  Ip = diag(rep(1, p))
  
  if (is.null(zero)) {
    zero = matrix(0, p, p)
  }
  
  if (is.null(one)) {
    one = matrix(0, p, p)
  }
  
  if (sum(one*zero) > 0){
    stop("Information on 0's and 1's overlaps!")
  }  
  
  ## To get the empirical correlation matrix 
  X = scale(X, center=TRUE, scale=TRUE)
  for (i in 2:p) {
    Y = matrix(X[, i], ncol=1)
    Xmat = matrix(X[, 1:(i-1)], ncol=i-1)
    
    ## Get the zero and one indices 
    infoInd = one[i, 1:(i-1)] - zero[i, 1:(i-1)]
    beta = matrix(0, i-1, length(lambda))
    
    if (sum((infoInd == 0)) == 0) {
      if (verbose) {
        cat("Complete information known!  ")
      }
      beta[which(infoInd == 1), ] = 1
      beta[which(infoInd == -1), ] = 0
    } else {
      if (verbose) {
        cat("Incomplete information known!  ")
      }
      
      if (sum((infoInd == -1)) == 0) {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: no 0's and no 1's!  ")
          }
          beta = glmnet.soft(x=Xmat, y=Y, lambda=lambda)	#glmnet if ncol(X) > 1, o.w. soft thresholding
        } else {
          if (verbose) {
            cat("Incomplete information: no 0's, but with 1's!  ")
          }
          if (sum(infoInd==1)>=n) {#if there are fewer samples than the number of 1's, simply adopt the one's and pass
            beta[which(infoInd == 1), ] = 1 
            beta[which(infoInd == 0), ] = 0 
          } else {
            Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
            Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
            fit.lm = lm(Y ~ 0 + Xmat1)
            beta[(infoInd == 1), ] = as.numeric(coef(fit.lm))
            res.lm = as.numeric(residuals(fit.lm)) ##get the residuals           
            beta[(infoInd == 0), ] = glmnet.soft(x=Xmat2, y=res.lm, lambda=lambda)
          }
        }
      } else {
        if (sum((infoInd == 1)) == 0) {
          if (verbose) {
            cat("Incomplete information: with 0's and no 1's!  ")
          }
          beta[(infoInd == -1), ] = 0
          Xnew = Xmat[, (infoInd != -1)]
          beta[(infoInd != -1), ] = glmnet.soft(x=Xnew, y=Y, lambda=lambda)
        } else {
          if (verbose) {
            cat("Incomplete information: with both 0's and 1's!  ")
          }
          beta[(infoInd == -1), ] = 0 #known non-edges 
          Xmat1 = matrix(Xmat[, (infoInd == 1)], ncol=sum(infoInd == 1)) ##known edges 
          Xmat2 = matrix(Xmat[, (infoInd == 0)], ncol=sum(infoInd == 0)) ##unknown edges 
          fit.lm = lm(Y ~ 0 + Xmat1)
          beta[(infoInd == 1), ] = as.numeric(coef(fit.lm))
          res.lm = as.numeric(residuals(fit.lm)) ##get the residuals           
          beta[(infoInd == 0), ] = glmnet.soft(x=Xmat2, y=res.lm, lambda=lambda)
        }
      }
    }
    beta[(abs(beta) < eps)] = 0
    Adj[i, 1:(i-1)] = as.vector(beta);#if (sum(is.na(beta))>0){print(i)}; 
  }
  
  infmat = solve(Ip - Adj)
  infmat[abs(infmat) < eps] <- 0
  
  return(list(Adj=Adj, infmat=infmat))
}


get.contrast <- function(InfMat, b){
  ##InfMat: a list of ncond information matrices
  ##b: the indicator for one pathway (a row vector of 0 and 1's)
  ncond = length(InfMat)
  p = nrow(InfMat[[1]])
  LC = matrix(0, nrow=ncond - 1, ncol = ncond * p)
  for (j in 1:(ncond-1)){
    L1 = (b %*% InfMat[[j]]) * b
    L2 = (b %*% InfMat[[j+1]]) * b
    tmp = c(-L1, L2)
    LC[j, ((j-1)*p + 1): ((j+1)*p)] = tmp 
  }
  return(LC)
}


netgsa.ftest <- function(s2g, s2e, D, DtD, DDInv, n_vec, B, beta_hat){
  ncond = length(D)
  npath = dim(B)[1] 
  p = dim(B)[2]
  Ip = diag(rep(1, p))
  
  teststat2 = matrix(0, npath, 1)
  df2 = matrix(0, npath, 1)
  
  ####C = inverse of Psi ' W^{-1} Psi is calculated in the main code
  Cmat = lapply(1:ncond, function(k) (s2g*Ip + s2e * DDInv[[k]])/n_vec[k])
  Cmat = as.matrix(bdiag(Cmat))
  
  for (rr in 1:npath){  	   
    ##obtain the contrast matrix L for a given pathway
    L = get.contrast(D, B[rr,]) 
    
    ##get the test statistic
    LCL = L %*% Cmat %*% t(L)
    LCL_inv = solve(LCL)
    q = rankMatrix(L)		
    teststat2[rr] = (t(unlist(beta_hat)) %*% t(L) %*% LCL_inv %*% L %*% unlist(beta_hat))/q
    
    ##find projection P and diagonal D such that LCL' = P'DP
    tmp <- eigen(LCL)
    D_diag <- diag(tmp$values)		
    
    ## Calculate the first-order derivatives of C wrt to parameters s2g and s2e	
    g1Mat = as.matrix(bdiag(lapply(1:ncond, function(ix) Ip/n_vec[ix])))
    g2Mat = as.matrix(bdiag(lapply(1:ncond, function(ix) DDInv[[ix]]/n_vec[ix])))
    
    ##calculate the empirical covariance matrix Kmat
    Sigma = lapply(1:ncond, function(ix) s2e * Ip + s2g * DtD[[ix]] )
    SigmaInv = lapply(1:ncond, function(ix) chol2inv(chol(Sigma[[ix]])) )
    SigmaInvD = lapply(1:ncond, function(ix) SigmaInv[[ix]] %*% DtD[[ix]] )
    SinvSinv = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInv[[ix]]))
    SinvDSinvD = lapply(1:ncond, function(ix) matTr(SigmaInvD[[ix]] %*% SigmaInvD[[ix]]))
    SinvSinvD = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInvD[[ix]]))
    EH11 = (1/2) * Reduce("+", SinvDSinvD)
    EH12 = (1/2) * Reduce("+", SinvSinvD)
    EH22 = (1/2) * Reduce("+", SinvSinv)
    
    Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
    KmatInv = solve(Kmat)
    
    Em <- 0				
    for(m in 1:q){
      lm <- matrix(L[m,], nrow=1)		## lm is the mth row of L
      gm1 <- lm %*% g1Mat %*% t(lm)
      gm2 <- lm %*% g2Mat %*% t(lm)
      gm <- c(gm1, gm2)  ## gm is the gradient of lm C lm' wrt theta
      vm <- (2*D_diag[m,m]^2) / (t(gm) %*% KmatInv %*% gm )
      if(vm > 2){
        Em <- Em + (vm/(vm-2))
      }
    }
    
    ##The first degree of freedom is q. We need to calculate the 2nd degree of freedom df2. 
    if(Em > q){
      df2[rr] = (2*Em) / (Em - q)
    }
    
    ## NOTE: matlab only accepts positive integers as df's for F-dist.
    ## Therefore, I had set the denom df to 1 if it is zero, otherwise round it to integer values
    df2[rr] <- (df2[rr] >= 1) * ceiling(df2[rr]) + (df2[rr] < 1) * 1
  }
  
  pvals = 1 - pf(abs(teststat2), q, df2) + pf(-abs(teststat2), q, df2)  
  
  return(list(teststat = teststat2, df = df2, p.value = pvals))
}

netgsa.ttest <- function(s2g, s2e, D, DtD, DDInv, n_vec, B, beta) {
  ncond = length(D)
  npath = dim(B)[1]
  p = dim(B)[2]
  Ip = matrix(0, p, p)
  diag(Ip) = 1
  
  n1 = n_vec[1]
  n2 = n_vec[2]
  teststat = matrix(0, npath, 1)
  num.tstat = matrix(0, npath, 1)
  
  ##Initializing the vector for degrees of freedom for the test statistics. 
  df = matrix(0, npath, 1)
  
  ##Building the "contrast" matrix L, see Result in the paper 
  L1 = (B %*% D[[1]]) * B
  L2 = (B %*% D[[2]]) * B
  LN = cbind(-L1, L2)
  
  ##----------------- 
  ##CALCULATING DEGREES OF FREEDOM & TEST STATS 
  ##----------------- 
  #matrices needed in calculatoin of degrees of freedom
  Sigma = lapply(1:ncond, function(ix) s2e * Ip + s2g * DtD[[ix]])
  SigmaInv = lapply(1:ncond, function(ix) chol2inv(chol(Sigma[[ix]])))
  SigmaInvD = lapply(1:ncond, function(ix) SigmaInv[[ix]] %*% DtD[[ix]])
  SinvSinv = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInv[[ix]]))
  SinvDSinvD = lapply(1:ncond, function(ix) matTr(SigmaInvD[[ix]] %*% SigmaInvD[[ix]]))
  SinvSinvD = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInvD[[ix]]))
  EH11 = (1/2) * Reduce("+",SinvDSinvD)
  EH12 = (1/2) * Reduce("+",SinvSinvD)
  EH22 = (1/2) * Reduce("+",SinvSinv)
  
  #In this version of the code, K matrix is calculated directly! 
  ## Kmat will be the expected information matrix. Here we need its inverse  
  ## in calculating the degrees of freedom. 
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  # print(Kmat)
  KmatInv = solve(Kmat)
  
  #These matrices are needed in the calculation of test statistics 
  mctildi = s2e * DDInv[[1]] + s2g * Ip
  mttildi = s2e * DDInv[[2]] + s2g * Ip
  
  ##----------------- 
  for (rr in 1:npath) {
    Lrow = t(as.matrix(LN[rr, ])) #single row of L3 
    
    Lrow1 = t(as.matrix(Lrow[, 1:p]))
    Lrow2 = t(as.matrix(Lrow[, (p + 1):(2 * p)]))
    
    LC11Lprime = (1/n2) * Lrow2 %*% mttildi %*% t(Lrow2) + (1/n1) * Lrow1 %*% mctildi %*% t(Lrow1)
    
    g1 = (1/n2) * (Lrow2 %*% t(Lrow2)) + (1/n1) * (Lrow1 %*% t(Lrow1))
    g2 = (1/n2) * Lrow2 %*% DDInv[[2]] %*% t(Lrow2) + (1/n1) * Lrow1 %*% DDInv[[1]] %*% t(Lrow1)
    g = matrix(c(g1, g2), 2, 1)
    
    #test statistic 
    num.tstat[rr] = Lrow2 %*% beta[[2]] + Lrow1 %*% beta[[1]]
    teststat[rr] = num.tstat[rr]/sqrt(LC11Lprime)
    
    #calculating df based on the Satterthwaite approximation method  
    #using the formula nu=(2*LCL'^2)/g'Kg with K being the empirical covariance matrix.
    #NOTE: If df2<2, it is set to 2 
    
    df[rr] = 2 * (LC11Lprime)^2/(t(g) %*% KmatInv %*% g)
    if (df[rr] < 2) 
      df[rr] = 2
    
  }
  pvals = 1 - pt(abs(teststat), df) + pt(-abs(teststat), df)
  
  return(list(teststat = teststat, df = df, p.value = pvals))
}

newton <-
  function(x0, lb, ub, f, g, h, alpha = 0.25, beta = 0.5, max.iter = 100, tol = 1e-2){
    count = 0  
    x = x0
    
    repeat {
      count = count + 1
      
      ## newton's step
      delta = - g(x)/h(x)
      
      ## line search to pick the stepsize 
      size = 1
      while ( (x + size * delta <=0) || (log(f(x + size * delta)) > log(f(x) + alpha * size * g(x) * delta)) ) {
        size = beta * size
      }
      ## Update
      x.new = x + size * delta
      
      if ( count >= max.iter || abs(x-x.new)< tol || x.new > ub || x.new < lb ){
        
        if(count == max.iter) warning("Maximum number of iterations reached!")
        break
      }
      x = x.new
      
    }  
    
    return(list(solution = x, iter = count, stepSize = size))
  }


normalizeAdj <-
  function(Amat, alpha = 1) {
    if (is.null(dim(Amat))) {
      ncond = length(Amat)
      p = nrow(Amat[[1]])
    } else {
      ncond <- 1
      p <- dim(Amat)[1]
      Amat <- list(Amat)
    }
    
    LapMat = Amat 
    Lmat = Amat
    normA = Amat
    InfMat = Amat
    
    Ip = diag(rep(1, p), p, p)
    
    for (i in 1:ncond) {
      LapResults = graphlaplacian(Amat[[i]])
      LapMat[[i]] = LapResults$Lnorm
      Lmat[[i]] = LapResults$Lunnorm
      normA[[i]] = alpha * LapMat[[i]]
      InfMat[[i]] = t(chol(pseudoinverse(Ip - normA[[i]]))) 
    }
    
    return(list(normA = normA, Lmat = Lmat, LapMat = LapMat, InfMat = InfMat))
  }

matTr <- function(z) sum(diag(z))

profile.newton.se <-
  function(x0, D, r, control = NULL) {
    ##Note D and r are both a list of length K
    if(is.null(control)){   
      lklMethod = 'REML'      
      lb = 0.01         	
      ub = 100 		       	
      s2profile = 'se' 		
      tol = 0.01          
    } else{    
      lklMethod = control$lklMethod
      lb = control$lb
      ub = control$ub
      s2profile = control$s2profile
      tol = control$tol
    }
    
    ncond = length(D)
    if (ncond==1){stop("Only one condition detected!")}
    
    p = dim(D[[1]])[1]
    Ip = diag(rep(1, p), p, p)
    n = sapply(r, ncol)
    
    ##This is based on the notation in Lindstrom and Bates (1988)
    N = sum(n) * p 
    
    ##residual matrices
    R = vector("list", ncond)
    for (k in 1:ncond){
      R[[k]] = matrix(0, p, p)
      tmp = r[[k]]
      for ( i in 1:n[k]){
        R[[k]] = R[[k]] + tmp[, i] %o% tmp[, i]      
      }
    }
    
    DD = lapply(D, function(m) m %*% t(m))
    
    ##Calculate the eigendecomposion of DD, record the eigenvalues and the eigenvectors
    eigDD = lapply(DD, eigen)
    eigVec = lapply(eigDD, function(m) m$vectors)
    eigval = lapply(eigDD, function(m) m$val)
    
    matTr <- function(z) sum(diag(z))
    
    #obj fn
    f <- function(tau) {
      
      V = lapply(DD, function(m) tau*m + Ip)
      Vinv = lapply(1:ncond, function(k) eigVec[[k]] %*% diag(1/(1+tau*eigval[[k]])) %*% t(eigVec[[k]]))
      # Vinv = lapply(V, function(m) chol2inv(chol(m)))
      
      tmp = ifelse((lklMethod == "REML"), N - ncond*p, N)
      # val = sum(sapply(1:ncond, function(k) n[k] * as.numeric(determinant(V[[k]])$modulus) )) 
      # + tmp * log(sum(sapply(1:ncond, function(k) matTr(Vinv[[k]] %*% R[[k]]) )))
      val = sum(sapply(1:ncond, function(k) n[k] * sum(1+tau*eigval[[k]]) )) 
      + tmp * log(sum(sapply(1:ncond, function(k) matTr(Vinv[[k]] %*% R[[k]]) )))
      
      ##for REML
      if (lklMethod == "REML"){
        # First sum the matrices, then take the determinant and finally take the logarithm
        m_list = lapply(1:ncond, function(k) n[k] * (t(D[[k]]) %*% Vinv[[k]] %*% D[[k]]) )
        tmp = Reduce("+", m_list)
        val = val + as.numeric(determinant(tmp)$modulus)                        
      }
      
      return(val)
    }
    
    #gradient fn
    g <- function(tau) {
      
      V = lapply(DD, function(m) tau*m + Ip)
      Vinv = lapply(1:ncond, function(k) eigVec[[k]] %*% diag(1/(1+tau*eigval[[k]])) %*% t(eigVec[[k]]))
      C = lapply(1:ncond, function(k) t(D[[k]]) %*% Vinv[[k]] %*% D[[k]])
      trace.C = sapply(C, matTr)
      
      tmp = ifelse((lklMethod == "REML"), N - ncond*p, N)
      
      # a = sapply(1:ncond, function(k) matTr(Vinv[[k]] %*% DD[[k]] %*% Vinv[[k]] %*% R[[k]]))
      a = sapply(1:ncond, function(k) matTr(eigVec[[k]] %*% diag(eigval[[k]]/(1+tau*eigval[[k]])^2) %*% eigVec[[k]] %*% R[[k]]))
      b = sapply(1:ncond, function(k) matTr(Vinv[[k]] %*% R[[k]])) 
      
      val = sum(sapply(1:ncond, function(k) n[k]*trace.C[k] )) - tmp * sum(a)/sum(b)
      
      ##for REML
      if (lklMethod == "REML") {
        C2 = lapply(C, function(m) m %*% m)
        tmp1 = Reduce("+", lapply(1:ncond, function(k) n[k] * C[[k]]))
        tmp2 = Reduce("+", lapply(1:ncond, function(k) n[k] * C2[[k]]))
        val = val - matTr(solve(tmp1) %*% tmp2)
      }
      
      return(val)
    }
    
    #hessian fn
    
    h <- function(tau) {
      V = lapply(DD, function(m) tau*m + Ip)
      Vinv = lapply(1:ncond, function(k) eigVec[[k]] %*% diag(1/(1+tau*eigval[[k]])) %*% t(eigVec[[k]]))
      C = lapply(1:ncond, function(k) t(D[[k]]) %*% Vinv[[k]] %*% D[[k]])
      C2 = lapply(C, function(m) m %*% m)
      
      tmp = ifelse((lklMethod == "REML"), N - ncond*p, N)
      
      a = sapply(1:ncond, function(k) matTr( eigVec[[k]] %*% diag(eigval[[k]]^2/(1+tau*eigval[[k]])^3) %*% eigVec[[k]] %*% R[[k]]))
      Vinv_R = lapply(1:ncond, function(k) Vinv[[k]] %*% R[[k]])
      hes1 = sum(a)/sum(sapply(Vinv_R, matTr))
      
      a = sapply(1:ncond, function(k) matTr(eigVec[[k]] %*% diag(eigval[[k]]/(1+tau*eigval[[k]])^2) %*% eigVec[[k]] %*% R[[k]]))
      hes2 = sum(a) / sum(sapply(Vinv_R, matTr))
      
      tmp2 = Reduce("+", lapply(1:ncond, function(k) n[k] * C2[[k]]))
      val = - matTr(tmp2) + tmp * 2 * hes1 - tmp * hes2^2
      
      ##for REML
      if (lklMethod == "REML") {
        C3 = lapply(C, function(m) m %*% m %*% m)
        tmp1 = Reduce("+", lapply(1:ncond, function(k) n[k] * C[[k]]))
        tmp3 = Reduce("+", lapply(1:ncond, function(k) n[k] * C3[[k]]))
        val = val + matTr(solve(tmp1) %*% (2 * tmp3 - tmp2))
      }
      return(val)
    }
    
    tau = newton(x0, lb, ub, f, g, h, tol=tol)$solution
    V = lapply(DD, function(m) tau*m + Ip)
    m_list = lapply(1:ncond, function(k) chol2inv(chol(V[[k]])) %*% R[[k]])
    
    tmp = sum(sapply(m_list, matTr))
    
    se = ifelse((lklMethod == "REML"), (1/(N - ncond*p)) * tmp, (1/N) * tmp)
    sg = tau * se
    
    return(list(s2e = se, s2g = sg, tau = tau))
  }



approxVarEst <-  function(se0, sg0, D, r, n_vec, control=NULL){
    #D and r are both lists. 
    if (is.null(control)) {
      tol = 0.01
      s2profile = "se"
      lklMethod = "REML"
    } else {
      tol = control$tol
      s2profile = control$s2profile
      lklMethod = control$lklMethod
    }
    
    matTr <- function(z) sum(diag(z))
    
    p = nrow(D[[1]])
    Ip = diag(rep(1, p))
    ncond = length(D)
    
    n = sum(n_vec)
    N = n * p
    
    ##residual matrices
    Rm = vector("list", ncond)
    for (loop_cond in 1:ncond){
      Rm[[loop_cond]] = matrix(0, p, p)
      tmp = r[[loop_cond]]
      for (i in 1:n_vec[loop_cond]){
        Rm[[loop_cond]] = Rm[[loop_cond]] +  tmp[, i] %o% tmp[, i]
      }
    }
    DtD = lapply(1:ncond, function(k) D[[k]] %*% t(D[[k]]))
    
    gap = 1
    sg = sg0
    se = se0
    cnt = 0
    ## Whether it's ML or REML
    lklConst = ifelse(lklMethod == "REML", N - ncond*p, N)
    
    while (gap > tol) {
      sg0 = sg
      se0 = se
      
      # print(cnt)
      cnt = cnt + 1
      
      SS = lapply(1:ncond, function(k) solve(sg * DtD[[k]] + se*Ip) %*% Rm[[k]] )
      tmp0 = lapply(SS, matTr)
      tmp0 = Reduce("+", tmp0)
      
      if (s2profile == "se") {     
        se =  tmp0 / lklConst
        tmp2 = vector("list", ncond)
        for (loop_cond in 1:ncond){
          tmp1 = cov(t(r[[loop_cond]]))
          tmp = min(diag(tmp1))
          tmp = min(se, tmp)
          tmp1 = tmp1 - diag(rep(tmp, p))
          tmp2[[loop_cond]] = solve(D[[loop_cond]]) %*% tmp1  
        }
        
        tmp3 = lapply(tmp2, function(A) mean(diag(A)))      
        sg = Reduce("+", tmp3)/ncond
        tau = sg/se     
      } else {     
        sg = tmp0/lklConst
        tmp2 = vector("list", ncond)
        for (loop_cond in 1:ncond){
          tmp1 = cov(t(r[[loop_cond]]))
          tmp = min(diag(tmp1))
          tmp = min(sg, tmp)
          tmp1 = tmp1 - diag(rep(tmp, p))
          tmp2[[loop_cond]] = solve(D[[loop_cond]]) %*% tmp1  
        }
        
        tmp3 = lapply(tmp2, function(A) mean(diag(A)))  
        se = Reduce("+", tmp3)/ncond
        tau = se/sg
      }
      
      gap = abs(sg - sg0) + abs(se - se0)
    }
    
    return(list(s2e = se, s2g = sg, tau = tau, finalgap = gap, niter = cnt))
  }


zeroInd <-
  function(Amat, r){
    if (sum(t(Amat)!=Amat)>0){
      stop("This method only works for symmetric matrix!")
    }
    p <- dim(Amat)[1]
    oneMat <- matrix(0, p, p)
    zeroMat <- matrix(0, p, p)
    
    one.pos <- which(Amat!=0, arr.ind = TRUE)
    zero.pos <- which(Amat==0, arr.ind = TRUE)
    
    zero.pos <- zero.pos[which(zero.pos[,1] > zero.pos[,2]) ,]
    sel.zero <- sample(seq(1, dim(zero.pos)[1]), r * dim(zero.pos)[1], replace = FALSE) 
    zeroMat[zero.pos[sel.zero, ]] <- 1
    zeroMat <- zeroMat + t(zeroMat)  
    zeroArr <- zero.pos[sel.zero, ]
    
    out <- list()
    out$zeroArr = zeroArr
    out$zeroMat = zeroMat
    
    if (dim(one.pos)[1] == 0){
      warning("The matrix is zero!")
      out$oneMat = matrix(0, p, p)
    } else 
    {
      one.pos <- one.pos[which(one.pos[,1] > one.pos[,2]) ,]
      if (is.null(dim(one.pos))){
        one.pos = matrix(one.pos, nrow = 1)
      }
      
      sel.one <- sample(seq(1, dim(one.pos)[1]), r * dim(one.pos)[1], replace = FALSE) 
      oneMat[one.pos[sel.one, ]] <- 1
      oneMat <- oneMat + t(oneMat)
      diag(oneMat) <- 0
      
      out$oneMat = oneMat  
    }
    
    return(out)  
  }


