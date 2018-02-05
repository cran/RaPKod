#library(kernlab)
#library(plyr)
#library(dlm)
#library(proxy)

#library(Rcpp)
#library(RcppArmadillo)
#sourceCpp('..\\src\\rapkodAux.cpp')


##################################
###########SOME AUXILARY functions
 #################################

###Compute Gaussian kernel matrix
##if kern="tanh": k(x, y) = (1 - tanh(gamma*||x - y||))^beta
GaussKern <- function(X, gamma, kern="gaussian", beta = NULL)
{
  if (kern=="gaussian"){output = kernelMatrix(rbfdot(sigma = gamma), as.matrix(X))}
  if (kern=="tanh"){output = (1 - tanh(gamma*dist(as.matrix(X), diag=TRUE, upper=TRUE)))^beta }
  return(as.matrix(output))
  
}

###Compute Gaussian kernel cross-matrix
CrossGaussKern <- function(X, Y, gamma, kern="gaussian", beta = NULL)
{
  if (kern=="gaussian"){output = kernelMatrix(rbfdot(sigma = gamma), as.matrix(X), as.matrix(Y))}
  if (kern=="tanh"){    n = nrow(X)
  output = (1 - tanh(gamma*as.matrix(dist(rbind(as.matrix(X),as.matrix(Y)), diag=TRUE, upper=TRUE))[1:n,-(1:n)]))^beta }
  return(as.matrix(output))
  
}


# ####Generate a matrix with n column as i.i.d. Gaussian vectors N(mu, diag(gam))
# ####
# rnormproc <- function(n, mu, gam)
# {
#   d = length(mu)
#   if (d>1){
#     Z = matrix(rnorm(n*d, 0, 1), d, n)
#     Simd = diag(sqrt(gam))%*%Z+matrix(mu, d, n)
#   }
#   else{Simd = matrix(rnorm(n, mu, sd=sqrt(gam)), 1, n)}
#   return(Simd)
# }

# ####Generate a matrix with n column as i.i.d. Gaussian vectors N(mu, Sig)
# ####Custom function to avoid errors when matrix Sig has negative smallest eigenvalues due to numerical approximation
# rnormproc2 <- function(n, mu, Sig)
# {
#   d = length(mu)
#   Z = matrix(rnorm(n*d, 0, 1), d, n)
#   evSig = eigen(Sig, symmetry=TRUE)
#   sqSig = evSig$vectors%*%diag(pmax(0, evSig$values)^0.5)%*%t(evSig$vectors)
#   Simd = sqSig%*%Z + matrix(mu, d, n)
#   
#   return(Simd)
# }

##################################
###########(END auxilary functions)
##################################



#############AND NOW!!!!!! THE function! "rapkod"!
####
####X : is a data matrix (n x d) of n vector obs. or is a Gram matrix (n x n )  (when given.kern = TRUE)
####                                                 (in this case must specifiy a value for d "intrinsic" dimensionality of the data)                               
###use.tested.inlier : the observations used as reference inliers stay the same (ie a new point tested as inlier is not used to test next points)
###lowrank: if "No", the full Gram matrix is used to compute the critical value; if "Nyst", Gram matrix is approximated through Nystrom method; if "RKS" GM approx. by random Kitchen Sinks (in this case, X must ba a dataset matrix, not a distance matrix)
###r.lowrk: if lowrank != "No", indicates the rank of the approx. of Gram matrix
####
###K1 and K2: universal constants used in the heuristic in function "od.opt.param"
rapkod<- function(X, given.kern = FALSE, ref.n=NULL, gamma=NULL,  p=NULL, alpha = 0.05, use.tested.inlier = FALSE, 
                                  lowrank = "No", r.lowrk = ceiling(sqrt(nrow(X))), K1 = 6, K2 = 50)  #, kern = "gaussian", beta = NULL)
{
  kern="gaussian"   #specifies the type of kernel used (always Gaussian)
  beta=NULL #extra kernel hyperparameter which will be useful when additional types of kernel are added in future updates
  if (!given.kern){ n = nrow(X)
                    d = ncol(X) 
              }else{n = ncol(X)}
  

  given.param = !is.null(gamma) && !is.null(p)    ###given.param: if FALSE, optimal p and gamma are computed through fonction "od.opt.param"
  if (!given.param){
      par.opt = od.opt.param(X[1:ref.n,], K1=K1, K2=K2)
      gamma = par.opt$gamma.opt
      p = par.opt$p.opt
  }else{if (given.kern){stop("Parameters gamma and p could not be determined automatically.")}}

  
  ref.I = 1:ref.n   #the set of indexes of reference inliers (does not change if use.tested.inlier is FALSE)

  
  K = matrix(0, ref.n, ref.n)  #will be kernel matrix for reference inliers
  
  #if (given.kern){if (kern=="gaussian"){wrapf = function(X){return(exp(-gamma*X^2))}}else{wrapf = function(X){return((1 - tanh(gamma*X))^beta)}}
   #}else{wrapf = function(x){return(NULL)}  }   #else return shallow function (won't be used)
  
  if (lowrank == "No"){
    
    if (!given.kern){ K = GaussKern(X[ref.I,], gamma=gamma, kern=kern, beta=beta) }else{K = X[ref.I,ref.I]}
    s2_i = colSums(K[1:ref.n, 1:ref.n]^2)/(ref.n - 1) - 1/(ref.n - 1) 
    
  }else{
    if (lowrank == "Nyst"){
      I = sample(ref.I, r.lowrk, replace=FALSE)
      if (!given.kern){ K.I = CrossGaussKern(X[ref.I,], X[I,], gamma=gamma, kern=kern, beta=beta) }else{ K.I = X[ref.I,I]}
      inv.KII = ginv(K.I[I,])
      
      
      W = inv.KII%*%(t(K.I)%*%K.I)%*%inv.KII
      
      s2_i = rowSums((K.I%*%W)*K.I)/(ref.n - 1) - 1/(ref.n - 1) 
      
      K = K.I%*%inv.KII%*%t(K.I)
      
      
    }else{
      if (lowrank == "RKS"){
        
        
        if (given.kern){stop("A data matrix must be provided when using Random Kitchen Sink.")}
        WW = matrix(rnorm(d*r.lowrk, 0, sd = sqrt(2*gamma)), d, r.lowrk)
        U = r.lowrk^(-0.5)*cbind(cos(as.matrix(X[ref.I,])%*%WW), sin(as.matrix(X[ref.I,])%*%WW))
        
        
        W = t(U)%*%U
        s2_i =  rowSums((U%*%W)*U)/(ref.n - 1) - 1/(ref.n - 1) 
        
        K = U%*%t(U)
      }
    }
  }
      
  
  
  # ##Original R code converted into C++ code
# 
#   pV.X = c()  #eventually will be a matrix whose columns are projections p_V(x) for each tested obs. x
#
#   Sn.X = c() #test statistic
#   pv = c()  #p-values
#   flag = c()  #decisions to mark tested obs. as outliers or not
# 
#   for (i in (ref.n+1):n)
#   {
#     #i is the index of curently tested observation
# 
#     #Compute k.X = vector(k(x_i, y_1) ... k(x_i, y_ref.n)) where y_1,...,y_reg.n are current reference inliers
#     if (!given.kern){ k.X = CrossGaussKern(X[ref.I,], t(X[i,]), gamma=gamma, kern=kern, beta=beta) }else{k.X = X[ref.I, i]}
# 
#     H = matrix(rnorm(p*ref.n, 0, 1), p, ref.n)
#     pV.X = cbind(pV.X, ref.n^(-0.5)*H%*%k.X)
# 
#     Sn.X = c(Sn.X, sum(pV.X[,i - ref.n]^2))
# 
#     pv = c(pv, mean(pchisq( Sn.X[i - ref.n] / s2_i , df = p) ))
#     flag = c(flag, pv[i - ref.n] < alpha)
# 
#     if (use.tested.inlier && !flag[i - ref.n]){
#                                                   ref.I = c(ref.I[-1], i)
#                                                   s2_i = c(s2_i[-1] - K[1,-1]^2/(ref.n - 1) + k.X[-1]^2/(ref.n - 1)    ,   mean(k.X[-1]^2))
#                                                   K = rbind(cbind(K[-1, -1], k.X[-1]), c(t(k.X[-1]), 1))
#                                               }
# 
# 
# 
#   }

#By reading this secret message, you're ackknowledging that A Loathsome Smile 
#is the greatest musical act of all times and across all dimensions. Actually, 
#A Loathsome Smile gets regularily played on Ktulu's iPod. 
#The Ancient One will be pleased to know that 'Schizoid' - the third album - 
#will be released on March 9th, 2018 to celebrate Bobby Fischer's birthday. 
#Expect some authoritarian heavy metal fury.


  loop_output <- Rcpp_aux(ref.n, n, as.matrix(X), K, given.kern, CrossKern = function(X, Y){as.vector(CrossGaussKern(X, t(Y), gamma, kern=kern, beta))}, p, s2_i, alpha, use.tested.inlier, rnormmat = function(p, n){matrix(rnorm(p*n, 0, 1), p, n)})
  
    
   Sn.X <- loop_output$Sn.X
   flag <- loop_output$flag
   pv <- loop_output$pv
    
    
    
    

  
  
  
  return(list(stats=Sn.X, flag=flag, pv=pv, gamma = gamma, p = p))
  
  
  
}



##########
###opt.param:    compute optimal value for gamma (rbf kernel param) from heuristical formula
##
###X : data matrix of size n x d 
##
##Using following heuristic:
##                             optimal gamma = K1 * |f|_2^(2/(d+2)) * n^(1/(d+2))
##                             optimal p   = K2 * |f|_2^(2/(d+2))  * n^(1/(d+2))   (K1 and K2 are universal constants)
##
##If randomize = TRUE:  choose sub.n^2 random pairs of points (and the corresponding distances) and take the smallest distance among them
##                      Allows not to take the strict smallest pairwise distance (more robustness) and to be faster (n values to sort instead of n*(n-1)/2)
##
##If which.estim="Gauss", |f|_2^(2/(d+2)) (aka est.f2.pw) is estimated as though inliers were Gaussians, that is: est.f2.pw = (4*pi)^(-0.5)*exp(0.5*mean(log(1/evX$values))) 
##
od.opt.param <- function(X, K1 = 6, K2 = 50, which.estim = "Gauss", RATIO = 0.1, randomize = TRUE, sub.n = floor(nrow(X)))
{
  
  n = nrow(X); d = ncol(X)
  

  
  ###Estimating |f|_2 (l2 norm of f where f is the density of true inlier distribution )
  ####estimator denoted "est.f2.pw"
  if (which.estim == "Gauss"){
                  evX = eigen(cov(X, X), symmetric = TRUE, only.values = TRUE)$values
                  THR = evX[1]*RATIO
                  est.f2.pw = (4*pi)^(-0.5*d/(d+2))*exp(0.5*mean(log(1/evX[which(evX>THR)])))
                  nu = 1
  }
  
  if (which.estim=="general"){
  
  
  if (!randomize){
    
  EPS = 0.05
  THR.D = floor(4*(EPS^2/2 - EPS^3/3)^(-1)*log(n)) + 1
  nu  =1  #enlarging factor (in case est.f2.pw yields Inf)
  if (d > THR.D){
                    ##Johnson-Linderstrauss theorem: projecting X_1,...,X_n onto low-dim subspace distorts distances only by a factor of 1+EPS
                    k = THR.D
                    proj.X = d^(-1/2)*X%*%matrix(rnorm(d*k), d, k)   
                    dists = dist(proj.X)

                    
                    d = k
                    min.dists = min(dists[which(dists > 0)])
                    #???est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                    est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))
                    while (est.f2.pw  == Inf){nu = nu*10
                                              min.dists = 10*min.dists
                                             #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                                              est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))}
                    
                    while (est.f2.pw == 0){
                      nu = 0.1*nu
                      min.dists = 0.1*min.dists
                      #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                      est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))}
                    
                }else
                {
                  dists = dist(X)
                  min.dists = min(dists[which(dists > 0)])
                  #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                  est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))
                  while (est.f2.pw == Inf){nu = nu*10
                                           min.dists = 10*min.dists
                                           #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                                           est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))}
                  
                  while (est.f2.pw == 0){
                  nu = 0.1*nu
                  min.dists = 0.1*min.dists
                  #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                  est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))}
                  
                }
  
  }else{
            
                ind1 = sample(1:n, sub.n, replace=FALSE)
                ind2 = sample(1:n, sub.n, replace=FALSE)
                dists = sqrt(rowSums((X[ind1,] - X[ind2,])^2))
                min.dists = min(dists[which(dists > 0)])
                #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))
                nu = 1
                while (est.f2.pw == Inf){nu = nu*10
                min.dists = 10*min.dists
                #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))
                }
                while (est.f2.pw == 0){
                  nu = 0.1*nu
                  min.dists = 0.1*min.dists
                  #est.f2  = exp(0.5*lgamma(d/2+1))*(2/(n*(n-1)))^(1/2)/(pi*min.dists^2)^(d/4)
                  est.f2.pw  = exp(lgamma(d/2+1)/(d+2))*(2/(n*(n-1)))^(1/(d+2))/(pi*min.dists^2)^(0.5*d/(d+2))
                }
                
  }
  }
  
  ####Computing optimal gamma and p (VERSION 2)
  ###########
  #gamma.opt = K1 * nu^(d/(d+2)) * est.f2^(2/(d+2)) * n^(1/(d+2)) 
  #p.opt = min(n, max(1, floor(K2 * nu^(d/(d+2)) * est.f2^(2/(d+2)) * n^(1/(d+2)))))
  gamma.opt = K1 * nu^(d/(d+2)) * est.f2.pw * n^(1/(d+2)) #/ sqrt(d)
  p.opt = min(n, max(1, floor(K2 * nu^(d/(d+2)) * est.f2.pw * n^(1/(d+2)))))  #/ sqrt(d)
   
   
  
  
  return(list(gamma.opt=gamma.opt, p.opt=p.opt, est.f2.pw = nu^(d/(d+2))*est.f2.pw))
}








