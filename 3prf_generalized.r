tprf <- function(x, ...) UseMethod("tprf")

tprf.default <- function(x, y, ...){
  
  x <- as.matrix()
  y <- as.numeric()
  
  est <- estimate.tprf(x, y, ...)
  
  est$call <- match.call()
  
  class(est) <- "tprf"
  est
}

autoProxy = function(y, X, L=1) {
  proxies = data.frame(matrix(NA, nrow = length(y), ncol = L))
  proxies[,1] = y 
  if (L>1) {
    for (i in 2:L) {
      testProxies = proxies[!is.na(proxies)]
      yk = estimate.3PRF(y, X, testProxies, 0) #lag??  
      r =  resid(yk$model.3PRF)
      proxies[,i] = r
    }
  }
  return(proxies)
}

estimate.3PRF = function(x, y, ...) {
  N = ncol(X)
  G = nrow(X) - lag
  
  #variance standardize to unit SD if not already standardized
  #standardize = sd(as.vector(as.matrix(X)))
  #X = X/standardize
  
  #chops off from data based on lag
#   y = y[(lag+1):elements(y)]
#   if(lag != 0) {
#     for (i in 1:lag) {  
#       X = X[-nrow(X),]
#     }
#   }
#   X = data.frame(X)
#   
  #generates autoproxies if user enters a scalar for L
  if(elements(L)==1) {
    Z = autoProxy(y, X, L)
    phi = data.frame(matrix(NA, nrow = 0, ncol = L)) 
    B = data.frame(matrix(NA, nrow = 0, ncol = L)) 
  }
  else {
    L = data.frame(L)
    Z = data.frame(L[((lag+1):nrow(L)),])
    phi = data.frame(matrix(NA, nrow = 0, ncol = ncol(L))) 
    B = data.frame(matrix(NA, nrow = 0, ncol = ncol(L))) 
  }
  
  #Step 1: Run time series regression of Xi on Z for each i = 1, ... ,N
  for (i in 1:N) {
    step1 = lm(X[,i] ~ . , data = Z)
    placeholder = data.frame(step1$coefficients)[-1,]
    phi = data.frame(rbind(phi,placeholder))  #these might not work
  }

  phi <- do.call("rbind", lapply(1:ncol(x), function(idx, X = x, Proxies = proxies) {
    (lm(X[,idx] ~ Proxies))$coef[-1]
  }))

  f <- do.call("rbind", lapply(1:nrow(x), function(idx, X = x, Phi = phi){
    lm(X[idx,] ~ Phi)$coefficients[-1]
  }))

  model <- lm(y ~ . , data = f)

  #Step 2: Run cross section regressions of Xt on phi for t = 1, ... , G
  for (i in 1:G) {
    step2 = lm(t(x[i,]) ~ . , data = data.frame(phi)) ######
    placeholder = data.frame(step2$coefficients)[-1,]
    B = data.frame(rbind(B, placeholder))
  }



  #Step 3: Run time series regression of yt+lag on predictive factors B
  threePRF = lm(y ~ . , data = B)
  return(list(model.3PRF = threePRF, step2 = B, step1 = phi, standardize = standardize))
}