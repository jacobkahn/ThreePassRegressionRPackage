#Functions for 3PRF
tprf <- function(x, ...) UseMethod("tprf")

tprf.default <- function(x, y, ...){
  
  x <- as.matrix()
  y <- as.numeric()
  
  est <- estimate.tprf(x, y, ...)
  
  est$call <- match.call()
  
  class(est) <- "tprf"
  est
}

#Generates automatic proxies. y is a numeric vector of your dependent variables
#X is a matrix/dataframe of your regressor variables
#L is an integer describing the number of automatic proxies you would like to generate
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

#Estimates 3PRF model
#y is a vector of your RHS variables, X is a matrix/dataframe of your X variables
#L is a dataframe of proxies or an integer if you would like to generate L automatic proxies
#lag is an integer describing how many periods you would like to forecast ahead
estimate.3PRF = function(x, y, L, lag = 0) {
  N = ncol(X)
  G = nrow(X) - lag
  
  #variance standardize to unit SD if not already standardized
  standardize = sd(as.vector(as.matrix(X)))
  X = X/standardize
  
  #chops off from data based on lag
   y = y[(lag+1):elements(y)]
    if(lag != 0) {
    for (i in 1:lag) {  
      X = X[-nrow(X),]
    }
}
  X = data.frame(X)
  
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

#model is the list object from the estimate.3PRF function
#X is a dataframe of regressors you'll use to forecast 
#obviously, should be the same variables as you used to estimate the model
forecast.3PRF = function(model, X) {
  X = X/model$standardize #this might be an issue do you include new OB?? I don't think so but I'm not sure
  G = nrow(X)
  H = data.frame(matrix(NA, nrow = 0, ncol = length(model$model.3PRF$coefficients)- 1))
  for (i in 1:T) {
    step2 = lm(t(X[i,]) ~ . , data=model$step1) ######
    placeholder = data.frame(step2$coefficients)[-1,]
    H = data.frame(rbind(F, placeholder))
  }
  #   forecast = predict(model$model.3PRF, F, interval = "prediction") predict() is a POS function
  intercept = rep(1,nrow(H))
  newF = data.frame(cbind(intercept,H))
  coefficients = model$model.3PRF$coefficients
  prediction = as.matrix(newF) %*% as.matrix(coefficients)
  se = summary(model$model.3PRF)$sigma
  lower = prediction - 1.96 * se
  upper = prediction + 1.96 * se
  forecast = data.frame(cbind(prediction, se, lower, upper))
  colnames(forecast) = c("forecast", "S.E.", "Lower", "Upper")
  return(forecast)
}

#re-estimates model and computes a new forecast each period, uses a rolling window       
rollcast.3PRF = function(y, X, L, lag, window) {
  forecast = data.frame(matrix(NA, nrow = 0, ncol = 1))
  
  y = y[(lag+1):length(y)]
  nolagX = X
  for (i in 0:(lag-1)) {
    X = X[-nrow(X),]
  }
  
  if(elements(L) > 1) {
    L = data.frame(L)
    Z = data.frame(L[((lag+1):length(y)),])
  }
  
  for(i in 1:(length(y) - window + 1)) {
    y = data.frame(y)
    X = data.frame(X)
    ywin = data.frame(y[(i:(i+window-1)),])
    Xwin = data.frame(X[(i:(i+window-1)),])
    ywin = t(t(ywin))
    if(elements(L) == 1) {
      #Zwin = autoProxy(ywin, Xwin, L)
      estimate = estimate.3PRF(ywin, Xwin, L, 0)
    }
    else {
      Z = data.frame(Z)
      Zwin = data.frame(Z[(i:(i+window-1)),])
      estimate = estimate.3PRF(ywin, Xwin, Zwin, 0)
    }
    #Zwin = data.frame(Zwin)
    #Xwin = data.frame(Xwin)
    passthru = nolagX[(i+window+lag-1),]
    periodforecast = forecast.3PRF(estimate, passthru)
    forecast = data.frame(rbind(forecast,periodforecast))
  } 
  return(forecast)
}

#re-estimates model and computes a new forecast each period, uses recursive estimation
#initial training period is specified by integer variable "window"
recursivecast.3PRF = function(y, X, L, lag, window) {
  forecast = data.frame(matrix(NA, nrow = 0, ncol = 1))
  
  y = y[(lag+1):length(y)]
  nolagX = X
  for (i in 0:(lag-1)) {
    X = X[-nrow(X),]
  }
  
  if(elements(L) > 1) {
    L = data.frame(L)
    Z = data.frame(L[((lag+1):length(y)),])
  }
  
  for(i in 1:(length(y) - window + 1)) {
    y = data.frame(y)
    X = data.frame(X)
    ywin = data.frame(y[(1:(i+window-1)),])
    Xwin = data.frame(X[(1:(i+window-1)),])
    ywin = t(t(ywin))
    if(elements(L) == 1) {
      #Zwin = autoProxy(ywin, Xwin, L)
      estimate = estimate.3PRF(ywin, Xwin, L, 0)
    }
    else {
      Z = data.frame(Z)
      Zwin = data.frame(Z[(1:(i+window-1)),])
      estimate = estimate.3PRF(ywin, Xwin, Zwin, 0)
    }
    #Zwin = data.frame(Zwin)
    #Xwin = data.frame(Xwin)
    passthru = nolagX[(i+window+lag-1),]
    periodforecast = forecast.3PRF(estimate, passthru)
    forecast = data.frame(rbind(forecast,periodforecast))
  } 
  return(forecast)
}

elements <- function(x) {
  if(is.list(x)) {
    do.call(sum,lapply(x, elements))
  } else {
    length(x)
  }
}
