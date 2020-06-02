lbm.NR <- function(x,y, vce="oim",mu=NULL,
                  start=NULL, noconstant=TRUE,
                  control,...)
{
  if(!is.null(mu)) {
    res.glm<-suppressWarnings(stats::glm.fit(x=x,y=y, family = binomial(link = log),
                              mustart=mu,control=control,...))
  }
  if(!is.null(start)) {
    res.glm<-suppressWarnings(stats::glm.fit(x=x,y=y, family = binomial(link = log),
                              start= start,control=control,...))
  }
  if(noconstant==F) {
    res.temp<-suppressWarnings(stats::glm.fit(x=x,y=y, family = binomial(link = log),
                              control=control,start= c(-1,rep(0,length(start)-1)),...))
    # print(res.temp$coefficients)
    if(res.temp$deviance<res.glm$deviance) res.glm<-res.temp
    # if(res.temp$deviance<res.glm$deviance) print(1) else print(2)
  }
  beta_old<-res.glm$coefficients
  # print(res.glm$coefficients)
  dev_old<-res.glm$deviance
  # name<-names(beta_old)
  u <- res.glm$fitted.values
  for(i in 1:25)
  {
    res_NR<-tryCatch({
      lbm.NR.fit(y=y,x=x,u=u,beta=beta_old)
      },
      error=function(err) {
        return(0)
      })
    if(is.list(res_NR)) {
      beta_new<-res_NR$beta
      u <- res_NR$u
      if(validmu(u)) dev_new<-lbm.NR.dev(y=y,u=u)
    } else {
      beta_new<-beta_old
      u <- u
      dev_new<-dev_old
    }
    if(!validmu(mu=u)) {
      res_retr<-tryCatch({
        lbm.retrieve(y=y,x=x,beta_new = beta_new,beta_old = beta_old)
      },
        error=function(err) {
          return(err<-0)
        })
      if(is.list(res_retr)) {
        beta_new<-res_retr$beta
        u<-res_retr$u
        dev_new<-lbm.NR.dev(y=y,u=u)
      } else {
        beta_new<-beta_old
        u <- res.glm$fitted.values
        dev_new<-dev_old
      }
    }
    if(abs(dev_new - dev_old)/(0.1 + dev_new) < control$epsilon) {
      coef<-as.vector(beta_new)
      dev<-dev_new
      fitted.values<-u
      names(fitted.values)<-row.names(x)
      break
    } else {
      beta_old<-beta_new
      u <- as.vector(exp(as.matrix(x) %*% beta_new))
      dev_old<-dev_new
    }
  }
  names(coef)<-colnames(x)
  if(dev_new!=dev_old) {

    vcov<-lbm.NR.vcov(x=x,y=y,fitted.values = fitted.values,
                      vce=vce, noconstant = noconstant)
  } else {
      p <- res.glm$rank
      p1 <- 1L:p
      Qr <- res.glm$qr
      vcov <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
  }
  dimnames(vcov)<-list(colnames(x),colnames(x))
  result <- list(coefficients=coef,deviance=dev,
                 vcov=vcov, fitted.values=fitted.values,
                 linear.predictors=log(fitted.values),
                 null.deviance=res.glm$null.deviance)
  return(result)
}

lbm.NR.fit<-function(y=NULL, x=NULL, u=NULL, beta=NULL)
{
  w<-ifelse((1-y)==0, 0, -(u*(1-y))/((1-u)^2))
  # Establish gradian vactor.
  g_temp<-ifelse((y-u)==0, 0, as.vector((y - u) / (1 - u)))
  g <-t(x) %*% g_temp
  # Observed Hessian matrix.
  H <- t(x) %*% diag(w) %*% as.matrix(x)
  # Variance-covariance matrix.
  vm <- -solve(H)
  H_temp <- as.vector(vm %*% g)
  # Update betas by Newton-Raphson.
  beta_new <- beta + H_temp
  u <- as.vector(exp(as.matrix(x) %*% beta_new))
  res<-list(beta=beta_new, u=u)
  return(res)
}

# Checking the validation of estimators.
validmu <- function(mu=NULL)
{
  if(any(mu > 1)) {
    ind <- FALSE
  } else if(any(mu < 0)) {
    ind <- FALSE
  } else {
    ind <- TRUE
  }
  return(ind)
}

# Deviance calculation.
lbm.NR.dev<-function(y=NULL,u=NULL)
{
  likelihood<-ifelse((1-y)==0, log(u), log(1-u))
  dev<--2*sum(likelihood)
  return(dev)
}

lbm.retrieve<-function(y=NULL,x=NULL,beta_new=NULL,beta_old=NULL)
{
  for(m in 1:20)
  {
    for(j in 1:10)
    {
      beta_new<-step_halving(beta_old = beta_old,beta_new=beta_new)
      u <- as.vector(exp(as.matrix(x) %*% beta_new))
      if(validmu(mu=u)) break
    }
    if(!(validmu(mu=u))) {
      u[which(u>=1)]<-0.9999
      res_NR<-lbm.NR.fit(y=y,x=x,u=u,beta=beta_new)
      beta_new<-res_NR$beta
      u <- res_NR$u
    }
    if(validmu(mu=u)) break
    if(m==20 & !(validmu(mu=u))) {
      stop("Fail to converge.")
    }
  }
  res<-list(u=u,beta=beta_new)
  return(res)
}

# Step halving method of R : take average between
# current estimator values and previous admissible
# estimator values.
step_halving <- function(beta_old = NULL, beta_new = NULL)
{
  # Check for fitted values outside domain.
  beta_new <- (beta_old+beta_new )/ 2
  return(beta_new)
}

lbm.NR.vcov<-function(x=x,y=y, fitted.values=NULL,
                      vce="oim", noconstant=TRUE)
{
  u <- fitted.values
  if(tolower(vce)=="eim") {
    w <- ifelse(1-u==0, 0, u / (1 - u))
  } else {
    w <- ifelse(1-y==0, 0, u * (1 - y) / (1 - u)^2)
  }
  H <- t(x) %*% diag(w) %*% as.matrix(x)
  vcov <- -solve(-H,tol = 1e-20)
  return(vcov)
}
