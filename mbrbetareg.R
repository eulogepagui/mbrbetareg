

########################################################################
########---------------  beta regression  script  ######################
########################################################################


##  the function make.dmu.deta is implemented in the "R package betareg" version 3.1-0 
## by Ioannis Kosmidis  and others
## Euloge C. Kenne Pagui <kenne@stat.unipd.it> [03/02/2017]

dyn.load("modification.so")
 

p.range = function(p, eps=10*.Machine$double.eps) 
{
   out = p;
   out[p<eps] = eps;
   out[p>(1-eps)] = (1-eps);
   return(out);
}

make.dmu.deta <- function(linkstr) switch(linkstr,
  "logit" = {
    logit_link <- make.link("logit")
    function(eta) logit_link$mu.eta(eta) * (1 - 2 * logit_link$linkinv(eta))
  },
  "probit" = function(eta) -eta * pmax(dnorm(eta), .Machine$double.eps),
  "cauchit" = function(eta) -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps),
  "cloglog" = function(eta) pmax((1 - exp(eta)) * exp(eta - exp(eta)), .Machine$double.eps),
  "loglog" = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
  "identity" = function(eta) rep.int(0, length(eta)),
  "log" = function(eta) pmax(exp(eta), .Machine$double.eps),
  "sqrt" = function(eta) rep.int(2, length(eta)),
  "1/mu^2" = function(eta) 3/(4 * eta^2.5),
  "inverse" = function(eta) 2/(eta^3)
)


########################################################################
##-- negative log-likelihood       ----------------------------------###
########################################################################

##      par:  vector parameter 
##        y:  response variable
##        X:  design matrix related to the mean
##        Z:  design matrix related to the precision parameter
##     link:  link function for the mean
## phi_link:  link function for the precision parameter

nllik = function(par,y,X,Z,link="logit",phi_link="log")
{
   linkobj <- make.link(link)
   linkinv <- linkobj$linkinv
   phi_linkobj <- make.link(phi_link)
   phi_linkinv <- phi_linkobj$linkinv
   p <- ncol(X)
   q <- ncol(Z)
   beta <- par[1:p]
   delta <- par[(p+1):(p+q)]
   eta <- drop(X%*%beta)
   zeta <- drop(Z%*%delta)
   mu <- linkinv(eta)
   phi <- phi_linkinv(zeta)
   TT <- log(y)
   U <- log(1-y)
   res= sum(phi*mu*(TT-U)+phi*U-lgamma(phi*mu)-lgamma(phi*(1-mu))+lgamma(phi))
   return(-res)
}   


######################################################
##########------score function-----------------#######
######################################################

score <- function(par,y,X,Z,link="logit",phi_link="log")
{
   linkobj <- make.link(link)
   linkinv <- linkobj$linkinv
   mu.eta <- linkobj$mu.eta
   phi_linkobj <- make.link(phi_link)
   phi_linkinv <- phi_linkobj$linkinv
   phi_mu.zeta <- phi_linkobj$mu.eta
   p <- ncol(X)
   q <- ncol(Z)
   beta <- par[1:p]
   delta <- par[(p+1):(p+q)]
   eta <- drop(X%*%beta)
   zeta <- drop(Z%*%delta)
   mu <- linkinv(eta)
   phi <- phi_linkinv(zeta)
   D1 <- mu.eta(eta)
   D2 <- phi_mu.zeta(zeta)
   lambda <- digamma(phi*mu)-digamma(phi)
   Tbar <- log(y)-lambda
   epsilon <- digamma(phi*(1-mu))-digamma(phi)
   Ubar <- log(1-y)-epsilon
   c(crossprod(X,phi*D1*(Tbar-Ubar)),crossprod(Z,D2*(mu*(Tbar-Ubar)+Ubar)))
} 

######################################################
#####--------Fisher information-----------------######
######################################################

info = function(par,X,Z,link="logit",phi_link="log")
{
   linkobj <- make.link(link)
   linkinv <- linkobj$linkinv
   mu.eta <- linkobj$mu.eta
   phi_linkobj <- make.link(phi_link)
   phi_linkinv <- phi_linkobj$linkinv
   phi_mu.zeta <- phi_linkobj$mu.eta  
   p <- ncol(X)
   q <- ncol(Z)
   n <- nrow(X)
   beta <- par[1:p]
   delta <- par[(p+1):(p+q)]
   eta <- drop(X%*%beta)
   zeta <- drop(Z%*%delta)
   mu <- linkinv(eta)
   # M <- diag(mu)
   phi <- phi_linkinv(zeta)
   # PHI <- diag(phi)
   D1 <- mu.eta(eta)
   D2 <- phi_mu.zeta(zeta)
   K2 <- trigamma(phi*mu)+trigamma(phi*(1-mu))
   OMEGA1 <- trigamma(phi)
   PSI1 <- trigamma(phi*(1-mu))
   one_n = rep(1,n)
   # ans <- matrix(NA,(p+q),(p+q))
   bb<- crossprod(X,D1^2*phi^2*K2*X)  
   bg <- crossprod(X,phi*D1*D2*(mu*K2-PSI1)*Z)
   gg <- crossprod(Z,D2^2*(mu^2*K2+(1-2*mu)*PSI1-OMEGA1)*Z)
   rbind(cbind(bb,bg),cbind(t(bg),gg))
}


nu.stu2 = function(par,X,Z,link="logit",phi_link="log")
{
   linkobj <- make.link(link)
   linkinv <- linkobj$linkinv
   mu.eta <- linkobj$mu.eta
   mu.eta.eta <-make.dmu.deta(link)
   phi_linkobj <- make.link(phi_link)
   phi_linkinv <- phi_linkobj$linkinv
   phi_mu.zeta <- phi_linkobj$mu.eta  
   phi_mu.zeta.zeta <-make.dmu.deta(phi_link)
   p <- ncol(X)
   q <- ncol(Z)
   n <- nrow(X)
   beta <- par[1:p]
   delta <- par[(p+1):(p+q)]
   eta <- drop(X%*%beta)
   zeta <- drop(Z%*%delta)
   mu <- linkinv(eta)
   # M <- diag(mu)
   phi <- phi_linkinv(zeta)
   # PHI <- diag(phi)
   D1 <- mu.eta(eta)
   Dp1 <- mu.eta.eta(eta)
   D2 <- phi_mu.zeta(zeta)
   Dp2 <- phi_mu.zeta.zeta(zeta)
   K2 <- trigamma(phi*mu)+trigamma(phi*(1-mu))
   K3 <- psigamma(phi*mu, 2)-psigamma(phi*(1-mu),2)
   OMEGA1 <- trigamma(phi)
   OMEGA2 <- psigamma(phi,  2)
   PSI1 <- trigamma(phi*(1-mu))
   PSI2 <- psigamma(phi*(1-mu), 2)
   
   ans <- array(NA,c((p+q),(p+q),c(p+q)))
   for(t in 1:p)
   {
      Xt<-X[,t]
      bb<- crossprod(X,phi^3*D1^3*K3*Xt*X)
      bg <- crossprod(X,phi^2*D1^2*D2*(mu*K3+PSI2)*Xt*Z)
      # ans[(p+1):(p+q),1:p,t] <- t(ans[1:p,(p+1):(p+q),t])
      gg <- crossprod(Z,phi*D1*D2^2*(mu^2*K3+2*mu*PSI2-PSI2)*Xt*Z)
      ans[t,,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
   }  
   for(s in 1:q)
   {
      Zs<-Z[,s]
      bb <- crossprod(X,phi^2*D1^2*D2*(mu*K3+PSI2)*Zs*X)
      bg <- crossprod(X,phi*D1*D2^2*(mu^2*K3+2*mu*PSI2-PSI2)*Zs*Z)
      gg<- crossprod(Z,D2^3*(mu^3*K3+(3*mu^2-3*mu+1)*PSI2-OMEGA2)*Zs*Z)
      ans[(p+s),,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
   }   
   ans
}





nu.s.tu2 = function(par,X,Z,link="logit",phi_link="log")
{
   linkobj <- make.link(link)
   linkinv <- linkobj$linkinv
   mu.eta <- linkobj$mu.eta
   mu.eta.eta <-make.dmu.deta(link)
   phi_linkobj <- make.link(phi_link)
   phi_linkinv <- phi_linkobj$linkinv
   phi_mu.zeta <- phi_linkobj$mu.eta  
   phi_mu.zeta.zeta <-make.dmu.deta(phi_link)
   p <- ncol(X)
   q <- ncol(Z)
   n <- nrow(X)
   beta <- par[1:p]
   delta <- par[(p+1):(p+q)]
   eta <- drop(X%*%beta)
   zeta <- drop(Z%*%delta)
   mu <- linkinv(eta)
   # M <- diag(mu)
   phi <- phi_linkinv(zeta)
   # PHI <- diag(phi)
   D1 <- mu.eta(eta)
   Dp1 <- mu.eta.eta(eta)
   D2 <- phi_mu.zeta(zeta)
   Dp2 <- phi_mu.zeta.zeta(zeta)
   K2 <- (trigamma(phi*mu)+trigamma(phi*(1-mu)))
   K3 <- (psigamma(phi*mu,2)-psigamma(phi*(1-mu),2))
   OMEGA1 <- trigamma(phi)
   OMEGA2 <- psigamma(phi,  2)
   PSI1 <- trigamma(phi*(1-mu))
   PSI2 <- psigamma(phi*(1-mu), 2)
   
   ans <- array(NA,c((p+q),(p+q),c(p+q)))
   for(t in 1:p)
   {
      Xt<-X[,t]
      bb<- crossprod(X,phi^2*D1*Dp1*K2*Xt*X)
      bg <- crossprod(X,phi*D1^2*D2*K2*Xt*Z)
      gg <- crossprod(Z,phi*D1*Dp2*(mu*K2-PSI1)*Xt*Z)
      ans[t,,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
   }  
   for(s in 1:q)
   {
      Zs<-Z[,s]
      bb <- crossprod(X,phi*D2*Dp1*(mu*K2-PSI1)*Zs*X)
      bg <- crossprod(X,D1*D2^2*(mu*K2-PSI1)*Zs*Z)
      gg <- crossprod(Z,D2*Dp2*(mu^2*K2+PSI1-2*mu*PSI1-OMEGA1)*Zs*Z)
      ans[(p+s),,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
   }   
   ans
}





######################################################
##---------- modification -------------------------###
######################################################



modification.c = function(par,X,Z,link="logit",phi_link="log")
{
   p <- ncol(X)
   q <- ncol(Z)
   info = info(par,X,Z,link,phi_link)
   InfoInv <- try(chol2inv(chol(info)),TRUE)
   if(failedInv <- inherits(InfoInv, "try-error")) 
   {
     warning("failed to invert the information matrix: iteration stopped prematurely")
     break
   }
   nu_r_s_t=nu.stu2(par,X,Z,link,phi_link)
   nu_r_st=nu.s.tu2(par,X,Z,link,phi_link)
   mod <- rep(0,(p+q))
   out <- .C('modificationc7',
       as.integer(p+q),
       as.double(InfoInv),
       as.double( nu_r_s_t),
       as.double(nu_r_st),
       mod=as.double(mod))
   out$mod
}


######################################################
##-------------Modified score function-------------###
##-------------------  AU+M -----------------------###
######################################################

scoremc=function(par,y,X,Z,link="logit",phi_link="log")
{
  linkobj <- make.link(link)
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  mu.eta.eta <-make.dmu.deta(link)
  phi_linkobj <- make.link(phi_link)
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.zeta <- phi_linkobj$mu.eta  
  phi_mu.zeta.zeta <-make.dmu.deta(phi_link)
  p <- ncol(X)
  q <- ncol(Z)
  n <- nrow(X)
  beta <- par[1:p]
  delta <- par[(p+1):(p+q)]
  eta <- drop(X%*%beta)
  zeta <- drop(Z%*%delta)
  mu <- linkinv(eta)
  # M <- diag(mu)
  phi <- phi_linkinv(zeta)
  # PHI <- diag(phi)
  D1 <- mu.eta(eta)
  Dp1 <- mu.eta.eta(eta)
  D2 <- phi_mu.zeta(zeta)
  Dp2 <- phi_mu.zeta.zeta(zeta)
  K2 <- (trigamma(phi*mu)+trigamma(phi*(1-mu)))
  K3 <- (psigamma(phi*mu,2)-psigamma(phi*(1-mu),2))
  OMEGA1 <- trigamma(phi)
  OMEGA2 <- psigamma(phi,  2)
  PSI1 <- trigamma(phi*(1-mu))
  PSI2 <- psigamma(phi*(1-mu), 2)
  
  nu.stu <- nu.s.tu <- array(NA,c((p+q),(p+q),c(p+q)))
  for(t in 1:p)
  {
    Xt<-X[,t]
    bb<- crossprod(X,phi^2*D1*Dp1*K2*Xt*X)
    bg <- crossprod(X,phi*D1^2*D2*K2*Xt*Z)
    gg <- crossprod(Z,phi*D1*Dp2*(mu*K2-PSI1)*Xt*Z)
    nu.s.tu[t,,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
    bb<- crossprod(X,phi^3*D1^3*K3*Xt*X)
    bg <- crossprod(X,phi^2*D1^2*D2*(mu*K3+PSI2)*Xt*Z)
    # ans[(p+1):(p+q),1:p,t] <- t(ans[1:p,(p+1):(p+q),t])
    gg <- crossprod(Z,phi*D1*D2^2*(mu^2*K3+2*mu*PSI2-PSI2)*Xt*Z)
    nu.stu[t,,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
  }  
  for(s in 1:q)
  {
    Zs<-Z[,s]
    bb <- crossprod(X,phi*D2*Dp1*(mu*K2-PSI1)*Zs*X)
    bg <- crossprod(X,D1*D2^2*(mu*K2-PSI1)*Zs*Z)
    gg <- crossprod(Z,D2*Dp2*(mu^2*K2+PSI1-2*mu*PSI1-OMEGA1)*Zs*Z)
    nu.s.tu[(p+s),,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
    bb <- crossprod(X,phi^2*D1^2*D2*(mu*K3+PSI2)*Zs*X)
    bg <- crossprod(X,phi*D1*D2^2*(mu^2*K3+2*mu*PSI2-PSI2)*Zs*Z)
    gg<- crossprod(Z,D2^3*(mu^3*K3+(3*mu^2-3*mu+1)*PSI2-OMEGA2)*Zs*Z)
    nu.stu[(p+s),,]<- rbind(cbind(bb,bg),cbind(t(bg),gg))
  }
  ##----score-----
  lambda <- digamma(phi*mu)-digamma(phi)
  Tbar <- log(y)-lambda
  epsilon <- digamma(phi*(1-mu))-digamma(phi)
  Ubar <- log(1-y)-epsilon
  score<-c(crossprod(X,phi*D1*(Tbar-Ubar)),crossprod(Z,D2*(mu*(Tbar-Ubar)+Ubar)))
  ##---- Fisher info and inverse----
  bb<- crossprod(X,D1^2*phi^2*K2*X)  
  bg <- crossprod(X,phi*D1*D2*(mu*K2-PSI1)*Z)
  gg <- crossprod(Z,D2^2*(mu^2*K2+(1-2*mu)*PSI1-OMEGA1)*Z)
  info<-rbind(cbind(bb,bg),cbind(t(bg),gg))
  InfoInv <- try(chol2inv(chol(info)),TRUE)
  if(failedInv <- inherits(InfoInv, "try-error")) 
  {
    warning("failed to invert the information matrix: iteration stopped prematurely")
    break
  }
  modification <- .C('modificationc5',
                     as.integer(p+q),
                     as.double(InfoInv),
                     as.double( nu.stu),
                     as.double(nu.s.tu),
                     mod=as.double(rep(0,(p+q))))$mod
  ##----- modified score function------##
  (1/diag(InfoInv))*InfoInv%*%score + modification
}





mbrbetareg.fit=function(par,y,X,Z,eps=1e-06,maxit=500,link="logit",phi_link="log")
{
   step <- .Machine$integer.max
      nIter <- 0
      test <- TRUE
      while ( test & (nIter < maxit))
      {
            nIter <- nIter + 1  
            info <-info(par,X,Z,link,phi_link)
            score <- score(par,y,X,Z,link,phi_link)
            InfoInv <- try(chol2inv(chol(info)),TRUE)
            # print(InfoInv)
            if(failedInv <- inherits(InfoInv, "try-error")) 
            {
               warning("failed to invert the information matrix: iteration stopped prematurely")
               break
            }
            
            mod <- modification.c(par,X,Z,link,phi_link)
            modscore <- score + info%*%mod
            par <- par+mod + InfoInv%*%score
            
            test <- sqrt(crossprod(drop(modscore))) > eps
            
      }
converged <- nIter < maxit
list(par=drop(par),convergence=converged,nIter=nIter,InfoInv=InfoInv,se=sqrt(diag(InfoInv)))
}



