## This file contains functions that have been copied
## from package 'betareg' version 3.1-0.
## The authors of betareg are:
## Achim Zeileis <Achim.Zeileis@R-project.org>
## Francisco Cribari-Neto <cribari@ufpe.br>
## Ioannis Kosmidis <ioannis@stats.ucl.ac.uk>
## Alexandre B. Simas 
## Andrea V. Rocha 

## In particular, 'mbrbetareg' and 'mbrbetareg.fit' were written using as basis the code 
## of 'betareg' and 'betareg.fit', respectively.

##  'print.mbrbetareg' is a copy of 'print.betareg'
##  'summary.mbrbetareg' is a copy of 'summary.betareg'
##  'print.summary.mbrbetareg' is a copy of 'print.summary.betareg'

## Euloge Clovis Kenne Pagui <kenne@stat.unipd.it> [01/09/2017]

library(Formula)
library(betareg)

mbrbetareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    link.phi = NULL, type = c("ML", "medianBR"),
                    control = betareg.control(...),
                    model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
   
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
 
  ## sanity checks on response variable (Y)
  if(length(Y) < 1) stop("empty model")
  if(!(min(Y) > 0 & max(Y) < 1)) stop("invalid dependent variable, all observations must be in (0, 1)")

  ## convenience variables
  n <- length(Y)

  ## type of estimator
  type <- match.arg(type)

  ## check links
  if(is.character(link)) link <- match.arg(link)
  if(is.null(link.phi)) link.phi <- if(simple_formula) "identity" else "log"
  if(is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## offsets in mean part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## offsets in precision part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for mean)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(mean = offsetX, precision = offsetZ)

  ## call the actual workhorse: betareg.fit()
  fitmodel <- mbrbetareg.fit(X, Y, Z, weights, offset, link, link.phi, type, control)

  ## further model information
  fitmodel$call <- cl
  fitmodel$formula <- oformula
  fitmodel$terms <- list(mean = mtX, precision = mtZ, full = mt)
  fitmodel$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  fitmodel$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
  if(model) fitmodel$model <- mf
  if(y) fitmodel$y <- Y
  if(x) fitmodel$x <- list(mean = X, precision = Z)

  class(fitmodel) <- "mbrbetareg"
  return(fitmodel)
}


mbrbetareg.control <- function(phi = TRUE, method = "BFGS",maxit = 5000,
                             hessian = FALSE, trace = FALSE, start = NULL,
                             fsmaxit = 200, fstol = 1e-8, ...) {
  res <- list(phi = phi,  hessian = hessian, method = method, trace = trace, 
              start = start, fsmaxit = fsmaxit, fstol = fstol)
  res <- c(res, list(...))
  if(!is.null(res$fnscale)) warning("fnscale must not be modified")
  res$fnscale <- -1
  if(is.null(res$reltol)) res$reltol <- .Machine$double.eps^(1/1.2)
  if (!is.numeric(res$fstol) || res$fstol <= 0)
    stop("value of 'fstol' must be > 0")
  res
}

mbrbetareg.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
                           link = "logit", link.phi = "log", type = "ML",
                           control = mbrbetareg.control()) {
  dyn.load("mod1.so")
  #######################################
  ## response (y),  design matrix (x, z)#
  #######################################
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(mean = offset, precision = rep.int(0, n))
  if(is.null(z)) {
    m <- 1L
    z <- matrix(1, ncol = m, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
    phi_const <- TRUE
  } else {
    m <- NCOL(z)
    if(m < 1L) stop("dispersion regression needs to have at least one parameter")
    phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[, 1L]), rep.int(1, n)))
  }
  
  ###########################################
  ## set links for the mean and dispersion ##
  ###########################################
  if(is.character(link)) {
    linkstr <- link
    if(linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      ## Enrich the family object with the required derivatives
      linkobj <- enrichwith::enrich(linkobj , with = "d2mu.deta")
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) -log(-log(mu)),
        linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
        mu.eta = function(eta) {
          eta <- pmin(eta, 700)
          pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
        },
        dmu.deta = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
        valideta = function(eta) TRUE,
        name = "loglog"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link
    linkstr <- link$name
    if(type != "ML" && is.null(linkobj$dmu.deta)) warning("link needs to provide dmu.deta component for medianBR")
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$d2mu.deta
  if(is.character(link.phi)) {
    phi_linkstr <- link.phi
    phi_linkobj <- make.link(phi_linkstr)
    phi_linkobj <- enrichwith::enrich(phi_linkobj , with = "d2mu.deta")
  } else {
    phi_linkobj <- link.phi
    phi_linkstr <- link.phi$name
    if(type != "ML" && is.null(phi_linkobj$dmu.deta)) warning("link.phi needs to provide dmu.deta component for  medianBR")
  }
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.eta <- phi_linkobj$mu.eta
  phi_dmu.deta <- phi_linkobj$d2mu.deta
  
  ## y* transformation ##
  ystar <- qlogis(y)
  
  ############################
  ## control unput parameters#
  ############################
  ocontrol <- control
  phi_full <- control$phi
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  control$phi <- control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL
  
  ###########################################
  ## starting values as suggested in  #######
  ## Ferrari & Cribari-Neto (2004, Section2)#
  ###########################################
  if(is.null(start)) {
    auxreg <- lm.wfit(x, linkfun(y), weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    sigma2 <- sum(weights * res^2)/((sum(weights) - k) * (dlink)^2)
    phi_y <- weights * yhat * (1 - yhat)/(sum(weights) * sigma2) - 1/n
    phi <- rep.int(0, ncol(z))
    phi[1L] <- suppressWarnings(phi_linkfun(sum(phi_y)))
    if(!isTRUE(phi_linkinv(phi[1L]) > 0)) {
      warning("no valid starting value for precision parameter found, using 1 instead")
      phi[1L] <- 1
    }
    start <- list(mean = beta, precision = phi)
  }
  if(is.list(start)) start <- do.call("c", start)
  
  ###################################################
  ## various fitted quantities useful to calculate ##
  ## the necessary expressions                     ##
  ###################################################
  fitfun <- function(par, deriv = 0L) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    # mean quantities
    eta <- as.vector(x %*% beta + offset[[1L]])
    mu <- linkinv(eta)

    # dispersion quantities
    phi_eta <- as.vector(z %*% gamma + offset[[2L]])
    phi <- phi_linkinv(phi_eta)

    mustar <- psi1 <- psi2 <-  omega1 <-  omega2 <- kappa2 <- kappa3 <- d1mu <- d2mu <- d1phi <- d2phi <- NULL
    if(deriv >= 1L) {
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      d1mu <- mu.eta(eta)
      d1phi <- phi_mu.eta(phi_eta)
    } 
    if(deriv >= 2L) {
      omega1 <- psigamma(phi, deriv = 1)
      omega2 <- psigamma(phi,  deriv = 2)
      psi1 <- psigamma(phi*(1-mu), deriv = 1)
      psi2 <- psigamma(phi*(1-mu), deriv = 2)
      kappa2 <- psigamma(phi*mu, deriv = 1) + psi1
      kappa3 <- psigamma(phi*mu, deriv = 2) - psi2
      d2mu <- dmu.deta(eta)
      d2phi <- phi_dmu.deta(phi_eta)
    } 
    list(
      beta = beta,
      gamma = gamma,
      eta = eta,
      phi_eta = phi_eta,
      mu = mu,
      phi = phi,
      d1mu = d1mu,
      d1phi = d1phi,
      d2mu = d2mu,
      d2phi =d2phi,
      mustar = mustar,
      omega1 = omega1,
      omega2 = omega2,
      psi1 = psi1,
      psi2 = psi2,
      kappa2 = kappa2,
      kappa3 = kappa3
    )
  }
  
  ###################################
  ## log-likelihood function ########
  ###################################
  loglikfun <- function(par, fit = NULL) {
    ## extract fitted quantities involved in the likelihood
    if(is.null(fit)) {
      fit <- fitfun(par)
    }
    with(fit, {
        alpha <- mu * phi
        beta <- (1 - mu) * phi
        ## compute log-likelihood
        if(any(!is.finite(phi)) | any(alpha > 1e300) | any(beta > 1e300)) NaN else { ## catch extreme cases without warning
          ll <- suppressWarnings(dbeta(y, alpha, beta, log = TRUE))
          if(any(!is.finite(ll))) NaN else sum(weights * ll) ## again: catch extreme cases without warning
          }  
    })
  }
  
  #############################################
  ## score function (gradient) by default    ## 
  ## or gradient contributions (sum = FALSE) ##
  #############################################
  gradfun <- function(par, sum = TRUE, fit = NULL) {
    ## extract fitted quantities involved in the score function
    if(is.null(fit)) {
      fit <- fitfun(par, deriv = 1L)
    }
    with(fit, {
        mu <- mu
        phi <- phi
        eta <- eta
        phi_eta <- phi_eta
        d1mu <- d1mu
        d1phi <- d1phi
        mustar <- mustar
        ## compute gradient contributions
        res <- cbind(
        phi * (ystar - mustar) * d1mu * weights * x,
        (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) * d1phi * weights * z
        )
        if(sum) colSums(res) else res
    })
  }
  
  #########################################
  ## Fisher information (default)  ########
  ## or covariance matrix (inverse=TRUE) ##
  #########################################
  information <- function(par, inverse = FALSE, fit = NULL) {
    ## extract fitted means/precisions
    if(is.null(fit)) {
      fit <- fitfun(par, deriv = 2L)
    }
    with(fit, {
        mu <- mu
        phi <- phi
        eta <- eta
        phi_eta <- phi_eta
        d1mu <- d1mu
        d1phi <- d1phi
        psi1 <- psi1
        omega1 <- omega1
        kappa2 <- kappa2
        ## compute elements of W (notation in Ospina et al., 2006, Section 2)
        wbb <- phi^2 * kappa2 * d1mu^2 
        wpp <- (mu^2*kappa2 + (1-2*mu)*psi1 - omega1) * d1phi^2    
        wbp <- phi * (mu * kappa2 - psi1) * d1mu * d1phi 
        ## compute elements of K
        kbb <- if(k > 0L) crossprod(sqrt(weights) * sqrt(wbb) * x) else crossprod(x)
        kpp <- if(m > 0L) crossprod(sqrt(weights) * sqrt(wpp) * z) else crossprod(z)
        kbp <- if(k > 0L & m > 0L) crossprod(weights * wbp * x, z) else crossprod(x, z)
        ## Fisher information
        K <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
        if(!inverse) K else chol2inv(chol(K))
    })  
  }
  
  ########################################################
  ## third order cumulants involved in the adjustment  ###
  ########################################################
  
  nuQuantities <- function(par, fit = NULL) {
    nu.stu <- nu.s.tu <- array(NA_real_,dim = c((k+m),(k+m),c(k+m))) 
    ## extract fitted means/precisions
    if(is.null(fit)) {
      fit <- fitfun(par, deriv = 2L)
    }
    with(fit, {
      mu <- mu
      phi <- phi
      eta <- eta
      phi_eta <- phi_eta
      d1mu <- d1mu
      d2mu <- d2mu
      d1phi <- d1phi
      d2phi <- d2phi
      psi1 <- psi1
      psi2 <- psi2
      omega1 <- omega1
      omega2 <- omega2
      kappa2 <- kappa2
      kappa3 <- kappa3
      for (t in 1:k ) {
        xt <- x[,t]
        if (k > 0L) {
          bb1 <- crossprod(x,weights*phi^2*d1mu*d2mu*kappa2*xt*x)
          bb2 <- crossprod(x,weights*phi^3*d1mu^3*kappa3*xt*x)
        } else bb1 <- bb2 <- crossprod(x)
        
        if ((k > 0L) & (m > 0L)) {
          bg1 <- crossprod(x,weights*phi*d1mu^2*d1phi*kappa2*xt*z)
          bg2 <- crossprod(x,weights*phi^2*d1mu^2*d1phi*(mu*kappa3+psi2)*xt*z)
        } else  bg1 <- bg2 <- crossprod(x ,z)
        
        if (m > 0L) {
          gg1 <- crossprod(z,weights*phi*d1mu*d2phi*(mu*kappa2-psi1)*xt*z)
          gg2 <- crossprod(z,weights*phi*d1mu*d1phi^2*(mu^2*kappa3+2*mu*psi2-psi2)*xt*z)
        }  else gg1 <- gg2 <- crossprod(z)
        
        nu.s.tu[t,,]<- rbind(cbind(bb1,bg1),cbind(t(bg1),gg1))
        nu.stu[t,,]<- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
      }
      
      for (s in 1:m) {
        zs<-z[,s]
        if (k > 0L) {
          bb1 <- crossprod(x,weights*phi*d1phi*d2mu*(mu*kappa2-psi1)*zs*x)
          bb2 <- crossprod(x,weights*phi^2*d1mu^2*d1phi*(mu*kappa3+psi2)*zs*x)
        } else bb1 <- bb2 <- crossprod(x)
        
        if ((k > 0L) & (m > 0L)) {
          bg1 <- crossprod(x,weights*d1mu*d1phi^2*(mu*kappa2-psi1)*zs*z)
          bg2 <- crossprod(x,weights*phi*d1mu*d1phi^2*(mu^2*kappa3+2*mu*psi2-psi2)*zs*z)
        } else  bg1 <- bg2 <- crossprod(x ,z)
        
        if (m > 0L) {
          gg1 <- crossprod(z,weights*d1phi*d2phi*(mu^2*kappa2+psi1-2*mu*psi1-omega1)*zs*z)
          gg2 <- crossprod(z,weights*d1phi^3*(mu^3*kappa3+(3*mu^2-3*mu+1)*psi2-omega2)*zs*z)
        }  else gg1 <- gg2 <- crossprod(z)
        
        nu.s.tu[(k+s),,]<- rbind(cbind(bb1,bg1),cbind(t(bg1),gg1))
        nu.stu[(k+s),,]<- rbind(cbind(bb2,bg2),cbind(t(bg2),gg2))
      }
      return(list(nu.stu=nu.stu, nu.s.tu=nu.s.tu))  
    })   
  }
  
  ###################################
  ## factor of the adjustment     ###
  ## adjustment = info%*%mod1     ###
  ###################################
  
  mod1 <- function(par, inverse_info, fit = NULL) {
    if (is.null(fit)) {
      fit <- fitfun(par, deriv = 2L)
    }
    nu <- nuQuantities(par, fit)
    nu.stu <- nu$nu.stu
    nu.s.tu <- nu$nu.s.tu
    neededQuant <- c(inverse_info, nu.stu, nu.s.tu)
    adj <- .C("mod1",  as.integer((k+m)), as.double(neededQuant), as.double(rep(0.0,(k+m))))[[3]]
    return(adj)
  }
  
  ## get the starting point
  # par <- start
  
  ########################################
  ## optimize the log likelihood       ###
  ########################################
  opt <- optim(par = start, fn = loglikfun, gr = gradfun,
               method = method, hessian = hessian, control = control)
  par <- opt$par
  
  
  # ## get first the MLE from betareg R  #
  # temp.fit <- betareg::betareg.fit(x = x, y = y, z = z, weights = weights, offset = offset,
  #                                  link = link, link.phi = link.phi, type = "ML",
  #                                  control = ocontrol)
  # 
  # ## get the BC estimates              ##
  # if (type == "BC")
  # temp.fit <- betareg::betareg.fit(x = x, y = y, z = z, weights = weights, offset = offset,
  #                                  link = link, link.phi = link.phi, type = "BC",
  #                                  control = ocontrol)
  # 
  # ## get the BR estimates             ##
  # if (type == "BR")
  #   temp.fit <- betareg::betareg.fit(x = x, y = y, z = z, weights = weights, offset = offset,
  #                                    link = link, link.phi = link.phi, type = "BR",
  #                                    control = ocontrol)
  
  ## get the median BR estimates     ##
  if(type == "medianBR" & fsmaxit <=0) 
    warning("BR cannot be performed with fsmaxit <= 0")
  if(type =="medianBR" & opt$convergence > 0) par <- start
  step <- .Machine$integer.max
  iter <- 0
  if(type == "medianBR" & fsmaxit > 0) {
    if( fsmaxit > 0) {
      for (iter in seq.int(fsmaxit)) {
        stepPrev <- step
        stepFactor <- 0
        testhalf <- TRUE
        while (testhalf & stepFactor < 11) {
          fit <- fitfun(par, deriv = 2L)
          score <- gradfun(par, sum = TRUE,  fit = fit)
          inverse_info <- try(information(par, inverse = TRUE, fit = fit))
          if (failed <- inherits(inverse_info, "try-error")) {
            warning("failed to invert the information matrix: iteration stopped prematurely")
            break
          }
          mod <- mod1(par, inverse_info, fit)
          par <- par + 2^(-stepFactor) * (step <- mod + inverse_info%*%score)
          stepFactor <- stepFactor + 1
          testhalf <- drop(crossprod(stepPrev) < crossprod(step))  
        }
        if (failed | (all(abs(step) < fstol))) {
          break
        }
      }
    }
  }
  
  if((fsmaxit == 0 & opt$convergence > 0) | iter >= fsmaxit) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }
  
  ## extract fitted values/parameters
  fit <- fitfun(par, deriv = 2L)
  beta <- fit$beta
  gamma <- fit$gamma
  eta <- fit$eta
  mu <- fit$mu
  phi <- fit$phi
  ## log-likelihood
  ll <- loglikfun(par, fit = fit)
  ## gradient
  grad <- gradfun(par,  sum = TRUE, fit = fit)
  ## Fisher information
  info <- information(fit = fit, inverse = FALSE)
  ## covariance matrix at optimized parameters
  vcov <- if (hessian & (type == "ML")) solve(-as.matrix(temp.fit$hessian)) else information(fit = fit, inverse = TRUE)
  ## R-squared
  pseudor2 <- if(var(eta) * var(ystar) <= 0) NA else cor(eta, linkfun(y))^2
  ## names
  names(beta) <- colnames(x)
  names(gamma) <- if(phi_const & phi_linkstr == "identity") "(phi)" else colnames(z)
  rownames(vcov) <- colnames(vcov)  <- c(colnames(x),
                                         if(phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"))
  
  res <- list(
    coefficients = list(mean = beta, precision = gamma),
    converged = converged,
    residuals = y - mu,
    fitted.values = structure(mu, .Names = names(y)),
    type = type,
    optim=opt,
    method = method,
    control = ocontrol,
    scoring = iter,
    start = start,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(mean = if(identical(offset[[1L]], rep.int(0, n))) NULL else offset[[1L]],
                            precision = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    df.null = nobs - 2,
    df.residual = nobs - k - m,
    phi = phi_full,
    loglik = ll,
    vcov = vcov,
    grad = if (type=="ML") grad else grad+mod%*%info,
    pseudo.r.squared = pseudor2,
    link = list(mean = linkobj, precision = phi_linkobj)
  )
  res
}


print.mbrbetareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients$mean)) {
      cat(paste("Coefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in mean model)\n\n")
    if(x$phi) {
      if(length(x$coefficients$precision)) {
        cat(paste("Phi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
        print.default(format(x$coefficients$precision, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      } else cat("No coefficients (in precision model)\n\n")
    }
  }

  invisible(x)
}

summary.mbrbetareg <- function(object, phi = NULL, type = "response", ...)
{
  ## treat phi as full model parameter?
  if(!is.null(phi)) object$phi <- phi

  ## residuals
  type <- match.arg(type, c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type

  ## extend coefficient table
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$precision)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], precision = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$precision) <- names(object$coefficients$precision)
  object$coefficients <- cf

  ## number of iterations
  mytail <- function(x) x[length(x)]
  object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$count))), "scoring" = as.vector(object$scoring))

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.mbrbetareg"
  object
}

print.summary.mbrbetareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    types <- c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
    Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals",
      "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$mean)) {
      cat(paste("\nCoefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in mean model)\n")

    if(x$phi) {
      if(NROW(x$coefficients$precision)) {
        cat(paste("\nPhi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
        printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
      } else cat("\nNo coefficients (in precision model)\n")
    }

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat("\nType of estimator:", x$type, switch(x$type,
      "ML" = "(maximum likelihood)",
      "medianBR" = "(median bias-reduced)"))
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    if(x$iterations[2L] > 0) {
      cat(paste("\nNumber of iterations:", x$iterations[1L],
        sprintf("(%s) +", x$method), x$iterations[2L], "(Fisher scoring) \n"))
    } else {
      cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
    }
  }

  invisible(x)
}







