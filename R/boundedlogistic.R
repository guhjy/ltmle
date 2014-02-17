
# boundedlogistic <- function(...) {
#   UseMethod("boundedlogistic", ...)
# }

boundedlogistic <- function(formula, data, subset, weights, na.action, contrasts = NULL, offset, a=0, b=1, method="L-BFGS-B") {
  #From lm()
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
      "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- model.response(mf, "numeric")

  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
           length(offset), NROW(y)), domain = NA)
  }  
  x <- model.matrix(mt, mf, contrasts)

  z <- if (is.null(w)) {
    boundedlogistic.fit(x, y, a=a, b=b, offset = offset, method=method)
  }
  else {
    boundedlogistic.fit(x, y, a=a, b=b, offset = offset, method=method, w=w)
  }

  class(z) <- "boundedlogistic"
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt

  z
}


predict.boundedlogistic <- function(obj, newdata=NULL, type="response", na.action=na.pass){

  if (is.null(newdata)) {
    xB <- obj$xB
  }  else {
    tt <- obj$terms
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = obj$xlevels)
    x <- model.matrix(Terms, m, obj$contrasts)
    offset <- rep(0, nrow(x))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
        "variables")[[i + 1]], newdata)
    if (!is.null(obj$call$offset)) 
      offset <- offset + eval(obj$call$offset, newdata)

    xB <- x %*% cbind(obj$coef) + offset
  }

  if (type=="link") {
    res <- xB
  }
  else if (type=="response") {
    res <- plogis(xB) * (obj$b-obj$a) + obj$a
  }
  else {
    stop("in bounded logistic, type ", type, " not supported")
  }
  res
}


print.boundedlogistic <- function(obj) print(coef(obj))



boundedlogistic.fit <- function(x, y, a, b, offset, method="L-BFGS-B", w=NULL) {
  ## Uses the transformed Y loss function
  ystar <- (y-a)/(b-a)
  ff <- function(beta) {
    xB <- if (!is.null(offset)) {
      x %*% cbind(beta) + offset
      } else {
      x %*% cbind(beta)
    }

    unweighted_loglik <- ystar * plogis(xB, log.p=TRUE) + (1-ystar) * plogis(xB, log.p=TRUE, lower.tail=FALSE)

    if (is.null(w)) {
      -2 * sum(unweighted_loglik)
    } else {
      -2 * sum(w * unweighted_loglik)
    }
  }

  gg <- function(beta) {
    xB <- if (!is.null(offset)) {
      x %*% cbind(beta) + offset
      } else {
      x %*% cbind(beta)
    }
    p <- plogis(xB)

    if (is.null(w)) 2 * crossprod(x, cbind(p - ystar))
    else 2 * crossprod(x, w*cbind(p - ystar))
  }

  hh <- function(beta) {
    xB <- if (!is.null(offset)) {
      x %*% cbind(beta) + offset
      } else {
      x %*% cbind(beta)
    }
    p <- plogis(xB)
    S <- c(p * (1-p))
    if (is.null(w)) 2* crossprod(x, S *  x)
    else 2* crossprod(x, (S*w) *  x)
  }

  nlminb.res <- nlminb(rep(0, ncol(x)), ff, gg, hh)
  # if (nlminb.res$convergence != 0)
  #   stop(nlminb.res$message)

  beta <- nlminb.res$par
  names(beta) <- colnames(x)

  xB <- if (!is.null(offset)) {
    x %*% cbind(beta) + offset
    } else {
    x %*% cbind(beta) 
    }



  res <- list(coefficients=beta, a=a, b=b, xB=xB, nlminb.res=nlminb.res)

  class(res) <- "boundedlogistic"

  res
}