##%######################################################%##
#                                                          #
####                   GRaF functions                   ####
#                                                          #
##%######################################################%##


#Functions adapted from the package GRaF (https://github.com/goldingn/GRaF)
#to fit models usign a Gausian-Bayesian approach (Golding & Purse, 2016. doi:10.1111/2041-210X.12523) 


#Fitting GRaF Models----

graf <-
  
  function (y, x, error = NULL, weights = NULL, prior = NULL, l = NULL, opt.l = FALSE,
            
            theta.prior.pars = c(log(10), 1), hessian = FALSE, opt.control = list(),
            
            verbose = FALSE, method = c('Laplace', 'EP')) {
    
    
    
    method <- match.arg(method)
    
    
    
    # optionally optimise graf (by recursively calling this function)
    
    if (opt.l) {
      
      
      
      # get all visible object as a list
      
      args <- capture.all()
      
      
      
      # get the expected objects
      
      expected_args <- names(formals(graf))
      
      
      
      # remove any unexpected arguments
      
      args <- args[names(args) %in% expected_args]
      
      
      
      # pass this to optimiser
      
      fit <- optimise.graf(args)
      
      
      
      # skip out of this function and return
      
      return (fit)
      
      
      
    }
    
    
    
    
    
    if (!is.data.frame(x)) stop ("x must be a dataframe")
    
    
    
    # convert any ints to numerics
    
    for(i in 1:ncol(x)) if (is.integer(x[, i])) x[, i] <- as.numeric(x[, i])
    
    
    
    obsx <- x
    
    k <- ncol(x)
    
    n <- length(y)
    
    
    
    if (is.null(weights)) {
      
      # if weights aren't provided
      
      weights <- rep(1, n)
      
    } else {
      
      # if they are, run some checks
      
      # throw an error if weights are specified with EP
      
      if (method == 'EP') {
        
        stop ('weights are not implemented for the EP algorithm (yet)')
        
      }
      
      # or if any are negative
      
      if (any(weights < 0)) {
        
        stop ('weights must be positive or zero')
        
      }
      
    }
    
    
    
    # find factors and convert them to numerics
    
    notfacs <- 1:k
    
    facs <- which(unlist(lapply(x, is.factor)))
    
    if (length(facs) > 0) notfacs <- notfacs[-facs]
    
    for (fac in facs) {
      
      x[, fac] <- as.numeric(x[, fac])
      
    }
    
    x <- as.matrix(x)
    
    
    
    # scale the matrix, retaining scaling
    
    scaling <- apply(as.matrix(x[, notfacs]), 2, function(x) c(mean(x), sd(x)))
    
    for (i in 1:length(notfacs)) {
      
      x[, notfacs[i]] <- (x[, notfacs[i]] - scaling[1, i]) / scaling[2, i]
      
    }
    
    
    
    # set up the default prior, if not specified
    
    exp.prev <- sum(weights[y == 1]) / sum(weights)
    
    if (is.null(prior))  mnfun <- function(x) rep(exp.prev, nrow(x))
    
    else mnfun <- prior
    
    
    
    # give an approximation to l, if not specified (or optimised)
    
    if (is.null(l)) {
      
      l <- rep(0.01, k)
      
      l[notfacs] <- apply(x[y == 1, notfacs, drop = FALSE], 2, sd) * 8
      
      l[l == 0] <- 1
      
    }
    
    
    
    # calculate mean (on unscaled data and probability scale)
    
    mn <- mnfun(obsx)
    
    
    
    # fit model
    
    if (method == 'Laplace') {
      
      # by Laplace approximation
      
      fit <- graf.fit.laplace(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = verbose)
      
    } else {
      
      # or using the expectation-propagation algorithm
      
      fit <- graf.fit.ep(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = FALSE)
      
    }
    
    
    
    fit$mnfun = mnfun
    
    fit$obsx <- obsx
    
    fit$facs <- facs
    
    fit$hessian <- hessian
    
    fit$scaling <- scaling
    
    fit$peak = obsx[which(fit$MAP == max(fit$MAP))[1], ]
    
    class(fit) <- "graf"
    
    fit
    
  }



capture.all <- function() {
  
  # capture all visible objects in the parent environment and pass to a list
  
  env <- parent.frame()
  
  object_names <- objects(env)
  
  objects <- lapply(object_names,
                    
                    get,
                    
                    envir = env)
  
  names(objects) <- object_names
  
  return (objects)
  
}

#Fit Laplace----
graf.fit.laplace <-
  
  function (y, x, mn, l, wt, e = NULL, tol  = 10 ^ -6, itmax = 50,
            
            verbose = FALSE) {
    
    
    
    if (is.vector(x)) x <- as.matrix(x)
    
    mn <- qnorm(mn)
    
    n <- length(y)
    
    
    
    # create the covariance matrix
    
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    
    
    
    # an identity matrix for the calculations
    
    eye <- diag(n) 
    
    
    
    # initialise
    
    a <- rep(0, n)
    
    f <- mn
    
    obj.old <- Inf
    
    obj <- -sum(wt * d0(f, y))
    
    it <- 0
    
    
    
    # start newton iterations
    
    while (obj.old - obj > tol & it < itmax) {
      
      it <- it + 1
      
      obj.old <- obj
      
      W <- -(wt * d2(f, y))
      
      rW <- sqrt(W)
      
      cf <- f - mn
      
      mat1 <- rW %*% t(rW) * K + eye
      
      L <- tryCatch(chol(mat1),
                    
                    error = function(x) return(NULL))
      
      b <- W * cf + wt * d1(f, y)
      
      mat2 <- rW * (K %*% b)
      
      adiff <- b - rW * backsolve(L, forwardsolve(t(L), mat2)) - a 
      
      dim(adiff) <- NULL
      
      
      
      # find optimum step size using Brent's method
      
      res <- optimise(psiline, c(0, 2), adiff, a, as.matrix(K), y, d0, mn, wt)
      
      a <- a + res$minimum * adiff
      
      f <- K %*% a + mn
      
      obj <- psi(a, f, mn, y, d0, wt)
      
      
      
    }
    
    
    
    # recompute key components
    
    cf <- f - mn
    
    W <- -(wt * d2(f, y))
    
    rW <- sqrt(W)
    
    mat1 <- rW %*% t(rW) * K + eye
    
    L <- tryCatch(chol(mat1),
                  
                  error = function(x) return(NULL))
    
    
    
    # return marginal negative log-likelihood
    
    mnll <- (a %*% cf)[1, 1] / 2 + sum(log(diag(L)) - (wt * d0(f, y)))
    
    
    
    # get partial gradients of the objective wrt l
    
    
    
    # gradient components
    
    W12 <- matrix(rep(rW, n), n)
    
    R <- W12 * backsolve(L, forwardsolve(t(L), diag(rW)))
    
    C <- forwardsolve(t(L), (W12 * K))
    
    
    
    # partial gradients of the kernel
    
    dK <- cov.SE.d1(x, e, l)
    
    
    
    # rate of change of likelihood w.r.t. the mode
    
    s2 <- (diag(K) - colSums(C ^ 2)) / 2 * d3(f, y)
    
    
    
    # vector to store gradients
    
    l_grads <- rep(NA, length(l))
    
    
    
    for (i in 1:length(l)) {
      
      
      
      grad <- sum(R * dK[[i]]) / 2
      
      grad <- grad - (t(a) %*% dK[[i]] %*% a) / 2
      
      b <- dK[[i]]%*% d1(f, y)
      
      grad <- grad - t(s2) %*% (b - K %*% (R %*% b))
      
      l_grads[i] <- -grad
      
      
      
    }
    
    
    
    if(verbose ) cat(paste("  ", it, "Laplace iterations\n"))
    
    if(it == itmax) print("timed out, don't trust the inference!")
    
    return(list(y = y, x = x, MAP = f, ls = l, a = a, W = W, L = L, K = K,
                
                e = e, obsx = x, obsy = y, mnll = mnll, wt = wt,
                
                l_grads = l_grads))
    
  }

#Fit EP----
graf.fit.ep <-
  
  function(y, x, mn, l, wt, e = NULL, tol = 1e-6, itmax = 50, itmin = 2,
           
           verbose = FALSE) {
    
    # fit a GRaF model using expectation-propagation
    
    # as implemented in the gpml matlab library
    
    # If parallel  = TRUE, the EP uses the parallel update
    
    # as implemented in GPStuff
    
    
    
    # whether to use parallel EP (rather than sequential)
    
    parallel = FALSE
    
    
    
    if (is.vector(x)) { 
      
      x <- as.matrix(x)
      
    }
    
    
    
    # mn to probability scale
    
    mn <- qnorm(mn)
    
    n <- length(y)
    
    # covariance matrix
    
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    
    # identity matrix
    
    eye <- diag(n)
    
    
    
    # convert observations to +1, -1, save 0, 1 version
    
    oldy <- y
    
    y <- y * 2 - 1
    
    
    
    # initialise
    
    
    
    # \tilde{\mathbf{\nu}} = \mathbf{0}
    
    tnu <- rep(0, n)
    
    # \tilde{\mathbf{\tau}} = \mathbf{0}
    
    ttau <- rep(0, n)
    
    # \mathbf{\mu} = \mathbf{0}
    
    mu <- rep(0, n)
    
    # \Sigma = \mathbf{K}  (only used in sequential EP)
    
    Sigma <- K
    
    
    
    # calculate marginal negative log likelihood at ttau = tnu = mu = 0s
    
    z <- mu / sqrt(1 + diag(K))
    
    mnll <- -sum(pnorm(y * z), log.p = TRUE)
    
    
    
    # ~~~~~~~~~~~~~~~~~~~~~
    
    # set up for EP algorithm
    
    # set up damping factor (for parallel EP) as in GPStuff
    
    df <- 0.8
    
    
    
    # set old mnll to Inf & start iteration counter
    
    mnll_old <- Inf
    
    logM0_old <- logM0 <- rep(0, n)
    
    it <- 1
    
    converged <- FALSE
    
    
    
    while (!converged) {
      
      
      
      if (parallel) {
        
        
        
        # calculate key parameters
        
        dSigma <- diag(Sigma)
        
        tau <- 1 / dSigma - ttau
        
        nu <- 1 / dSigma * mn - tnu
        
        mu <- nu / tau
        
        sigma2 <- 1 / tau
        
        
        
        # get marginal moments of posterior
        
        lis <- ep.moments(y, sigma2, mu)
        
        
        
        # recalculate parameters from these moments
        
        # \delta\tilde{\tau} (change in \tilde{\tau})
        
        dttau <- 1 / lis$sigma2hat - tau - ttau
        
        
        
        # \tilde{\tau}
        
        ttau <- ttau + df * dttau
        
        
        
        # \delta\tilde{\nu} (change in \tilde{\nu})
        
        dtnu <- 1 / lis$sigma2hat * lis$muhat - nu - tnu
        
        
        
        # \tilde{\nu}
        
        tnu <- tnu + df * dtnu
        
        
        
      } else {
        
        # otherwise use sequential EP
        
        
        
        lis <- list(Sigma = Sigma,
                    
                    ttau = ttau,
                    
                    tnu = tnu,
                    
                    mu = mu)
        
        
        
        # cycle through in random order
        
        for (i in sample(1:n)) {
          
          lis <- update.ep(i, y, mn, lis)    
          
        } # end permuted for loop
        
        
        
        Sigma <- lis$Sigma
        
        ttau <- lis$ttau
        
        tnu <- lis$tnu
        
        mu <- lis$mu
        
        
        
        sigma2 <- diag(Sigma)
        
        
        
      }  # end parallel / sequential
      
      
      
      # recompute the approximate posterior parameters \Sigma and \mathbf{\mu}
      
      # using eq. 3.53 and eq. 3.68
      
      
      
      # sW = \tilde{S}^{\frac{1}{2}} = \sqrt{\mathbf{\tilde{\tau}}}
      
      sW <- sqrt(ttau)
      
      
      
      # L = cholesky(I_n  + \tilde{S}^{\frac{1}{2}} * K * \tilde{S}^{\frac{1}{2}})
      
      L <- chol(sW %*% t(sW) * K + eye)
      
      
      
      # V = L^T \\ \tilde{S}^{\frac{1}{2}} * K
      
      sWmat <- matrix(rep(sW, n), n) # byrow = TRUE?
      
      V <- backsolve(L, sWmat * K, transpose = TRUE)
      
      
      
      # \Sigma = \mathbf{K} - \mathbf{V}^T \mathbf{V}
      
      Sigma <- K - t(V) %*% V
      
      
      
      # \mathbf{\mu} = \Sigma \tilde{\mathbf{\nu}}
      
      mu <- Sigma %*% tnu  # + mn
      
      
      
      # calculate new mnll and assess convergence
      
      # compute logZ_{EP} using eq. 3.65, 3.73 and 3.74 and the existing L
      
      # \mathbf{\sigma}^2 = diag(\Sigma)
      
      sigma2 <- diag(Sigma)
      
      tau_n <- 1 / sigma2 - ttau
      
      nu_n <- mu / sigma2 - tnu + mn * tau_n
      
      
      
      z <- nu_n / tau_n / (sqrt(1 + 1 / tau_n))
      
      lZ <- pnorm(y * z, log.p = TRUE)
      
      
      
      # split the final equation up into 5 bits...
      
      mnll.a <- sum(log(diag(L))) - sum(lZ) - t(tnu) %*% Sigma %*% tnu / 2
      
      mnll.b <- t(nu_n - mn * tau_n)
      
      mnll.c <- ((ttau / tau_n * (nu_n - mn * tau_n) - 2 * tnu) / (ttau + tau_n)) / 2
      
      mnll.d <- sum(tnu ^ 2 / (tau_n + ttau)) / 2
      
      mnll.e <- sum(log(1 + ttau / tau_n)) / 2
      
      
      
      mnll <- as.numeric(mnll.a - mnll.b %*% mnll.c + mnll.d - mnll.e)
      
      
      
      # improvement in negative log marginal likelihood
      
      dmnll <- abs(mnll - mnll_old)
      
      
      
      # improvement in log of the 0th moment
      
      dlogM0 <- max(abs(logM0 - logM0_old))
      
      
      
      # both under tolerance?
      
      sub_tol <- dmnll < tol & dlogM0 < tol
      
      
      
      if ((sub_tol & it >= itmin) | it >= itmax) {
        
        # stop if there was little improvement and there have been at least
        
        # itmin iterations or if there have been itmax or more
        
        # iterations
        
        converged <- TRUE
        
      } else {
        
        it <- it + 1
        
        mnll_old <- mnll
        
      }
      
      
      
    } # end while loop
    
    
    
    # throw an error if the iterations maxed out before convergence
    
    if (it >= itmax) {
      
      stop(paste0('maximum number of iterations (', itmax, ') reached
                  
                  without convergence'))
      
    }
    
    
    
    # calculate posterior parameters
    
    sW <- sqrt(ttau)
    
    alpha = tnu - sW * backsolve(L, backsolve(L, sW * (K %*% tnu), transpose = TRUE))
    
    f = crossprod(K, alpha) + mn
    
    
    
    # return relevant parameters
    
    return (list(y = oldy,
                 
                 x = x,
                 
                 MAP = f,
                 
                 ls = l,
                 
                 a = alpha,
                 
                 W = ttau,
                 
                 L = L,
                 
                 K = K,
                 
                 e = e,
                 
                 obsx = x,
                 
                 obsy = oldy,
                 
                 mnll = mnll,
                 
                 wt = wt,
                 
                 l_grads = NULL))
    
  }

#Cov.SE----
cov.SE <- function(x1, x2 = NULL, e1 = NULL, e2 = NULL, l) {
  
  n1 <- nrow(x1)
  
  n2 <- ifelse(is.null(x2), n1, nrow(x2))
  
  n3 <- ncol(x1)
  
  
  
  # distance matrices
  
  if(is.null(x2)) {
    
    e2 <- e1
    
    # if no second matrix do with distance matrices for speed up
    
    dists <- lapply(1:n3, function(i, x) dist(x[, i]) ^ 2, x1)
    
  } else {
    
    dists <- list()
    
    for (i in 1:n3) {
      
      dists[[i]] <-   x1[, i] ^ 2 %*% t(rep(1, n2)) +
        
        rep(1, n1) %*% t(x2[, i] ^ 2) - 2 * x1[, i] %*% t(x2[, i])
      
    }
    
  }
  
  
  
  # with error matrices
  
  if (!is.null(e1)) {
    
    E1 <- list()
    
    ones <- t(rep(1, n2))
    
    for (i in 1:n3) {
      
      E1[[i]] <- e1[, i] %*% ones
      
    }
    
    
    
    if (!is.null(e2)) {
      
      E2 <- list()
      
      ones <- t(rep(1, n1))
      
      for (i in 1:n3) {
        
        E2[[i]] <- t(e2[, i] %*% ones)
        
      }
      
    } else {
      
      E2 <- as.list(rep(0, n3))
      
    }
    
    
    
    # run through each covariate
    
    
    
    sumdiffs <- 0
    
    denom <- 1
    
    lower <- lower.tri(E1[[1]])
    
    for (i in 1:n3) {
      
      err <- E1[[i]] + E2[[i]]
      
      if (is.null(x2)) {
        
        err <- err[lower] # save only lower portion for speed up
        
      }
      
      sumdiffs <- sumdiffs + dists[[i]] / (err + l[i])
      
      denom <- denom * (1 + err / l[i])
      
    }
    
    # inverse kronecker delta
    
    ikds <- as.numeric(sumdiffs > 0)
    
    diag(ikds <- 1)
    
    denom <- sqrt(denom) * ikds
    
    K <- exp(-0.5 * sumdiffs) / denom
    
    
    
  } else {
    
    # without error matrices
    
    sumdiffs <- 0
    
    for (i in 1:n3) {
      
      sumdiffs <- sumdiffs + dists[[i]] / l[i]
      
    }
    
    K <- exp(-0.5 * sumdiffs)  # to matrix?
    
  }
  
  
  
  if(class(sumdiffs) == 'dist') {
    
    K <- as.matrix(K)
    
    diag(K) <- 1
    
  }
  
  K
  
}

#d0----
d0 <-
  
  function(z, y) {
    
    if(length(y) != length(z)) y <- rep(y, length(z))
    
    pr <- y > 0 & y < 1
    
    npr <- !pr
    
    ans <- vector('numeric', length(y))
    
    y[pr] <- qnorm(y[pr])
    
    y[npr] <- 2 * y[npr] - 1
    
    ans[pr] <- dnorm(y[pr], z[pr], log = TRUE)
    
    ans[npr] <- pnorm(y[npr] * z[npr], log.p = TRUE)
    
    ans
    
  }
#d1----

d1 <-
  
  function(z, y) {
    
    pr <- y > 0 & y < 1
    
    npr <- !pr
    
    ans <- vector('numeric', length(y))
    
    y[pr] <- qnorm(y[pr])
    
    y[npr] <- 2 * y[npr] - 1
    
    ans[pr] <- y[pr] - z[pr]
    
    ans[npr] <- y[npr] * dnorm(z[npr]) / pnorm(y[npr] * z[npr])
    
    ans
    
  }

#d2----

d2 <-
  
  function(z, y) {
    
    pr <- y > 0 & y < 1
    
    npr <- !pr
    
    ans <- vector('numeric', length(y))
    
    y[npr] <- 2 * y[npr] - 1
    
    ans[pr] <- -1
    
    a <- dnorm(z[npr]) ^ 2 / pnorm(y[npr] * z[npr]) ^ 2
    
    b <- y[npr] * z[npr] * dnorm(z[npr]) / pnorm(y[npr] * z[npr])
    
    ans[npr] <- -a - b
    
    ans
    
  }
#d3----

d3 <- function(z, y) {
  
  pr <- y > 0 & y < 1
  
  npr <- !pr
  
  ans <- vector('numeric', length(y))
  
  y[npr] <- 2 * y[npr] - 1
  
  ans[pr] <- 0
  
  n <- dnorm(z[npr])
  
  p <- pnorm(y[npr] * z[npr])
  
  f <- z[npr]
  
  a <- n / p ^ 3
  
  b <- f * (y[npr] ^ 2 + 2) * n * p
  
  c <- y[npr] * (f ^ 2 - 1) * p ^ 2 + 2 * y[npr] * n ^ 2
  
  ans[npr] <- a * (b + c)
  
  ans
  
}

#psiline----
psiline <-
  
  function(s, adiff, a, K, y, d0, mn = 0, wt) {
    
    a <- a + s * as.vector(adiff)
    
    f <- K %*% a + mn
    
    psi(a, f, mn, y, d0, wt)
    
  }

#psi----


psi <-
  
  function(a, f, mn, y, d0, wt) {
    
    0.5 * t(a) %*% (f - mn) - sum(wt * d0(f, y))
    
  }
#cov.SE.d1----
cov.SE.d1 <- function (x, e = NULL, l) {
  
  # get gradients (matrices) of the kernel wrt. the parameters
  
  # CURRENTLY IGNORES e!!
  
  
  
  # number of parameters
  
  n <- length(l)
  
  
  
  # assign vector for gradients
  
  grads <- list()
  
  
  
  # get full covariance matrix
  
  K <- cov.SE(x1 = x, e1 = e, l = l)
  
  
  
  # loop through them
  
  for (i in 1:n) {
    
    
    
    # squared distances
    
    d2_i <- as.matrix(dist(x[, i]) ^ 2)
    
    
    
    # gradient for each parameter
    
    grads[[i]] <- K * (1 / l[i] ^ 2) * d2_i / 2
    
    
    
  }
  
  
  
  # return as a list
  
  return (grads)
  
  
  
}
#predict----
predict.graf.raster <-
  
  function (object, x, type, CI, maxn, ...) {
    
    
    
    # the following adapted from dismo:::predict
    
    
    
    variables <- colnames(object$x)
    
    
    
    if (is.null(CI)) {
      
      
      
      # if CIs aren't required
      
      out <- raster(x)
      
      
      
    } else {
      
      
      
      if (CI == 'std') {
        
        # if the SD is needed (latent case)
        
        out <- brick(x,
                     
                     nl = 2)
        
        
        
      } else {
        
        # if proper CIs are required
        
        out <- brick(x,
                     
                     nl = 3)
        
        
        
      }
      
      
      
    }
    
    
    
    ncols <- ncol(out)
    
    
    
    firstrow <- 1
    
    firstcol <- 1
    
    
    
    if (!canProcessInMemory(out, 3)) {
      
      
      
      filename <- rasterTmpFile()
      
      out <- writeStart(out, filename = filename, ...)
      
      inMemory <- FALSE
      
      
      
    } else {
      
      
      
      v <- array(NA, dim = c(nrow(out), ncol(out), nlayers(out)))
      
      inMemory <- TRUE
      
      
      
    }
    
    
    
    tr <- blockSize(out, n = nlayers(x) + 2)
    
    pb <- pbCreate(tr$n)
    
    
    
    for (i in 1:tr$n) {
      
      
      
      rr <- firstrow + tr$row[i] - 1
      
      
      
      rowvals <- getValuesBlock(x,
                                
                                row = rr,
                                
                                nrows = tr$nrows[i], 
                                
                                firstcol,
                                
                                ncols)
      
      
      
      rowvals <- rowvals[, variables, drop = FALSE]
      
      
      
      res <- matrix(NA,
                    
                    nrow = nrow(rowvals),
                    
                    ncol = nlayers(out))
      
      
      
      rowvals <- na.omit(rowvals)
      
      
      
      if (length(rowvals) > 0) {
        
        
        
        newdata <- data.frame(rowvals)
        
        
        
        # loop through and convert factors to factors
        
        for (col in ncol(newdata)) {
          
          
          
          if (is.factor(object$obsx[, col])) {
            
            
            
            newdata[, col] <- factor(newdata[, col])
            
            
            
          }
          
          
          
        }
        
        
        
        # predict to this batch of pixels
        
        p <- predict.graf(object,
                          
                          newdata = newdata,
                          
                          type = type,
                          
                          CI = CI,
                          
                          maxn = maxn,
                          
                          ...)
        
        
        
        naind <- as.vector(attr(rowvals, "na.action"))
        
        
        
        if (!is.null(naind)) {
          
          
          
          res[-naind, ] <- p
          
          
          
        } else {
          
          
          
          res <- p
          
          
          
        }
        
        
        
      }
      
      
      
      if (inMemory) {
        
        
        
        res <- array(res, dim = c(ncol(out), tr$nrows[i], nlayers(out)))
        
        res <- aperm(res, c(2, 1, 3))
        
        rows <- tr$row[i]:(tr$row[i] + tr$nrows[i] - 1)
        
        v[rows, , ] <- res
        
        
        
      } else {
        
        
        
        out <- writeValues(out, res, tr$row[i])
        
        
        
      }
      
      
      
      pbStep(pb, i)
      
      
      
    }
    
    
    
    pbClose(pb)
    
    
    
    if (inMemory) {
      
      
      
      # permute v so that it gets vectorized in the right direction
      
      v <- aperm(v, c(2, 1, 3))
      
      
      
      for (i in 1:nlayers(out)) {
        
        
        
        out <- setValues(out,
                         
                         as.vector(v[, , i]),
                         
                         layer = i)
        
        
        
      }
      
      
      
    } else {
      
      
      
      out <- writeStop(out)
      
      
      
    }
    
    
    
    names(out) <- colnames(p)
    
    
    
    return (out)
    
    
    
  }

predict.graf <-
  
  function(object, newdata = NULL, type = c('response', 'latent'),
           
           CI = 0.95, maxn = NULL, ...) {
    
    
    
    if (class(newdata) %in% c('Raster', 'RasterBrick', 'RasterStack')) {
      
      
      
      # predict to a raster
      
      ans <- predict.graf.raster(object = object,
                                 
                                 x = newdata,
                                 
                                 type = type,
                                 
                                 CI = CI,
                                 
                                 maxn = maxn,
                                 
                                 ...)
      
      
      
      return (ans)
      
      
      
    }
    
    
    
    type = match.arg(type)
    
    if (is.null(maxn)) maxn <- round(nrow(object$x) / 10)
    
    # set up data
    
    if (is.null(newdata)) {
      
      # use already set up inference data if not specified
      
      newdata <- object$x
      
      # get mean on raw data
      
      mn <- object$mnfun(object$obsx)
      
      
      
    } else {
      
      # convert any ints to numerics
      
      for(i in 1:ncol(newdata)) if (is.integer(newdata[, i])) newdata[, i] <- as.numeric(newdata[, i])
      
      
      
      if (is.data.frame(newdata) & all(sapply(object$obsx, class) == sapply(newdata, class))) {
        
        
        
        # get mean on raw data
        
        mn <- object$mnfun(newdata)
        
        
        
        k <- ncol(newdata)
        
        # numericize factors
        
        for (fac in object$facs) {
          
          newdata[, fac] <- as.numeric(newdata[, fac])
          
        }
        
        # convert to a matrix
        
        newdata <- as.matrix(newdata)
        
        # scale, if needed
        
        if (!is.null(object$scaling)) {
          
          notfacs <- (1:k)
          
          if(length(object$facs) > 0) notfacs <- notfacs[-object$facs]
          
          for(i in 1:length(notfacs)) {
            
            newdata[, notfacs[i]] <- (newdata[, notfacs[i]] - object$scaling[1, i]) / object$	scaling[2, i]
            
          }
          
        }
        
      } else {
        
        stop('newdata must be either a dataframe with the same elements as used for inference, or NULL')
        
      }
      
    }
    
    
    
    # check CI
    
    if (!is.null(CI)) {
      
      if (!(CI == 'std' & type == 'latent')) {
        
        if (CI >= 1 | CI <= 0) {
          
          stop("CI must be a number between 0 and 1, or NULL")
          
        }
        
        err <- qnorm( 1 - (1 - CI) / 2 )
        
      }
      
    }
    
    # latent case
    
    if (type == 'latent') {
      
      
      
      if (is.null(CI)) {
        
        # if CIs aren't wanted
        
        ans <- pred(newdata, object, mn, std = FALSE, maxn = maxn)
        
        colnames(ans) <- "posterior mean"
        
      } else if (CI == 'std') { # if standard deviations are wanted instead
        
        ans <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        
        colnames(ans) <- c("posterior mean", "posterior std")
        
      } else {
        
        # if they are
        
        pred <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        
        upper <- pred[, 1] + err * pred[, 2]
        
        lower <- pred[, 1] - err * pred[, 2]
        
        ans <- cbind(pred[, 1], lower, upper)
        
        colnames(ans) <- c("posterior mean", paste("lower ", round(100 * CI), "% CI", sep = ""),
                           
                           paste("upper ", round(100 * CI), "% CI", sep = ""))
        
      }
      
    } else {
      
      # response case
      
      if (is.null(CI)) {
        
        # if CIs aren't required
        
        ans <- pnorm(pred(newdata, object, mn, std = FALSE, maxn = maxn))
        
        colnames(ans) <- "posterior mode"
        
      } else {
        
        # if CIs are required
        
        pred <- pred(newdata, object, mn, std = TRUE, maxn = maxn)
        
        upper <- pred[, 1] + err * pred[, 2]
        
        lower <- pred[, 1] - err * pred[, 2]
        
        ans <- pnorm(cbind(pred[, 1], lower, upper))
        
        colnames(ans) <- c("posterior mode", paste("lower ", round(100 * CI), "% CI", sep = ""),
                           
                           paste("upper ", round(100 * CI), "% CI", sep = ""))
        
      }
      
    }
    
    ans
    
  }

pred <-
  
  function(predx, fit, mn, std = TRUE, maxn = 250, same = FALSE) {
    
    predx <- as.matrix(predx)
    
    n <- nrow(predx)  
    
    if (n > maxn & !same) {
      
      inds <- split(1:n, ceiling((1:n) / maxn))
      
      fun <- function(ind, X, fit, std, maxn) {
        
        pred(X[ind, , drop = FALSE], fit, mn[ind], std, maxn, same)
        
      }    
      
      prediction <- lapply(inds, fun, predx, fit, std, maxn)
      
      prediction <- do.call('rbind', prediction)
      
    } else {
      
      if(same) {
        
        # if predicting back to input data re-use covariance matrix
        
        Kx <- fit$K
        
        prediction <- fit$MAP
        
      } else {
        
        Kx <- cov.SE(fit$x, predx, e1 = fit$e, l = fit$ls)
        
        mpred <- qnorm(mn)
        
        prediction <- crossprod(Kx, fit$a) + mpred
        
      }
      
      
      
      if (std) {
        
        v <- backsolve(fit$L, sqrt(as.vector(fit$W)) * Kx, transpose = T)
        
        # using correlation matrix, so diag(kxx) is all 1s, no need to compute kxx
        
        predvar <- 1 - crossprod(v)
        
        prediction <- cbind(prediction, sqrt(diag(predvar)))
        
        colnames(prediction) <- c("MAP", "std")
        
      }
      
    }
    
    prediction
    
  }

#Plot
plot.graf <-
  function(x, vars = NULL, resolution = 50, CI = 0.95, prior = FALSE, data = TRUE, jitter = 1,
           peak = FALSE, ...) {
    # find factors
    facs <- x$facs
    # if not plotting at peak, find column menas, accounting for factors
    if (!peak) {
      if(length(facs) > 0) {
        x$peak[(1:ncol(x$obsx))[-facs]] <- colMeans(x$obsx[, -facs, drop = FALSE])
        x$peak[, facs] <- sapply(data.frame(x$obsx[, facs, drop = FALSE]),
                                 function(x) names(sort(table(x), decreasing = TRUE))[1])
        for(i in facs) x$peak[, i] <- factor(x$peak[, i], levels = levels(x$obsx[, i]))
      } else {
        x$peak <- colMeans(x$obsx)
      }
    }
    pars <- x$peak
    k <- ncol(x$obsx)
    if (is.null(vars)) vars <- 1:k
    X <- as.data.frame(lapply(pars, rep, resolution))
    
    # calculate scaled data for mean function
    scaledX <- X
    if (!is.null(x$scaling)) {
      notfacs <- (1:k)
      if(length(facs) > 0) notfacs <- notfacs[-facs]
      for(i in 1:length(notfacs)) {
        #      if (i %in% notfacs) {
        scaledX[, notfacs[i]] <- (scaledX[, notfacs[i]] - x$scaling[1, i]) / x$scaling[2, i]
        #      }
      }
    }
    for (i in vars) {
      predX <- X
      #    predscaledX <- scaledX
      if(i %in% facs) {
        # for factors
        seqx <- levels(X[, i])
        
        m <- length(seqx)
        predX[1:m, i] <- seqx
        predX <- predX[1:m, ]
        
        #      predscaledX[1:m, i] <- seqx
        #      predscaledX <- predscaledX[1:m, ]
        
        p <- predict(x, newdata = predX, CI = CI)
        mn <- p[, 1]
        up <- p[, 3]
        low <- p[, 2]
        width = 0.2
        
        x.plot <- as.numeric(predX[, i])
        
        nam <- ifelse(is.null (colnames(x$x)[i]), "covariate", colnames(x$x)[i]) 
        # set up plot
        plot(mn ~ predX[, i], ylim = c(0, 1), type = 'p', xlab = nam, border ='white',
             ylab = 'probability of presence', ...)
        
        # plot CIs
        rect(xleft = x.plot - width, ybottom = low, xright = x.plot + width, ytop = up, col = 'light grey',
             border = NA)
        # plot mean
        for(j in 1:m) {
          lines(x = c(x.plot[j] - width, x.plot[j] + width), y = rep(mn[j], 2), col = rgb(0.4, 0.4, 0.4), lwd = 2)
        }
        
        if (prior) {
          xs <- sort(x.plot)
          xs[c(1, m)] <- c(-1, m + 1)        
          lines(x$mnfun(predX) ~ xs, lty = 2, col = rgb(0.4, 0.4, 0.4))
        }
        
        if (data) {
          rug(jitter(x$x[x$obsy == 1, i], jitter), side = 3, col = 'dark grey')
          rug(jitter(x$x[x$obsy == 0, i], jitter), side = 1, col = 'dark grey')
          box()
        }
        
      } else {
        # for continuous variables
        seqx <- seq(min(x$obsx[, i]), max(x$obsx[, i]), length.out = resolution)
        
        predX[, i] <- seqx
        
        #     predscaledX[, i] <- (seqx - x$scaling[1, i]) / x$scaling[2, i]
        
        
        p <- predict(x, newdata = predX, CI = CI)
        mn <- p[, 1]
        up <- p[, 3]
        low <- p[, 2]
        nam <- ifelse(is.null (colnames(x$x)[i]), "covariate", colnames(x$x)[i]) 
        plot(mn ~ seqx, ylim = c(0, 1), type = 'n', xlab = nam,
             ylab = 'probability of presence', ...)
        polygon(x = c(seqx, rev(seqx)), y = c(up, rev(low)), col = 'light grey',
                border = NA)
        lines(mn ~ seqx, col = rgb(0.4, 0.4, 0.4), lwd = 2)
        if (prior) {
          lines(x$mnfun(predX) ~ seqx, lty = 2, col = 'dark grey')
        }
        if (data) {
          rug(jitter(x$obsx[x$obsy == 1, i], jitter), side = 3, col = 'dark grey')
          rug(jitter(x$obsx[x$obsy == 0, i], jitter), side = 1, col = 'dark grey')
          box()
        }
      }
    }
  }