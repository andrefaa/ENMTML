maxnet2 <- function(p,
                    data,
                    f = maxnet.formula(p, data),
                    regmult = 1.0,
                    regfun = maxnet::maxnet.default.regularization,
                    ...)
  
{
  mm <- stats::model.matrix(f, data)
  reg <- regfun(p, mm) * regmult
  weights <- p + (1 - p) * 100
  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  model <-
    glmnet::glmnet(
      x = mm,
      y = as.factor(p),
      family = "binomial",
      standardize = F,
      penalty.factor = reg,
      lambda = 10 ^ (seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
        sum(p) / sum(weights),
      weights = weights
    )
  class(model) <- c("maxnet", class(model))
  if (length(model$beta) < 200)
    stop("Error: glmnet failed to complete regularization path")
  bb <-
    model$beta[, ncol(model$beta)]#Was originally [,200], changed to [,ncol(model$beta)]
  model$betas <- bb[bb != 0]
  model$alpha <- 0
  rr <-
    predict.maxnet(model, data[p == 0, , drop = TRUE], type = "exponent", clamp =
                     F)
  raw <- as.matrix(rr / sum(rr))
  raw <- raw[raw != 0]
  model$entropy <- -sum(raw * log(raw))
  model$alpha <- -log(sum(rr))
  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)
  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)
  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv],
                       function(n)
                         which.max(table(data[p == 1, n, drop = FALSE])))
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- c(means, majorities)
  model$levels <- lapply(data, levels)
  model
}

predict.maxnet <-
  function(object,
           newdata,
           clamp = T,
           type = c("link", "exponential", "cloglog", "logistic"),
           ...)
  {
    if (clamp) {
      for (v in intersect(names(object$varmax), names(newdata))) {
        newdata[, v] <-
          pmin(pmax(newdata[, v], object$varmin[v]), object$varmax[v])
      }
    }
    terms <-
      sub("hinge\\((.*)\\):(.*):(.*)$",
          "hingeval(\\1,\\2,\\3)",
          names(object$betas))
    terms <-
      sub("categorical\\((.*)\\):(.*)$",
          "categoricalval(\\1,\\2)",
          terms)
    terms <-
      sub("thresholds\\((.*)\\):(.*)$",
          "thresholdval(\\1,\\2)",
          terms)
    f <- stats::formula(paste("~", paste(terms, collapse = " + "), "-1"))
    mm <- stats::model.matrix(f, data.frame(newdata))
    if (clamp)
      mm <-
      t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
             object$featuremaxs[names(object$betas)]))
    link <- (mm %*% object$betas) + object$alpha
    type <- match.arg(type)
    if (type == "link")
      return(link)
    if (type == "exponential")
      return(exp(link))
    if (type == "cloglog")
      return(1 - exp(0 - exp(object$entropy + link)))
    if (type == "logistic")
      return(1 / (1 + exp(-object$entropy - link)))
  }
