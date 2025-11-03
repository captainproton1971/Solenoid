library("numDeriv")
## ======================== York-style poly regression (ML + sandwich + PI) ========================

## ---------- helpers ----------
.york_basis <- function(x, degree, intercept = TRUE) {
  stopifnot(degree %in% c(1, 2))
  if (degree == 1) {
    if (intercept) cbind(1, x) else cbind(x)
  } else {
    if (intercept) cbind(1, x, x^2) else cbind(x, x^2)
  }
}

.york_dfdx <- function(x, theta, intercept = TRUE) {
  # theta ordered to match .york_basis
  if (length(theta) == 1L) {
    # linear, no intercept: f'(x) = a1
    rep(theta[1L], length(x))
  } else if (length(theta) == 2L) {
    # either linear+intercept or quadratic no-intercept
    # linear+int: (a0,a1) -> f' = a1
    # quad no-int: (a1,a2) -> f' = a1 + 2 a2 x
    if (intercept) rep(theta[2L], length(x)) else theta[1L] + 2*theta[2L]*x
  } else {
    # quadratic + intercept: (a0,a1,a2) -> f' = a1 + 2 a2 x
    theta[2L] + 2*theta[3L]*x
  }
}

# Gradient of f'(x; theta) w.r.t. theta (row vector per x)
.york_dfdx_grad_row <- function(xi, p, degree, intercept) {
  # returns length-p vector: d f'(xi) / d theta
  if (degree == 1) {
    if (intercept) {
      # theta=(a0,a1), f' = a1
      c(0, 1)
    } else {
      # theta=(a1)
      c(1)
    }
  } else { # degree==2
    if (intercept) {
      # theta=(a0,a1,a2), f' = a1 + 2 a2 x
      c(0, 1, 2*xi)
    } else {
      # theta=(a1,a2)
      c(1, 2*xi)
    }
  }
}

# collapse accidental "name.name" duplicates; also strip stray backticks
.clean_coef_names <- function(nm) {
  nm <- as.character(nm)
  nm <- sub("^(.+)\\.?\\1$", "\\1", nm, perl = TRUE)  # "(Intercept).(Intercept)" -> "(Intercept)", "x.x"->"x"
  nm <- sub("^`(.+)`$", "\\1", nm, perl = TRUE)
  nm
}

# Build Jacobian J = d beta_orig / d alpha_std and back-transform mapping.
# Returns list(beta, J) where alpha_std are standardized-space coefficients.
.york_backtransform <- function(alpha, degree, intercept,
                                mx, sx, my, sy,
                                mode = c("z","scale")) {
  mode <- match.arg(mode)
  if (degree == 1) {
    if (intercept) {
      # y' = a0 + a1 x', x'=(x-mx)/sx, y = sy*y' + my
      a0 <- alpha[1L]; a1 <- alpha[2L]
      beta0 <- sy * (a0 - a1 * mx / sx) + my
      beta1 <- sy * a1 / sx
      beta  <- c("(Intercept)" = beta0, "x" = beta1)
      J <- rbind(
        c(sy,         -sy*mx/sx),
        c(0,           sy/sx)
      )
    } else {
      # No-intercept linear. Enforce scale-only to respect zero-intercept structure.
      a1 <- alpha[1L]
      beta1 <- sy * a1 / sx
      beta <- c("x" = beta1)
      J <- matrix(sy/sx, nrow = 1L, ncol = 1L)
    }
  } else { # degree == 2
    if (intercept) {
      # y' = a0 + a1 x' + a2 x'^2, x'=(x-mx)/sx, y = sy*y' + my
      a0 <- alpha[1L]; a1 <- alpha[2L]; a2 <- alpha[3L]
      beta2 <- sy * a2 / sx^2
      beta1 <- sy * (a1 / sx - 2 * a2 * mx / sx^2)
      beta0 <- sy * (a0 - a1 * mx / sx + a2 * mx^2 / sx^2) + my
      beta <- c("(Intercept)" = beta0, "x" = beta1, "I(x^2)" = beta2)
      J <- rbind(
        c(sy,            -sy*mx/sx,         sy*mx^2/sx^2),
        c(0,              sy/sx,           -2*sy*mx/sx^2),
        c(0,              0,                sy/sx^2)
      )
    } else {
      # Quad without intercept: y' = a1 x' + a2 x'^2 => y = (sy a1/sx) x + (sy a2/sx^2) x^2
      a1 <- alpha[1L]; a2 <- alpha[2L]
      beta1 <- sy * a1 / sx
      beta2 <- sy * a2 / sx^2
      beta <- c("x" = beta1, "I(x^2)" = beta2)
      J <- rbind(
        c(sy/sx,     0),
        c(0,         sy/sx^2)
      )
    }
  }
  list(beta = beta, J = J)
}

# Standardize (center/scale) according to mode.
# Returns x_s, y_s, sx_s, sy_s and the used mx, sx, my, sy.
.york_standardize <- function(x, y, sx, sy,
                              intercept, mode = c("auto","z","scale")) {
  mode <- match.arg(mode)
  if (mode == "auto") mode <- if (intercept) "z" else "scale"

  mx <- if (mode == "z") mean(x, na.rm = TRUE) else 0
  my <- if (mode == "z") mean(y, na.rm = TRUE) else 0

  sxdat <- sd(x, na.rm = TRUE)
  sydat <- sd(y, na.rm = TRUE)
  if (!is.finite(sxdat) || sxdat <= 0) stop("Cannot standardize: sd(x) must be positive.")
  if (!is.finite(sydat) || sydat <= 0) stop("Cannot standardize: sd(y) must be positive.")

  x_s  <- (x - mx) / sxdat
  y_s  <- (y - my) / sydat
  sx_s <- sx / sxdat
  sy_s <- sy / sydat

  list(
    x = x_s, y = y_s, sx = sx_s, sy = sy_s,
    mx = mx, sxdat = sxdat, my = my, sydat = sydat,
    mode = mode
  )
}

.safe_solve <- function(A, b = NULL) {
  out <- tryCatch({
    if (is.null(b)) solve(A) else solve(A, b)
  }, error = function(e) {
    if (is.null(b)) chol2inv(chol(A)) else qr.solve(A, b)
  })
  out
}

# negative log-likelihood (sum over i) for given theta and s2 (in working space)
.york_nll <- function(theta, log_s, x, y, sx, sy, X, intercept, degree) {
  mu <- as.numeric(X %*% theta)
  fp <- .york_dfdx(x, theta, intercept)
  s2 <- exp(2*log_s)
  v  <- sy^2 + (fp^2) * sx^2 + s2
  v  <- pmax(v, .Machine$double.eps)
  r  <- y - mu
  0.5 * sum(log(v) + (r^2)/v)
}

# Per-observation score wrt theta (holding s2 fixed) for sandwich; returns n x p matrix
.york_scores_theta <- function(theta, s2, x, y, sx, sy, X, intercept, degree) {
  mu <- as.numeric(X %*% theta)
  r  <- y - mu
  fp <- .york_dfdx(x, theta, intercept)
  v  <- sy^2 + (fp^2) * sx^2 + s2
  v  <- pmax(v, .Machine$double.eps)

  p <- length(theta)
  n <- length(y)
  S <- matrix(0, nrow = n, ncol = p)
  # dmu/dtheta = X row; dv/dtheta = 2 fp * (dfp/dtheta) * sx^2
  for (i in seq_len(n)) {
    Xi <- X[i, , drop = TRUE]
    dmu <- Xi
    dfp <- .york_dfdx_grad_row(x[i], p, degree, intercept)
    dv  <- 2 * fp[i] * dfp * sx[i]^2
    # score_i = (r/v) dmu + 0.5*(r^2/v^2 - 1/v) dv
    S[i, ] <- (r[i]/v[i]) * dmu + 0.5 * ( (r[i]^2)/(v[i]^2) - 1/v[i] ) * dv
  }
  S
}

# Sandwich covariance for theta via numerical Hessian and analytic scores; returns p x p
.york_cov_sandwich <- function(theta, s2, x, y, sx, sy, X, intercept, degree) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) return(NULL)
  nll_theta_only <- function(th) .york_nll(th, log_s = 0.5*log(pmax(s2, 1e-16)),
                                           x, y, sx, sy, X, intercept, degree)
  H <- try(numDeriv::hessian(nll_theta_only, theta), silent = TRUE)
  if (inherits(H, "try-error")) return(NULL)
  S <- .york_scores_theta(theta, s2, x, y, sx, sy, X, intercept, degree)
  J <- crossprod(S)  # sum_i s_i s_i^T
  Hinv <- .safe_solve(H)
  Hinv %*% J %*% Hinv
}

## ---------- core fitter ----------
york_poly_fit <- function(x, y, sx, sy,
                          degree = 1,
                          intercept = TRUE,
                          estimate_sigma_int = TRUE,
                          maxit = 200, tol = 1e-10,
                          na.action = na.omit,
                          ridge_factor = 0,
                          standardize = FALSE,
                          standardize_mode = c("auto","z","scale"),
                          covariance = c("sandwich","naive"),
                          track_objective = TRUE) {
  covariance <- match.arg(covariance)
  standardize_mode <- match.arg(standardize_mode)
  stopifnot(length(x)==length(y), length(x)==length(sx), length(x)==length(sy))

  # sanitize & validate
  sx <- as.numeric(sx); sy <- as.numeric(sy)
  if (any(!is.finite(sx) | !is.finite(sy))) {
    stop("All sx and sy must be finite numbers.")
  }
  if (any(sx < 0 | sy < 0, na.rm = TRUE)) {
    warning("Negative sx or sy detected; replacing with absolute values.")
    sx <- abs(sx); sy <- abs(sy)
  }
  # clip tiny variances to avoid degeneracy
  sx[sx == 0] <- 0
  sy[sy == 0] <- 0

  degree <- as.integer(degree); stopifnot(degree %in% c(1L,2L))

  dat <- data.frame(x=x, y=y, sx=sx, sy=sy)
  d  <- na.action(dat)
  na <- attr(d, "na.action")
  x <- d$x; y <- d$y; sx <- d$sx; sy <- d$sy

  n <- length(y)

  # ====== (1) Optional standardization ======
  if (isTRUE(standardize)) {
    std <- .york_standardize(x, y, sx, sy, intercept, mode = standardize_mode)
    x_fit  <- std$x;  y_fit  <- std$y
    sx_fit <- std$sx; sy_fit <- std$sy
  } else {
    std <- list(mx = 0, my = 0, sxdat = 1, sydat = 1, mode = "none")
    x_fit <- x; y_fit <- y; sx_fit <- sx; sy_fit <- sy
  }

  X <- .york_basis(x_fit, degree, intercept)
  p <- ncol(X)

  # OLS init (ignores errors) in working space
  theta <- as.numeric(coef(lm.fit(x = X, y = y_fit)))
  if (length(theta) < p) theta <- c(theta, rep(0, p - length(theta)))

  s2_int <- 0
  conv <- FALSE
  iters <- 0L
  obj_old <- Inf

  # IRLS skeleton (for a solid starting point)
  for (iter in 1:maxit) {
    iters <- iter

    # weights from current theta
    fprime <- .york_dfdx(x_fit, theta, intercept)
    v <- sy_fit^2 + (fprime^2) * sx_fit^2 + s2_int
    w <- 1 / pmax(v, .Machine$double.eps)

    XtW  <- t(X) * w
    XtWX <- XtW %*% X
    if (ridge_factor > 0) {
      XtWX <- XtWX + diag(ridge_factor * mean(diag(XtWX)), ncol(XtWX))
    }
    XtWy <- XtW %*% y_fit

    theta_new <- as.numeric(.safe_solve(XtWX, XtWy))
    mu <- as.numeric(X %*% theta_new)
    r  <- y_fit - mu

    # recompute fprime for theta_new BEFORE estimating s2_int
    if (estimate_sigma_int) {
      fprime_new <- .york_dfdx(x_fit, theta_new, intercept)
      target <- max(1, n - p)
      f <- function(s2) {
        vtmp <- sy_fit^2 + (fprime_new^2) * sx_fit^2 + pmax(s2, 0)
        sum((r^2) / vtmp) - target
      }
      lo <- 0
      hi <- max(1e6, 10*mean(sy_fit^2, na.rm=TRUE))
      tries <- 0L
      while (is.finite(f(hi)) && f(hi) > 0 && tries < 6L) {
        hi <- hi * 10
        tries <- tries + 1L
      }
      s2_int_new <- s2_int
      flo <- f(lo); fhi <- f(hi)
      if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
        root <- try(uniroot(f, c(lo, hi))$root, silent=TRUE)
        if (!inherits(root, "try-error")) s2_int_new <- max(0, root)
      } else {
        s2_int_new <- if (is.finite(flo) && flo <= 0) 0 else hi
      }
      s2_int <- s2_int_new
    }

    # optional objective tracking for convergence: full NLL
    if (isTRUE(track_objective)) {
      nll <- .york_nll(theta_new, log_s = 0.5*log(pmax(s2_int, 1e-16)),
                       x_fit, y_fit, sx_fit, sy_fit, X, intercept, degree)
      if (abs(nll - obj_old) < tol * (1 + abs(obj_old))) { theta <- theta_new; conv <- TRUE; break }
      obj_old <- nll
    } else {
      if (max(abs(theta_new - theta)) < tol * (1 + max(1, max(abs(theta))))) {
        theta <- theta_new; conv <- TRUE; break
      }
    }
    theta <- theta_new
  }
  if (!conv) warning("york_poly_fit: IRLS stage did not fully converge (continuing to ML polish).")

  # ====== (1b) ML polish over (theta, log s) in working space ======
  if (isTRUE(estimate_sigma_int)) {
    par0 <- c(theta, log_s = 0.5*log(pmax(s2_int, 1e-16)))
  } else {
    # lock sigma_int=0
    par0 <- c(theta, log_s = -30)
  }

  npar <- length(par0)
  nll_all <- function(par) {
    th <- par[seq_len(npar-1)]
    ls <- par[npar]
    .york_nll(th, ls, x_fit, y_fit, sx_fit, sy_fit, X, intercept, degree)
  }
  opt <- try(optim(par0, nll_all, method = "BFGS", control = list(maxit = 400)), silent = TRUE)
  if (!inherits(opt, "try-error") && opt$convergence == 0) {
    theta <- opt$par[seq_len(npar-1)]
    s2_int <- if (isTRUE(estimate_sigma_int)) exp(2*opt$par[npar]) else 0
    conv <- TRUE
  } else {
    warning("ML polish did not converge; keeping IRLS estimates.")
  }

  # ====== (2) Final summaries in working space ======
  fprime <- .york_dfdx(x_fit, theta, intercept)
  v <- sy_fit^2 + (fprime^2) * sx_fit^2 + s2_int
  w <- 1 / pmax(v, .Machine$double.eps)

  XtW <- t(X) * w; XtWX <- XtW %*% X; XtWy <- XtW %*% y_fit
  if (ridge_factor > 0) {
    XtWX <- XtWX + diag(ridge_factor * mean(diag(XtWX)), ncol(XtWX))
  }

  alpha_hat <- as.numeric(.safe_solve(XtWX, XtWy)) # standardized-space coefficients
  names(alpha_hat) <- if (degree==1) c(if (intercept) "(Intercept)", "x") else
    c(if (intercept) "(Intercept)", "x", "I(x^2)")

  # Covariance in standardized space
  cov_unscaled_std_naive <- tryCatch({
    .safe_solve(XtWX)
  }, error = function(e) {
    .safe_solve(XtWX)
  })

  cov_unscaled_std <- cov_unscaled_std_naive
  if (covariance == "sandwich") {
    cov_sand <- .york_cov_sandwich(alpha_hat, s2_int, x_fit, y_fit, sx_fit, sy_fit, X, intercept, degree)
    if (!is.null(cov_sand)) cov_unscaled_std <- cov_sand
  }

  # ====== (3) Back-transform coefficients & covariance to original scale ======
  if (isTRUE(standardize)) {
    bt <- .york_backtransform(
      alpha = alpha_hat,
      degree = degree,
      intercept = intercept,
      mx = std$mx, sx = std$sxdat,
      my = std$my, sy = std$sydat,
      mode = std$mode
    )
    beta <- bt$beta
    names(beta) <- .clean_coef_names(names(beta))
    J    <- bt$J
    cov_unscaled <- J %*% cov_unscaled_std %*% t(J)
    sigma_int_out <- sqrt(max(std$sydat^2 * s2_int, 0)) # back to original y units
  } else {
    beta <- alpha_hat
    names(beta) <- .clean_coef_names(names(beta))
    cov_unscaled <- cov_unscaled_std
    sigma_int_out <- sqrt(max(s2_int, 0))
  }

  # Fitted, residuals, weights recomputed on ORIGINAL scale so summaries match user scale
  X_orig <- .york_basis(x, degree, intercept)
  mu <- as.numeric(X_orig %*% beta)
  r  <- y - mu

  # For weights on original scale, recompute v using original (x, sx, sy, sigma_int_out) and original derivative
  fprime_orig <- .york_dfdx(x, beta, intercept)
  v_orig <- sy^2 + (fprime_orig^2) * sx^2 + sigma_int_out^2
  w_orig <- 1 / pmax(v_orig, .Machine$double.eps)
  wr <- sqrt(w_orig) * r

  df_resid <- max(0, n - length(beta))
  sigma2_hat <- if (df_resid > 0) sum(w_orig * r^2) / df_resid else NA_real_
  sigma_hat  <- sqrt(sigma2_hat)

  se <- sqrt(pmax(diag(cov_unscaled), .Machine$double.eps))
  tvals <- beta / se
  pvals <- 2 * pt(-abs(tvals), df = df_resid)

  # weighted pseudo-R^2 (on original scale)
  if (intercept) {
    ybar_w <- sum(w_orig*y)/sum(w_orig); TSS_w <- sum(w_orig*(y - ybar_w)^2)
  } else {
    TSS_w <- sum(w_orig*y^2)
  }
  RSS_w <- sum(w_orig*r^2)
  r2 <- if (TSS_w > 0) 1 - RSS_w/TSS_w else NA_real_
  adj_r2 <- if (!is.na(r2) && df_resid > 0) 1 - (1 - r2)*(n - as.integer(intercept))/df_resid else NA_real_

  # F-stat
  p_eff <- length(beta)
  df_num <- p_eff - as.integer(intercept)
  fstat_vec <- if (df_num > 0 && df_resid > 0 && is.finite(RSS_w) && is.finite(TSS_w)) {
    fval <- ((TSS_w - RSS_w)/df_num) / (RSS_w/df_resid)
    c(value = fval, numdf = df_num, dendf = df_resid)
  } else NA

  out <- list(
    call = match.call(),
    terms = if (degree==1)
      as.formula(if (intercept) "y ~ x" else "y ~ 0 + x") else
        as.formula(if (intercept) "y ~ x + I(x^2)" else "y ~ 0 + x + I(x^2)"),
    degree = degree,
    intercept = intercept,
    coefficients = beta,
    fitted.values = mu,
    residuals = r,
    weighted.residuals = wr,
    weights = w_orig,
    v = v_orig,
    sx = sx, sy = sy,
    sigma_int = sigma_int_out,
    cov.unscaled = cov_unscaled,
    df = c(p_eff, df_resid, p_eff),
    r.squared = r2,
    adj.r.squared = adj_r2,
    fstatistic = fstat_vec,
    na.action = na,
    converged = conv,
    iterations = iters,
    model = data.frame(x = x, y = y),
    # standardization metadata
    standardize = isTRUE(standardize),
    standardize_mode = if (isTRUE(standardize)) std$mode else "none",
    scaling = if (isTRUE(standardize)) list(
      x_mean = std$mx, y_mean = std$my,
      x_sd = std$sxdat, y_sd = std$sydat
    ) else NULL,
    covariance = covariance
  )
  class(out) <- c("yorklm","list")
  out
}

## ---------- slope CI (delta method) ----------
slope_ci <- function(object, level = 0.95, at = 0) {
  stopifnot(inherits(object, "yorklm") || inherits(object, "summary.yorklm"))

  if (inherits(object, "summary.yorklm")) {
    coef <- as.numeric(object$coefficients[, "Estimate"])
    V    <- object$cov.unscaled
    df   <- object$df[2L]
    degree <- object$degree
    intercept <- object$intercept
  } else {
    coef <- object$coefficients
    V    <- object$cov.unscaled
    df   <- object$df[2L]
    degree <- object$degree
    intercept <- object$intercept
  }

  key <- paste0(degree, "_", as.integer(intercept))
  g <- switch(key,
              "1_1" = c(0, 1),           # (a0, a1) -> slope = a1
              "1_0" = c(1),              # (a1)
              "2_1" = c(0, 1, 2*at),     # (a0, a1, a2)
              "2_0" = c(1, 2*at)         # (a1, a2)
  )
  slope <- switch(key,
                  "1_1" = coef[2L],
                  "1_0" = coef[1L],
                  "2_1" = coef[2L] + 2*coef[3L]*at,
                  "2_0" = coef[1L] + 2*coef[2L]*at
  )
  se2 <- as.numeric(t(g) %*% V %*% g)
  se  <- sqrt(pmax(se2, .Machine$double.eps))
  crit <- if (is.finite(df) && df > 0) qt(1 - (1 - level)/2, df) else qnorm(1 - (1 - level)/2)

  list(est = slope, se = se, lwr = slope - crit*se, upr = slope + crit*se)
}

## ---------- summary / print ----------
summary.yorklm <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  p <- length(object$coefficients)
  df_resid <- object$df[2L]

  se    <- sqrt(pmax(diag(object$cov.unscaled), .Machine$double.eps))
  tvals <- object$coefficients / se
  pvals <- 2 * pt(-abs(tvals), df = df_resid)

  crit <- if (is.finite(df_resid) && df_resid > 0) qt(0.975, df_resid) else qnorm(0.975)

  lwr <- object$coefficients - crit * se
  upr <- object$coefficients + crit * se

  coef_tab <- cbind(
    Estimate     = object$coefficients,
    `Std. Error` = se,
    `t value`    = tvals,
    `2.5 %`      = lwr,
    `97.5 %`     = upr,
    `Pr(>|t|)`   = pvals
  )
  coef_tab <- as.matrix(coef_tab)
  rownames(coef_tab) <- .clean_coef_names(names(object$coefficients))

  ans <- list(
    call = object$call,
    residuals = object$weighted.residuals,
    coefficients = coef_tab,
    aliased = setNames(rep(FALSE, p), rownames(coef_tab)),
    sigma = sqrt(if (is.finite(df_resid))
      sum(object$weights * object$residuals^2) / max(1, df_resid)
      else NA_real_),
    df = object$df,
    fstatistic = object$fstatistic,
    r.squared = object$r.squared,
    adj.r.squared = object$adj.r.squared,
    cov.unscaled = object$cov.unscaled,
    correlation = if (isTRUE(correlation)) {
      csd <- sqrt(diag(object$cov.unscaled))
      object$cov.unscaled / outer(csd, csd)
    } else NULL,
    symbolic.cor = isTRUE(symbolic.cor),
    na.action = object$na.action,
    # carry through for slope_ci
    degree = object$degree,
    intercept = object$intercept
  )
  class(ans) <- "summary.yorklm"
  ans
}

print.yorklm <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nConverged:", x$converged, "in", x$iterations, "iterations\n")
  cat(sprintf("Estimated sigma_int (original units): %.6g\n", x$sigma_int))
  if (isTRUE(x$standardize)) {
    s <- x$scaling
    cat(sprintf("Standardization: %s; x_mean=%.6g, x_sd=%.6g; y_mean=%.6g, y_sd=%.6g\n\n",
                x$standardize_mode, s$x_mean, s$x_sd, s$y_mean, s$y_sd))
  } else {
    cat("\n")
  }
  tmp <- x$coefficients
  names(tmp) <- .clean_coef_names(names(tmp))
  cat("Coefficients (original scale):\n"); print(tmp)
  invisible(x)
}

print.summary.yorklm <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nWeighted Residuals:\n")
  rq <- stats::quantile(x$residuals, na.rm = TRUE)
  print(structure(rq, names = c("Min","1Q","Median","3Q","Max")))
  cat("\nCoefficients (original scale):\n"); printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat(sprintf("\nResidual standard error: %.5g on %d degrees of freedom\n",
              x$sigma, x$df[2]))
  if (all(is.finite(x$fstatistic))) {
    cat(sprintf("Multiple R-squared (weighted): %.5g, Adjusted R-squared: %.5g\n",
                x$r.squared, x$adj.r.squared))
    cat(sprintf("F-statistic: %.5g on %d and %d DF\n",
                x$fstatistic[1], x$fstatistic[2], x$fstatistic[3]))
  } else {
    cat(sprintf("Multiple R-squared (weighted): %.5g, Adjusted R-squared: %.5g\n",
                x$r.squared, x$adj.r.squared))
  }
  ci0 <- try(slope_ci(x, level = 0.95, at = 0), silent = TRUE)
  if (!inherits(ci0, "try-error")) {
    cat(sprintf("Slope at x=0 (95%% CI): [%g, %g]\n", ci0$lwr, ci0$upr))
  }
  invisible(x)
}

## ---------- predict (adds prediction intervals) ----------
predict.yorklm <- function(object, newdata = NULL, se.fit = FALSE,
                           interval = c("none","confidence","prediction"),
                           level = 0.95,
                           sx_new = NULL, sy_new = NULL, ...) {
  interval <- match.arg(interval)
  newx <- if (is.null(newdata)) object$model$x else if (!is.null(newdata$x)) newdata$x else as.numeric(newdata)
  Xg <- .york_basis(newx, object$degree, object$intercept)
  fit <- as.numeric(Xg %*% object$coefficients)
  if (!se.fit && interval == "none") return(fit)

  V <- object$cov.unscaled
  var_fit <- rowSums((Xg %*% V) * Xg)   # coeff-uncertainty

  if (interval == "prediction") {
    # include measurement + intrinsic scatter at newx
    if (is.null(sx_new)) sx_new <- rep(median(object$sx, na.rm=TRUE), length(newx))
    if (is.null(sy_new)) sy_new <- rep(median(object$sy, na.rm=TRUE), length(newx))
    fp_new <- .york_dfdx(newx, object$coefficients, object$intercept)
    v_new  <- sy_new^2 + (fp_new^2) * sx_new^2 + object$sigma_int^2
    var_fit <- var_fit + v_new
  }

  se <- sqrt(pmax(var_fit, .Machine$double.eps))
  df <- object$df[2L]
  crit <- if (is.finite(df) && df > 0) qt(1 - (1 - level)/2, df) else qnorm(1 - (1 - level)/2)
  list(fit = fit, se.fit = se, lwr = fit - crit*se, upr = fit + crit*se, x = newx)
}

## ---------- ggplot2 method wrapper ----------
# ggplot2 calls predictdf(model, xseq, se, level)
predictdf.yorklm <- function(model, xseq, se, level) {
  Xg  <- .york_basis(xseq, model$degree, model$intercept)
  fit <- as.numeric(Xg %*% model$coefficients)

  if (isTRUE(se)) {
    V      <- model$cov.unscaled
    varfit <- rowSums((Xg %*% V) * Xg)
    sefit  <- sqrt(pmax(varfit, .Machine$double.eps))
    df     <- model$df[2L]
    mult   <- if (is.finite(df) && df > 0) qt(1 - (1 - level)/2, df) else qnorm(1 - (1 - level)/2)
    data.frame(x = xseq, y = fit, ymin = fit - mult*sefit, ymax = fit + mult*sefit)
  } else {
    data.frame(x = xseq, y = fit, ymin = fit, ymax = fit)
  }
}

york_smooth <- function(formula, data,
                        degree = 1,
                        intercept = TRUE,
                        sx_col = "sx", sy_col = "sy",
                        sx = NULL, sy = NULL,
                        estimate_sigma_int = TRUE,
                        standardize = FALSE,
                        standardize_mode = c("auto","z","scale"),
                        covariance = c("sandwich","naive"),
                        ...) {
  covariance <- match.arg(covariance)
  standardize_mode <- match.arg(standardize_mode)
  tt <- terms(formula, data = data)
  has_intercept <- attr(tt, "intercept") == 1
  intercept <- intercept && has_intercept

  mf <- model.frame(tt, data = data)
  y <- model.response(mf)
  X <- model.matrix(tt, mf)

  if (ncol(X) < 1)
    stop("Formula must have a single predictor, e.g. y ~ x or y ~ 0 + x")

  x <- as.numeric(X[, ncol(X)])

  # pull sx/sy either from explicit args or from data columns
  if (!is.null(sx) && !is.null(sy)) {
    sxv <- as.numeric(sx); syv <- as.numeric(sy)
  } else {
    if (is.null(data[[sx_col]]) || is.null(data[[sy_col]]))
      stop(sprintf("Columns '%s' and '%s' not found in data.", sx_col, sy_col))
    sxv <- as.numeric(data[[sx_col]]); syv <- as.numeric(data[[sy_col]])
  }

  york_poly_fit(x, y, sxv, syv,
                degree = degree, intercept = intercept,
                estimate_sigma_int = estimate_sigma_int,
                standardize = standardize,
                standardize_mode = standardize_mode,
                covariance = covariance,
                ...)
}

## ---------- S3 helpers ----------
vcov.yorklm <- function(object, ...) object$cov.unscaled
coef.yorklm <- function(object, ...) object$coefficients
confint.yorklm <- function(object, parm = NULL, level = 0.95, ...) {
  b <- coef(object); V <- vcov(object)
  if (is.null(parm)) parm <- names(b)
  idx <- match(parm, names(b))
  se <- sqrt(pmax(diag(V)[idx], .Machine$double.eps))
  df <- object$df[2L]
  crit <- if (is.finite(df) && df > 0) qt(1 - (1 - level)/2, df) else qnorm(1 - (1 - level)/2)
  out <- cbind(lwr = b[idx] - crit*se, upr = b[idx] + crit*se)
  rownames(out) <- names(b)[idx]
  out
}

## ======================== end file ========================
