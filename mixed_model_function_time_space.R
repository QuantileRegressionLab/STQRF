# B-spline basis construction for spatio-temporal quantile regression (space, time, interaction terms)

acomb <- function(...) abind(..., along = 3)

wcrossprod <- function(x, y, w) {
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  if (!missing(y)) {
    if (is.vector(y)) {
      y <- as.matrix(y)
    }
    if (nrow(x) != nrow(y)) stop("x and y not conformable")
  }
  if (missing(w)) {
    if (missing(y)) return(crossprod(x)) else return(crossprod(x, y))
  } else if (
    length(w) == 1 || (is.vector(w) && sd(w) < sqrt(.Machine$double.eps))
  ) {
    if (missing(y)) {
      return(w[1] * crossprod(x))
    } else {
      return(w[1] * crossprod(x, y))
    }
  } else {
    if (is.vector(w)) {
      if (length(w) != nrow(x)) {
        stop("w is the wrong length")
      }
      if (missing(y)) {
        return(crossprod(x, w * x))
      } else {
        return(crossprod(x, w * y))
      }
    } else {
      if (nrow(w) != ncol(w) || nrow(w) != nrow(x)) {
        stop("w is the wrong dimension")
      }
      if (missing(y)) {
        return(crossprod(x, w %*% x))
      } else {
        return(crossprod(x, w %*% y))
      }
    }
  }
}

# Quantile loss
quantile_loss <- function(epsilon, tau) {
  epsilon * (tau - (epsilon < 0))
}

# EVALUATION
qloss <- function(y_true, y_pred, tau) {
  residuals <- y_true - y_pred
  loss <- residuals * (tau - (residuals < 0))
  return(mean(loss))
}

# ALD likelihood
ald_likelihood <- function(sigma, y, predicted_qrf, predicted_space_time, tau) {
  ald_density <- dal(
    y,
    mu = predicted_qrf + predicted_space_time,
    sigma,
    tau,
    log = T
  )
  # ald_density[is.na(ald_density)]=0
  likelihood <- sum(ald_density)
  return(likelihood)
}

# power function
tpower <- function(x, t, p) {
  (x - t)^p * (x >= t)
}

# Construct B-spline basis
bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  dx <- (xr - xl) / nseg
  knot <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knot, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)

  # Make B-splines exactly zero beyond their end knots
  nb <- ncol(B)
  sk <- knot[(1:nb) + deg + 1]
  Mask <- outer(x, sk, '<')
  B <- B * Mask
  return(B)
}

# Box product
box.prod <- function(coord, k) {
  x1.basis <- bbase(coord[, 1], nseg = k - 3)
  x2.basis <- bbase(coord[, 2], nseg = k - 3)
  B1. <- kronecker(x1.basis, t(rep(1, k)))
  B2. <- kronecker(t(rep(1, k)), x2.basis)
  Q <- B1. * B2.
  return(Q)
}

# Box product interaction space-time
box.prod.inter <- function(x1.basis, x2.basis, k1, k2) {
  B1. <- kronecker(x1.basis, t(rep(1, k2^2)))
  B2. <- kronecker(t(rep(1, k1)), x2.basis)
  Q <- B1. * B2.
  return(Q)
}


##################

QuantileRandomForest_BSplines <- function(
  data,
  tau,
  df = 5,
  lambda,
  lambda.ridge,
  theta,
  difo,
  ntree = 1e3,
  tol = 10e-04,
  max_iter = 1000,
  seed = NULL,
  nthreads = nthreads
) {
  if (missing(theta)) {
    theta = 1
  }
  if (missing(difo)) {
    difo = 2
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  y = c(data[[which(names(data) == "y")]])
  x = data[, which(
    !names(data) %in%
      c(
        "y",
        "longitude",
        "latitude",
        "state",
        "geometry",
        "time",
        "location",
        "date",
        "polygon_id"
      )
  )] #ho aggiunto date e polygon_id !!!!
  x = as.data.frame(x)
  longitude = data$longitude
  latitude = data$latitude
  time = data$time
  #Space
  basis_spatial <- box.prod(cbind(longitude, latitude), k = df[1])
  #Time
  basis_time <- bbase(time, nseg = df[2] - 3)
  #Interaction space-time
  basis_interaction <- box.prod.inter(basis_time, basis_spatial, df[2], df[1])

  #Create model matrix
  basis_combined <- as.matrix(cbind(
    basis_spatial,
    basis_time,
    basis_interaction
  ))

  ###Create penalty submatrices
  # Create penalty
  #Space
  D.mat <- diff(diag(df[1]), differences = difo)
  I1.mat <- diag(df[1])
  I2.mat <- diag(df[1])
  P1.mat <- kronecker(D.mat, I2.mat)
  P2.mat <- kronecker(I1.mat, D.mat)
  S1.mat <- crossprod(P1.mat)
  S2.mat <- crossprod(P2.mat)
  penalty.space <- lambda[1] * S1.mat + lambda[2] * S2.mat
  #Time
  P.time <- diff(diag(df[2]), differences = difo)
  S.time <- t(P.time) %*% P.time
  penalty.time <- S.time * lambda[3]
  #Interaction space year
  D1.inter <- diff(diag(df[2]), differences = difo)
  D2.inter <- diff(diag(df[1]^2), differences = difo)
  I1.inter <- diag(df[2])
  I2.inter <- diag(df[1]^2)
  P1.inter <- kronecker(D1.inter, I2.inter)
  P2.inter <- kronecker(I1.inter, D2.inter)
  S1.inter <- crossprod(P1.inter)
  S2.inter <- crossprod(P2.inter)
  penalty.inter <- lambda[4] * S1.inter + lambda[5] * S2.inter

  ##Combine penalty matrices together
  penalty <- matrix(
    0,
    nrow = df[1]^2 + df[2] + (df[1]^2 * df[2]),
    ncol = df[1]^2 + df[2] + (df[1]^2 * df[2])
  )
  penalty[1:df[1]^2, 1:df[1]^2] <- penalty.space
  penalty[
    (df[1]^2 + 1):(df[1]^2 + df[2]),
    (df[1]^2 + 1):(df[1]^2 + df[2])
  ] <- penalty.time
  penalty[
    (df[1]^2 + df[2] + 1):(df[1]^2 + df[2] + (df[1]^2 * df[2])),
    (df[1]^2 + df[2] + 1):(df[1]^2 + df[2] + (df[1]^2 * df[2]))
  ] <- penalty.inter
  penalty <- as.matrix(penalty)

  #Include the ridge penalty
  ridge.penalty <- matrix(
    0,
    nrow = df[1]^2 + df[2] + (df[1]^2 * df[2]),
    ncol = df[1]^2 + df[2] + (df[1]^2 * df[2])
  )
  ridge.penalty[1:df[1]^2, 1:df[1]^2] <- diag(colSums(basis_spatial))
  ridge.penalty[
    (df[1]^2 + 1):(df[1]^2 + df[2]),
    (df[1]^2 + 1):(df[1]^2 + df[2])
  ] <- diag(colSums(basis_time))
  ridge.penalty[
    (df[1]^2 + df[2] + 1):(df[1]^2 + df[2] + (df[1]^2 * df[2])),
    (df[1]^2 + df[2] + 1):(df[1]^2 + df[2] + (df[1]^2 * df[2]))
  ] <- diag(colSums(basis_interaction))
  ridge.penalty <- as.matrix(lambda.ridge * ridge.penalty)

  n = length(y)

  # Set initial estimates
  # quant_model <- lm.fit(x=cbind(1,x), y=y)
  # predicted_qrf <- quant_model$fitted.values

  #### andare in parallelo in caso sulla rf
  qrf_model <- quantregForest(
    x = x,
    y = y,
    ntree = ntree,
    keep.inbag = T,
    nthreads = nthreads
  )
  predicted_qrf <- predict(qrf_model, what = tau)
  # predicted_qrf <- rep(0, n)
  alpha = rep(0, ncol(basis_combined))
  y_star = y - (predicted_qrf + basis_combined %*% alpha)
  # %%

  sigma = mean(quantile_loss(y_star, tau = tau), na.rm = T)

  diff <- Inf
  iter <- 0
  llkold <- -10^250
  start_time = Sys.time()

  while (diff > tol & iter < max_iter) {
    y_star <- y - predicted_qrf

    # space
    # space_model<-rq(y_star ~ basis_spatial - 1, tau = tau, method = "fn")
    # alpha.new=space_model$coefficients
    # space_time_model<-pirls(y_star = y_star, basis_combined = basis_combined, tau = tau, penalty = penalty, ridge.penalty=ridge.penalty,
    #                         verbose = FALSE)
    space_time_model <- pirls(
      y_star = y_star,
      basis_combined = basis_combined,
      tau = tau,
      penalty = penalty,
      ridge.penalty = ridge.penalty,
      # tol = 10e-03
      tol = 10e-05
    )
    alpha.new = space_time_model$alpha
    predicted_space_time = basis_combined %*% alpha.new

    # random forest
    y_star = as.vector(y - predicted_space_time)
    qrf_model <- quantregForest(
      x = x,
      y = y_star,
      ntree = ntree,
      keep.inbag = T,
      nthreads = nthreads
    )
    # qrf_model <- quantregForest(x=x, y=y_star, ntree = ntree, keep.inbag = T, nthreads = parallel::detectCores() - 1)
    predicted_qrf <- predict(qrf_model, what = tau) #T for out of bag predictions, gives you information about which information is in and out of bag

    # scale
    sigma.new = mean(
      quantile_loss(y - predicted_qrf - predicted_space_time, tau = tau),
      na.rm = T
    )

    llk <- ald_likelihood(
      sigma = sigma.new,
      y = y,
      predicted_qrf = predicted_qrf,
      predicted_space_time = predicted_space_time,
      tau = tau
    )

    diff <- abs(llk - llkold) / abs(llkold)

    cat(
      "\niteration:",
      iter,
      "\nLog likelihood:",
      llk,
      "\nLikelihood difference:",
      diff
    )

    alpha = alpha.new
    sigma = sigma.new
    llkold <- llk
    iter = iter + 1
  }

  end_time = Sys.time()
  timetot = end_time - start_time

  output = list()
  output$diff = diff
  output$iter = iter
  output$timetot = timetot
  output$alpha = alpha
  output$lambda = lambda
  output$lambda.ridge = lambda.ridge
  output$qrf_model = qrf_model
  output$space_model = space_time_model
  output$predicted_qrf = predicted_qrf
  output$predicted_space = predicted_space_time # basis_combined%*%alpha.new
  output$sigma = sigma
  output$crit = space_time_model$crit

  return(output)
}

surfaceplotSpace <- function(
  longitude,
  latitude,
  nu = 1e2,
  nv = 1e2,
  k,
  modelcoef,
  main = "",
  cx.axis = 1.5,
  cx.lab = 1.5,
  cx.main = 1.5
) {
  coord = cbind(longitude, latitude)
  u <- seq(min(coord[, 1]), max(coord[, 1]), length = nu)
  v <- seq(min(coord[, 2]), max(coord[, 2]), length = nv)
  U. <- outer(rep(1, nv), u)
  V. <- outer(v, rep(1, nu))
  U <- as.vector(U.)
  V <- as.vector(V.)
  Bgrid <- box.prod(cbind(U, V), k)

  zgrid <- Bgrid %*% modelcoef
  Fitgrid <- matrix(as.vector(zgrid), nu, nv, byrow = T)
  d2lim = range(Fitgrid)

  # col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
  col.br <- rev(RColorBrewer::brewer.pal(10, "Spectral"))
  # col.br <- RColorBrewer::brewer.pal(10, "Spectral")
  ### Plot the surface in 2D with contour lines
  # image(U.[1,  ], V.[, 1], Fitgrid, xlab = "longitude", ylab ="latitude",
  #       zlim=d2lim,col=col.br(25),main=main,cex.axis=cx.axis,cex.lab=cx.lab,cex.main=cx.main)
  image(
    U.[1, ],
    V.[, 1],
    Fitgrid,
    xlab = "longitude",
    ylab = "latitude",
    zlim = d2lim,
    col = col.br,
    main = main,
    cex.axis = cx.axis,
    cex.lab = cx.lab,
    cex.main = cx.main
  )
  contour(
    U.[1, ],
    V.[, 1],
    Fitgrid,
    add = T,
    zlim = d2lim,
    lty = 4,
    # nlevels = 20
    nlevels = 50
  )
}
