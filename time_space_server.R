# Server-optimized spatio-temporal model fitting script; includes simulation result post-processing

rm(list = ls())
gc()
graphics.off()

lib = getwd()
repos = "http://cran.uk.r-project.org"
.libPaths(c(.libPaths(), lib))

Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OMP_NUM_THREADS = 1)

library(splines)
library(quantreg)
library(quantregForest)
library(qgam)
library(mboost)
library(lqmm)
library(MASS)
library(abind)
library(parallel)
library(doParallel)

source("mixed_model_function_time_space.R")
source("PIRLS_time_space.R")

# Create cluster
# ncores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
ncores <- parallel::detectCores()
ncores
cl <- makeCluster(ncores)
registerDoParallel(cl)
invisible(clusterEvalQ(
  cl = cl,
  c(
    source("mixed_model_function_time_space.R"),
    source("PIRLS_time_space.R"),
    library(splines),
    library(quantreg),
    library(qgam),
    library(mboost),
    library(MASS),
    library(quantregForest),
    library(lqmm)
  )
))


# Data simulation
nspace = 50 # number of locations
ntime = 10 # number of years or days -> 365*nyr
n = nspace * ntime
tau_vec <- c(0.05, 0.5, 0.95)
scenario = "non_linear" # "linear" or "non_linear"

df <- c(5, 5) # B-spline degrees of freedom space and time
ntree = 300 # number of trees
max_iter = 300 # maximum number of iteration
tol = 10e-04 # tolerance of the algorithm
B = 250 # number of Monte Carlo replicates
#Lambda
lambda_values <- 1 / (n^seq(0.01, 0.5, length.out = 10))
lambda_grid <- matrix(lambda_values, nrow = length(lambda_values), ncol = 5)
colnames(lambda_grid) = paste("lambda", 1:5, sep = "")
lambda.ridge = 0.0001
models = c("QR", "STQR", "QGAM", "QBOOST", "QRF", "STQRF")


t0 = Sys.time()
res.out = foreach(b = 1:B) %dopar%
  {
    print(b)
    set.seed(b)

    longitude <- runif(nspace, 0, 1) # x-coordinates
    latitude <- runif(nspace, 0, 1) # y-coordinates

    dat <- do.call(
      rbind,
      lapply(1:ntime, function(t) {
        data.frame(
          time = t,
          longitude = longitude,
          latitude = latitude,
          x1 = rnorm(nspace), # covariate 1
          x2 = rnorm(nspace) # covariate 2
        )
      })
    )

    # Create a spatial effect (three examples)
    spatial_effect <- 0.2 *
      dat$longitude +
      0.25 * dat$latitude +
      3 * dat$longitude * dat$latitude # stable
    # spatial_effect = sin(10 * longitude) * cos(10 * latitude) + 0.5 * sin(5 * longitude * latitude) # wiggly
    # spatial_effect <- sin(longitude / 2) + cos(latitude / 3)

    #Create a temporal effect
    temporal_effect <- 1 + 0.3 * scale(dat$time) + 0.7 * scale(dat$time)^2

    # Draw the surface plot over a grid of 100 points
    #uv = seq(0, 1, length = 1e2)
    #z <- outer(uv, uv, function(longitude, latitude) spatial_effect <- 0.2 * longitude + 0.25 * latitude + 3 * longitude * latitude)
    #image(uv, uv, z, col = colorRampPalette(c("blue", "cyan", "yellow", "red"))(25))
    #contour(uv, uv, z,add=T,lty=4)

    # Create the outcome variable (signal + spatial + temporal + noise)
    Alltrend = as.vector(temporal_effect + spatial_effect)
    Alltrend.mat = matrix(Alltrend, ncol = ntime, nrow = nspace, byrow = T)

    if (scenario == "linear") {
      y <- 3 + 2 * dat$x1 + 0.5 * dat$x2 + Alltrend + rnorm(n, mean = 0, sd = 1) # rchisq(n, df = 3)/sqrt(6)
    } else if (scenario == "non_linear") {
      y <- 3 +
        2 * dat$x1 +
        0.5 * dat$x2^2 +
        dat$x1 * dat$x2 +
        Alltrend +
        rnorm(n, mean = 0, sd = 1) # rchisq(n, df = 3)/sqrt(6)
    }

    y <- as.vector(y)

    # Organize the data
    data.s <- data.frame(
      y = y,
      x1 = dat$x1,
      x2 = dat$x2,
      longitude = dat$longitude,
      latitude = dat$latitude,
      time = dat$time,
      location = rep(1:nspace, ntime)
    )

    train_idx <- sample(1:nspace, floor(nspace * 0.7), replace = F)
    train <- data.s[data.s$location %in% train_idx, ]
    test <- data.s[!(data.s$location %in% train_idx), ]

    # MODEL COMPARISON
    # Quantile regression model
    qr_model <- rq(y ~ x1 + x2, data = train, tau = tau_vec)
    pred_qr <- predict(qr_model, newdata = test[c("x1", "x2")])

    # Spatial temporal quantile regression model
    sqr_model <- rq(
      y ~ x1 +
        x2 +
        bs(longitude, df = df[1], Boundary.knots = c(0, 1)) +
        bs(latitude, df = df[1], Boundary.knots = c(0, 1)) +
        bs(time, df = df[2]),
      data = train,
      tau = tau_vec
    )
    pred_sqr <- predict(
      sqr_model,
      newdata = test[c("x1", "x2", "longitude", "latitude", "time")]
    )

    # QGAM cubic regression spline with shrinkage (tensor product)
    qgam_model <- mqgam(
      y ~ x1 +
        x2 +
        ti(longitude, latitude, bs = "cr", k = df[1]) +
        s(time, bs = "cr", k = df[2]),
      data = train,
      qu = tau_vec
    )
    pred_qgam <- qdo(
      qgam_model,
      qu = tau_vec,
      fun = predict,
      newdata = test[c("x1", "x2", "longitude", "latitude", "time")]
    )
    names(pred_qgam) <- tau_vec
    pred_qgam <- do.call(cbind, pred_qgam)

    # Qboost regression model
    pred_qboost <- matrix(NA, nrow(test), length(tau_vec))
    qboost_model <- gamboost(
      y ~ bols(x1, intercept = F) +
        bols(x2, intercept = F) +
        bspatial(longitude, latitude, df = df[1], boundary.knots = c(0, 1)) +
        bbs(time, df = df[2]),
      data = train,
      family = QuantReg(tau = 0.5),
      control = boost_control(mstop = 300, nu = 0.4)
    )
    cv05f <- cv(model.weights(qboost_model), type = "kfold", B = 5)
    cvmsmo <- cvrisk(qboost_model, folds = cv05f)

    for (t in seq_along(tau_vec)) {
      tau <- tau_vec[t]
      qboost_model <- gamboost(
        y ~ bols(x1, intercept = F) +
          bols(x2, intercept = F) +
          bspatial(longitude, latitude, df = df[1], boundary.knots = c(0, 1)) +
          bbs(time, df = df[2]),
        data = train,
        family = QuantReg(tau = tau),
        control = boost_control(mstop = mstop(cvmsmo), nu = 0.4)
      )
      pred_qboost[, t] <- predict(
        qboost_model[mstop(cvmsmo), ],
        newdata = test[c("x1", "x2", "longitude", "latitude", "time")]
      )
    }

    # Quantile Random Forest
    qrf_model <- quantregForest(
      x = train[c("x1", "x2")],
      y = train$y,
      ntree = ntree,
      keep.inbag = T,
      nthreads = 1
    )
    pred_qrf <- predict(
      qrf_model,
      what = tau_vec,
      newdata = test[c("x1", "x2", "latitude", "longitude", "time")]
    )

    # Spatio-temporal quantile Random Forest
    info_sqrf <- matrix(NA, 3, length(tau_vec))
    rownames(info_sqrf) = c("lambda.idx", "iter", "timetot")
    colnames(info_sqrf) = tau_vec
    # Loop to choose the optimal lambda
    tau_lambda = 0.5
    criteria_results <- sapply(1:nrow(lambda_grid), function(l) {
      QuantileRandomForest_BSplines(
        train,
        tau_lambda,
        df = df,
        lambda = lambda_grid[l, ],
        lambda.ridge = lambda.ridge,
        ntree = 300,
        max_iter = max_iter,
        tol = tol,
        seed = b,
        nthreads = 1
      )$crit
    })
    criteria_results <- as.data.frame(t(criteria_results))
    # Select best lambda
    lambda <- as.numeric(lambda_grid[which.min(criteria_results$SIC), ])
    info_sqrf[1, ] <- rep(which.min(criteria_results$SIC), length(tau_vec))

    # Initialize list of predictions for each tau
    pred_sqrf <- matrix(NA, nrow(test), length(tau_vec))
    mixed_model <- vector("list", length(tau_vec))

    # Loop for each tau
    for (t in seq_along(tau_vec)) {
      tau <- tau_vec[t]

      mixed_model[[t]] <- QuantileRandomForest_BSplines(
        train,
        tau,
        df = df,
        lambda = lambda,
        lambda.ridge = lambda.ridge,
        ntree = ntree,
        max_iter = max_iter,
        tol = tol,
        seed = b,
        nthreads = 1
      )
      alpha = mixed_model[[t]][["alpha"]]
      basis_spatial = box.prod(test[c("longitude", "latitude")], k = df[1])
      basis_time <- bbase(test$time, nseg = df[2] - 3)
      basis_interaction <- box.prod.inter(
        basis_time,
        basis_spatial,
        df[2],
        df[1]
      )
      basis_combined <- cbind(basis_spatial, basis_time, basis_interaction)
      pred_sqrf[, t] = c(
        predict(mixed_model[[t]][["qrf_model"]], what = tau, newdata = test) +
          basis_combined %*% alpha
      )
      info_sqrf[2, t] = mixed_model[[t]][["iter"]]
      info_sqrf[3, t] = mixed_model[[t]][["timetot"]]
    }

    # Results
    res.tau <- matrix(NA, nrow = length(models), ncol = length(tau_vec))
    rownames(res.tau) = models
    colnames(res.tau) = tau_vec
    for (t in seq_along(tau_vec)) {
      tau <- tau_vec[t]
      res.tau["QR", t] <- qloss(test$y, pred_qr[, t], tau)
      res.tau["STQR", t] <- qloss(test$y, pred_sqr[, t], tau)
      res.tau["QGAM", t] <- qloss(test$y, pred_qgam[, t], tau)
      res.tau["QBOOST", t] <- qloss(test$y, pred_qboost[, t], tau)
      res.tau["QRF", t] <- qloss(test$y, pred_qrf[, t], tau)
      res.tau["STQRF", t] <- qloss(test$y, pred_sqrf[, t], tau)
    }

    out = list()
    out$res.tau = res.tau
    out$info_sqrf = info_sqrf

    return(out)
  }
t1 = Sys.time() - t0

res.tau = sapply(1:B, function(b) res.out[[b]]$res.tau, simplify = "array")
info_sqrf = sapply(1:B, function(b) res.out[[b]]$info_sqrf, simplify = "array")

apply(res.tau, 1:2, mean)
apply(res.tau, 1:2, median)

apply(info_sqrf, 1:2, mean)
apply(info_sqrf, 1:2, median)

boxplot(t(res.tau[, 1, ]), pch = 19, main = "Qloss per tau = 0.05")
boxplot(t(res.tau[, 2, ]), pch = 19, main = "Qloss per tau = 0.50")
boxplot(t(res.tau[, 3, ]), pch = 19, main = "Qloss per tau = 0.95")

#################################
#################################
#################################

ts_tag <- format(Sys.time(), "%Y%m%d")
fname <- sprintf(paste0("res_norm_", scenario, "%s_B", B, ".RData"), ts_tag)
save.image(fname)


stopCluster(cl)
