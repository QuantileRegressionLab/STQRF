# Server-optimized spatio-temporal model fitting script; includes simulation result post-processing

rm(list = ls())
gc()
graphics.off()

lib = getwd()
repos = "http://cran.uk.r-project.org"
.libPaths(c(.libPaths(), lib))
# install.packages(c("mboost"), lib = lib, repos = repos)

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
# library(car)
library(abind)
library(parallel)
library(doParallel)

source("../functions/mixed_model_function_time_space.R")
source("../functions/PIRLS_time_space.R")

# Create cluster
# ncores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")) # detectCores()
ncores <- parallel::detectCores() #aggiunto
ncores
cl <- makeCluster(ncores)
registerDoParallel(cl)
invisible(clusterEvalQ(
  cl = cl,
  c(
    .libPaths(getwd()),
    source("../functions/mixed_model_function_time_space.R"), #aggiunto
    source("../functions/PIRLS_time_space.R"), #aggiunto
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
# tau_vec <- c(0.1, 0.5, 0.9)
# tau_vec <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95) #nella prova con 250 simulazioni
tau_vec <- c(0.05, 0.5, 0.95) #aggiunta per tau estremi
scenario = "non_linear" # "linear" or "non_linear"

df <- c(5, 5) # B-spline degrees of freedom space and time
ntree = 300 # number of trees
max_iter = 300 # maximum number of iteration
tol = 10e-04 # tolerance of the algorithm
# B = 100 # 16 * 4 # number of Monte Carlo replicates
B = 250 # 16 * 4 # number of Monte Carlo replicates
#Lambda
lambda_values <- 1 / (n^seq(0.01, 0.5, length.out = 10)) # 1 / (n^seq(0.1, 1, length.out = 10)) # 1 / (n^seq(0.8, 10, length.out = 10))
# lambda_values <- seq(1e-9, 1e-5, length.out = 5)
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
      # y <- 3 + 2 * dat$x1 + 2 * dat$x2 + Alltrend + rnorm(n, mean = 0, sd = 1)
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

    # train_idx <- sample(1:n, floor(n * 0.7), replace = F)
    # train <- data.s[train_idx, ]
    # test <- data.s[-train_idx, ]

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
    # pred_qboost <- matrix(NA, nrow(test), length(tau_vec))

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

# boxplot(t(res.tau[, 1, ]), main = "Qloss per tau = 0.1")
# boxplot(t(res.tau[, 2, ]), main = "Qloss per tau = 0.5")
# boxplot(t(res.tau[, 3, ]), main = "Qloss per tau = 0.9")

#################################
#################################
#################################

ts_tag <- format(Sys.time(), "%Y%m%d")
fname <- sprintf(paste0("res_norm_", scenario, "%s_B", B, ".RData"), ts_tag)
save.image(fname)


stopCluster(cl) #aggiunto

### Post-processing and plots for non-linear scenario
load("res_norm_nonlin_20260112_B250.RData")
# library(gridExtra)
# library(grid)
# ts_tag <- format(Sys.time(), "%Y%m%d")
# name <- sprintf("Boxplot_Sim_B250.pdf", ts_tag)
# CairoPDF(name, width = 12, height = 10)
# for (t in 1:length(tau_vec)) {
#   boxplot(t(res.tau[, t, ]), main = paste("Qloss per tau =", tau_vec[t]))
# }
res.tau.nonlin <- res.tau
# # Nuova pagina per la tabella
# grid.newpage()
tab_nonlinear <- round(apply(res.tau, 1:2, mean), 4)
tab_nonlinear <- tab_nonlinear[, c(1, 4, 7)]
tab_nonlinear_sd <- round(apply(res.tau, 1:2, sd), 4)
tab_nonlinear_sd <- tab_nonlinear_sd[, c(1, 4, 7)]


# #Tabella per frequenze di lambda
# lambda_chosen <- info_sqrf["lambda.idx", , ]
# # Poiché lambda è lo stesso per tutti i tau (scelto a tau=0.5),
# # basta prendere una riga qualsiasi, es. tau=0.5
# lambda_vec <- lambda_chosen["0.5", ] # vettore di 250 valori

# freq_table <- table(lambda_vec)
# freq_rel <- round(prop.table(freq_table), 5)
# freq_perc <- round(prop.table(freq_table) * 100, 2) # in percentuale

# # tab iterazioni
# iter_mat <- info_sqrf["iter", , ] # matrice 7 x 250

# iter_summary <- t(apply(iter_mat, 1, function(x) {
#   c(
#     Mean = round(mean(x), 1),
#     # Median = round(median(x), 1),
#     Min = min(x),
#     Max = max(x)
#   )
# }))

# library(xtable)
# xtable(iter_summary[c(1, 4, 7), ])

### Post-processing and plots for linear scenario
load("res_norm_linear20260314_B250.RData")
library(gridExtra)
library(grid)
# ts_tag <- format(Sys.time(), "%Y%m%d")
# name <- sprintf(paste0("Boxplot_Sim_B250_", scenario, ".pdf"), ts_tag)
# # CairoPDF(name, width = 12, height = 10)
# for (t in 1:length(tau_vec)) {
#   boxplot(t(res.tau[, t, ]), main = paste("Qloss per tau =", tau_vec[t]))
# }
res.tau.lin <- res.tau
# # Nuova pagina per la tabella
# grid.newpage()
tab_linear <- round(apply(res.tau.lin, 1:2, mean), 4)
tab_linear_sd <- round(apply(res.tau.lin, 1:2, sd), 4)
# tab_linear <- tab_linear[, c(2, 3, 4)]
grid.table(tab)
# dev.off()
# xtable(cbind(tab_linear, tab_nonlinear), digits = 4)
# #Tabella per frequenze di lambda
lambda_chosen <- info_sqrf["lambda.idx", , ]
# # Poiché lambda è lo stesso per tutti i tau (scelto a tau=0.5),
# # basta prendere una riga qualsiasi, es. tau=0.5
lambda_vec <- lambda_chosen["0.5", ] # vettore di 250 valori

freq_table <- table(lambda_vec)
freq_rel <- round(prop.table(freq_table), 5)
# freq_perc <- round(prop.table(freq_table) * 100, 2) # in percentuale

# # tab iterazioni
iter_mat <- info_sqrf["iter", , ] # matrice 7 x 250

iter_summary <- t(apply(iter_mat, 1, function(x) {
  c(
    Mean = round(mean(x), 1),
    # Median = round(median(x), 1),
    Min = min(x),
    Max = max(x)
  )
}))

# library(xtable)
xtable(iter_summary[])

#### Boxplot with ggplot ##
library(ggplot2)
# per tau=0.05
df_0.05 <- df_0.5 <- df_0.95 <- matrix(
  data = NA,
  nrow = 2 * B,
  ncol = 6,
  dimnames = list(NULL, c(names(res.tau.lin[, 1, 1])))
)
for (b in 1:B) {
  df_0.05[b, ] <- res.tau.lin[, 1, b] #tau=0.05
  df_0.05[b + 250, ] <- res.tau.nonlin[, 1, b] #tau=0.05

  df_0.5[b, ] <- res.tau.lin[, 2, b] #tau=0.5
  df_0.5[b + 250, ] <- res.tau.nonlin[, 4, b] #tau=0.5

  df_0.95[b, ] <- res.tau.lin[, 3, b]
  df_0.95[b + 250, ] <- res.tau.nonlin[, 7, b] #tau=0.95
}
df_0.05 <- as.data.frame(df_0.05)
df_0.5 <- as.data.frame(df_0.5)
df_0.95 <- as.data.frame(df_0.95)
v1 <- rep("linear", 250)
v2 <- rep("non-linear", 250)
df_0.05$Scenario <- df_0.5$Scenario <- df_0.95$Scenario <- c(v1, v2)


df_0.05_long <- df_0.05 %>%
  pivot_longer(
    cols = c("QR", "STQR", "QGAM", "QBOOST", "QRF", "STQRF"),
    names_to = "Model",
    values_to = "Value"
  )

df_0.5_long <- df_0.5 %>%
  pivot_longer(
    cols = c("QR", "STQR", "QGAM", "QBOOST", "QRF", "STQRF"),
    names_to = "Model",
    values_to = "Value"
  )

df_0.95_long <- df_0.95 %>%
  pivot_longer(
    cols = c("QR", "STQR", "QGAM", "QBOOST", "QRF", "STQRF"),
    names_to = "Model",
    values_to = "Value"
  )


#Boxplot:
bp_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 15, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(1.0, 'cm'),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    plot.margin = margin(5, 10, 5, 5)
  )

gg_0.05 <- ggplot(df_0.05_long, aes(x = Model, y = Value, fill = Scenario)) +
  geom_boxplot(linewidth = 0.6, width = 0.6, outlier.size = 0.8) +
  bp_theme +
  labs(x = "Model at τ = 0.05")

gg_0.5 <- ggplot(df_0.5_long, aes(x = Model, y = Value, fill = Scenario)) +
  geom_boxplot(linewidth = 0.6, width = 0.6, outlier.size = 0.8) +
  bp_theme +
  labs(x = "Model at τ = 0.50")

gg_0.95 <- ggplot(df_0.95_long, aes(x = Model, y = Value, fill = Scenario)) +
  geom_boxplot(linewidth = 0.6, width = 0.6, outlier.size = 0.8) +
  bp_theme +
  labs(x = "Model at τ = 0.95")


library(ggpubr)
ts_tag <- format(Sys.time(), "%Y%m%d")
name <- sprintf("Boxplot_Sim_B250_%s.pdf", ts_tag)
# Dimensioni consigliate per A4 con margini standard:
# Cairo::CairoPDF(name, width = 13, height = 4.5)
Cairo::CairoPDF(name, width = 13, height = 5)
ggarrange(gg_0.05, gg_0.5, gg_0.95, nrow = 1, common.legend = TRUE)
dev.off()
