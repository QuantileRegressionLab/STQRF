# Penalized Iteratively Reweighted Least Squares solver for spatio-temporal quantile regression

pirls <- function(
  y_star,
  basis_combined,
  tau,
  penalty,
  ridge.penalty,
  theta = 1,
  tol = 10^-3,
  difo = 2
) {
  n = nrow(basis_combined)

  ### Step 1: get initial estimates
  residuo <- y_star

  obj.0 <- 0
  i <- 0
  criterio <- 1

  #################################################
  while (criterio > 0) {
    i <- i + 1
    print(paste("Iteration", i))
    #Weights
    w <- (tau - as.numeric(residuo < 0)) / (2 * residuo)
    wlimit <- 20
    w[w >= wlimit] <- wlimit

    #Calculate the new estimates
    alpha.new <- solve(
      wcrossprod(basis_combined, basis_combined, w) + penalty + ridge.penalty
    ) %*%
      wcrossprod(basis_combined, y_star, w)
    # alpha.new <- ginv(wcrossprod(basis_combined, basis_combined, w) + penalty + ridge.penalty)%*%wcrossprod(basis_combined, y_star, w)
    #Calculate difference between estimates based on residuals
    alpha <- alpha.new
    fits <- as.vector(basis_combined %*% alpha)
    residuo <- y_star - fits
    obj.1 <- sum(quantile_loss(epsilon = residuo, tau = tau))
    diferencia <- obj.1 - obj.0
    print(abs(diferencia))

    criterio <- as.numeric((abs(diferencia)) >= tol)
    obj.0 <- obj.1

    #Clean up to save memory
    gc()
  }

  #Calculate "hat" matrix
  # Maux0 <- t(sweep(basis_combined,MARGIN=1,w,"*"))
  # Maux1 <- tcrossprod(Maux0,t(basis_combined))+penalty + ridge.penalty
  # Maux2 <- solve(Maux1)
  # # Maux2 <- ginv(Maux1)
  # Hat <- as.matrix(crossprod(t(Maux2),tcrossprod(Maux0,t(basis_combined))))

  # Hat.diag <- rowSums(basis_combined * t(solve(wcrossprod(basis_combined, basis_combined, w) + penalty + ridge.penalty)%*%t(basis_combined)%*%diag(w)))
  # Hat.diag <- rowSums(basis_combined * t(ginv(wcrossprod(basis_combined, basis_combined, w) + penalty + ridge.penalty)%*%t(basis_combined)%*%diag(w)))
  Hat.diag <- rowSums(
    basis_combined *
      t(
        solve(
          wcrossprod(basis_combined, basis_combined, w) +
            penalty +
            ridge.penalty
        ) %*%
          t(basis_combined * w)
      )
  )

  #Calculate fit criteria: SIC, GCV
  obj <- sum(quantile_loss(epsilon = residuo, tau = tau))
  dof <- sum(Hat.diag) # sum(diag(Hat))
  SIC <- log(obj / n) + (log(n) * dof) / (2 * n)
  GCV <- mean(residuo^2) / (1 - theta * dof / n)^2
  GCVq <- obj / (n - theta * dof)
  crit = c(SIC, GCV, GCVq)
  names(crit) = c("SIC", "GCV", "GCVq")

  out = list()
  out$w = w
  out$alpha = alpha
  out$fits = fits
  out$crit = crit
  out$obj = obj
  out$dof = dof

  return(out)
}
