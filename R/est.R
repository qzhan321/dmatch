callmtCov_t <- function(X, nu){
  require(MASS)
  res <- cov.trob(y1, nu = nu)$cov
  return(res)
}


run_simple_est <- function(X, Y, x1, y1, steps = 500){

  n <- nrow(X)
  m <- nrow(Y)
  mu_x <-  apply(x1, 2, mean)
  mu_y <-  apply(y1, 2, mean)
  D1 <- cov(y1)
  D2 <- cov(x1)
  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  D1_inv <- solve(D1)
  D2_inv <- solve(D2)
  D12_inv_sum <- D1_inv + D2_inv

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- (beta1+beta2)/2

  for(i in 1:steps){

    A_t <- estA

    omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
    omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)

    C1 <- YY%*%A_t%*%D1_inv/m + XX%*%A_t%*%D2_inv/n
    C2 <- t(Y)%*%omega1%*%D1_inv/m + t(X)%*%omega2%*%D2_inv/n

    A_gradient <- C1 - C2

    errorX <- apply(X%*%estA, 1, function(z){z-mu_x})
    errorY <- apply(Y%*%estA, 1, function(z){z-mu_y})

    ave_errorX <- apply(errorX, 1, mean)
    ave_errorY <- apply(errorY, 1, mean)

    dmat_part <- ave_errorX%*%D2_inv + ave_errorY%*%D1_inv

    d_gradient <- estd%*%D12_inv_sum - dmat_part

    pre <- ceiling(log10(abs(A_gradient)))
    pre2 <- ceiling(log10(abs(d_gradient)))
    estA <- estA - 10^(-(pre+2))*A_gradient
    estd <- estd - 10^(-(pre2+2))*d_gradient
  }

  return(list(A = estA, d = estd))
}




rotation_matrix <- function(theta){
  a1 <- sin(theta*pi)
  a2 <- cos(theta*pi)
  return(matrix(c(a2, -a1, a1, a2), ncol = 2))
}


run_simple_angle_est <- function(X, Y, x1, y1, steps = 1000){

  n <- nrow(X)
  m <- nrow(Y)
  mu_x <-  apply(x1, 2, mean)
  mu_y <-  apply(y1, 2, mean)
  D1 <- cov(y1)
  D2 <- cov(x1)
  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  D1_inv <- solve(D1)
  D2_inv <- solve(D2)
  D12_inv_sum <- D1_inv + D2_inv

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  for(kk in 1:steps){

    A_t <- estA

    omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
    omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
    C1 <- YY%*%A_t%*%D1_inv + XX%*%A_t%*%D2_inv
    C2 <- t(Y)%*%omega1%*%D1_inv + t(X)%*%omega2%*%D2_inv
    A_gradient <- C1 - C2

    angles <- c(-1/180, -1/360, 1/360, 1/180)
    changes <- NULL

    for(i in 1:length(angles)){
      theta <- angles[i]
      A_t_prime <- rotation_matrix(theta)
      C1 <- YY%*%A_t_prime%*%D1_inv/m + XX%*%A_t_prime%*%D2_inv/n
      C2 <- t(Y)%*%omega1%*%D1_inv/m + t(X)%*%omega2%*%D2_inv/n
      changes[[i]] <- abs((C1 - C2)/A_gradient)
    }

    changes_median <- sapply(changes, median)
    aa <- which(changes_median<1)
    errorX <- apply(X%*%estA, 1, function(z){z-mu_x})
    errorY <- apply(Y%*%estA, 1, function(z){z-mu_y})

    ave_errorX <- apply(errorX, 1, mean)
    ave_errorY <- apply(errorY, 1, mean)

    dmat_part <- ave_errorX%*%D2_inv + ave_errorY%*%D1_inv

    d_gradient <- estd%*%D12_inv_sum - dmat_part

    if(length(aa)>=1){
      selected.angle <- angles[aa[1]]
      estA <- rotation_matrix(selected.angle)
    } else {
      estA <- matrix(c(1, 0, 0, 1), ncol = 2)
    }


    #  pre <- ceiling(log10(abs(A_gradient)))
    #  estA <- estA - 10^(-(pre+2))*A_gradient
    #  pre <- ceiling(log10(abs(d_gradient)))
    pre2 <- ceiling(log10(abs(d_gradient)))
    estd <- estd - 10^(-(pre2+2))*d_gradient

  }

  return(list(A = estA, d = estd))
}


run_simple_angle_likelihood_est <- function(X, Y, x1, y1, steps = 100){

  n <- nrow(X)
  m <- nrow(Y)
  mu_x <- apply(x1, 2, mean)
  mu_y <- apply(y1, 2, mean)
  D1 <- cov(y1)
  D2 <- cov(x1)
  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  D1_inv <- solve(D1)
  D2_inv <- solve(D2)
  D12_inv_sum <- D1_inv + D2_inv

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  for(i in 1:steps){

    A_t <- estA

    omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
    omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
    C1 <- YY%*%A_t%*%D1_inv/n + XX%*%A_t%*%D2_inv
    C2 <- t(Y)%*%omega1%*%D1_inv/n + t(X)%*%omega2%*%D2_inv
    A_gradient <- C1 - C2

    angles <- seq(-5/180, 5/180, by = 1/360)
    changes <- rep(0, length(angles))
    for(i in 1:length(angles)){
      theta <- angles[i]
      A_t_prime <- rotation_matrix(theta)
      changes[i] <- likelihood_para(X, Y, mu_x, mu_y, A_t_prime, estd, D1, D2)
    }

    errorX <- apply(X%*%estA, 1, function(z){z-mu_x})
    errorY <- apply(Y%*%estA, 1, function(z){z-mu_y})

    ave_errorX <- apply(errorX, 1, mean)
    ave_errorY <- apply(errorY, 1, mean)

    dmat_part <- ave_errorX%*%D2_inv + ave_errorY%*%D1_inv

    d_gradient <- estd%*%D12_inv_sum - dmat_part

    selected.angle <- angles[which.max(changes)]
    estA <- rotation_matrix(selected.angle)

    #  pre <- ceiling(log10(abs(A_gradient)))
    #  estA <- estA - 10^(-(pre+2))*A_gradient
    #  pre <- ceiling(log10(abs(d_gradient)))
    pre2 <- ceiling(log10(abs(d_gradient)))
    estd <- estd - 10^(-(pre2+2))*d_gradient

  }

  return(list(A = estA, d = estd))
}


run_simple_angle_t_distribution_est <- function(X, Y, x1, y1, steps = 20, nu = 5){

  require(MASS)

  nu_factor <- nu/(nu-2)
  n <- nrow(X)
  m <- nrow(Y)
  mu_x <- apply(x1, 2, mean)
  mu_y <- apply(y1, 2, mean)


  Sigma1 <- cov.trob(y1, nu = nu)$cov/nu_factor
  Sigma2 <- cov.trob(x1, nu = nu)$cov/nu_factor

  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  Sigma1_inv <- solve(Sigma1)
  Sigma2_inv <- solve(Sigma2)

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  for(tt in 1:steps){

    A_t <- estA
    A <- t(estA)
    omega1 <- mu_y + estd
    omega2 <- mu_x + estd

    A_gradient_part1 <- matrix(c(0, 0, 0, 0), ncol = 2)
    A_gradient_part2 <- matrix(c(0, 0, 0, 0), ncol = 2)

    d_gradient_part1 <- c(0, 0)
    d_gradient_part2 <- c(0, 0)

    for(j in 1:m){
      ay_d_u <- A%*%Y[j, ] - omega1
      S_ay_d_u <- Sigma1_inv%*%ay_d_u
      C_1u <- Y[j, ]%*%t(Y[j, ])%*%A_t%*%Sigma1_inv - Y[j, ]%*%t(omega1)%*%Sigma1_inv
      C_1d <- t(ay_d_u)%*%S_ay_d_u
      A_gradient_part1 <- A_gradient_part1 + C_1u/(1+C_1d[1])
      d_gradient_part1 <- d_gradient_part1 - S_ay_d_u/(1+C_1d[1])
    }

    for(i in 1:n){
      ax_d_u <- A%*%X[i, ] - omega2
      S_ax_d_u <- Sigma2_inv%*%ax_d_u
      C_2u <- X[i, ]%*%t(X[i, ])%*%A_t%*%Sigma2_inv - X[i, ]%*%t(omega2)%*%Sigma2_inv
      C_2d <- t(ax_d_u)%*%S_ax_d_u
      A_gradient_part2 <- A_gradient_part2 + C_2u/(1+C_2d[1])
      d_gradient_part2 <- d_gradient_part2 - S_ax_d_u/(1+C_2d[1])
    }

    A_gradient <- A_gradient_part1/m + A_gradient_part2/n
    d_gradient <- d_gradient_part1/m + d_gradient_part2/n

    angles <- seq(-5/180, 5/180, by = 1/360)
    changes <- NULL
    for(kk in 1:length(angles)){
      theta <- angles[kk]
      A_t_prime <- rotation_matrix(theta)
      A_prime <- t(A_t_prime)

      A_gradient_prime_part1 <- matrix(c(0, 0, 0, 0), ncol = 2)
      A_gradient_prime_part2 <- matrix(c(0, 0, 0, 0), ncol = 2)

      for(j in 1:m){
        ay_d_u <- A_prime%*%Y[j, ] - omega1
        S_ay_d_u <- Sigma1_inv%*%ay_d_u
        C_1u <- Y[j, ]%*%t(Y[j, ])%*%A_t_prime%*%Sigma1_inv - Y[j, ]%*%t(omega1)%*%Sigma1_inv
        C_1d <- t(ay_d_u)%*%S_ay_d_u
        A_gradient_prime_part1 <- C_1u/(1+C_1d[1])
      }

      for(i in 1:n){
        ax_d_u <- A_prime%*%X[i, ] - omega2
        S_ax_d_u <- Sigma2_inv%*%ax_d_u
        C_2u <- X[i, ]%*%t(X[i, ])%*%A_t_prime%*%Sigma2_inv - X[i, ]%*%t(omega2)%*%Sigma2_inv
        C_2d <- t(ax_d_u)%*%S_ax_d_u
        A_gradient_prime_part2 <- C_2u/(1+C_2d[1])
      }

      changes[[kk]]<- abs(A_gradient_prime_part1+A_gradient_prime_part2)/abs(A_gradient)
    }

    changes_median <- sapply(changes, median)
    aa <- which(changes_median<1)

    if(length(aa)>=1){
      selected.angle <- angles[aa[1]]
      estA <- rotation_matrix(selected.angle)
    } else {
      estA <- matrix(c(1, 0, 0, 1), ncol = 2)
    }

  }
  pre2 <- ceiling(log10(abs(d_gradient)))
  estd <- estd - 10^(-(pre2+2))*d_gradient
  return(list(A = estA, d = estd))

}


run_simple_angle_t_distribution_likelihood_est <- function(X, Y, x1, y1, steps = 20, nu = 5){

  require(MASS)
  nu_factor <- nu/(nu-2)
  n <- nrow(X)
  m <- nrow(Y)
  mu_x <- apply(x1, 2, mean)
  mu_y <- apply(y1, 2, mean)


  Sigma1 <- cov.trob(y1, nu = nu)$cov/nu_factor
  Sigma2 <- cov.trob(x1, nu = nu)$cov/nu_factor

  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  Sigma1_inv <- solve(Sigma1)
  Sigma2_inv <- solve(Sigma2)
  Sigma12_inv_sum <- Sigma1_inv + Sigma2_inv

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  for(tt in 1:steps){

    A_t <- estA
    A <- t(estA)
    omega1 <- mu_y + estd
    omega2 <- mu_x + estd

    A_gradient_part1 <- matrix(c(0, 0, 0, 0), ncol = 2)
    A_gradient_part2 <- matrix(c(0, 0, 0, 0), ncol = 2)

    d_gradient_part1 <- c(0, 0)
    d_gradient_part2 <- c(0, 0)

    for(j in 1:m){
      ay_d_u <- A%*%Y[j, ] - omega1
      S_ay_d_u <- Sigma1_inv%*%ay_d_u
      C_1u <- Y[j, ]%*%t(Y[j, ])%*%A_t%*%Sigma1_inv - Y[j, ]%*%t(omega1)%*%Sigma1_inv
      C_1d <- t(ay_d_u)%*%S_ay_d_u
      A_gradient_part1 <- A_gradient_part1 + C_1u/(1+C_1d[1])
      d_gradient_part1 <- d_gradient_part1 - S_ay_d_u/(1+C_1d[1])
    }

    for(i in 1:n){
      ax_d_u <- A%*%X[i, ] - omega2
      S_ax_d_u <- Sigma2_inv%*%ax_d_u
      C_2u <- X[i, ]%*%t(X[i, ])%*%A_t%*%Sigma2_inv - X[i, ]%*%t(omega2)%*%Sigma2_inv
      C_2d <- t(ax_d_u)%*%S_ax_d_u
      A_gradient_part2 <- A_gradient_part2 + C_2u/(1+C_2d[1])
      d_gradient_part2 <- d_gradient_part2 - S_ax_d_u/(1+C_2d[1])
    }

    A_gradient <- A_gradient_part1/m + A_gradient_part2/n
    d_gradient <- d_gradient_part1/m + d_gradient_part2/n

    angles <- seq(-5/180, 5/180, by = 1/360)
    changes <- rep(0, length(angles))
    for(kk in 1:length(angles)){
      theta <- angles[kk]
      A_t_prime <- rotation_matrix(theta)
      changes[kk] <- likelihood_t_para(X, Y, mu_x, mu_y, A_t_prime, estd, Sigma1, Sigma2, nu)
    }

    selected.angle <- angles[which.max(changes)]
    estA <- rotation_matrix(selected.angle)

    pre2 <- ceiling(log10(abs(d_gradient)))
    estd <- estd - 10^(-(pre2+2))*d_gradient

  }

  return(list(A = estA, d = estd))

}



run_simple_t_distribution_est <-  function(X, Y, x1, y1, steps = 1000, nu = 5){

  require(MASS)
  nu_factor <- nu/(nu-2)
  n <- nrow(X)
  m <- nrow(Y)
  mu_x <- apply(x1, 2, mean)
  mu_y <- apply(y1, 2, mean)


  Sigma1 <- cov.trob(y1, nu = nu)$cov/nu_factor
  Sigma2 <- cov.trob(x1, nu = nu)$cov/nu_factor

  p <- ncol(X)

  estd <- (apply(X, 2, mean) - mu_x)/2 + (apply(Y, 2, mean) - mu_y)/2

  XX <- t(X)%*%X
  YY <- t(Y)%*%Y

  Sigma1_inv <- solve(Sigma1)
  Sigma2_inv <- solve(Sigma2)
  Sigma12_inv_sum <- Sigma1_inv + Sigma2_inv

  omega1 <- matrix(rep(mu_y + estd, m), ncol=p, byrow=T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol=p, byrow=T)
  beta1 <- solve(YY)%*%t(Y)%*%omega1
  beta2 <- solve(XX)%*%t(X)%*%omega2

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  for(tt in 1:steps){

    A_t <- estA
    A <- t(estA)
    omega1 <- mu_y + estd
    omega2 <- mu_x + estd

    A_gradient_part1 <- matrix(c(0, 0, 0, 0), ncol = 2)
    A_gradient_part2 <- matrix(c(0, 0, 0, 0), ncol = 2)

    d_gradient_part1 <- c(0, 0)
    d_gradient_part2 <- c(0, 0)

    for(j in 1:m){
      ay_d_u <- A%*%Y[j, ] - omega1
      S_ay_d_u <- Sigma1_inv%*%ay_d_u
      C_1u <- Y[j, ]%*%t(Y[j, ])%*%A_t%*%Sigma1_inv - Y[j, ]%*%t(omega1)%*%Sigma1_inv
      C_1d <- t(ay_d_u)%*%S_ay_d_u
      A_gradient_part1 <- A_gradient_part1 + C_1u/(1+C_1d[1])
      d_gradient_part1 <- d_gradient_part1 - S_ay_d_u/(1+C_1d[1])
    }

    for(i in 1:n){
      ax_d_u <- A%*%X[i, ] - omega2
      S_ax_d_u <- Sigma2_inv%*%ax_d_u
      C_2u <- X[i, ]%*%t(X[i, ])%*%A_t%*%Sigma2_inv - X[i, ]%*%t(omega2)%*%Sigma2_inv
      C_2d <- t(ax_d_u)%*%S_ax_d_u
      A_gradient_part2 <- A_gradient_part2 + C_2u/(1+C_2d[1])
      d_gradient_part2 <- d_gradient_part2 - S_ax_d_u/(1+C_2d[1])
    }

    A_gradient <- A_gradient_part1/m + A_gradient_part2/n
    d_gradient <- d_gradient_part1/m + d_gradient_part2/n


    pre <- ceiling(log10(abs(A_gradient)))
    estA <- estA - 10^(-(pre+2))*A_gradient
    pre2 <- ceiling(log10(abs(d_gradient)))
    estd <- estd - 10^(-(pre2+2))*d_gradient

  }

  return(list(A = estA, d = estd))
}




