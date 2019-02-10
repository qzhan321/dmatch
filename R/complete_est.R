

likelihood_para <- function(X, Y, mu_x, mu_y, A, d, D1, D2){

  likelihood <- 0
  moveX <- t(apply(X%*%A, 1, function(z){z - d}))
  moveY <- t(apply(Y%*%A, 1, function(z){z - d}))

  for(i in 1:nrow(moveX)){
    likelihood <- likelihood + dmvnorm(moveX[i, ], mu_x, D2, log = TRUE)
  }

  for(i in 1:nrow(moveY)){
    likelihood <- likelihood + dmvnorm(moveY[i, ], mu_y, D1, log = TRUE)
  }

  return(likelihood)

}

likelihood_t_para <- function(X, Y, mu_x, mu_y, A, d, Sigma1, Sigma2, df){
  library(mvtnorm)
  likelihood <- 0
  moveX <- t(apply(X%*%A, 1, function(z){z - d}))
  moveY <- t(apply(Y%*%A, 1, function(z){z - d}))

  for(i in 1:nrow(moveX)){
    likelihood <- likelihood + dmvt(moveX[i, ], mu_x, Sigma2, df, log = TRUE)
  }

  for(i in 1:nrow(moveY)){
    likelihood <- likelihood + dmvt(moveY[i, ], mu_y, Sigma1, df, log = TRUE)
  }

  return(likelihood)

}



select_stepsize_for_A_one <- function(Y, X, A, rowID, colID, gra_A, ll, gamma, mu_x, mu_y, d, D1, D2, down){

  A_one = A[rowID, colID]
  gra_one = gra_A[rowID, colID]
  gra_one2 = gra_one*gra_one*gamma
  start = sqrt(abs(A_one/gra_one))/2

  aa = start
  selected = A_one
  while(aa > 0){
    aa2 = aa*aa
    A_one_prime = A_one - aa2*gra_one
    A[rowID, colID] = A_one_prime
    lA_one_prime = likelihood_para(X, Y, mu_x, mu_y, A, d, D1, D2)
    if(lA_one_prime - ll - aa2*gra_one2 > 0) {
      selected = A_one_prime
      break
    }
    aa = aa - start*down
  }
  return(selected)

}


select_stepsize_for_d_one <- function(Y, X, d, ID, gra_d, ll, gamma, A, mu_x, mu_y, D1, D2, down){

  d_one = d[ID]
  gra_one = gra_d[ID]
  gra_one2 = gra_one*gra_one*gamma
  start = sqrt(abs(d_one/gra_one))/2

  if(start > 0.1) start = 0.1
  aa = start
  selected = d_one
  while(aa > 0){
    aa2 = aa*aa
    d_one_prime = d_one - aa2*gra_one
    d[ID] = d_one_prime
    ld_one_prime = likelihood_para(X, Y, mu_x, mu_y, A, d, D1, D2)
    if(ld_one_prime - ll + aa2*gra_one2 < 0 ) {
      selected = d_one_prime
      break
    }
    aa = aa - start*down
  }
  return(selected)

}



gradient_d_A <- function(X, XX, Y, YY, estA, estd, D1_inv, D2_inv, D12_inv_sum, mu_x, mu_y, m, n, p){

  omega1 <- matrix(rep(mu_y + estd, m), ncol = p, byrow = T)
  omega2 <- matrix(rep(mu_x + estd, n), ncol = p, byrow = T)

  C1 <- YY%*%estA%*%D1_inv/m + XX%*%estA%*%D2_inv/n
  C2 <- t(Y)%*%omega1%*%D1_inv/m + t(X)%*%omega2%*%D2_inv/n

  A_gradient <- C1 - C2
  errorX <- apply(X%*%estA, 1, function(z){z - mu_x})
  errorY <- apply(Y%*%estA, 1, function(z){z - mu_y})

  ave_errorX <- apply(errorX, 1, mean)
  ave_errorY <- apply(errorY, 1, mean)

  dmat_part <- ave_errorX%*%D2_inv + ave_errorY%*%D1_inv

  d_gradient <- - dmat_part + estd%*%D12_inv_sum

  return(list(A_gradient = A_gradient, d_gradient = d_gradient))

}


gamma = 0.75
steps = 50
down = 0.1

gradient_descent_d_A <- function(X, XX, Y, YY, estA, estd, D1, D2, D1_inv, D2_inv, D12_inv_sum, mu_x, mu_y, m, n, p, gamma, down, steps){


  gradient = gradient_d_A(X, XX, Y, YY, estA, estd, D1_inv, D2_inv, D12_inv_sum, mu_x, mu_y, m, n, p)
  gra_A = gradient$A_gradient
  gra_d = gradient$d_gradient
  ll = likelihood_para(X, Y, mu_x, mu_y, estA, estd, D1, D2)

  A_prime = matrix(0, p, p)
  d_prime = rep(0, p)

  for(i in 1:steps){

    for(k in 1:p){
      for(j in 1:p){
        if(abs(gra_A[k, j]) >= 0.001){
          A_prime[k, j] = select_stepsize_for_A_one(Y, X, estA, k, j, gra_A, ll, gamma, mu_x, mu_y, estd, D1, D2, down)
        } else {
          A_prime[k, j] = estA[k, j]
        }
      }
    }

    for(k in 1:p){
      if(abs(gra_d[k]) >= 0.001){
        d_prime[k] = select_stepsize_for_d_one(Y, X, estd, k, gra_d, ll, gamma, estA, mu_x, mu_y, D1, D2, down)
      } else {
        d_prime[k] = estd[k]
      }
    }

    estA = A_prime
    estd = d_prime

    gradient = gradient_d_A(X, XX, Y, YY, estA, estd, D1_inv, D2_inv, D12_inv_sum, mu_x, mu_y, m, n, p)
    gra_A = gradient$A_gradient
    gra_d = gradient$d_gradient
    ll = likelihood_para(X, Y, mu_x, mu_y, estA, estd, D1, D2)
    print(ll)

  }

  return(list(A = estA, d = estd, ll = ll))
}


run_est <- function(X, Y, x1, y1, gamma = 0.75, down = 0.1, steps = 25){

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

  estA <- matrix(c(1, 0, 0, 1), ncol = 2)

  res <- gradient_descent_d_A(X, XX, Y, YY, estA, estd, D1, D2, D1_inv, D2_inv, D12_inv_sum, mu_x, mu_y, m, n, p, gamma, down, steps)

  return(res)

}
