#' @export
make_matern_prec_matrix_2d_M <- function(a, dim, nu) {
  Q <- Matrix::bandSparse(
    n = dim^2,
    k = c(-1 - dim, -1, 0, 1, 1 + dim),
    diagonals = list(
      rep(-1, dim^2),
      rep(-1, dim^2),
      rep(4 + a, dim^2),
      rep(-1, dim^2),
      rep(-1, dim^2)
    )
  )
  
  for (i in seq_len(nu)) {
    Q <- Q %*% Q
  }
  
  Q 
}

#' @export
cholesky_bandsparse_M <- function(Q) {
  n <- nrow(Q)
  bandwidth <- sqrt(n) + 1
  # Create a matrix to store the result
  L <- Matrix::Matrix(0, n, n)

  for (i in 1:n) {
    for (j in max(1, i - bandwidth):i) {
      if (i == j) {
        # Diagonal elements
        L[i, j] <- sqrt(Q[i, i] - sum(L[i, max(1, i - bandwidth):(j - 1)]^2))
      } else {
        # Off-diagonal elements
        L[i, j] <- (Q[i, j] - sum(L[i, max(1, i - bandwidth):(j - 1)] * L[j, max(1, i - bandwidth):(j - 1)])) / L[j, j]
      }
    }
  }

  return(L)
}


#' @export
compute_marginal_variances_M <- function(L) {
  n <- nrow(L)
  bandwidth <- sqrt(n) + 1
  Sigma <- Matrix::Matrix(0, n, n)

  for (i in n:1) {
    I_i <- seq(i, min(n, i + bandwidth))
    I_i <- I_i[I_i != i]
    for (j in c(I_i, i)) {
      value <- (i == j) / L[i, i]^2
      value <- value - 1 / L[i, i] * (L[i, I_i] %*% Sigma[I_i, j])
      Sigma[i, j] <- Sigma[j, i] <- as.numeric(value)
    }
  }

  sigmas <- Matrix::diag(Sigma)
  

  return(Matrix::Diagonal(n, sigmas))
}


#' @export
make_standardized_matern_M <- function(a, dim, nu) {
  Q <- make_matern_prec_matrix_2d_M(a, dim, nu)
  L <- cholesky_bandsparse_M(Q) |> Matrix::t()
  sigmas <- compute_marginal_variances_M(L) |> sqrt()


  sigmas %*% Q %*% sigmas
}
