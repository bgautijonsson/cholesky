#' @export
make_AR_prec_matrix <- function(dim, rho = 0.5) {

  scaling <- 1 / (1 - rho^2)
  off_diags <- -rho * scaling
  diag <- (1 + rho^2) * scaling

  Q <- Matrix::bandSparse(
    n = dim,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(off_diags, dim),
      c(scaling, rep(diag, dim - 2), scaling),
      rep(off_diags, dim)
    )
  )
  Q |> as.matrix()
}

#' @export
make_matern_prec_matrix <- function(dim, rho) {
  Q <- make_AR_prec_matrix(dim, rho)

  I <- Matrix::Diagonal(n = dim, x = 1)

  out <- Matrix::kronecker(Q, I) + Matrix::kronecker(I, Q)
  
  out |> as.matrix()
}

#' @export
cholesky_bandsparse <- function(Q) {
  n <- nrow(Q)
  bandwidth <- sqrt(n)
  # Create a matrix to store the result
  L <- matrix(0, n, n)
  
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
compute_marginal_variances <- function(L) {
  n <- nrow(L)
  bandwidth <- sqrt(n)
  Sigma <- matrix(0, bandwidth + 1, bandwidth + 1)
  out <- numeric(n)
  
  for (i in n:1) {
    
    I_i <- seq(i, min(n, i + bandwidth)) |> rev()
    I_i <- I_i[I_i != i]
    index <- rev(seq_along(I_i)) + 1
    
    
    Sigma[1, index] <- Sigma[index, 1] <- -1 / L[i, i] * L[I_i, i] %*% Sigma[index, index]
    Sigma[1, 1] <- 1 / L[i, i]^2 - 1 / L[i, i] * (L[I_i, i] %*% Sigma[index, 1])
    
    out[i] <- Sigma[1, 1]
    
    end_seq <- seq(2, bandwidth + 1)
    start_seq <- seq(1, bandwidth)
    Sigma[end_seq, end_seq] <- Sigma[start_seq, start_seq]
  }
  
  
  return(out)
}



#' @export
make_standardized_matern <- function(dim, rho) {
  Q <- make_matern_prec_matrix(dim, rho)
  L <- cholesky_bandsparse(Q) 
  sigmas <- compute_marginal_variances(L) |> sqrt() |> diag()
  
  
  sigmas %*% Q %*% sigmas
}
