library(cholesky)
library(Matrix)


f1a <- function(dim) {
  make_matern_prec_matrix_2d(dim = dim) |> cholesky_bandsparse(dim + 1)
  return(1)
}

f1b <- function(dim) {
  make_matern_prec_matrix_2d(dim = dim) |> cholesky_bandsparse_Matrix()
  return(1)
}

f2a <- function(dim) {
  make_matern_prec_matrix_2d(dim = dim) |> as.matrix() |> cholesky_bandsparse_cpp(dim + 1)
  return(1)
}

f2b <- function(dim) {
  make_matern_prec_matrix_2d(dim = dim) |> cholesky_bandsparse_cpp_eigen()
  return(1)
}


f3 <- function(dim) {
  matern_cholesky_cpp(a = 0, dim = dim, nu = 0)
  return(1)
}

f4 <- function(dim) {
  matern_cholesky_cpp_eigen(a = 0, dim = dim, nu = 0)
  return(1)
}

dim <- 10

bench::mark(
  "All in R" = f1a(dim),
  "All Sparse in R" = f1b(dim),
  "Matrix in R" = f2a(dim),
  "Sparse Matrix in R" = f2b(dim),
  "All in Rcpp" = f3(dim),
  "All Sparse in Rcpp::Eigen" = f4(dim),
  min_iterations = 10
)




dim <- 200
bench::mark(
  "All Sparse in R" = f1b(dim),
  "Sparse Matrix in R" = f2b(dim),
  "All Sparse in Rcpp::Eigen" = f4(dim), 
  min_iterations = 10
)


dim <- 3

Q <- make_matern_prec_matrix_2d_cpp_eigen(
  a = 0,
  dim = dim,
  nu = 0
)

L <- cholesky_bandsparse_cpp_eigen(Q)

margvars <- compute_marginal_variances_cpp(L)
sigmas <- sqrt(margvars)

diag(solve(sigmas %*% Q %*% sigmas))

library(cholesky)
library(Matrix)

matern_workflow_cpp(
  a = 0,
  dim = 3,
  nu = 0
) |> 
  solve() |> 
  diag()
