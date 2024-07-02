# Custom sparse matrix representation
create_sparse_matrix <- function(n_rows, n_cols) {
  list(
    values = numeric(0),
    row_indices = integer(0),
    col_pointers = integer(n_cols + 1),
    n_rows = n_rows,
    n_cols = n_cols
  )
}

add_element <- function(sparse_mat, row, col, value) {
  sparse_mat$values <- c(sparse_mat$values, value)
  sparse_mat$row_indices <- c(sparse_mat$row_indices, row - 1) # 0-based indexing
  sparse_mat$col_pointers[(col + 1):(sparse_mat$n_cols + 1)] <-
    sparse_mat$col_pointers[(col + 1):(sparse_mat$n_cols + 1)] + 1
  sparse_mat
}


# Function to create a sparse Matérn precision matrix
make_sparse_matern_prec_matrix_2d <- function(a, dim) {
  n <- dim * dim
  Q <- create_sparse_matrix(n, n)
  
  for (i in 1:n) {
    Q <- add_element(Q, i, i, 4 + a)
    
    if (i %% dim != 1) Q <- add_element(Q, i, i - 1, -1)
    if (i %% dim != 0) Q <- add_element(Q, i, i + 1, -1)
    if (i > dim) Q <- add_element(Q, i, i - dim, -1)
    if (i <= n - dim) Q <- add_element(Q, i, i + dim, -1)
  }
  
  Q
}
