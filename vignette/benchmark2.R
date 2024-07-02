library(cholesky)
library(Matrix)
library(purrr)
library(tidyverse)
library(glue)

dim <- 20
results <- bench::mark(
  "C++" = make_standardized_matern_sparse(dim, 0.5),
  filter_gc = FALSE,
  iterations = 10,
  check = FALSE
)

my_fun <- function(dim) {
  bench::mark(
    "C++" = make_standardized_matern_sparse(dim, 0.5),
    filter_gc = FALSE,
    iterations = 10,
    check = FALSE
  ) |> 
    mutate(
      dim = dim
    )
}


results <- map(c(10, 20, 40, 80), my_fun)

results |> 
  list_rbind() |> 
  select(Q_size = dim, time = median, memory = mem_alloc) |> 
  mutate(
    Q_size = glue("{Q_size^2}x{Q_size^2}")
  )



  my_fun <- function(rho) {
    bench::mark(
      "C++" = make_standardized_matern_sparse(20, rho),
      filter_gc = FALSE,
      iterations = 10,
      check = FALSE
    ) |> 
      mutate(
        rho = rho
      )
  }
  
  
  results <- map(c(0, 0.25, 0.5, 0.75, 1), my_fun)
  
  results |> 
    list_rbind() |> 
    select(rho, time = median, memory = mem_alloc) 
