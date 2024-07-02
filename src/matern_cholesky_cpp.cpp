#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Function to create a band-sparse Matérn precision matrix
// [[Rcpp::export]]
NumericMatrix make_matern_prec_matrix_2d_cpp(double a, int dim) {
  int n = dim * dim;
  NumericMatrix Q(n, n);
  
  Q(0, 0) = 4 + a;
  Q(n - 1, n - 1) = 4 + a;
  Q(0, 1) = -1;
  Q(0, dim + 1) = -1;
  Q(n - 1, n - 2) = -1;
  Q(n - 1, n - dim - 2) = -1;
  
  for (int i = 1; i < (n - 1); ++i) {
    Q(i, i) = 4 + a;
    
    Q(i, i -1) = -1;
    Q(i, i + 1) = -1;
    
    if (i >= dim + 1) { // Top neighbor
      Q(i, i - dim - 1) = -1;
    }
    if (i < n - dim - 1) { // Bottom neighbor
      Q(i, i + dim + 1) = -1;
    }
  }
  
  return Q;
}


// [[Rcpp::export]]
NumericMatrix cholesky_bandsparse_cpp(NumericMatrix A) {
  int n = A.nrow();
  NumericMatrix L(n, n);
  int bandwidth = std::sqrt(n) + 1;
  
  for (int j = 0; j < n; ++j) {
    for (int i = std::max(0, j - bandwidth); i <= j; ++i) {
      double sum = 0.0;
      
      if (i == j) { // Diagonal elements
        for (int k = std::max(0, j - bandwidth); k < j; ++k) {
          sum += L(k, j) * L(k, j);
        }
        L(j, j) = std::sqrt(A(j, j) - sum);
      } else { // Off-diagonal elements
        for (int k = std::max(0, j - bandwidth); k < i; ++k) {
          sum += L(k, i) * L(k, j);
        }
        L(i, j) = (A(i, j) - sum) / L(i, i);
      }
    }
  }
  
  return L;
}

// [[Rcpp::export]]
NumericVector compute_marginal_variances_cpp(NumericMatrix L) {
  int n = L.nrow();
  int bandwidth = sqrt(n) + 1;
  NumericMatrix Sigma(bandwidth + 1, bandwidth + 1);
  NumericVector out(n);
  
  for (int i = n - 1; i >= 0; i--) {
    IntegerVector I_i = seq(i, std::min(n, i + bandwidth));
    std::reverse(I_i.begin(), I_i.end());
    I_i = I_i[I_i != i];
    
    IntegerVector index = seq_len(I_i.size());
    std::reverse(index.begin(), index.end());
    for (int j = 0; j < index.size(); j++) {
      double sum_result = 0;
      for (int k = 0; k < I_i.size(); k++) {
        sum_result += L(i, I_i[k]) * Sigma(index[j], index[k]);
      }
      Sigma(0, index[j]) = Sigma(index[j], 0) = -1 / L(i, i) * sum_result;
    }


    
    double sum_result_2 = 0;
    for (int k = 0; k < I_i.size(); k++) {
      sum_result_2 += L(i, I_i[k]) * Sigma(index[k], 0);
    }
    Sigma(0, 0) = 1 / (L(i, i) * L(i, i)) - 1 / L(i, i) * sum_result_2;
    
    out[i] = Sigma(0, 0);
    
    for (int j = bandwidth - 1; j >= 0; j--) {
      for (int k = bandwidth - 1; k >= 0; k--) {
        Sigma(j + 1, k + 1) = Sigma(j, k);
      }
    }
    
    
  }
  
  return out;
}


// Function to create the standardized Matérn precision matrix
// [[Rcpp::export]]
NumericMatrix make_standardized_matern_cpp(double a, int dim) {
    // Create the Matérn precision matrix
    NumericMatrix Q = make_matern_prec_matrix_2d_cpp(a, dim);
    
    // Perform Cholesky decomposition
    NumericMatrix L = cholesky_bandsparse_cpp(Q);
    
    // Compute marginal variances
    NumericVector sigmas = compute_marginal_variances_cpp(L);
    
    
    // Compute standardized precision matrix
    int n = Q.nrow();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q(i, j) = sigmas[i] * Q(i, j) * sigmas[j];
        }
    }
    
    return Q;
}

