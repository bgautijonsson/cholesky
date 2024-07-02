#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// Function to create a band-sparse Matérn precision matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> make_matern_prec_matrix_2d_cpp_eigen(double a, int dim) {
  int n = dim * dim;
  std::vector<Eigen::Triplet<double>> tripletList;
  
  
  
  // Fill the matrix with appropriate values
  tripletList.push_back(Eigen::Triplet<double>(0, 0, 4 + a));
  tripletList.push_back(Eigen::Triplet<double>(n - 1, n - 1, 4 + a));
  tripletList.push_back(Eigen::Triplet<double>(0, 1, -1));
  tripletList.push_back(Eigen::Triplet<double>(0, dim + 1, -1));
  tripletList.push_back(Eigen::Triplet<double>(n - 1, n - 2, -1));
  tripletList.push_back(Eigen::Triplet<double>(n - 1, n - dim - 2, -1));
  
  for (int i = 1; i < n - 1; ++i) {
    tripletList.push_back(Eigen::Triplet<double>(i, i, 4 + a));
    tripletList.push_back(Eigen::Triplet<double>(i, i - 1, -1));
    tripletList.push_back(Eigen::Triplet<double>(i, i + 1, -1));
    
    if (i >= dim + 1) {
      tripletList.push_back(Eigen::Triplet<double>(i, i - dim - 1, -1));
    }
    if (i < n - dim - 1) {
      tripletList.push_back(Eigen::Triplet<double>(i, i + dim + 1, -1));
    }
  }
  
  Eigen::SparseMatrix<double> Q(n, n);
  Q.setFromTriplets(tripletList.begin(), tripletList.end());
  
  return Q;
}

// Function to perform Cholesky factorization on a band-sparse matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> cholesky_bandsparse_cpp_eigen(Eigen::SparseMatrix<double> A) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        Rcpp::stop("Cholesky factorization failed. The matrix might not be positive definite.");
    }
    return solver.matrixU();
}




// [[Rcpp::export]]
Eigen::VectorXd compute_marginal_variances_cpp_eigen(const Eigen::SparseMatrix<double>& L) {
  int n = L.rows();
  int bandwidth = std::sqrt(n) + 1;
  Eigen::SparseMatrix<double> Sigma(n, n);
  Sigma.reserve(VectorXi::Constant(n, bandwidth));
  Eigen::VectorXd sigmas(n);

  for (int i = n - 1; i >= 0; --i) {
    std::vector<int> I_i;
    for (int j = i + 1; j < std::min(n, i + bandwidth); ++j) {
      I_i.push_back(j);
    }

    for (int j : I_i) {
      double value = 0.0;
      for (int k : I_i) {
        value += L.coeff(i, k) * Sigma.coeff(k, j);
      }
      value = -value / L.coeff(i, i);
      Sigma.coeffRef(i, j) = Sigma.coeffRef(j, i) = value;
    }

    // Diagonal element
    double diag_value = 1.0 / (L.coeff(i, i) * L.coeff(i, i));
    for (int k : I_i) {
      diag_value -= L.coeff(i, k) * Sigma.coeff(k, i) / L.coeff(i, i);
    }
    Sigma.coeffRef(i, i) = diag_value;
    sigmas[i] = std::sqrt(diag_value);
  }

  return sigmas;
}


// [[Rcpp::export]]
Eigen::MatrixXd make_standardized_matern_cpp_eigen(double a, int dim) {
    // Create the Matérn precision matrix
    Eigen::SparseMatrix<double> Q = make_matern_prec_matrix_2d_cpp_eigen(a, dim);
    
    // Perform Cholesky decomposition using the provided function
    Eigen::SparseMatrix<double> L = cholesky_bandsparse_cpp_eigen(Q);
    
    // Compute marginal variances (now returns a diagonal sparse matrix)
    Eigen::VectorXd  sigmas = compute_marginal_variances_cpp_eigen(L);
    
    
    // Convert to dense matrix for output
    return Q;
}