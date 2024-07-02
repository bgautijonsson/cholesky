
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cholesky

This directory is for development of fast and memory-efficient code that
creates Matérn precision matrices that have been standardized so that
their inverse is a correlation matrix.

The package can be installed with

``` r
pak::pak("bgautijonsson/cholesky")
```

``` r
library(cholesky)
#> Loading required package: Matrix
```

``` r
make_standardized_matern(dim = 3, rho = 0.5)
#> 9 x 9 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,]  1.155093 -0.2808520  .        -0.2808520  .          .          .       
#>  [2,] -0.280852  1.2291667 -0.280852  .         -0.2661131  .          .       
#>  [3,]  .        -0.2808520  1.155093  .          .         -0.2808520  .       
#>  [4,] -0.280852  .          .         1.2291667 -0.2661131  .         -0.280852
#>  [5,]  .        -0.2661131  .        -0.2661131  1.2962963 -0.2661131  .       
#>  [6,]  .         .         -0.280852  .         -0.2661131  1.2291667  .       
#>  [7,]  .         .          .        -0.2808520  .          .          1.155093
#>  [8,]  .         .          .         .         -0.2661131  .         -0.280852
#>  [9,]  .         .          .         .          .         -0.2808520  .       
#>                           
#>  [1,]  .          .       
#>  [2,]  .          .       
#>  [3,]  .          .       
#>  [4,]  .          .       
#>  [5,] -0.2661131  .       
#>  [6,]  .         -0.280852
#>  [7,] -0.2808520  .       
#>  [8,]  1.2291667 -0.280852
#>  [9,] -0.2808520  1.155093
```
