# Here we implment the couplings between fields of different spin

kg <- function(s1, s2, s3, d = 5, N = 1){
  # This function gives the value of the coupling between three fields of spin s1, s2 and s3 in d dimensions.
  # N must be used as a fitting parameter
  result <- pi^((d-3)/4) * 2^((3 * d - 1 + s1 + s2 + s3)/2) / (sqrt(N) * gamma(d + s1 +s2 +s3 -3))
  result <- result * sqrt(gamma(s1 + (d-1)/2) * gamma(s2 + (d-1)/2) * gamma(s3 + (d-1)/2)/(gamma(s1 + 1) * gamma(s2 + 1) * gamma(s3 + 1)))
  result
}

#' @export
kU1 <- function(J, N = 1) kg(s1 = 1 , s2 = 1, s3 = J, d = 5, N = N)
#' @export
kPhi <- function(J, N = 1) kg(s1 = 0 , s2 = 0, s3 = J, d = 5, N = N)
