# Compute u as a function of z
#' @export
ufun <- Vectorize(function(r) integrate(Gfun,z[1],r)$value)

#' @export
VectorMesonSpectra <- function(nEigenvalues = 5)
{
  u <- ufun(z)
  vf0 <- (1+W1*lambda+W2*lambda^2)
  VMPotential <- Ader1^2 / (4.0 * G^2) - Ader1 * dG / (2.0 * G^3) - 4 * a1 * Ader1 * lambdader1 / (G^2 * (4+3*a1*lambda))+W1 * Ader1 * lambdader1 / (2* G^2 * vf0) + W2 * lambda * Ader1 * lambdader1 / (G^2 * vf0) + 4 * a1 * dG * lambdader1 / (G^3 * (4+3*a1*lambda)) - W1 * dG * lambdader1 / (2*G^3 * vf0)- W2 * lambda * dG * lambdader1 / (G^3 * vf0) + 28 * a1^2 * lambdader1^2 / (G^2 * (4+3*a1*lambda)^2)- W1^2 * lambdader1^2/(4*G^2*vf0^2) - W1 * W2 * lambda * lambdader1^2/(G^2 * vf0^2)-W2^2 * lambda^2 * lambdader1^2 /(G^2 * vf0^2) + W2 * lambdader1^2 / (G^2 * vf0^2)- 4 * a1 * W1 * lambdader1^2 / (G^2 * (4+3*a1*lambda)*vf0)- 8 * a1 * W2 * lambda * lambdader1^2 / (G^2 * (4+3*a1*lambda)*vf0) - a0 * tau * Ader1 * tauder1 / G^2+ a0 * tau * dG * tauder1 / G^3 + 8 * a0 * a1 * tau * lambdader1 * tauder1 / (G^2 * (4+3*a1*lambda))-a0 * W1 * tau * lambdader1 * tauder1 / (G^2 * vf0) - 2 * a0 * W2 * lambda * tau * lambdader1 * tauder1 / (G^2 * vf0)-a0 * tauder1^2 / G^2 + a0^2 * tau^2 * tauder1^2 / G^2 + Ader2 / (2*G^2)-4 * a1 * lambdader2 / (G^2 * (4+3*a1*lambda)) + W1 * lambdader2 / (2*G^2 * vf0)+W2 * lambda * lambdader2 / (G^2 * vf0) - a0 * tau * tauder2 / G^2
  # Now that we have the potential we can compute the specturm
  #plot(A, VMPotential, ylim = c(0,140000))
  s <- computeSpectrum(u, VMPotential, nEigen = nEigenvalues)
  s
}
