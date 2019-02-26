clearCache()
init()

massValues <- function(parameters)
{
  muG   <- parameters[1]
  mq    <- parameters[2]
  sigma <- parameters[3]
  print(paste('Parameters: muG =', muG, ', mq = ', mq, ', sigma = ', sigma))
  # Compute the background
  sol <-solve(MDSW(),muG = muG, mq = mq, sigma = sigma, h = 0.001)
  # Compute the spectrum
  Meson  <- computeSpectrum(z, urho, nEigen = 2)
  Scalar <- computeSpectrum(z, u0, nEigen = 1)
  Tensor <- computeSpectrum(z, u2, nEigen = 1)
  mesonMasses <- sqrt(Meson$energies)
  scalarMasses <- sqrt(Scalar$energies)
  tensorMasses <- sqrt(Tensor$energies)
  mesonExp <- c(0.775, 1.282)
  scalarExp <- c(1.475)
  tensorExp <- c(2.4)
  convFactor <- mesonExp[1] / mesonMasses[1]
  mesonMasses <- convFactor * mesonMasses
  scalarMasses <- convFactor * scalarMasses
  tensorMasses <- convFactor * tensorMasses
  print(paste('Predicted meson masses: rho =', mesonMasses[1], ', rho* = ', mesonMasses[2]))
  print(paste('Predicted 0^{++} mass: 0^{++} =', scalarMasses[1]))
  print(paste('Predicted 2^{++} mass: 2^{++} =', tensorMasses[1]))
  print(paste('Exp meson masses: rho =', mesonExp[1], ', rho* = ', mesonExp[2]))
  print(paste('Lattice 0^{++} mass: 0^{++} =', scalarExp[1]))
  print(paste('Lattice 2^{++} mass: 2^{++} =', tensorExp[1]))
}

J <- function(parameters)
{
  muG   <- parameters[1]
  mq    <- parameters[2]
  sigma <- parameters[3]
  flog.debug(paste('Parameters: muG =', muG, ', mq = ', mq, ', sigma = ', simga))
  # Compute the background
  sol <-solve(MDSW(),muG = muG, mq = mq, sigma = sigma, h = 0.001)
  # Compute the spectrum
  Meson  <- computeSpectrum(z, urho, nEigen = 2)
  Scalar <- computeSpectrum(z, u0, nEigen = 1)
  Tensor <- computeSpectrum(z, u2, nEigen = 1)
  mesonMasses <- sqrt(Meson$energies)
  scalarMasses <- sqrt(Scalar$energies)
  tensorMasses <- sqrt(Tensor$energies)
  mesonExp <- c(0.775, 1.282)
  scalarExp <- c(1.475)
  tensorExp <- c(2.4)
  j <- ( (mesonMasses[2]/mesonMasses[1] - mesonExp[2]/mesonExp[1])^2 +
           (scalarMasses[1]/mesonMasses[1] - scalarExp[1]/mesonExp[1])^2 +
           (tensorMasses[1]/mesonMasses[1] - tensorExp[1]/mesonExp[1])^2) / 3.0
  if ( (muG < 0) || (mq < 0) || (sigma < 0) ) j <- j + 100000
  j
}

#Parameters used in the paper arXiv 1303.6929 that predict the rho meson better
pars <- c(0.43, 0.005, 0.240)
res <- optim(par = pars, fn = J, method = "Nelder-Mead")
massValues(res$par)

#res$par = c(6.731078e-01 1.828313e-07 4.564821e-02)
