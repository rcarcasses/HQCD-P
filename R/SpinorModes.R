VL <- function(z, m5 = 2.5) {
  # This function computes the potential for a left-handed mode
  # Returns a float if z is a float or an array if z is an array
  return( exp(2*AsSpline(z)) * m5^2 - m5 * AsSpline(z, deriv = 1) * exp(AsSpline(z)) )
}

VR <- function(z, m5 = 2.5) {
  # This function computes the potential for a right-handed mode
  # Returns a float if z is a float or an array if z is an array
  return( exp(2*AsSpline(z)) * m5^2 + m5 * AsSpline(z, deriv = 1) * exp(AsSpline(z)) )
}

getSpinorLModeRaw <- function(m5 = 2.5)
{
  # Computes the lowest L eigenvalue wavefunction for a bulk field of mass m5
  # It returns gL^2 spline function
  # m5 = 5/2 is the value for a spinor field dual to an operator of Delta = 9/2
  # Compute the potential for the left-handed mode
  V  <- VL(z, m5)
  # We identify the proton as the lowest energy state
  data <- computeSpectrum(z, V, nEigen = 1)$wfs[[1]]
  wf <- splinefun(data$x, data$y^2)
  return(wf)
}

getSpinorRModeRaw <- function(m5 = 2.5)
{
  # Computes the lowest R eigenvalue wavefunction for a bulk field of mass m5
  # It returns gR^2 spline function
  # m5 = 2.5 is the value for a spinor field dual to an operator of Delta = 9/2
  # Compute the potential for the left-handed mode
  V  <- VR(z, m5)
  # We identify the proton as the lowest energy state
  data <- computeSpectrum(z, V, nEigen = 1)$wfs[[1]]
  wf <- splinefun(data$x, data$y^2)
  return(wf)
}

getExternalProtonFactorRaw <- function(m5 = 2.5)
{
  # Returns gL^2 + gR^2 that will be usefull later
  gL2fun <- getSpinorLModeRaw(m5)
  gR2fun <- getSpinorRModeRaw(m5)
  protonExtFactorFun <- splinefun(z, gL2fun(z) + gR2fun(z))
  return(protonExtFactorFun)
}

#' @export
getSpinorLMode <- cache(getSpinorLModeRaw)

#' @export
getSpinorRMode <- cache(getSpinorRModeRaw)

#' @export
getExternalProtonFactor <- cache(getExternalProtonFactorRaw)



