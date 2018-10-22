# Load necessary libraries
library(schrodinger)
library(HQCDP)

init()
s <- iHQCD()$solve()

# Load z, A, Phi, dA, dPhi, d2A, d2Phi
z     <- s$z
zmin  <- z[1]
zmax  <- 10
z     <- z[z <= zmax]
A     <- s$A[s$z <= zmax]
Phi   <- s$Phi[s$z <= zmax]
dA    <- s$Ader1[s$z <= zmax]
d2A   <- s$Ader2[s$z <= zmax]
dPhi  <- s$Phider1[s$z <= zmax]
d2Phi <- s$Phider2[s$z <= zmax]

# Define the spline functions
zfun     <- splinefun(A, z)
Phifun   <- splinefun(A, Phi)
dAfun    <- splinefun(A, dA)
d2Afun   <- splinefun(A, d2A)
dPhifun  <- splinefun(A, dPhi)
d2Phifun <- splinefun(A, d2Phi)
# Now that the background fields are defined we can compute the normalizable modes
# Definition of the system of ODEs we want to solve
fun <- function(t, y, eig){
  A <- (dAfun(t)^2 * (1 - dPhifun(t)) + d2Afun(t)) / dAfun(t)^2
  B <- eig / dAfun(t)^2
  dy1 <- y[2]
  dy2 <- - A * y[2] - B * y[1]
  return(list(c(dy1, dy2)))
}
# AdS boundary conditions
mass     <- sqrt(11.170932) # 2.906454  6.900528 11.170932 15.741255 20.456861
y0AdS    <- mass * z[1] * besselJ(mass * z[1], 1)
dy0AdS   <- mass * mass * z[1] * besselJ(mass * z[1], 0)

init <- c(y0AdS,  dy0AdS)
end  <- c(NA, 0)

sol <- bvpshoot( yini = init, yend = end,
                 x = A, func = fun,
                 extra = 6.9)

attributes(sol)$roots

numSol <- data.frame(sol)
x    <- numSol$x  # Array of independent variable
yST  <- numSol$X1 # ST solution for the first eigenvalue
ySTfun <- splinefun(z, yST)
plot(z, ySTfun(z))
# Now i want to see how compatible are these eigenfunctions to the wfs
# of the equivalent Schrodinger problem
# The ST solution and Sch solution are related by ST = exp((Phi-A)/2) Sch

V <- 0.25 * ( dA * dA - 2 * dA * dPhi + dPhi * dPhi + 2 * d2A - 2 * d2Phi )

# Now let's compute the spectrum
spec        <- computeSpectrum(z, V, 5)
eigenvalues <- spec$energies
eigenfun    <- splinefun(spec$wfs[[1]]$x, spec$wfs[[1]]$y)
yScfun <- splinefun(x, exp(0.5 * (Phifun(x)-Afun(x)))*eigenfun(x))

# Plot to compare the two functions
plot(x, yScfun(x), col = "red", type = "l")
lines(x, ySTfun(x) , col = "blue")

# There is a difference in the shape of the curves. I checked that if
# I tune the first derivative i can make them match. i wanna check if this
# is some normalization problem.
f1 <- splinefun(x, exp(Afun(x) - Phifun(x)) * ySTfun(x) * ySTfun(x))
f2 <- splinefun(x, exp(Afun(x) - Phifun(x)) * yScfun(x) * yScfun(x))
# Compute the normalization constant
C1 <- integrate(f = f1, lower = x[1], upper = x[length(x)])$value
C2 <- integrate(f = f2, lower = x[1], upper = x[length(x)])$value
# plot the functions now
ySTfun <- splinefun(x, ySTfun(x) / sqrt(C1))
yScfun <- splinefun(x, yScfun(x) / sqrt(C2))
plot(x,  yScfun(x) , col = "red", type = "l")
lines(x, ySTfun(x) , col = "blue")
plot(x, ySTfun(x)/yScfun(x), ylim = c(0,2))

# Now the differences between two functions are purely numerical
# It turns out to be a good idea to use the equivalent Sch problem
# to guess the eigenvalues
