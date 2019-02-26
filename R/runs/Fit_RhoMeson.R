library(HQCDP)
init()

flog.threshold(TRACE)

# Create VMP objects
vmpdsigma <- VMPDSigma("Phi")
vmpsigma  <- VMPSigma("Phi")

attr(vmpdsigma, 'IzNBar') <- function(dsigma, row, spec, zstar, hpars) 1
attr(vmpsigma, 'IzNBar')  <- function(dsigma, row, spec, zstar, hpars) 1

kernel_pars <- c(invls = 1 / 0.153, a = -4.35, b = 1.41, c = 0.626, d = - 0.117)
kpars <- c(1.0, 0.0, 0.0)

fixed <- list(
  # zstar <- 1.0
  pars = kernel_pars
  # hpars = kpars
)

HQCDP(rsslog = TRUE, fixed = fixed) %>%
  addKernel(numReg = 4) %>%
  addProcessObservable(vmpdsigma) %>%
  addProcessObservable(vmpsigma) -> p
attr(p, 'cacheSpectra') <- FALSE
attr(p, 'addSPconstraint') <- 1e4
ts <- seq(-1.2, 0.0, len = 12)
attr(p, 'useTVals') <- ts

fit(p, pars = c(invls = NA, a = NA, b = NA, c = NA, d = NA),
    zstar = 1.0,
    hpars = kpars)

warnings()
