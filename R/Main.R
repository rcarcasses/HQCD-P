#' @import futile.logger
#' @import lubridate
#' @import schrodinger
#' @import bvpSolve
#' @import nloptr
#' @import parallel
#' @import magrittr
#' @import rootSolve
#' @import plotrix
#' @import rredis
#' @import cmaes
#' @import scales
#' @import tictoc
#' @import optimx
#' @import testthat
#' @import latex2exp
#' @import minpack.lm

#' @export
init <- function(chebPoints = 400, useRedis = TRUE) {
  flog.debug("[HQCD-P] Initializing")
  if(useRedis)
    startRedis()

  ihqcd <- iHQCD()
  #populate the ihqcd environment
  ihqcd$solve()
  # set the method we want to use to compute the eigenvalues
  schrodinger::setSchroMethod('cheb', chebPoints);
}

# we have to define the value to something otherwise the <<- won't work
cl <- NULL
.onLoad <- function(lib, pkg) {
  # Calculate the number of cores, left one for the system
  cores <- detectCores() - 1
  cat('Number of cores available:', cores, '\n')
  # Initiate cluster, all the variables available to all nodes
  cl <<- makeCluster(cores)
}

#' @export
dumpList <- function(l) {
  paste(mapply(function(n, v) paste(n,'=', v), names(l), l), collapse = ', ')
}

on.exit({
  cat('Stopping cluster\n')
  if(!is.null(cl))
  parallel::stopCluster(cl)
})

#' @export
copyEnv <- function(from, to, names=ls(from, all.names=TRUE)) {
  mapply(assign, names, mget(names, from), list(to),
         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  invisible(NULL)
}

#' @export
catList <- function(pars) {
  mapply(function(n, v) cat(paste0(n,' = ', v, ', ')), names(pars), pars)
  cat('\n')
}

#' @export
catListStart <- function(pars) {
  mapply(function(n, v) cat(paste0(n,' = c(0, ', v, ', 0), ')), names(pars), pars)
  cat('\n')
}

#' @export
exportMathematica <- function(x, y, name = 'mathematica.txt') {
  str <- '{'
  mapply(function(xval, yval)
          str <<- paste0(str, '{', xval, ', ', yval, '},')
    , x, y)

  str <- substr(str, 1, nchar(str)-1)
  str <- paste0(str, '}')

  write(str, file = name)
}

# this is a nice helper function to produce nice labels inside plots
#' @export
getParExpression <- function(par, val, index = NULL)  {
  eqStr <- paste0('==', val)
  str <- 'expression("'
  str <- paste0(str, substr(par, 1, 1), '"')
  if(nchar(par) > 1)
    str <- paste0(str, '[', substr(par, 2, nchar(par)), ']')
  str <- paste0(str, eqStr,')')
  eval(parse(text = str))
}

# some helpful function to remove all the things already loaded
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
# detachAllPackages()

# Some test function (actually the one being used in the NMC paper)
#' @export
UJgTest <- function(J, invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117) {
  # the deformation of the potential of the graviton
  u2 + (J - 2) * ((2 * e2As * invls^2) * (1 + d / l1_2) + e2A * (J + 2)
                  + (a * aF + b * bF + c * cF))
}
