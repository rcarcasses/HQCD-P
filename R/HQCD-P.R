#' @import futile.logger
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
#' @import latex2exp
#' @import minpack.lm

#' @export
init <- function(chebPoints = 400, userRedis = TRUE) {
  flog.debug("[HQCD-P] Initializing")
  if(userRedis)
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

# TODO, not working when the name of a parameter is longer than the value, like when b0 = 1 case
catBeautifulNamedVector <- function(v, namTag = '', valTag = '', noNames = FALSE, justNames = FALSE) {
  formatedName <- c()
  mapply(function(n, val) {
                      valLen <- length(strsplit(paste0(val), split = '')[[1]])
                      if(valLen < 4)
                        valLen <- length(strsplit(format(val, digits = 4), split = '')[[1]])

                      nameLen      <- length(strsplit(n, split = '')[[1]])
                      offsetX      <- as.integer((valLen - nameLen) / 2)
                      leftSpace    <- paste0(rep(' ', len = offsetX), collapse = '')
                      rightSpace   <- paste0(rep(' ', len = valLen - nameLen - offsetX), collapse = '')
                      formatedName <<- c(formatedName, paste0(leftSpace, n, rightSpace))
                    }, names(v), v)

  if(!noNames)
    cat(paste0('\n', namTag, paste0(formatedName, collapse = ' | ')))
  if(!justNames)
    cat(paste0('\n', valTag, paste0(v, collapse = ' | ')))
}

on.exit({
  cat('Stopping cluster\n')
  if(!is.null(cl))
  parallel::stopCluster(cl)
})

#' This is freaking helpful!
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

detachAllPackages()
