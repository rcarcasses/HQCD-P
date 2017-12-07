
#' @export
potential <- function(model = 'bcde', addAlfonsoTerm = TRUE) {
  model <- strsplit(model, split = '')[[1]]
  terms <- c()
  pars <- c()
  lastU <- NULL
  potentialSet <- NULL

  addWeakCouplingTerm <- function(){
    concatTerms('exp(2 * As) * (2 + J + (lambda / (8 * pi^2)) * (1 - 4 / (J *(J - 1)) - 4 / ((J + 1) * (J + 2)) + 4 * log(J))')
    #cat('\nPotential is now:\n', u2jBody())
    lastU <- u()
  }

  createTerm <- function(letter) paste0('(', letter, '0 + ', letter, '1 * (J - 2)) * ', letter,'factor')
  createVecPars <- function(letter) c(paste0(letter, '0 = 0'), paste0(letter, '1 = 0'))

  concatTerms <- function(term, params = NULL) {
    terms <<- c(terms, term)
    if(!is.null(params))
      pars <<- c(pars, params)
  }

  # add the terms required by the specified model:
  # ex. model = c('b', 'c') will generate a potential with the b0, b1, c0, c1 terms, etc
  lapply(model, function(letter) concatTerms(createTerm(letter), createVecPars(letter)))

  u2jBody <- function() {
    header <- pars %>% paste(., collapse = ', ') %>% paste('function(J, ', ., ') {')    # get the string of the header of the function
    body <- ''
    body <- paste('\n  if(ls < 0.01) ',
                    '    ls <- 0.01',
                    '  ', sep = '\n')
    body <- paste(body, '\n  u2 + (J - 2) * (', paste(terms, collapse = ' + '), ')\n}')
    paste(header, body)
  }

  z          <- NULL
  A          <- NULL
  As         <- NULL
  lambda     <- NULL
  u2         <- NULL
  dress      <- NULL
  lambdader1 <- NULL
  lambdader2 <- NULL
  Asder1     <- NULL
  Asder2     <- NULL
  afactor    <- NULL
  bfactor    <- NULL
  cfactor    <- NULL
  dfactor    <- NULL
  efactor    <- NULL

  loadAll <- function() {
    A          <<- get('A', envir = ihqcdEnv)
    As         <<- get('As', envir = ihqcdEnv)
    lambda     <<- get('lambda', envir = ihqcdEnv)
    u2         <<- get('u2', envir = ihqcdEnv)
    z          <<- get('z', envir = ihqcdEnv)
    dress      <<- get('dress', envir = ihqcdEnv)
    lambdader1 <<- get('lambdader1', envir = ihqcdEnv)
    lambdader2 <<- get('lambdader2', envir = ihqcdEnv)
    Asder1     <<- get('Asder1', envir = ihqcdEnv)
    Asder2     <<- get('Asder2', envir = ihqcdEnv)
    afactor    <<- 1
    bfactor    <<- exp(2 * A) * lambda^(5/6)
    cfactor    <<- (lambdader2 / lambda) - (lambdader1 / lambda)^2
    dfactor    <<- Asder2 - Asder1^2
    efactor    <<- (lambdader1 / lambda)^2
  }

  u <- function() {
    if(is.null(potentialSet)) {
      # load all the possible terms
      loadAll()
      body <- u2jBody()
      lastU <<- eval(parse(text = body));
      lastU
    } else {
      # the potential has been set, use it
      cat('Using custom potential\n')
      potentialSet
    }
  }

  rho <- function(b0 = 0, c0 = 0, d0 = 0, e0 = 0) {
    loadAll()
    list(x = z, y = b0 * bfactor + c0 * cfactor + d0 * dfactor + e0 * efactor)
  }

  plotRho0 <- function(b0 = 0, c0 = 0, d0 = 0, e0 = 0) {
    rho0 <- rho(b0, c0, d0, e0)
    plot(z, rho0)
  }

  addTerm <- function(term = c('(2/ls^2) * exp(2 * As)'), params = c('ls = 0.136')) {
    concatTerms(term, params)
    lastU <- u()
  }

  # add by default the ls term
  addTerm()

  addAlfonsoTermFun <- function(){
    concatTerms('exp(2 * A) * (J + 2)')
    # cat('\nPotential is now:\n', u2jBody())
    lastU <- u()
  }

  if(addAlfonsoTerm)
    addAlfonsoTermFun()

  #cat('Potential created:\n', u2jBody())
  getLastU <- function(...) {
    if(is.null(lastU))
      u()

    lastU(...)
  }

  getPars <- function() {
    if(is.null(potentialSet)) {
      # pars is a vector of strings
      eval(parse(text = paste0('c(', paste(pars, collapse = ','), ')')))
    } else {
      # pars is already a vector of parameters
      pars
    }
  }
  # this block allows to parse the parameters of the potential and those which are extra options
  breakParams <- function(...) {
    parsName <- c('J', names(p$getPars()))
    opts <- list()
    pars     <- list()
    mapply(function(n, v) {
            if(is.element(n, parsName))
              pars[[n]] <<- v
            else
              opts[[n]] <<- v
          }, names(list(...)), list(...))
    list(pars = pars, opts = opts)
  }

  isPotentialParameter <- function(str) {
    is.element(str, names(getPars()))
  }

  extractPotentialParameters <- function(...) {
    pars <- list(...)
    goodPars <- NULL
    mapply(function(n, p) {
                             if(isPotentialParameter(n))
                                goodPars[[n]] <<- p
                          }, names(pars), pars)

    goodPars
  }

  Mz <- function(...) {
    potPars <- extractPotentialParameters(...)
    zVal <- list(...)$z
    if(is.null(zVal)) {
      cat('Please enter a value for z to plot the effective mass term\n')
      return()
    }

    u2 <- get('u2', envir = ihqcdEnv)
    As <- get('As', envir = ihqcdEnv)
    u2fun <- splinefun(z, u2)
    Asfun <- splinefun(z, As)
    Js <- seq(0.5, 4, len = 20)

    Mz <- sapply(Js, function(j) {
      args <- as.list(c(J = j, potPars))
      u2j  <- do.call(lastU, args)
      Ufun <- splinefun(z, u2j)
      (Ufun(zVal) - u2fun(zVal)) * exp(-2 * Asfun(zVal))
    })

    list(Js = Js, Mz = Mz)
  }

  plotlsVsNew <- function(zVal = NULL) {
    if(is.null(zVal)) {
      cat('Please enter a value for z to plot the effective mass term\n')
      return()
    }

    A          <- get('A', envir = ihqcdEnv)
    As         <- get('As', envir = ihqcdEnv)
    z          <- get('z', envir = ihqcdEnv)

    Asfun <- splinefun(z, As)
    Afun  <- splinefun(z, A)

    oldTerm <- function(J) {
      ls <- 0.15
      (2/ls^2) * exp(2 * Asfun(zVal)) * (J - 2) + exp(2 * Afun(zVal)) * (J^2 - 4)
    }

    newTerm <- function(J) {
      exp(2 * Afun(zVal)) * (J^2 - J - 4)
    }

    Js <- seq(0.5, 4, len = 200)
    plot(Js, oldTerm(Js), type = 'l', main = paste0('z=', zVal))
    lines(Js, newTerm(Js), col = 'blue')
  }

  plotMz <- function(...) {
    mz <- Mz(...)
    zVal <- list(...)$z
    if(is.null(zVal)) {
      cat('Please enter a value for z to plot the effective mass term\n')
      return()
    }

    if(is.null(mz))
      return()

    plot.new()
    sapply(1:length(mz$Mz[,1]), function(i) {
      if(i == 1)
        plot(x = mz$Js, y = mz$Mz[1,], ylab = TeX('$m_{J}^2(z)$'), xlab = 'J', type = 'l')
      else
        lines(mz$Js, mz$Mz[i,], type = 'l')

      label <- paste0('z = ', zVal[i])
      Mzfun <- splinefun(mz$Js, mz$Mz[i,])
      xLab  <- mz$Js[as.integer(0.9 * length(mz$Js))]
      boxed.labels(xLab, Mzfun(xLab), labels = label, cex = 0.7, xpad = 1, ypad = 1, border = FALSE, bg = 'white')
    })
    abline(v = 2, h = 0, lty = 3, col = 'gray30')
  }

  # set the potential function externally
  # parameters has to be a vector
  setNewPotential <- function(f, parameters) {
    # set the potential
    potentialSet <<- f
    # update the private pars list
    pars <<- parameters
    # reset the last potential
    lastU <<- u()
  }

  p <- list(u = u,
            z = z, addTerm = addTerm,
            setNewPotential = setNewPotential,
            addWeakCouplingTerm = addWeakCouplingTerm,
            isPotentialParameter = isPotentialParameter,
            extractPotentialParameters = extractPotentialParameters,
            u2jBody = u2jBody,
            breakParams = breakParams,
            getPars = getPars,
            getLastU = getLastU,
            getTerms = function() terms,
            Mz = Mz,
            plotMz = plotMz,
            rho = rho,
            plotlsVsNew = plotlsVsNew)

  class(p) <- append(class(p), 'potential')
  p
}

#' @export
print.potential <- function(p, ...) print(paste('potential is:', p$u2jBody()))
#' @export
show.potential <- function(p) print(p)
#' @export
cat.potential <- function(p) print(p)

#' @export
plot.potential <- function(p, ...) {
  brokenPars <- p$breakParams(...) # distingish between potential and plot parameters
  pars <- brokenPars$pars          # assign the potential parameters
  u <- do.call(p$getLastU, pars)
  defaultOpts <- list(p$z, u, ylab = expression('U'[2]), xlab = 'z', ylim = c(-10, 50), xlim = c(0, 6), lwd = 3, type = 'l')
  opts <- modifyList(defaultOpts, brokenPars$opts)
  par(xpd=FALSE) # clip
  plot.new()
  do.call(plot, opts)
  parsNoJ <- rep(pars)
  parsNoJ$J <- NULL   # remove the J from this list
  #if(length(parsNoJ) > 0)
    #legend(x = 1.5, y = 50, legend = mapply(function(n, v) getParExpression(n, v),
     #                                      names(parsNoJ), parsNoJ))

  abline(h = 0, v = 0, col = 'gray60')
  if(is.null(opts$col))  cl <- 'blue' else cl <- opts$col
  labels <- getParExpression('J', pars$J)
  i <- as.integer(pars$J * 0.2 * length(p$z))
  boxed.labels(p$z[[i]], u[[i]], labels = labels, col = cl, cex = 0.7, xpad = 1, ypad = 1, border = FALSE, bg = 'white')
}

#' @export
lines.potential <- function(p,...) {
  brokenPars <- p$breakParams(...) # distingish between potential and plot parameters
  pars <- brokenPars$pars          # assign the potential parameters
  u <- do.call(p$getLastU, pars)
  defaultOpts <- list(p$z, u, lwd = 3, type = 'l')
  opts <- modifyList(defaultOpts, brokenPars$opts)
  do.call(lines, opts)
  i <- as.integer(pars$J * 0.2 * length(p$z))
  labels <- getParExpression('J', pars$J)
  if(is.null(opts$col))  cl <- 'blue' else cl <- opts$col
  boxed.labels(p$z[[i]], u[[i]], labels = labels, col = cl, cex = 0.7, xpad = 1, ypad = 1, border = FALSE, bg = 'white')
}

