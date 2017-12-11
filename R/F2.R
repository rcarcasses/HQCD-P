#' @export
F2 <- function(potPars = 'bcde',
                 numReg = 3,
                 data = NULL,
                 kernelName = '',
                 considerNonMinimalCoupling = TRUE) {
  if(is.null(data))
    data <- getHERAF2()

  # at initialization get the data that will be needed later
  As      <- get('As', envir = ihqcdEnv)
  Asder1  <- get('Asder1', envir = ihqcdEnv)
  Asder2  <- get('Asder2', envir = ihqcdEnv)
  Phider1 <- get('Phider1', envir = ihqcdEnv)
  Phider2 <- get('Phider2', envir = ihqcdEnv)
  z       <- get('z', envir = ihqcdEnv)
  A0      <- get('A0', envir = ihqcdEnv)
  h       <- get('h', envir = ihqcdEnv)

  # create a kernel
  kernel <- kernelBCDQ(potPars, numReg, kernelName = kernelName)

  # create a regge object
  r <- regge(potPars, numReg = numReg)
  # this variable stores the value of the last computation of the intercepts
  Js  <- NULL
  # this variables store the values for the masses for J=2,4,6 for the first and second trajectories
  tr1 <- NULL
  tr2 <- NULL

  # this function set the potential for the kernel
  setNewPotential <- function(f, pars) {
    kernel$setNewPotential(f, pars)
  }

  plotEffectiveExponent <- function(...) {
    # plot the experimental results
    eff <- data$eff
    plot(eff$Q2, eff$ep, log = 'x', ylim = c(0, 0.5), pch = 20, xlab = TeX('Q^2'), ylab = TeX('$\\epsilon'))
    abline(h = 0.1 * (0:5), v = c(0.2, 2, 20, 100), col = "lightgray", lty = 3)
    # plot the error bars
    arrows(eff$Q2, eff$ep - eff$epErr, eff$Q2, eff$ep + eff$epErr,
           length = 0.02, angle = 90, code = 3)

    # plot the model estimation
    Q2L <- length(data$Q2s)
    effTheo <- data.frame(Q2 = numeric(Q2L), ep = numeric(Q2L), epErr = numeric(Q2L))
    F2fun <- function(Q2, x) do.call(F2, c(list(Q2 = Q2, x = x), list(...)))
    for(i in 1:Q2L) {
      Q2 <- data$Q2s[i]
      cat('Q2 =', Q2, ' x =', unlist(data$xs[[i]]), '\n')
      if(length(data$xs[[i]]) > 2) {
        F2s <- unlist(lapply(data$xs[[i]], function(xx) {
          val <- F2fun(Q2, xx)
          val
        }))
        # now we do the same thing we did as for the experimental data
        f2xFit <- data.frame(x = data$xs[[i]], F2 = F2s)
        fit <- lm( log(F2) ~ log(x),
                   weights = data$weights[[i]],
                   data = f2xFit)
        s <- summary(fit)$coefficients
        effTheo$Q2[i] <- Q2
        effTheo$ep[i] <- -s['log(x)', 'Estimate']
      }
    }

    lines(effTheo$Q2, effTheo$ep, col = 'red', type = 'o', pch = 21)
  }

  extractCoeffsAndPars <- function(...) {
    parArgs <- list(...)
    coeff   <- NULL
    pars    <- NULL
    mapply(function(n, p) {
                            if(n == 'A0') { # this is a special case
                              # solve IHQCD with this A0, get some mode to trigger the whole computation
                              p <- min(3.4, max(2.2, p))
                              A0 <<- p
                              m <- getMode(0.15, A0 = A0, h)
                              As <<- get('As', envir = ihqcdEnv)
                              z  <<- get('z', envir = ihqcdEnv)
                            } else {
                              # here we split the coefficients, which are the constants multiplying the
                              # overall fn, of the parameters of the potential
                              if(!kernel$potential$isPotentialParameter(n))
                                coeff[[n]] <<- p
                              else
                                pars[[n]] <<- p
                            }
                          }, names(parArgs), parArgs)
    list(coeff = coeff, pars = pars)
  }

  # Returns the value of F2 for the given parameters
  # coeff are the s0, s1, s2, etc
  # pars are ls, b0, b1, c0, etc
  value <- function(Q2, x, ...) {
    parArgs <- extractCoeffsAndPars(...)
    coeff   <- parArgs$coeff
    pars    <- parArgs$pars

    # if no coefficients are passed, then assume that we want to use the optimal ones
    if(is.null(coeff)) {
      coeff <- findBestCoeff(...)$coefficients
      # cat('\n [F2 eval] Using optimal coefficients', coeff)
    }
    # compute the kernel
    spect <- force(do.call(kernel$findKernel, as.list(pars)))
    tr1   <<- spect$tr1
    tr2   <<- spect$tr2
    spect$tr1 <- NULL
    spect$tr2 <- NULL

    Js <<- lapply(spect, function(spec) spec$js)
    # cat('coeff =', unlist(coeff), '\n')
    # if there is any 0 just return NA
    if(length(Filter(function(js) js == 0, Js)) > 0)
      return(NA)

    # the [1:numReg] is because sometimes the kernel was computed with 4 reggeons and stored in memory
    # but in the current computation we may need only say 2.
    tfn <- lapply(spect[1:numReg], function(spec) Fn(Q2, spec$js, spec$wf) * x^(1 - spec$js))

    # here we consider the contribution to F2 from the non-minimal coupling case if require
    if(considerNonMinimalCoupling)
      tfn <- c(tfn, lapply(spect[1:numReg],
                           function(spec)
                             FnNMC(Q2, spec$js, spec$wf) * x^(1 - spec$js)))

    ts <- coeff * unlist(tfn)
    # sum all the reggeon contributions
    sum(ts)
  }

  # this return a dataframe of Fn times the x factor part
  # for each value of x and Q2
  # fn can be Fn or FnNMC currently
  getGenFns <- function(fn, ...) {
    # compute the kernel if necessary
    spect <- force(do.call(kernel$findKernel, list(...)))[1:numReg]
    # update the Js values
    Js <<- lapply(spect, function(spec) spec$js)

    dfCmdString <- paste0('data.frame(', paste0(lapply(0:(numReg-1), function(i) paste0('f', i, ' = numeric()')), collapse = ', '), ')')
    df <- eval(parse(text = dfCmdString))

    mapply(function(Q2, x){
              r <- unname(unlist(lapply(spect,
                             function(spec) {
                               fn(Q2, spec$js, spec$wf) * x^(1 - spec$js)
                            })))
              df <<- rbind(df, r)
    }, data$Q2, data$x)
    df
  }

  getFns    <- function(...) getGenFns(Fn, ...)
  getFnsNMC <- function(...) getGenFns(FnNMC, ...)
  getAllFns <- function(...) {
    cbind(getGenFns(Fn, ...), getGenFns(FnNMC, ...))
  }

  getStoredFn <- function(index = 1) {
    fns <- rredis::redisSUnion('fns')
    # get the size of the list of fns
    sizeFns <- length(fns[[1]])
    if(sizeFns < index + 1){
      cat('fns data found have only', sizeFns, ' functions, you are asking for', index, ' please check your computation\n')
      return()
    }
    q <- c()
    f <- c()
    lapply(fns, function(fn) {
      q <<- c(q, fn[[1]])
      f <<- c(f, fn[[index + 1]])
    })

    d <- data.frame(list(Q2 = q, fn = f))
    # return the ordered version respect to the first column Q2
    d[order(d[,1]), ]
  }

  plotStoredFn <- function(index = 1) {
    fn <- getStoredFn(index)
    # here the new convention is used: reggeons label start at 1 (before we use 0)
    plot(fn$Q2, fn$fn, log = 'xy', xlab = expression(Q^2),
         ylab = eval(parse(text = paste0('expression(f[', index, '])'))), type = 'o', pch = 16)

    abline(v = c(1, 5, 20, 100), h = c(0.2, 1, 2), col = 'lightgray', lty = 3)
  }

  plotAllStoredFn <- function() {
    plotStoredFn(1)
    fns <- rredis::redisSUnion('fns')
    lapply(2:length(fns[[1]]), function(index) {
      fn <- getStoredFn(index)
      lines(fn$Q2, fn$fn, type = 'o', pch = 16 - index)
    })
  }

  rssFunBase <- function (fn = NULL, dof = NULL, ...) {
    if(is.null(fn))
      F2fun <- function(Q2, x) do.call(F2, c(list(Q2 = Q2, x = x), list(...)))
    else
      F2fun <- function(Q2, x) do.call(fn, c(list(Q2 = Q2, x = x), list(...)))

    if(is.null(dof))
      dof <- length(list(...))

    chi2vec <- mapply(function(F2, Q2, x, err) {
                                        ((F2fun(Q2, x) - F2) / err)^2
                                      }
                                    , data$F2, data$Q2, data$x, data$err)
    list(dof = dof, chi2vec = chi2vec)
  }

  rssFun <- function(fn = NULL, dof = NULL, ...) {
    # compute the fit degrees of freedom
    dof <- length(data$F2) - length(list(...)) - numReg
    if(considerNonMinimalCoupling)
      dof <- dof - numReg

    fit <- findBestCoeff(...)

    # now return the best data found
    s    <- sum((resid(fit)/data$err)^2)
    chi2 <- s / dof
    list(rss = s, realRss = s, dof = dof, chi2 = chi2, coeff = fit$coefficients,
         size = length(data$F2), Js = Js, tr1 = tr1, tr2 = tr2, kernel = kernel)
  }

  findBestCoeff <- cache(findBestCoeffFun, kernelName)
  findBestCoeffFun <- function(...) {
    # first get the columns of the fns
    fs <- getFns(...)
    if(considerNonMinimalCoupling)
      fs <- cbind(fs, getFnsNMC(...))

    # now get the best linear fit for the coefficients
    fsString <- paste0(lapply(1:length(fs), function(i) paste0('fs[,', i, ']')), collapse = ' + ')
    # this returns a fit model: the output of the lm function
    eval(parse(text = paste0('lm(data$F2~', fsString,' - 1, weights = 1/data$err^2)')))
  }

  rssFunOld <- function(fn = NULL, dof = NULL, ...) {
    base    <- rssFunBase(fn, dof, ...)
    chi2vec <- base$chi2vec
    dof     <- base$dof

    s       <- sum(chi2vec)
    # let's force a better match for the lowest Q2
    ss      <- s#sum(chi2vec /(data$Q2^(0.2)))
    chi2    <- s / (length(data$F2) - dof)
    list(rss = ss, realRss = s, dof = dof, chi2 = chi2, size = length(data$F2), Js = Js, tr1 = tr1, tr2 = tr2, kernel = kernel)
  }

  rssFunLambda <-  function(fn = NULL, dof = NULL, ...) {
    base    <- rssFunBase(fn, dof, ...)
    chi2vec <- base$chi2vec
    dof     <- base$dof

    sLambda <- sum(unlist(lapply(chi2vec, function(c2) log(1 + 0.18 * c2))))
    s       <- sum(chi2vec)
    chi2    <- s / (length(data$F2) - dof)
    list(rss = sLambda, dof = dof, chi2 = chi2, size = length(data$F2), Js = Js, tr1 = tr1, tr2 = tr2, kernel = kernel)
  }

  Fn    <- function(Q2, J, wf) Q2^J * P13bar(Q2, J, wf)
  FnNMC <- function(Q2, J, wf) Q2^J * P13barNMC(Q2, J, wf)

  P13bar <- function(Q2, J, wf) {
    # Alfonso's note eq 2.48
    t1fun <- splinefun(z, exp((-J + 1.5) * As))
    t2fun <- getMode(Q2, A0, h)$factor
    t3fun <- splinefun(wf$x, wf$y)

    # I think one of the biggest problems I had was because I was just summing
    # instead of doing a proper integration, so here we go.
    integral <- integrate(function(x) t1fun(x) * t2fun(x) * t3fun(x), z[1], z[length(z)], stop.on.error = FALSE)
    integral$value
  }

  P13barNMC <- function(Q2, J, wf) {
    mode    <- getMode(Q2, A0, h)
    f2Q     <- mode$f2Q
    f2Qder1 <- mode$f2Qder1
    # here is where the non-minimal coupling changes are implemented
    t1fun   <- splinefun(z, exp((-J - 0.5) * As))
    wffun   <- splinefun(wf$x, wf$y)
    wfEval  <- wffun(z)
    wfder1  <- wffun(z, deriv = 1)
    wfder2  <- wffun(z, deriv = 2)
    # contribution from the transverse part
    t2fun   <- splinefun(z, f2Q(z) * ((Asder1 * Phider1 - 0.5 * Asder1^2) * wfEval + Asder1 * wfder1))
    # contribution from the longitudinal part
    t3fun   <- splinefun(z, (f2Qder1(z) / Q2) * (wfder2 + 2 * (Phider1 - Asder1) * wfder1 +
         (Phider2 + (J - 1.5) * Asder2 - Asder1 * (2 * Phider1 + J - 1) + Phider1^2 + 0.75 * Asder1^2) * wfEval))

    integral <- integrate(function(x) t1fun(x) * (t2fun(x) + t3fun(x)), z[1], z[length(z)], stop.on.error = FALSE)
    # cat('int ', integral$value, '\n')
    integral$value
  }

  P13barVecQ2    <- Vectorize(P13bar, 'Q2')
  P13barVecQ2NMC <- Vectorize(P13barNMC, 'Q2')
  FnVecQ2        <- Vectorize(Fn, 'Q2')

  P13Bar <- function(n, Q2 = NULL, ...) {
    if(is.null(Q2))
      Q2 <- Q2s
    spect <- force(do.call(kernel$findKernel, list(...)))
    Js <<- lapply(spect, function(spec) spec$js)
    P13barVecQ2(Q2, Js[[n]], spect[[n]]$wf)
  }

  P13BarNMC <- function(n, Q2 = NULL, ...) {
    if(is.null(Q2))
      Q2 <- Q2s
    spect <- force(do.call(kernel$findKernel, list(...)))
    Js <<- lapply(spect, function(spec) spec$js)
    P13barVecQ2NMC(Q2, Js[[n]], spect[[n]]$wf)
  }
  spectrum <- function(...) {
    spectrumBCDQ(do.call(kernel$findKernel, list(...)))
  }

  plotSpectrum <- function(...) plot.spectrumBCDQ(spectrum(...))

  plotRegge <- function(...) {
    do.call(r$plot, list(...))
  }

  plotGlueballMasses <- function() {
    r$plotGlueballMasses()
  }

  plotMesonsMasses <- function() {
    r$plotMesonsMasses()
  }

  fi <- function(n, Q2 = NULL, ...) {
    if(is.null(Q2))
      Q2 <- Q2s
    spect <- force(do.call(kernel$findKernel, list(...)))
    Js <<- lapply(spect, function(spec) spec$js)
    FnVecQ2(Q2, Js[[n]], spect[[n]]$wf)
  }

  # a parameter for F2 is a coefficient starting by s or is a parameter for the potential
  isF2parameter <- function(str) is.element(str, sapply(0:(numReg - 1), function(i) paste0('s',i))) || kernel$potential$isPotentialParameter(str)

  plotF2 <- function(...) {
    parArgs <- list(...)

    rate <- 1
    if(!is.null(parArgs$rate)) {
      rate <- parArgs$rate
      parArgs$rate <- NULL   #remove it
    }

    # figure out which parameters are for plotting and which for calling F2
    f2Pars   <- NULL
    plotPars <- NULL
    mapply(function(n, p)
      if(isF2parameter(n))
        f2Pars[[n]] <<- p
      else if(n != 'F2')
        plotPars[[n]] <<- p
      , names(parArgs), parArgs)

    # we would like eventually pass a different F2 function for testing
    F2in <- function(Q2, x) do.call(F2, c(list(Q2 = Q2s[ii], x = x), f2Pars))
    if(!is.null(parArgs$F2))
      F2in <- parArgs$F2

    data <- getHERAF2()
    xLimMax <- 2 * max(data$x)
  	plotPars <- c(list(x = data$x, y = data$F2, type = "n", log = 'x',
         xlim = c(1e-6, xLimMax), ylim = c(0.6 * min(data$F2), 1.2 * max(data$F2)),
         xlab = expression("x"), ylab = expression("F"[2])), xaxt="n", plotPars)
  	plot.new()
    do.call(plot, plotPars)
    hTicks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
    vTicks <- seq(0.5, 2, len = 4)
  	abline(v = hTicks, h = vTicks, col = 'gray90', lty = 3)
  	axis(1, at = hTicks)
  	axis(2, at = vTicks)

    Q2s <- unique(data$Q2)

  	# first let's create a bunch of colors to differentiate the graphs
  	cl <- sapply(1:2, function(i) paste0('grey', i * 15))
  	cl <- c('blue', 'red', 'green')
  	cl <- rep_len(cl, as.integer(length(Q2s)))
  	pchs <- rep_len(0:6, as.integer(length(Q2s)))

  	for(i in 1:as.integer(length(Q2s) / rate))
  	{
  	  ii <- i * rate
  		# now we need to extract the columns that we are interested for a given value of Q2
  		f2x <- NULL
  	  mapply(function(F2, Q2, x, err)
        	      if(Q2 == Q2s[ii]) {
        	        f2x$x   <<- c(f2x$x, x)
        	        f2x$F2  <<- c(f2x$F2, F2)
        	        f2x$err <<- c(f2x$err, err)
        	      }
        	    , data$F2, data$Q2, data$x, data$err)
  	  # Draw the fit on the plot by getting the prediction from the fit at 200 x-coordinates across the range of xdata
  	  x <- c(0.5 * min(f2x$x))
  	  # put 10 theoretical points between two experimental
  	  lapply(f2x$x, function(xVal) x <<- c(x, seq(x[length(x)], xVal, len = 10)))
  	  x <- c(x, 2 * max(f2x$x))
  	  #x <- seq(1e-6, xLimMax, len = 200)
  		fitPlot <- sapply(x, function(x) do.call(F2in, list(Q2 = Q2s[ii], x = x)))
  		lines(x, fitPlot, col = cl[ii], type = 'l')
  		# plot the dots
  		lines(f2x$x, f2x$F2, type = "p", col= cl[ii], cex = 0.5, pch=pchs[ii])
  		# plot the error bars
  		arrows(f2x$x, f2x$F2 - f2x$err, f2x$x, f2x$F2 + f2x$err,
  		       length = 0.02, angle = 90, code = 3, col = cl[ii])

  		# put the Q2 value at the end of the lined
  		textPosX <- 0.35 * min(f2x$x)
  		#textPosX <- 10^(0.65 * (log(x, base = 10) + 5))
  		#cat('textX', textPosX,'\n')
  		textPosY <- do.call(F2in, list(Q2 = Q2s[ii], x = textPosX))
      boxed.labels(textPosX, textPosY, labels = paste(Q2s[ii]), col = cl[ii], cex = 0.7, xpad = 1, ypad = 1, border = FALSE, bg = 'white')
  	}
  }

  fDL <- c(function(Q2, x)  0.001461 * (1 + Q2/9.108)^(0.5*1.43-0.5) * (Q2/(1 + Q2/9.108))^(1.43)*x^(1 - 1.43), # Hard Pomeron
          function(Q2, x) 0.5954 * (Q2/(1 + Q2/0.5894))^(1.08)*x^(1 - 1.08), # Soft Pomeron
          function(Q2, x) 1.154 * (Q2/(1 + Q2/0.2605))^(0.55)*x^(1 - 0.55))  # Meson trajectories

  # the fit of DL
  F2DL <- function(Q2, x) sum(sapply(fDL, function(f) f(Q2, x)))

  # Marko's F2 proposal
  #' @export
  F2M <- function(Q2, x, rho = 0.7792, g02 = 103.14, Qp = 0.433, z0 = 4.96) {
    Q <- sqrt(Q2)
    Qratio <- Q / Qp
    ## if you want the so-called 'error function'
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    ## (see Abramowitz and Stegun 29.2.29)
    ## and the so-called 'complementary error function'
    erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
    tau <- log(rho * Qratio / (2 * x))
    t1 <- exp(-(log(Qratio)^2) / (rho * tau)) / sqrt(tau)
    eta <- (log(Q * Qp * z0^2) + (rho * tau)) / sqrt(rho * tau)
    F <- 1 - 2 * sqrt(rho * pi * tau) * exp(eta^2) * erfc(eta)
    t2 <- F *  exp(-(log(1/(Q * Qp * z0 ^2))^2) / (rho * tau)) / sqrt(tau)

    f2m <- (g02) * rho^(1.5) * Qratio * exp((1 - rho) * tau)/ (32 * pi^(2.5))
    f2m <- f2m * (t1 + t2)
    return(f2m)
  }

  # default values are the one obtained by Alfonso
  F2A <- function(Q2, x, s0 = 0.213, s1 = -0.129, j0 = 1.14, j1 = 1.08, lambda = 0.36) {
  #F2A <- function(Q2, x, s0 = 0.2819369, s1 = -0.1665779, j0 = 1.1090260, j1 = 1.0512662, lambda = 0.2953093) {
    tau <- Q2 * x^lambda
    s0 * Q2 * tau^((1 - j0) / lambda) + s1 * Q2 * tau^((1 - j1) / lambda)
  }

  F2A2 <- function(Q2, x, j0=1.16, j1=1.08 , j2 = 1, lambda = 0.334 , s0 = 0.157 , s1 = -0.064, s2 = 0) {
    f0 <- s0 * x^(1-j0) * (Q2^(1 - (1/lambda) * tanh(Q2) * (j0-1)))
    f1 <- s1 * x^(1-j1) * (Q2^(1 - (1/lambda) * tanh(Q2) * (j1-1)))
    f2 <- s2 * x^(1-j2) * (Q2^(1 - (1/lambda) * tanh(Q2) * (j2-1)))
    f0 + f1 + f2
  }

  extractPotentialParameters <- function(...) do.call(kernel$extractPotentialParameters, list(...))

  f2 <- list(value = value,
             kernel = kernel,
             setNewPotential = setNewPotential,
             extractPotentialParameters = extractPotentialParameters,
             rss = rssFun,
             rssLambda = rssFunLambda,
             fi = fi,
             ld1 = r$ld1,       #|
             ld2 = r$ld2,       #| experimental
             ra = r$ra,         #| data
             wf = r$wf,         #|
             getTrs = function() list(tr1 = tr1, tr2 = tr2),  # get the last values of the mass spectrum computed
             getJs = function() unlist(Js),
             getFns = getFns,
             getFnsNMC = getFnsNMC,
             getAllFns = getAllFns,
             P13Bar = P13Bar,
             spectrum = spectrum,
             plotRegge = plotRegge,
             plotGlueballMasses = plotGlueballMasses,
             plotMesonsMasses = plotMesonsMasses,
             plotSpectrum = plotSpectrum,
             plotEffectiveExponent = plotEffectiveExponent,
             fDL = fDL,
             F2DL = F2DL,
             F2A = F2A,
             F2A2 = F2A2,
             F2M = F2M,
             plot = plotF2)
  class(f2) <- append(class(f2), 'F2')
  f2
}

#' @export
plot.F2 <- function(f2, ...) f2$plot(...)
