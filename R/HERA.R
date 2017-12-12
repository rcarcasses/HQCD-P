
#' @export
heraEnv <- new.env()
# the following are the intercepts currently being obtained
#' @export
jn <- c(1.4, 1.096, 0.55)
#' @export
e0 <- 1 - jn

#' @export
Q2s <- NULL

#' @export
getHERAF2 <- function(maxX = 0.01, maxF2 = 5, maxQ2 = 1110, minQ2 = 0.1) {
  flog.debug(paste('Loading HERA data with maxX', maxX, ' and ', minQ2, ' <= Q2 <=', maxQ2))
  # read the HERA nce+p data
  nceppPath <- system.file('extdata', 'd09-158.nce+p.txt', package = 'HQCDP')
  flog.debug(paste('[HERA] Loading DIS HERA data from ', nceppPath))
  ncepp <- read.table(nceppPath, header = TRUE)

  # remove all the high x together with some points with "weird" F2
  data <- ncepp[ncepp$x < maxX & ncepp$F2 < maxF2 & ncepp$Q2 <= maxQ2 & ncepp$Q2 >= minQ2,]
  #data <- ncepp[ncepp$x < maxX & ncepp$F2 < maxF2 & ncepp$Q2 < maxQ2 & ncepp$Q2 > 7,]
  flog.debug(paste('[HERA] Q2 range [', min(data$Q2),',', max(data$Q2), '], number of data points', length(data$Q2)))
  f2x <- data[,c('F2', 'Q2', 'x', 's_r', 'tot')]
	# this list contains all the different Q2 entries
	Q2s <- unique(data[, c("Q2")])
  # let's also compute here what is the effective intercept by fitting the data using f(Q)x^e(Q)
	Q2L <- length(Q2s)
	W <- list()
	X <- list()
	eff <- data.frame(Q2 = numeric(Q2L), ep = numeric(Q2L), epErr = numeric(Q2L), minX = numeric(Q2L), maxX = numeric(Q2L))
	for(i in 1:Q2L) {
	  Q2 <- Q2s[i]
	  # now we need to extract the columns that we are interested for a given value of Q2
	  f2xFit <- data[data$Q2 == Q2,][,c("x","F2")]
	  s_r <- data[data$Q2 == Q2,][,c("s_r")]
	  tot <- data[data$Q2 == Q2,][,c("tot")]
	  err <- s_r * tot / 100
	  w <- 1 / (err^2)

	  # skip those data which are too small
	  if(length(f2xFit$x) > 2) {
  	  # now let's try to fit it
      fit <- lm( log(F2) ~ log(x),# p0 * x^(-ep),
                    data = f2xFit,
                    weights = w)#,
                    #start = list(p0 = 1, ep = 1))
  	  s <- summary(fit)$coefficients
  	  eff$Q2[i]    <- Q2
  	  eff$ep[i]    <- -s['log(x)', 'Estimate']
  	  # let's use 3 sigma as the error uncertainty
  	  eff$epErr[i] <- 3 * s['log(x)', 'Std. Error']
  	  #cat('e = ', eff$ep[i], 'de = ', eff$epErr[i], '\n')
  	  eff$minX[i] <- min(f2xFit$x)
  	  eff$maxX[i] <- max(f2xFit$x)
	  }
	  W[[i]]  <- w
	  X[[i]]  <- f2xFit$x
	}

	eff <- eff[eff$Q2 !=0,]

  list(F2 = f2x$F2, Q2 = f2x$Q2, x = f2x$x, err = f2x$s_r * f2x$tot / 100, eff = eff, weights = W, xs = X, Q2s = Q2s)
}

#' @export
plotHERARangeXvsQ <- function(maxX = 0.01, maxF2 = 5, maxQ2 = 1110, minQ2 = 0.1) {
  ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)
  dataAbove <- ncepp[ncepp$x < maxX & ncepp$F2 < maxF2 & ncepp$Q2 <= maxQ2 & ncepp$Q2 >= minQ2,]
  dataBelow <- ncepp[ncepp$x > maxX & ncepp$F2 < maxF2 & ncepp$Q2 <= maxQ2 & ncepp$Q2 >= minQ2,]
  cat('data below x <', maxX, length(dataAbove$x), ', total points ', length(ncepp$x),'\n')
  plot(dataAbove$Q2, 1/dataAbove$x, log = 'xy', type = 'p', pch = 16, cex = 0.5,
       xlab = expression(Q^2), ylab = expression(1/x), ylim = c(1, 1e6))
  abline(h = c(1, 1e2, 1e4, 1e6), v = c(0.5, 5, 50, 500), col = "lightgray", lty = 'dashed' )
  abline(h = c(1e2), col = "black", lwd = 2 )
  lines(dataBelow$Q2, 1/dataBelow$x, type = 'p', cex = 0.5)
}

#' @export
loadHERA <- function(useIHQCD = TRUE, js = NULL, plotF2 = TRUE) {
  A <- function(z) {return(-log(z / 1))}
  if(useIHQCD)
	  A <- splinefun(get('z', envir = ihqcdEnv), get('A', envir = ihqcdEnv))

  if(is.null(js))
    js = jn

	# read the HERA nce+p data
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

	# remove all the high x together with some points with "weird" F2
	data <- ncepp[ncepp$x < 0.01 & ncepp$F2 < 5 & ncepp$Q2 > 0.10 & ncepp$Q2 < 2200,]

	# this list contains all the different Q2 entries
	Q2s <- unique(data[, c("Q2")])

	f2x <- data[,c("x","F2")]
	if(plotF2) {
	  # initialize the plot
	  plot(log(f2x$x, 10), f2x$F2, type="n",
	       main=expression("F"[2]),
	       xlab = expression("log"[10]*"x"), ylab = expression("F"[2]),
	       xlim = c(-6,-2), ylim = c(0,1.5))
	}


	# first let's create a bunch of colors to differentiate the graphs
	cl <- rainbow(length(Q2s))

	paramVsQ2 <- data.frame( "Q2" = numeric(0),
              					   "p0" = numeric(0),
              					   "p1" = numeric(0),
              					   # "p2" = numeric(0),
              					   "z"  = numeric(0),
              					   stringsAsFactors=FALSE)
	rss <- 0

	for(i in 2:(length(Q2s)))
	{
		# now we need to extract the columns that we are interested for a given value of Q2
		f2x <- data[data$Q2 == Q2s[i],][,c("x","F2")]
		s_r <- data[data$Q2 == Q2s[i],][,c("s_r")]
		tot <- data[data$Q2 == Q2s[i],][,c("tot")]
	  err <- s_r * tot / 100
	  w <- 1 / (err^2)

		# skip those data which are too small
		if(length(f2x$F2) < 3)
      next

		# now let's try to fit it
		tryCatch({
			fit <- nlsLM( F2 ~ p0 * x^(1 - js[1]) +  p1 * x^(1 - js[2]),# + p2 * x^(1 - js[3]),
				data = f2x,
				weights = w,
				start = list(p0 = 1, p1 = 1))#, p2 = 1))

			# sum the residuals to have a control of the quality of the fit
			rss <- rss + sum(residuals(fit)^2)

			# get the parameters for the given value of Q2
			# to find z we need to invert Q = exp(A(z))
			qvsz <- function(z) { return(exp(A(z)) - sqrt(Q2s[i])) }
			zSol = uniroot(qvsz, c(0, 7), tol = 1e-9)
      # print(paste(zSol$root, " in AdS should be ", 1/sqrt(Q2s[i])))
			row = union(union(Q2s[i],fit$m$getAllPars()), zSol$root)
			paramVsQ2[nrow(paramVsQ2) + 1, ] <- row

			if(plotF2) {
			  # Draw the fit on the plot by getting the prediction from the fit at 200 x-coordinates across the range of xdata
  			fitPlot = data.frame(x = seq(min(f2x$x), max(f2x$x), len = 200))
  			lines(log(fitPlot$x, 10), predict(fit, newdata = fitPlot), col = cl[i])

  			# plot the dots
  			lines(log(f2x$x, 10), f2x$F2, type = "p", col= cl[i])
			}
		},
			error = function(e){
				print(paste("Unable to fit for Q2=", Q2s[i], " data ", f2x))
				print(paste("-> ERROR :",conditionMessage(e), "\n"))
		})
	}

	# now make these available through the file environment
	assign("paramVsQ2", paramVsQ2, envir = heraEnv)
	assign("js", js, envir = heraEnv)
	assign("A", A, envir = heraEnv)
	assign("rss", rss, envir = heraEnv)
	cat('rss for ', js, '=', rss, '\n')
	assign("data", data, envir = heraEnv)
	z          <- paramVsQ2$z    # beware this is not the z of IHQCD, this is the z for each Q2, so is a different list.
	Asfun      <- splinefun(get('z', envir = ihqcdEnv), get('As', envir = ihqcdEnv))
	lambdafun  <- splinefun(get('z', envir = ihqcdEnv), get('lambda', envir = ihqcdEnv))
	ff         <- lapply(js, function(J) z^(-2 * J) * exp((-J + 0.5) * Asfun(z)) * lambdafun(z))
  # reference: F2 ~ Q2^J P13 exp((-J + 0.5) * As) exp(phi)
	# for phi = cte ~ z^(-2J) z^(J - 0.5) ~ z^(-J - 0.5)
	# P13 ~ delta(z - 1/Q)
	assign("phi0z", paramVsQ2$p0 / ff[[1]], envir = heraEnv)
	assign("phi1z", paramVsQ2$p1 / ff[[2]], envir = heraEnv)
	#assign("phi2z", paramVsQ2$p2 / ff[[3]], envir = heraEnv)
	#computeU(A)
	return(list( z = paramVsQ2$z, paramVsQ2 = paramVsQ2, data = data, rss = rss))
}

#' @export
showBestJs <- function(nX = 10, nY = 10) {
  j0s <- seq(0.9, 1.2, len = nX)
  j1s <- seq(1.21, 1.6, len = nY)

  minRss <- list(j0 = 0, j1 = 0, rss = 1e10)
  rss <- matrix(nrow = nX, ncol = nY)
  for (i in 1:length(j0s)) {
    for (j in 1:length(j1s)) {
      hera <- loadHERA(js = c(j0s[i], j1s[j]), plotF2 = FALSE)
      rss[[i, j]] = hera$rss
      if(hera$rss < minRss$rss){
        minRss$rss <- hera$rss
        minRss$j0 <- j0s[i]
        minRss$j1 <- j1s[j]
      }
    }
  }

  plotData <- list( x = j0s, y = j1s, z = rss)
  contour(plotData, levels = c(0.03, 0.04, 0.05, 0.06, 0.07, 0.1, 0.2, 0.3, 1, 3),
          xlab = expression('j'[0]),
          ylab = expression('j'[1]),
          main = 'Residues squared sum')
  cat('Minimun rss =', minRss$rss,' found for j0 =', minRss$j0, ' j1 =', minRss$j1, '\n')
  points(x = minRss$j0, y = minRss$j1)
  return(minRss)
}

#' @export
computeU <- function(A = function(z) {return(-log(z / 1))}) {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z

  integrand <- function(z) { return(-exp(A(z))) }
  u <- lapply(paramVsQ2$z, function(z) { return (integrate(integrand, 10, z)$value) })

  phi1u = splinefun(u, exp(-0.5 * A(z) * (0.5 + j0[1])) * paramVsQ2$p1)
  phi2u = splinefun(u, exp(-0.5 * A(z) * (0.5 + j0[2])) * paramVsQ2$p2)

  # now make these available through the file environment
	assign("phi1u", phi1u, envir = heraEnv)
	assign("phi2u", phi2u, envir = heraEnv)
  assign("u", u, envir = heraEnv)
}

#' @export
reconstructVu <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z
  u = get("u", envir = heraEnv)
  A = get("A", envir = heraEnv)
  phi1 = get("phi1u", envir = heraEnv)
  phi2 = get("phi2u", envir = heraEnv)

  # first reconstruction
  Vu1 = lapply(u, function(u) {
      V = (phi1(u, deriv = 2) / phi1(u)) - j0[1]
      return(V)
    })

 	plot.new()
	plot(u, Vu1, type="p", main="Reconstructing with first wavefunction (in blue)", xlab = "u", ylab = "V(u)", ylim = c(-500, 500))
  lines(u, Vu1, type="o")
  lines(u, 3000 * phi1(u), type="l", col = "blue")

  # second reconstruction
  Vu2 = lapply(u, function(u) {
      V = (phi2(u, deriv = 2) / phi2(u)) - j0[2]
      return(V)
    })

 	plot.new()
	plot(u, Vu2, type="p", main="Reconstructing with second wavefunction (in blue)", xlab = "u", ylab = "V(u)")
  lines(u, Vu2, type="o")
  lines(u, 1000 * phi2(u), type="l", col = "blue")

}

derivative <- function(f, x, ..., order = 1, delta = 0.1, sig = 6) {
  # Numerically computes the specified order derivative of f at x
  vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
  grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
  vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
  for (i in 2:(order + 1)) {
    for (j in 1:(order - i + 2)) {
      stepsize <- grid[i + j - 1] - grid[i + j - 2]
      vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
    }
  }
  return(signif(vals[order + 1, 1], sig))
}

#' @export
plotUvsZ <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z
  u = get("u", envir = heraEnv)
 	plot.new()
	plot(z, u, type="p", main=expression(paste("u vs. z")), xlab = "z", ylab = "u")
	lines(z, -log(z / 10), type = "l", col = "blue")
	legend("topright", expression(paste("Line represents AdS")), col=c("black", "blue"))
}

#' @export
plotZvsQ <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z
	plot.new()
	plot(sqrt(paramVsQ2$Q2), paramVsQ2$z, type="p", main=expression(paste("z vs. Q")), xlab = "Q", ylab = "z")
	lines(sqrt(paramVsQ2$Q2), 1/sqrt(paramVsQ2$Q2), type = "l", col = "blue")
	legend("topright", expression(paste("Line represents AdS")), col=c("black", "blue"))
}

# plot of the parameters as a function of Q2
#' @export
plotPvsQ2 <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  Q2 = paramVsQ2$Q2
	# yLim <- c(0.9 * min(paramVsQ2$p0, paramVsQ2$p1), 1.1 * max(paramVsQ2$p0, paramVsQ2$p1))
	yLim <- c(-0.5, 1.1 * max(paramVsQ2$p0, paramVsQ2$p1))
	plot.new()
	plot(Q2, paramVsQ2$p0, type = "o",
			xlab = expression(paste(Q^2,(GeV^2))), log = 'x',
			ylab = expression(paste("f"[0],",  f"[1])),
			xlim = c(1e-1, max(Q2)),
			ylim = yLim)
	#lines(Q2, paramVsQ2$p0, type="o", col="red")
	lines(Q2, paramVsQ2$p1, type="o", col="blue")
	#lines(z, paramVsQ2$p2 / 5, type="o", col="green")
	#abline(h = 0, col = "gray60")
  abline(h = seq(yLim[1], yLim[2], len = 5), v = seq(0.2, 250, len = 5), col = "lightgray", lty = 3)
}

# let's try to show how the functions should look like in z
#' @export
plotPvsZ <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z
	plot.new()
	plot(z, type = "n",
			xlab = "z", log = 'x',
			ylab = expression(paste("f"[1],",  f"[2])),
			xlim = c(7e-2, 2.1),
			ylim = c(-0.5, 0.6))
	lines(z, paramVsQ2$p0, type="o", col="red")
	lines(z, paramVsQ2$p1, type="o", col="blue")
	#lines(z, paramVsQ2$p2 / 5, type="o", col="green")
	axis(3, at = z, labels = paramVsQ2$Q2, col.axis="blue", cex.axis=0.7, tck=-0.01)
	mtext(expression(paste(Q^2,(GeV^2))), side=3, at=c(2.2), col="blue", line=1.5)
}

# now the p1 and p2 functions are related non trivially to the wavefunctions,
# let's try to guess how the wavefunctions should actually look like
#' @export
plotPhivsZ <- function(maxValue = 1) {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  z = paramVsQ2$z
  phi0z = get("phi0z", envir = heraEnv)
  phi1z = get("phi1z", envir = heraEnv)
  # phi2z = get("phi2z", envir = heraEnv)
  if(maxValue > 0) {
    phi0z = (maxValue / max(phi0z)) * phi0z
    phi1z = (maxValue / max(phi1z)) * phi1z
    # phi2z = (maxValue / phi2z[2]) * phi2z
  }
  rangeY = 1#max(max(phi0z, phi1z, phi2z))
	plot.new()
	plot(z, type = "n",
		 xlab = "z", #log = 'x',
		 ylab = expression(paste(psi[0],", ", psi[1])),
		 xlim = c(0.05, 1.1 * max(z)),
		 ylim = c(-0.2 * rangeY, 1.1 * rangeY))

	phis = list(phi0z, phi1z)
	cols = c('red', 'blue')
	mapply(function(phi, color) {
	  lines(z, phi, type="p", col=color)
	  lines(z, phi, type="l", col=alpha(color, 0.3))
	}, phis, cols)

	axis(3, at = z, labels = paramVsQ2$Q2, cex.axis=0.7, tck=-0.01)
	mtext(expression(paste(Q^2,(GeV^2))), side=3, at=c(1.7), line=1.5)
	abline(h = 0, col = "gray60")
  abline(h = seq(-0.2, 1, by = 0.2), v = seq(0.5, 2.5, by = 0.5), col = "lightgray", lty = 3)
	# lines(z, phi2z /20, type="o", col="green")
}

# What about the last in the u variable?
#' @export
plotPhivsU <- function() {
  paramVsQ2 = get("paramVsQ2", envir = heraEnv)
  u = get("u", envir = heraEnv)
  A = get("A", envir = heraEnv)
  z = paramVsQ2$z
  j0 = get("j0", envir = heraEnv)
  phi1 = exp(-0.5 * A(z) * (0.5 + j0[1])) * paramVsQ2$p1
  phi2 = exp(-0.5 * A(z) * (0.5 + j0[2])) * paramVsQ2$p2
	plot.new()
	plot(z, type = "n",
		 main = expression(paste(phi[0]," and ", phi[1]," vs u")),
		 xlab = "u",
		 ylab = expression(paste(phi[0],", ", phi[1])),
		 xlim = c(0, 7),
		 ylim = c(0, 0.3))
	lines(u, phi1, type="o", col="red")
	lines(u, phi2, type="o", col="blue")
	#axis(3, at = z, labels = paramVsQ2$Q2, col.axis="blue", cex.axis=0.7, tck=-0.01)
}
