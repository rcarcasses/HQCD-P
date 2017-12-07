
#' @export
plotF2vsPredictor <- function(F2pred = F2DL, rate = 1) {
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)
	# remove all the high x together with some points with "weird" F2
	data <- ncepp[ncepp$x < 0.01 & ncepp$F2 < 5,]
  # this list contains all the different Q2 entries
	Q2s = unique(data[, c("Q2")])

	f2x <- data[,c("x","F2")]
  # initialize the plot
  #par(xpd=NA)
	plot.new()
  plot(f2x$x, f2x$F2, type="n", log = 'xy',
       main=expression("F"[2]),
       xlab = expression("x"), ylab = expression("F"[2]))

	# first let's create a bunch of colors to differentiate the graphs
	cl <- topo.colors(as.integer(length(Q2s) / 2))
	cl <- c(cl, cl)

	for(i in 1:as.integer(length(Q2s) / rate))
	{
	  ii <- i * rate
		# now we need to extract the columns that we are interested for a given value of Q2
		f2x <- data[data$Q2 == Q2s[ii],][,c("x","F2")]
		s_r <- data[data$Q2 == Q2s[ii],][,c("s_r")]
		tot <- data[data$Q2 == Q2s[ii],][,c("tot")]
	  err <- s_r * tot / 100
	  # Draw the fit on the plot by getting the prediction from the fit at 200 x-coordinates across the range of xdata
	  len <- 200
		x <- seq(0.5 * min(f2x$x), 5 * max(f2x$x), len = len)
		#x <- seq(min(f2x$x), max(f2x$x), len = 200)
		fitPlot <- F2pred(Q2s[ii], x)
		lines(x, fitPlot, col = cl[ii])

		# plot the dots
		lines(f2x$x, f2x$F2, type = "p", col= cl[ii], pch=19)
		# plot the error bars
		arrows(f2x$x, f2x$F2 - err, f2x$x, f2x$F2 + err,
		       length = 0.05, angle = 90, code = 3, col = cl[ii])

		# put the Q2 value at the end of the line
		text(x = 0.9 * x[1], y = 1.1 * fitPlot[1], labels = paste(Q2s[ii]), col = cl[ii])
	}
}


# This function just look for the best parameters that give
# the closest to the f1 of Donnachie and Landshoff (DL)
#' @export
fitfsToDL <- function(n = 1) {
  q <- Q2s[1:32]  # we want to fit the tail
  fn <- function(pars) {
    kernel <<- do.call(findKernel, pars)
    pars[3] * Fn(q, kernel$js[[n + 1]], kernel$wfs[[n + 1]])
  }

  fs <- fDL[[n + 1]]
  residFun <- function(pars) {
    val <- log((fs(q, x = 1) / fn(pars))^2)          # x = 1 gives the Q coefficient only
    val <- val^2 + 1000 * (kernel$js[[2]] - 1.09)^2  # force the intercept to be the SP one
    val <- sum(val)
    cat('\n', pars, ' val=',val, '| j1 =', kernel$js[[2]], ' | ')
    return(val)
  }

  op <- optim(list(b0 = -10, c0 = 4, s0 = 0.1),
        fn = residFun,
        #method = c('SANN'),
        upper = c(b0 = 0, c0 = 5, s0 = 1),
        lower = c(b0 = -100, c0 = -5, s0 = 0.001),
        control = list(trace = 10, REPORT = 500, maxit = 1e6))

  plotBestFitf1DL <- function() {
    plot.new()
    fDL <- fs(Q2 = q, x = 1)
    plot(q, fDL, type = 'o', log = 'xy', col = 'blue',
         ylab = expression('f'[n]),
         xlab = expression('Q'^2),
         ylim = c(min(fDL), 1.2 * max(fDL)))
    lines(q, fn(op$par), type = 'o')
  }

  list(optim = op, residFun = residFun, plot = plotBestFitf1DL)
}

# Find the best match for f0 and f1 of Donnachie and Landshoff (DL)
#' @export
fitToDL <- function() {
  q <- Q2s[1:32]  # we want to fit the tail
  residFunOptim <- function(pars) {
    kernel <- findKernel(b0 = pars[1], b1 = pars[2], c0 = pars[3], c1 = pars[4], ls = 0.136)
    j0  <- kernel$js[[1]]
    wf0 <- kernel$wfs[[1]]
    j1  <- kernel$js[[2]]
    wf1 <- kernel$wfs[[2]]
    s0 <- pars[5]
    s1 <- pars[6]
    fn0 <-  s0 * Fn(q, j0, wf0)
    fn1 <-  s1 * Fn(q, j1, wf1)
    val <- log((f0(q, x = 1)/fn0)^2)^2 + log((f1(q, x = 1)/fn1)^2)^2    # x = 1 gives the Q coefficient only
    val <- sum(val ^2)
    cat('\n', pars, ' val=',val)
    return(val)
  }

  optim(list(b0 = 0, b1 = 0, c0 = 0, c1 = 0, s0 = 0.01, s1 = -17),
        fn = residFunOptim,
        #method = c('SANN'),
        upper = c(b0 = 10, b1 = 10, c0 = 8, c1 = 8, s0 = 1, s1 = -1),
        lower = c(b0 = -10, b1 = -10, c0 = -8, c1 = -8, s0 = 0.00001, s1 = -100),
        control = list(trace = 10, REPORT = 500, maxit = 1e6))
}



#' @export
fitF2 <- function() {	# read the HERA nce+p data
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

	# remove all the high x together with some points with "weird" F2
	data <- ncepp[ncepp$x < 0.01 & ncepp$F2 < 5,]
	err <- data$s_r * data$tot / 100
	w <- 1 / (err^2)

	# residual function
	residFun <- function(pars, i) {
	  val <- data$F2[i] - max(F2(data$Q2[i], data$x[i],
	                         b0 = pars[1],
	                         b1 = pars[2],
	                         c0 = 0, #pars[3],
	                         c1 = 0, #pars[4],
	                         ls = pars[3],
	                         s0 = pars[4], s1 = pars[5], s2 = pars[6]), 1e-15)
	  val <- val / err[i]
	  val^2
	}

	residFunOptim <- function(pars) {
	  cat('\n', pars)
	  i <- c(1:length(data$F2))
	  sum(residFun(pars), i)
	}

	# result was 3.315313 -3.324192 2.661214 -2.662756 0.23
  optim(list(b0 = 0, b1 = 0, ls = 0.136, s0 = 0.017, s1 = -58, s2 = 1),
  #optim(list(ls = 0.136, s0 = 1.7, s1 = -5),
             fn = residFunOptim,
             #method = c('SANN'),
             upper = c(b0 = 20, b1 = 20, ls = 0.15, s0 = 10, s1 = -1, s2 = 100),
             lower = c(b0 = -20, b1 = -20, ls = 0.13, s0 = 0.001, s1 = -100, s2 = -100),
             control = list(trace = 10, REPORT = 500, maxit = 1e6))


	#return(nlsLM(F2 ~ F2(Q2, x, b0, b1, c0, c1, ls, s0, s1),
  #          data = data,
  #           weights = w,
             #start = list(b0 = 286 , b1 = 311 , c0 = -36.8 , c1 = -40 , ls = 0.136 , s0 = 0.01768549, s1 = -58.7),
  #           upper = c(b0 = 8, b1 = 8, c0 = 8, c1 = 8, ls = 0.3, s0 = 10, s1 = -1),
  #           start = c(b0 = 0, b1 = 0, c0 = 0, c1 = 0, ls = 0.136549, s0 = 0.01, s1 = -50),
  #           lower = c(b0 = -8, b1 = -8, c0 = -8, c1 = -8, ls = 0.117, s0 = 0.000000001, s1 = -100),
  #           trace = TRUE))
}

#' @export
fitLessF2 <- function() {	# read the HERA nce+p data
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

	# remove all the high x together with some points with "weird" F2
	data <- ncepp[ncepp$x < 0.01 & ncepp$F2 < 5,]
	err <- data$s_r * data$tot / 100
	w <- 1 / (err^2)
  return(nlsLM(F2 ~ F2(Q2, x, b0, b1, s0, s1),
             data = data,
             weights = w,
             #start = list(b0 = 286 , b1 = 311 , c0 = -36.8 , c1 = -40 , ls = 0.136 , s0 = 0.01768549, s1 = -58.7),
             start = list(b0 = 80 , b1 = 80 , s0 = 0.024, s1 = -50),
             upper = c(b0 = 180, b1 = 180, s0 = 0.1, s1 = -10),
             lower = c(b0 = 10, b1 = 10, s0 = 0.01, s1 = -100),
             trace = TRUE))
}
#' @export
fitBruteForceF2 <- function() {	# read the HERA nce+p data
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

	# remove all the high x together with some points with "weird" F2
	data <- ncepp[ncepp$x < 0.005 & ncepp$F2 < 5,]
	st <- expand.grid(b0 = seq(280, 290, len = 3),
	                  b1 = seq(300, 320, len = 3),
	                  c0 = seq(-40, -32, len = 4),
	                  c1 = seq(-45, -35, len = 4),
	                  s0 = seq(0.01, 0.02, len = 4),
	                  s1 = seq(-70, -50, len = 4),
	                  ls = seq(0.136, 0.14, len = 1))
  return(nls2(F2 ~ F2(Q2, x, b0, b1, c0, c1, ls, s0, s1),
             data = data,
             start = st,
             algorithm = 'brute-force',
             trace = TRUE))
}

#' @export
fitGeometricScalingF2 <- function(maxX = 0.01, maxF2 = 5, maxQ2 = 90, minQ2 = 0.25) {
  # read the HERA nce+p data
	ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

  data <- ncepp[ncepp$x < maxX & ncepp$F2 < maxF2 & ncepp$Q2 <= maxQ2 & ncepp$Q2 >= minQ2,]
  data <- data[,c('F2', 'Q2', 'x', 's_r', 'tot')]
	err <- data$s_r * data$tot / 100
	data$s_r <- NULL
	data$tot <- NULL
	data$err <- err

	f2 <- F2BCD('')

	sift.nls(fn = f2$F2A2, data, traces = FALSE,
	          start = list(j0 = 1.16, j1 = 1.08 , j2 = 1, lambda = 0.334 , s0 = 0.157 , s1 = -0.064, s2 = 0),
	          options = list(yName = 'F2', errName = 'err', normal = FALSE))
}

loadFitData <- function(maxX = 0.01, maxF2 = 5, maxQ2 = 990, minQ2 = -0.25) {
  # read the HERA nce+p data
  ncepp <- read.table("d09-158.nce+p.txt", header = TRUE)

  data <- ncepp[ncepp$x < maxX & ncepp$F2 < maxF2 & ncepp$Q2 <= maxQ2 & ncepp$Q2 >= minQ2,]
  data <- data[,c('F2', 'Q2', 'x', 's_r', 'tot')]
  err  <- data$s_r * data$tot / 100
  data$s_r <- NULL
  data$tot <- NULL
  data$err <- err
  data
}

#' @export
fitMarko <- function() {
  data <- loadFitData()
  f2 <- F2BCDQ('')

  sift.nls(fn = f2$F2M, data, traces = TRUE,
           start = list(rho = 0.7792, g02 = 103.14, Qp = 0.433, z0 = 4.96),
           options = list(yName = 'F2', errName = 'err', normal = TRUE))
}
