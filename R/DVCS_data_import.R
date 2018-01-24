# create a list with all the files inside the specific folder
filesAtFolder <- function(folder, n) unlist(lapply(n, function(n) paste0(folder, '/Table', n, '.csv')))


# load all the data and put it in the proper format
loadDVCS_sigma <- function() {
  # This function creates vectors with the values of W, Q2, sigma and dsigma
  # of the data used by Marko in his DVCS paper with Miguel
  # Getting data from 761990 folder
  # Getting data from Table1.csv file
  W_data      <- c(82,82,82,82) # W = sqrt(s) and |t| < 1
  Q2_data     <- c(8.75,15.5,25.0,55.0) # Q2 is in GeV**2
  sigma_data  <- c(3.59,1.38,0.58,0.13) # sigma is in NB
  # Summing the square of statistical and systematic errors assuming they are independent
  dsigma_data <- c(sqrt(0.21*0.21+0.41*0.41),sqrt(0.1*0.1+0.21+0.21),sqrt(0.09*0.09+0.09*0.09),sqrt(0.03*0.03+0.04*0.04))
  # Getting data from Table2.csv file
  W_data      <- c(W_data,c(45.0,70.0,90.0,110.0,130.0))
  Q2_data     <- c(Q2_data,c(8.0,8.0,8.0,8.0,8.0))
  sigma_data  <- c(sigma_data,c(2.91,3.96,4.78,5.55,6.56))
  dsigma_data <- c(dsigma_data,c(sqrt(0.2*0.2+0.25*0.25),sqrt(0.32*0.32+0.37*0.37),sqrt(0.41**2+0.57**2),sqrt(0.57**2+0.88**2),sqrt(1.17**2+1.77**2)))
  # Getting data from Table3.csv file
  W_data      <- c(W_data,c(45.0,70.0,90.0,110.0,130.0,45.0,70.0,90.0,110.0,130.0,45.0,70.0,90.0,110.0,130.0))
  Q2_data     <- c(Q2_data,c(8.0,8.0,8.0,8.0,8.0,15.5,15.5,15.5,15.5,15.5,25.0,25.0,25.0,25.0,25.0))
  sigma_data  <- c(sigma_data,c(2.6,3.15,5.25,5.11,5.88,0.94,1.54,0.95,1.69,2.06,0.35,0.36,0.83,0.9,0.9))
  dsigma_data <- c(dsigma_data,c(sqrt(0.24**2+0.24**2),sqrt(0.4**2+0.33**2),sqrt(0.55**2+0.55**2),sqrt(0.71**2+0.76**2),sqrt(1.89**2+1.26**2)))
  dsigma_data <- c(dsigma_data,c(sqrt(0.1**2+0.1**2),sqrt(0.17**2+0.14**2),sqrt(0.2**2+0.17**2),sqrt(0.31**2+0.33**2),sqrt(0.51**2+0.56**2)))
  dsigma_data <- c(dsigma_data,c(sqrt(0.13**2+0.07**2),sqrt(0.1**2+0.05**2),sqrt(0.18**2+0.09**2),sqrt(0.23**2+0.18**2),sqrt(0.36**2+0.32**2)))
  # Getting data from 805219 folder
  # Getting data from Table1.csv file
  W_data      <- c(W_data,c(104.0,104.0,104.0,104.0,104.0,104.0))
  Q2_data     <- c(Q2_data,c(3.25,7.5,12.5,20.0,32.5,70.0))
  sigma_data  <- c(sigma_data,c(21.28,5.87,3.27,1.23,0.55,0.16))
  dsigma_data <- c(dsigma_data,c(sqrt(0.92**2+max(1.02,1.34)**2),sqrt(0.42**+max(0.14,0.3)**2),sqrt(0.33**2+max(0.07,0.16)**2),sqrt(0.21**2+max(0.05,0.08)**2),sqrt(0.18**2+0.04**2),sqrt(0.07**2+0.02**2)))
  # Gettind data from Table2.csv file
  W_data      <- c(W_data,c(50.0,70.0,90.0,110.0,130.0,155.0))
  Q2_data     <- c(Q2_data,c(3.2,3.2,3.2,3.2,3.2,3.2))
  sigma_data  <- c(sigma_data,c(14.47,20.38,17.95,20.65,26.42,27.6))
  dsigma_data <- c(dsigma_data,c(sqrt(1.05**2+max(0.5,0.88)**2),sqrt(1.57**2+max(1.01,1.99)**2),sqrt(1.35**2+max(0.64,0.93)**2),sqrt(1.26**2+max(0.59,1.16)**2),sqrt(1.84**2+max(0.79,0.88)**2),sqrt(3.74**2+max(2.01,3.34)**2)))
  # Getting data from Table3.csv file
  W_data      <- c(W_data,c(52.5,77.5,102.5,127.5,155.0))
  Q2_data     <- c(Q2_data,c(2.4,2.4,2.4,2.4,2.4))
  sigma_data  <- c(sigma_data,c(27.06,22.36,26.49,35.94,35.72))
  dsigma_data <- c(dsigma_data,c(sqrt(3.44**2+max(4.37,4.42)**2),sqrt(3.11**2+max(3.4,1.73)**2),sqrt(1.89**2+max(0.89,1.44)**2),sqrt(2.63**2+max(1.81,1.89)**2),sqrt(9.47**2+max(3.01,2.94)**2)))
  W_data      <- c(W_data,c(155.0,155.0,155.0))
  Q2_data     <- c(Q2_data,c(6.2,9.9,18.0))
  sigma_data  <- c(sigma_data,c(16.93,6.15,2.21))
  dsigma_data <- c(dsigma_data,c(sqrt(2.43**2+max(1.37,1.4)**2),sqrt(1.67**2+0.51**2),sqrt(0.82**2+0.18**2)))
  # Getting data from 827347 folder
  # Getting data from Table1.csv file
  W_data      <- c(W_data,c(82.0,82.0,82.0,82.0))
  Q2_data     <- c(Q2_data,c(8.75,15.5,25.0,55.0))
  sigma_data  <- c(sigma_data,c(3.87,1.46,0.55,0.16))
  dsigma_data <- c(dsigma_data,c(sqrt(0.15**2+0.41**2),sqrt(0.07**2+0.18**2),sqrt(0.07**2+0.08**2),sqrt(0.02**2+0.03)))
  # Getting data from Table2.csv file
  W_data      <- c(W_data,c(45.0,70.0,90.0,110.0,130.0))
  Q2_data     <- c(Q2_data,c(10.0,10.,10.,10.,10.))
  sigma_data  <- c(sigma_data,c(2.23,2.92,3.63,3.71,4.37))
  dsigma_data <- c(dsigma_data,c(sqrt(0.11**2+0.19**2),sqrt(0.16**2+0.27**2),sqrt(0.22**2+0.4**2),sqrt(0.29**2+0.61**2),sqrt(0.6**2+1.16**2)))
  # Getting data from Table3.csv file
  W_data      <- c(W_data,c(45.0,70.0,90.0,110.0,130.0))
  Q2_data     <- c(Q2_data,c(8.0,8.0,8.0,8.0,8.0))
  sigma_data  <- c(sigma_data,c(3.06,3.54,4.93,5.16,5.62))
  dsigma_data <- c(dsigma_data,c(sqrt(0.18**2+0.25**2),sqrt(0.29**2+0.34**2),sqrt(0.39**2+0.52**2),sqrt(0.51**2+0.74**2),sqrt(1.34**2+1.19**2)))

  return (data.frame(W_data, Q2_data, sigma_data, dsigma_data))
  }
