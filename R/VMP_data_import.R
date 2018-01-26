
# Masses of the vector mesons. Taken from arxiv:1307.0009
m_Jpsi   <- 3.096916
m_rho0   <- 0.77549
m_phi    <- 1.019445
m_omega  <- 0.78265

#' @export
loadData.VMPSigma <- function(vmp) {
  # This function creates vectors with the values of m,  W,  Q2,  sigma and deltaSigma
  # of the data used by arxiv:1307.0009
  # Getting data from 694372 folder
  # Getting data from Table1.csv file
  m           <- c(m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi)
  W           <- c(90., 90., 90., 90., 90., 90., 90., 90.)
  Q2          <- c(0.05, 2.5, 4.0, 6.3, 10.0, 15.8, 25.0, 47.3)
  sigma       <- c(73.1, 37.3, 31.7, 21.8, 13.3, 7.53, 3.43, 0.6)
  deltaSigma <- c(sqrt(1.1**2+6.4**2), sqrt(3.9**2+3.6**2), sqrt(2.7**2+3.0**2), sqrt(2.4**2+2.1**2), sqrt(1.8**2+1.3**2), sqrt(1.24**2+0.72**2), sqrt(0.81**2+0.33**2), sqrt(0.24**2+0.06**2))
  # Getting data from Table2.csv file
  m           <- c(m, c(m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi))
  W           <- c(W, c(44.8, 54.8, 64.8, 74.8, 84.9, 94.9, 104.9, 119.5, 144.1, 144.9, 162.5, 177.3, 194.8, 219.6, 219.6, 244.8, 267.2, 292.3))
  Q2          <- c(Q2, c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))
  sigma       <- c(sigma, c(46.0, 48.5, 59.7, 62.7, 72.6, 78.6, 82.6, 91.5, 98.3, 98.6, 114.0, 126.0, 143.0, 187.0, 133.0, 171.0, 173.0, 194.0))
  deltaSigma <- c(deltaSigma, c(sqrt(2.4**2+4.0**2), sqrt(2.3**2+4.3**2), sqrt(2.8**2+5.3**2), sqrt(3.2**2+5.5**2), sqrt(3.4**2+6.4**2), sqrt(3.7**2+6.9**2), sqrt(4.0**2+7.3**2)))
  deltaSigma <- c(deltaSigma, c(sqrt(3.5**2+8.1**2), sqrt(4.4**2+8.7**2), sqrt(6.6**2+9.6**2), sqrt(8.0**2+11.0**2), sqrt(8.0**2+12.0**2), sqrt(10.0**2+15.0**2), sqrt(14.0**2+25.0**2)))
  deltaSigma <- c(deltaSigma, c(sqrt(10.0**2+18.0**2), sqrt(13**2+17**2), sqrt(13.0**2+18**2), sqrt(19.0**2+23.0**2)))
  # Getting data from Table3.csv file
  m           <- c(m, c(m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi))
  W           <- c(W, c(53.3, 83.9, 114.1, 144.2))
  Q2          <- c(Q2, c(3.2, 3.2, 3.2, 3.2))
  sigma       <- c(sigma, c(25.1, 30.0, 41.5, 45.0))
  deltaSigma <- c(deltaSigma, c(sqrt(2.9**2+2.4**2), sqrt(3.4**2+2.9**2), sqrt(5.1**2+4.0**2), sqrt(8.8**2+4.5**2)))
  # Getting data from Table4.csv file
  m           <- c(m, c(m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi))
  W           <- c(W, c(53.3, 83.9, 114.1, 144.2))
  Q2          <- c(Q2, c(7.0, 7.0, 7.0, 7.0))
  sigma       <- c(sigma, c(12.9, 14.5, 24.7, 24.1))
  deltaSigma <- c(deltaSigma, c(sqrt(2.5**2+1.2**2), sqrt(2.5**2+1.4**2), sqrt(4.1**2+2.4**2), sqrt(6.2**2+2.5**2)))
  # Getting data from Table5.csv file
  m           <- c(m, c(m_Jpsi, m_Jpsi, m_Jpsi, m_Jpsi))
  W           <- c(W, c(53.4, 83.9, 114.1, 144.2))
  Q2          <- c(Q2, c(22.4, 22.4, 22.4, 22.4))
  sigma       <- c(sigma, c(3.19, 4.04, 5.29, 6.1))
  deltaSigma <- c(deltaSigma, c(sqrt(0.69**2+0.31**2), sqrt(0.7**2+0.39**2), sqrt(1.0**2+0.5**2), sqrt(1.6**2+0.6**2)))
  # Getting data from 835609 folder
  # Getting data from Table1.csv file
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0))
  W           <- c(W, c(75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0))
  Q2          <- c(Q2, c(2.65, 2.95, 3.30, 3.70, 4.15, 4.65, 5.20, 5.85, 6.55, 7.35, 8.20, 9.20, 10.3, 11.5, 12.9))
  Q2          <- c(Q2, c(14.5, 16.5, 18.8, 21.7, 25.0, 29.3, 35.0, 46.0))
  sigma       <- c(sigma, c(563.0, 423.0, 383.0, 295.0, 232.0, 198.0, 154.0, 131.0, 102.0, 79.1, 56.5, 53.0, 39.3))
  sigma       <- c(sigma, c(30.8, 25.8, 15.7, 12.5, 9.22, 5.99, 3.54, 2.24, 1.68, 0.742))
  deltaSigma  <- c(deltaSigma, c(sqrt(32**2+59**2), sqrt(24**2+49**2), sqrt(20**2+43**2), sqrt(17**2+33**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(16**2+27**2), sqrt(13**2+20**2), sqrt(9**2+17**2), sqrt(9**2+12**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(8**2+8**2), sqrt(6.4**2+6.6**2), sqrt(4.3**2+5.5**2), sqrt(3.9**2+4.2**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(3.9**2+4.2**2), sqrt(2.9**2+2.8**2), sqrt(2.2**2+2.7**2), sqrt(1.3**2+1.7**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(1.0**2+1.3**2), sqrt(0.82**2+1.11**2), sqrt(0.58**2+0.67**2), sqrt(0.42**2+0.39**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(0.31**2+0.31**2), sqrt(0.27**2+0.24**2), sqrt(0.105**2+0.101**2)))
  # Getting data from Table3.csv file
  m           <- c(m, c(m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi, m_phi))
  W           <- c(W, c(75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0))
  Q2          <- c(Q2, c(2.71, 3.21, 3.82, 4.52, 5.36, 6.35, 7.60, 9.30, 12.0, 14.85, 19.20, 32.15))
  sigma       <- c(sigma, c(72.6, 64.5, 46.4, 35.0, 25.1, 18.2, 12.6, 7.04, 5.34, 2.25, 1.28, 0.371))
  deltaSigma  <- c(deltaSigma, c(sqrt(8.5**2+10.6**2), sqrt(6.2**2+8.3**2), sqrt(4.5**2+5.2**2), sqrt(4.0**2+4.1**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(2.7**2+2.6**2), sqrt(2.0**2+1.9**2), sqrt(1.5**2+1.2**2), sqrt(0.90**2+0.78**2)))
  deltaSigma  <- c(deltaSigma, c(sqrt(0.63**2+0.56**2), sqrt(0.38**2+0.3**2), sqrt(0.24**2+0.17**2), sqrt(0.076**2+0.049**2)))
  # Getting data from Table9.csv
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(41, 54, 67, 80, 93))
  Q2          <- c(Q2, c(3.3, 3.3, 3.3, 3.3, 3.3))
  sigma       <- c(sigma, c(308, 294, 346, 416, 397))
  deltaSigma  <- c(deltaSigma, c(sqrt(14**2+30**2), sqrt(16**2+26**2), sqrt(19**2+30**2), sqrt(23**2+40**2), sqrt(27**2+40**2)))
  # Getting data from Table10.csv
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(48, 64, 80, 96, 114))
  Q2          <- c(Q2, c(6.6, 6.6, 6.6, 6.6, 6.6))
  sigma       <- c(sigma, c(72.6, 97.2, 99.3, 120, 115))
  deltaSigma  <- c(deltaSigma, c(sqrt(4.3**2+7.8**2), sqrt(5.9**2+9.4**2), sqrt(6.2**2+7.2**2), sqrt(9**2+10**2), sqrt(10**2+10**2)))
  # Getting data from Table11.csv
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(59, 77, 95, 113, 131))
  Q2          <- c(Q2, c(11.9, 11.9, 11.9, 11.9, 11.9))
  sigma       <- c(sigma, c(25.3, 33.1, 32.1, 27.5, 34.9))
  deltaSigma  <- c(deltaSigma, c(sqrt(2.6**2+2.4**2), sqrt(2.7**2+2.4**2), sqrt(3.5**2+3.0**2), sqrt(4.1**2+3.5**2), sqrt(3.2**2+3.1**2)))
  # Getting data from Table12.csv
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(61, 83, 105, 127, 149))
  Q2          <- c(Q2, c(19.5, 19.5, 19.5, 19.5, 19.5))
  sigma       <- c(sigma, c(6.5, 9.6, 9.3, 9.8, 16.9))
  deltaSigma  <- c(deltaSigma, c(sqrt(0.6**2+0.6**2), sqrt(0.9**2+0.7**2), sqrt(1.0**2+0.8**2), sqrt(1.1**2+0.8**2), sqrt(1.8**2+1.4**2)))
  # Getting data from Table13.csv
  m           <- c(m, c(m_rho0, m_rho0, m_rho0, m_rho0, m_rho0))
  W           <- c(W, c(71, 97, 116, 139, 165))
  Q2          <- c(Q2, c(35.6, 35.6, 35.6, 35.6, 35.6))
  sigma       <- c(sigma, c(1.2, 2.0, 2.3, 3.4, 2.9))
  deltaSigma  <- c(deltaSigma, c(sqrt(0.2**2+0.1**2), sqrt(0.4**2+0.2**2), sqrt(0.5**2+0.2**2), sqrt(0.6**2+0.3**2), sqrt(0.7**2+0.3**2)))
  # Getting data from Table17.csv
  m           <- c(m, c(m_phi, m_phi, m_phi, m_phi, m_phi))
  W           <- c(W, c(41, 54, 67, 80, 93))
  Q2          <- c(Q2, c(3.3, 3.3, 3.3, 3.3, 3.3))
  sigma       <- c(sigma, c(41.2, 55.1, 49.2, 57.5, 69.6))
  deltaSigma  <- c(deltaSigma, c(sqrt(4.3**2+5.3**2), sqrt(5.7**2+5.4**2), sqrt(6.9**2+6.2**2), sqrt(7.4**2+6.8**2), sqrt(8.4**2+7.5**2)))
  # Getting data from Table18.csv
  m           <- c(m, c(m_phi, m_phi, m_phi, m_phi, m_phi))
  W           <- c(W, c(48, 64, 80, 96, 114))
  Q2          <- c(Q2, c(6.6, 6.6, 6.6, 6.6, 6.6))
  sigma       <- c(sigma, c(13.2, 13.0, 20.5, 14.7, 23.3))
  deltaSigma  <- c(deltaSigma, c(sqrt(1.5**2+1.5**2), sqrt(1.7**2+1.5**2), sqrt(2.3**2+1.9**2), sqrt(2.4**2+1.7**2), sqrt(4.0**2+2.2**2)))
  # Getting data from Table19.csv
  m           <- c(m, c(m_phi, m_phi, m_phi, m_phi, m_phi))
  W           <- c(W, c(71, 97, 116, 139, 165))
  Q2          <- c(Q2, c(15.8, 15.8, 15.8, 15.8, 15.8))
  sigma       <- c(sigma, c(2.3, 2.5, 3.9, 4.4, 7.6))
  deltaSigma  <- c(deltaSigma, c(sqrt(0.3**2+0.2**2), sqrt(0.4**2+0.2**2), sqrt(0.6**2+0.3**2), sqrt(0.8**2+0.5**2), sqrt(3.2**2+0.6**2)))

  data.frame(m, W, Q2, sigma, deltaSigma)
}

#' @export
loadData.VMPDSigma  <- function(vmp) {
  # This function creates vectors with the values of m,  W,  Q2,  t,  dsigma and deltaDSigma
  # of the data used by arxiv:1307.0009
  # Getting data from 694372 folder
  m            <- c()
  W            <- c()
  Q2           <- c()
  t            <- c()
  dsigma       <- c()
  deltaDSigma  <- c()
  data.frame(m, W, Q2, t, dsigma, deltaDSigma)
}
