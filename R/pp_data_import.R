#' @export
loadData.ppDSigma <- function(pp) {
  setwd(paste0(system.file(package = 'HQCDP'), '/extdata/pp'))
  file_list <- list.files()
  # Getting data from the first file
  my_data <- data.frame(read.csv(file_list[1]))
  l <- c()
  for(i in 1:length(my_data$minus_t)) {l <- append(l, 23)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  new_data1 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,
                    my_data$error_plus)
  names(new_data1) <- c("sqrt_s","minus_t","dsigma","error")
  #plot(new_data1$minus_t,new_data1$dsigma, col="blue", type="l",xlab = "-t",ylab = "dsigma/dt")
  data <- new_data1
  # Getting the data from the second file
  my_data <- data.frame(read.csv(file_list[2]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 23.46)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.05*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data2 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data2) <- c("sqrt_s","minus_t","dsigma","error")
  #plot(new_data2$minus_t,new_data2$dsigma, col="blue", type="l",xlab = "-t",ylab = "dsigma/dt")
  data <- rbind(new_data1,new_data2)
  # Getting the data from the third file
  my_data <- data.frame(read.csv(file_list[3]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 30.54)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.05*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data3 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data3) <- c("sqrt_s","minus_t","dsigma","error")
  #plot(new_data3$minus_t,new_data3$dsigma, col="blue", type="l",xlab = "-t",ylab = "dsigma/dt")
  data <- rbind(data,new_data3)
  # Getting the data from the fourth file
  my_data <- data.frame(read.csv(file_list[4]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 31.0)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.1*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data4 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data4) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data4)
  # Getting the data from the fifth file
  my_data <- data.frame(read.csv(file_list[5]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 44.64)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.05*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data5 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data5) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data5)
  # Getting the data from the sixth file
  my_data <- data.frame(read.csv(file_list[6]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 52.81)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.05*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data6 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data6) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data6)
  # Getting the data from the seventh file
  my_data <- data.frame(read.csv(file_list[7]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 53.00)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.1*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data7 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data7) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data7)
  # Getting the data from the eighth file
  my_data <- data.frame(read.csv(file_list[8]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 62.0)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  new_data8 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,my_data$error_plus)
  names(new_data8) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data8)
  # Getting the data from the ninth file
  my_data <- data.frame(read.csv(file_list[9]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 62.0)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.1*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data9 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data9) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data9)
  # Getting the data from the tenth file
  my_data <- data.frame(read.csv(file_list[10]))
  l <- c()
  for(i in 1:length(my_data$minus_t)){l <- append(l, 62.07)}
  l <- data.frame(l)
  names(l) <- "sqrt_s"
  my_data$sqrt_s <- l
  # Compute total error
  exp_error <- sqrt((0.05*my_data$dsigma)^2 + (my_data$error_plus)^2)
  new_data10 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data10) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data10)
  # Getting the data from the eleventh file
  my_data <- data.frame(read.csv(file_list[11]))
  my_data$sqrt_s <- rep(7000,length(my_data$minus_t))
  # Compute total error
  exp_error <- sqrt((my_data$syst_plus)^2 + (my_data$error_plus)^2)
  new_data11 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data11) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data11)
  # Getting the data from the twelveth file
  my_data <- data.frame(read.csv(file_list[12]))
  my_data$sqrt_s <- rep(7000,length(my_data$minus_t))
  # Compute total error
  exp_error <- sqrt((my_data$syst_plus)^2 + (my_data$error_plus)^2)
  new_data12 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data12) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data12)
  #write.csv(data, file = "pp_data.csv", row.names=FALSE)
  # Getting the data from the thirteenth file
  my_data <- data.frame(read.csv(file_list[13]))
  my_data$sqrt_s <- rep(8000,length(my_data$minus_t))
  # Compute total error
  exp_error <- sqrt((my_data$syst_plus)^2 + (my_data$error_plus)^2)
  new_data13 <- data.frame(my_data$sqrt_s,my_data$minus_t,my_data$dsigma,exp_error)
  names(new_data13) <- c("sqrt_s","minus_t","dsigma","error")
  data <- rbind(data,new_data13)
  full_data <- data.frame(Q2 = rep(NA, length(data$dsigma)), W = data$sqrt_s, t = -data$minus_t,
             dsigma = data$dsigma, deltaDSigma = data$error)
  full_data
}
