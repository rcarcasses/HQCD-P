# load all the data and put it in the proper format

#' @export
loadData.ppDSigma <- function(pp) {
  # This function imports data from the file pp_data.csv
  # and loads it in a pp object.
  # Data comes from the references 1, 10, 17 and 18 of 1309.1292
  data_loaded <- data.frame(read.csv("~/Desktop/HQCD-P/inst/extdata/pp/pp_data.csv"))
  names(data_loaded) <- c("sqrt_s","t", "dsigma", "deltaDSigma")
  data_loaded
}
