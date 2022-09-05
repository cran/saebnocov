#' @title Small Area Estimation method with Empirical Bayes and its RRMSE value by Naive Method
#'
#' @param data.dir direct estimator information from function direct.est
#' @param pcap pcap (the weighted sample mean), vardir (the weighted sample variance),yt (the total number of the "success" category from each area), and nt (the total number of sample from each area)
#' @param param Alpha and Beta estimator
#'
#' @return This function returns a list with following objects :
#' \item{eb.est}{EB estimator in each area}
#' \item{mse}{MSE of EB estimator obtained by Naive method}
#' \item{rrmse}{RRMSE of EB estimator obtained by Naive method}
#'
#' @export
#'
#' @import descr
#' @import dplyr
#' @importFrom stats aggregate
#' @importFrom stats weighted.mean
#' @importFrom rlang .data
#'
#' @examples
#' ## load dataset with no weight value
#' data(dataEB)
#' temp = pcapdir(dataEB[,-c(3)])
#'
#' ## estimates alpha and beta parameter
#' ## in EB estimate with Moment method by J.N.K.Rao
#' temp1 = alphabetaEB(data.dir = temp$direst ,pcap = temp$pcap,
#'                       method = "rao", opt = "moment",
#'                       maxiter = 100,tol = 0.00001)
#'
#' ## calculates EB estimator
#' ## and its MSE by naive method
#' estEBnaive(data.dir = temp$direst, pcap = temp$pcap, param = temp1)
#'
estEBnaive <- function(data.dir,pcap,param){
  m <- nrow(data.dir)
  gamma_i <- data.dir$ni/(data.dir$ni + param[,1] + param[,2])
  est_teta_EB <- (gamma_i * data.dir$p) + ((1 - gamma_i) * pcap[1])
  MSE_EB_Naive_u <- (data.dir$yi + param[,1]) * (data.dir$ni - data.dir$yi + param[,2])
  MSE_EB_Naive_l <- (data.dir$ni + param[,1] + param[,2] + 1) * (data.dir$ni + param[,1] + param[,2])^2
  MSE_EB_Naive <- MSE_EB_Naive_u / MSE_EB_Naive_l
  return(list(eb.est=est_teta_EB,mse = MSE_EB_Naive, rrmse=sqrt(MSE_EB_Naive)/est_teta_EB*100))
}
