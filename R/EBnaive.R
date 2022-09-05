#' @title Small Area Estimation method with Empirical Bayes and its RRMSE value by Naive Method
#'
#' @param data the data must contain two or three columns : code, y, and weight data if exist.
#' @param method Method to estimate alpha and beta parameter according to person(rao or claire)
#' @param opt Method to estimate alpha and beta parameter according to the way of calculation (moment or nr)
#' @param maxiter the Maximum iteration value with default 100
#' @param tol Tolerance error value at iteration with default 0.00001
#'
#' @return This function returns a list with following objects :
#' \item{finalres}{an information about direct estimatior and EB estimator in each area}
#' \item{estimation}{an information about EB estimator and its RRMSE value obtained by Naive method}
#' \item{parameter}{Alpha and beta estimator}
#' \item{pcap}{pcap (the weighted sample mean), vardir (the weighted sample variance),yt (the total number of the "success" category from each area), and nt (the total number of sample from each area)}
#' \item{dir.est}{an information about direct estimator}
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
#' ## Calculates EB estimator
#' ## with its RRMSE value by Naive method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by J.N.K.Rao
#' EBnaive(data = dataEB[,-c(3)],method = "rao",opt = "moment", maxiter = 100, tol = 1e-5)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' ## Calculates EB estimator
#' ## with its RRMSE value by Naive method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by Claire E.B.O.
#' EBnaive(data = dataEB, method = "claire",opt = "moment", maxiter = 100, tol = 1e-5)
#'
EBnaive <- function(data,method,opt,maxiter=100,tol=1e-5){
  data.dir <- direct.est(data)
  nt <- sum(data.dir$ni)
  yt <- sum(data.dir$yi)
  if(ncol(data)==3){ ## if weight exists
    prw <- prop.table(data$weight)
    wsum <- aggregate(prw~data$code,FUN=sum)[,2]
    pcap <- weighted.mean(data.dir$p,wsum)
    pcap <- c(pcap=pcap,vardir=sum(wsum*(data.dir$p-pcap)^2),yt=yt,nt=nt)
  }
  else if(ncol(data)==2){
    pcap <- sum(data.dir$ni/nt*data.dir$p)
    pcap <- c(pcap=pcap,vardir=sum(data.dir$ni/nt*(data.dir$p-pcap)^2),yt=yt,nt=nt)
  }
  m <- nrow(data.dir)
  param <- alphabetaEB(data.dir,pcap,method,opt,maxiter,tol)
  est <- estEBnaive(data.dir,pcap,param)
  est$mse <- replace(est$mse,est$mse < 0,0)
  est$rrmse <- sqrt(est$mse)/est$eb.est*100
  datafinal <- data.dir%>%mutate(eb.p=est$eb.est,mse.eb.naive=est$mse,
                                 rrmse.eb.naive=sqrt(est$mse)/est$eb.est*100)
  return(list(finalres = datafinal,estimation=est,parameter=param,pcap=pcap,dir.est=data.dir))
}
