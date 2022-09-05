#' @title Small Area Estimation method with Empirical Bayes and its RRMSE value by Jackknife Method
#'
#' @param data the data must contain two or three columns : code, y, and weight data if exist.
#' @param method Method to estimate alpha and beta parameter according to person(rao or claire)
#' @param opt Method to estimate alpha and beta parameter according to the way of calculation (moment or nr)
#' @param maxiter the Maximum iteration value with default 100
#' @param tol Tolerance error value at iteration with default 0.00001
#'
#' @return This function returns a list with following objects :
#' \item{finalres}{an information about direct estimator and EB estimator in each area with its RRMSE value obtained by jackknife method}
#' \item{eb.estimation}{an information about EB estimator in each area with its RRMSE value obtained by Naive method}
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
#' ## Calculates EB estimator with
#' ## its RRMSE value by Jackknife method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by J.N.K.Rao
#' jackknifeEB(data = dataEB[,-c(3)], method = "rao",
#'  opt = "moment", maxiter = 20, tol = 1e-5)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' ## Calculates EB estimator with
#' ## its RRMSE value by Jackknife method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by Claire E.B.O.
#' jackknifeEB(data = dataEB, method = "rao",
#'  opt = "moment", maxiter = 20, tol = 1e-5)
#'
jackknifeEB <- function(data,method,opt,maxiter=100,tol=1e-5){
  eb.est <- EBnaive(data,method,opt,maxiter,tol)
  data.dir <- eb.est$dir.est
  pcap <- eb.est$pcap
  eb.est$estimation$mse <- replace(eb.est$estimation$mse,eb.est$estimation$mse < 0,0)
  m <- nrow(data.dir)
  m1 <- m2 <- c()
  for(i in 1:m){
    datajack <- data.dir[-i,]
    prop.nj <- prop.table(datajack$ni)
    pcapjack <- sum(prop.nj*datajack$p)
    varjack <- sum(prop.nj*(datajack$p-pcapjack)^2)
    pcapjack <- c(pcap=pcapjack,vardir=varjack)
    param.jack <- alphabetaEB(datajack,pcapjack,method,opt,maxiter,tol)
    eb.jack <- estEBnaive(data.dir,pcap,param.jack)
    m1 <- cbind(m1, (eb.jack$mse-eb.est$estimation$mse))
    m2 <- cbind(m2,((eb.jack$eb.est-eb.est$estimation$eb.est)^2))
  }
  m1 <- eb.est$estimation$mse - (rowSums(m1))
  m2 <- c(m-1)/m * rowSums(m2)
  mse = m1+m2
  mse <- replace(mse,mse < 0,0)
  datafinal <- data.dir%>%mutate(eb.p=eb.est$estimation$eb.est,mse.eb.jack=mse,
                                 rrmse.eb.jack=sqrt(mse)/eb.est$estimation$eb.est*100)
  res <- list(finalres = datafinal,eb.estimation=eb.est)
  return(res)
}
