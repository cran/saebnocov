#' @title Small Area Estimation method with Empirical Bayes and its RRMSE value by Bootstrap Method
#'
#' @param data the data must contain two or three columns : code, y, and weight data if exist.
#' @param method Method to estimate alpha and beta parameter according to person(rao or claire)
#' @param opt Method to estimate alpha and beta parameter according to the way of calculation (moment or nr)
#' @param maxiter the Maximum iteration value with default 100
#' @param tol Tolerance error value at iteration with default 0.00001
#' @param B The number of iteration of bootstrap resampling with default 200
#' @param seed Setting a seed in set.seed() function to initialize a pseudorandom number generator with default number 0
#'
#' @return This function returns a list with following objects :
#' \item{finalres}{an information about direct estimator and EB estimator in each area with its RRMSE value obtained by bootstrap method}
#' \item{eb.estimation}{an information about EB estimator in each area with its RRMSE value obtained by Naive method}
#'
#' @export
#'
#' @import descr
#' @import dplyr
#' @importFrom stats aggregate
#' @importFrom stats weighted.mean
#' @importFrom rlang .data
#' @importFrom stats runif
#'
#' @references
#' Rao J, Peralta IM (2015). _Small Area Estimation Second Edition_. John Wiley & Sons, Inc.,Hoboken,
#' New Jersey, Canada. ISBN 978-1-118-73578-7.
#'
#' @examples
#' ## load dataset with no weight value
#' data(dataEB)
#' ## Calculates EB estimator with its
#' ## RRMSE value by Bootstrap method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by J.N.K.Rao
#' bootstrapEB(data = dataEB[,-c(3)], method = "rao",
#'  opt = "moment", maxiter = 20, tol = 1e-5,B=20,seed=0)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' ## Calculates EB estimator with its
#' ## RRMSE value by Bootstrap method.
#' ## Its alpha and beta estimator obtained
#' ## by Moment method by Claire E.B.O.
#' bootstrapEB(data = dataEB, method = "rao",
#'  opt = "moment", maxiter = 20, tol = 1e-5,B=20,seed=0)
#'


bootstrapEB <- function(data,method,opt,seed = NA ,maxiter=25,tol=1e-5,B=50){
  if (is.na(seed)) seed = round(runif(1,0,10000))
  eb.est <- EBnaive(data,method,opt,maxiter,tol)
  data.dir <- eb.est$dir.est
  pcap <- eb.est$pcap
  m1 <- 0
  m2 <- 0
  eb.est$estimation$mse <- replace(eb.est$estimation$mse,eb.est$estimation$mse < 0,0)
  param.boot0 <- data.frame(alpha_cap=1,beta_cap=1)
  for(i in 1:B){
    set.seed(seed+i)
    databoot <- data%>%group_by(code)%>%sample_frac(1,replace = T)
    if(ncol(data)==3){
      dtwboot <- data.frame(aggregate(weight~code,databoot,sum),sumw = data.dir$sumw)
      calib <- left_join(databoot,dtwboot,'code')
      colnames(calib)[c(3,4)] <- c('weight','sumweight')
      calib$weight <- (calib$sumw/calib$sumweight)*calib$weight
      databoot <- calib[,c(1:3)]
    }
    direst.boot <- pcapdir(databoot)
    param.boot <- alphabetaEB(direst.boot$direst,direst.boot$pcap,method,opt,maxiter,tol)
    if(any(unlist(param.boot)%in%c(0,-Inf,Inf))) param.boot <- param.boot0
    eb.boot <- estEBnaive(data.dir,pcap,param.boot)
    m1 <- m1+eb.boot$mse
    m2 <- m2+(eb.boot$eb.est-eb.est$estimation$eb.est)^2
    param.boot0 <- param.boot
  }
  m1 <- (2*eb.est$estimation$mse)-(m1/B)
  m2 <- m2/B
  mse <- m1+m2
  mse <- replace(mse,mse<0,0)
  datafinal <- data.dir%>%mutate(eb.p=eb.est$estimation$eb.est,mse.eb.boot=mse,
                                 rrmse.eb.boot=sqrt(mse)/eb.est$estimation$eb.est*100)
  res <- list(finalres = datafinal,eb.estimation=eb.est)
  return(res)
}
