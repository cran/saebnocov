#' @title Estimates alpha and beta parameter to obtain EB estimator
#'
#' @param data.dir Direct estimates of the data from function pcapdir
#' @param pcap weighted sample mean and variance from function pcapdir
#' @param method Method to estimate alpha and beta parameter according to person(rao or claire)
#' @param opt Method to estimate alpha and beta parameter according to the way of calculation (moment or nr)
#' @param maxiter the Maximum iteration value
#' @param tol Tolerance error value at iteration
#'
#' @return This function returns a data frame with following objects :
#' \item{alpha_cap}{an alpha estimator by user's choice method}
#' \item{beta_cap}{an beta estimator by user's choice method}
#'
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
#' ## estimates alpha and beta parameter
#' ## in EB estimate with Moment method by J.N.K.Rao
#' alphabetaEB(data.dir = temp$direst ,pcap = temp$pcap,
#' method = "rao", opt = "moment",maxiter = 100,tol = 0.00001)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' temp = pcapdir(dataEB)
#' ## estimates alpha and beta parameter
#' ## in EB estimate with Moment method by Claire E.B.O.
#' alphabetaEB(data.dir = temp$direst ,pcap = temp$pcap,
#' method = "claire", opt = "moment",maxiter = 100,tol = 0.00001)
#'
alphabetaEB <- function(data.dir,pcap,method,opt,maxiter,tol){
  if(method=='rao'){
    if(tolower(opt)=='moment'){
      return (momentRao(data.dir,pcap))
    }
    else if(tolower(opt)=='nr'){
      return (newtonRaphsonR(data.dir,pcap,maxiter,tol))
    }
    else stop("You should only choose between 'moment' and 'nr'")
  }
  else if(tolower(method)=='claire'){
    if(tolower(opt)=='moment'){
      return (momentClaire(data.dir,pcap))
    }
    else if(tolower(opt)=='nr'){
      return (newtonRaphsonC(data.dir,pcap,maxiter,tol))
    }
  }
  else stop("You should only choose between 'rao' and 'claire'")
}
