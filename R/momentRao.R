#' @title Estimates alpha and beta parameter with Moment method by J.N.K.Rao
#'
#' @param data.dir Direct estimates of the data from function pcapdir
#' @param pcap weighted sample mean and variance from function pcapdir
#'
#' @return This function returns a data frame with following objects :
#' \item{alpha_cap}{an alpha estimator by Moment method of Claire E.B.O.}
#' \item{beta_cap}{an beta estimator by Moment method of Claire E.B.O.}
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
#' momentRao(data.dir = temp$direst, pcap = temp$pcap)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' temp = pcapdir(dataEB[,-c(3)])
#' momentRao(data.dir = temp$direst, pcap = temp$pcap)
#'
momentRao <- function(data.dir,pcap){
  nt <- sum(data.dir$ni)
  yt <- sum(data.dir$yi)
  m <- nrow(data.dir)
  beta_alpha_cap_upper <- pcap[1]*(1-pcap[1])*(nt - sum((data.dir$ni)^2/nt) - (m - 1))
  beta_alpha_cap_lower <- (nt*pcap[2]) - (pcap[1] * (1 - pcap[1]) * (m - 1))
  alpha_cap <- pcap[1] * ((beta_alpha_cap_upper/beta_alpha_cap_lower)-1)
  beta_cap <- pcap[1] * ((beta_alpha_cap_upper/beta_alpha_cap_lower)-1) * ((1/pcap[1])-1)
  return(data.frame(alpha_cap=alpha_cap,beta_cap=beta_cap))
}
