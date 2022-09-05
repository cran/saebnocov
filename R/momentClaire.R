#' @title Estimates alpha and beta parameter with Moment method by Claire E.B.O.
#'
#' @param data.dir Direct estimates of the data from function pcapdir
#' @param pcap weighted sample mean and variance from function pcapdir
#'
#' @return This function returns a data frame with following objects :
#' \item{alpha_cap}{an alpha estimator by Moment method of Claire E.B.O.}
#' \item{beta_cap}{a beta estimator by Moment method of Claire E.B.O.}
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
#' momentClaire(data.dir = temp$direst, pcap = temp$pcap)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' temp = pcapdir(dataEB[,-c(3)])
#' momentClaire(data.dir = temp$direst, pcap = temp$pcap)
#'
momentClaire <- function(data.dir,pcap){
  nt <- sum(data.dir$ni)
  yt <- sum(data.dir$yi)
  m <- nrow(data.dir)
  pbar <- pcap[1]
  s2 <- pcap[2]
  alpha_cap <- pbar*((pbar*(1-pbar)/s2)-1)
  beta_cap <- (1-pbar)*((pbar*(1-pbar)/s2)-1)
  return(data.frame(alpha_cap=alpha_cap,beta_cap=beta_cap))
}
