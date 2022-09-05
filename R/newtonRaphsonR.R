#' @title Estimates alpha and beta parameter with Newton Raphson method by J.N.K. Rao
#'
#' @param data.dir Direct estimates of the data from function pcapdir
#' @param pcap weighted sample mean and variance from function pcapdir
#' @param maxiter the Maximum iteration value
#' @param tol Tolerance error value in iteration
#'
#' @return This function returns a data frame with following objects :
#' \item{alpha_cap}{an alpha estimator by Newton Raphson method of J.N.K.Rao}
#' \item{beta_cap}{an beta estimator by Newton Raphson method of J.N.K.Rao}
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
#' newtonRaphsonR(data.dir = temp$direst, pcap = temp$pcap,
#'  maxiter = 100, tol = 0.00001)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' temp = pcapdir(dataEB)
#' newtonRaphsonR(data.dir = temp$direst, pcap = temp$pcap,
#'  maxiter = 100, tol = 0.00001)
#'
newtonRaphsonR <- function(data.dir,pcap,maxiter,tol){
  nt <- sum(data.dir$ni)
  yt <- sum(data.dir$yi)
  m <- nrow(data.dir)
  beta_alpha_cap_upper <- pcap[1]*(1-pcap[1])*(nt - sum((data.dir$ni)^2/nt) - (m - 1))
  beta_alpha_cap_lower <- (nt*pcap[2]) - (pcap[1] * (1 - pcap[1]) * (m - 1))
  alpha_cap <- pcap[1] * ((beta_alpha_cap_upper/beta_alpha_cap_lower)-1)
  beta_cap <- pcap[1] * ((beta_alpha_cap_upper/beta_alpha_cap_lower)-1) * ((1/pcap[1])-1)
  diff <- tol+1
  iter <- 0
  theta_0 <- cbind(alpha_cap,beta_cap)
  while(max(diff) > tol & iter < maxiter){
    iter <- iter+1
    g <- vectorRao(theta_0[1,1],theta_0[1,2],data.dir$ni,data.dir$yi)
    G <- matrixRao(theta_0[1,1],theta_0[1,2],data.dir$ni,data.dir$yi)
    G1g <- t(solve(G)%*%matrix(g,ncol=1))
    theta_1 <- theta_0 - G1g
    diff <- theta_1-theta_0
    theta_0 <- theta_1
  }
  return(data.frame(alpha_cap=theta_0[,1],beta_cap=theta_0[,2]))
}
