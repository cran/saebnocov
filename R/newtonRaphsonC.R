#' @title Estimates alpha and beta parameter with Newton Raphson method by Claire E.B.O.
#'
#' @param data.dir Direct estimates of the data from function pcapdir
#' @param pcap weighted sample mean and variance from function pcapdir
#' @param maxiter the Maximum iteration value
#' @param tol Tolerance error value in iteration
#'
#' @return This function returns a data frame with following objects :
#' \item{alpha_cap}{an alpha estimator by Newton Raphson method of Claire E.B.O.}
#' \item{beta_cap}{an beta estimator by Newton Raphson method of Claire E.B.O.}
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
#' newtonRaphsonC(data.dir = temp$direst, pcap = temp$pcap,
#'  maxiter = 100, tol = 0.00001)
#'
#' ##load dataset with weight value
#' data(dataEB)
#' temp = pcapdir(dataEB[,-c(3)])
#' newtonRaphsonC(data.dir = temp$direst, pcap = temp$pcap,
#'  maxiter = 100, tol = 0.00001)
#'
newtonRaphsonC <- function(data.dir,pcap,maxiter,tol){
  nt <- pcap[4]
  yt <- pcap[3]
  m <- nrow(data.dir)
  pbar <- pcap[1]
  s2 <- pcap[2]
  alpha_cap <- pbar*((pbar*(1-pbar)/s2)-1)
  beta_cap <- (1-pbar)*((pbar*(1-pbar)/s2)-1)
  theta_0 <- cbind(alpha_cap,beta_cap)
  diff <- matrix(rep(tol+1,2*m),ncol=2)
  iter <- 0
  while(max(diff) > tol & iter < maxiter){
    iter <- iter+1
    g <- vectorClaire(theta_0[1,1],theta_0[1,2],data.dir$p)
    G <- matrixClaire(theta_0[1,1],theta_0[1,2])
    G1g <- t(solve(G)%*%matrix(g,ncol=1))
    theta_1 <- theta_0 - G1g
    diff <- theta_1-theta_0
    theta_0 <- theta_1
  }
  return(data.frame(alpha_cap=theta_0[,1],beta_cap=theta_0[,2]))
}
