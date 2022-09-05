#' @title Matrix G in Newton Raphson method by Claire E.B.O.
#'
#' @param alpha An alpha estimate value on iterating process
#' @param beta A beta estimate value on iterating process
#'
#' @return This function returns a value of matrix G.
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
#' ##calculates matrix G
#' matrixClaire(alpha = temp1$alpha_cap, beta = temp1$beta_cap)
#'
matrixClaire <- function(alpha,beta){
  G <- matrix(0,ncol=2,nrow=2)
  G[1,1] <- trigamma(alpha)-trigamma(alpha+beta)
  G[1,2] <- G[2,1] <- -trigamma(alpha+beta)
  G[2,2] <- trigamma(beta)-trigamma(alpha+beta)
  return(G)
}
