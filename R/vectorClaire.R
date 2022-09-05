#' @title Vector g in Newton Raphson Method by Claire E.B.O.
#'
#' @param alpha An alpha estimate value on iterating process
#' @param beta A beta estimate value on iterating process
#' @param p direct estimator or proportion value
#'
#' @return This function returns a value of vector g.
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
#' ##calculates vector g
#' vectorClaire(alpha = temp1$alpha_cap, beta = temp1$beta_cap, p = temp$direst$p)
#'
vectorClaire <- function(alpha,beta,p){
  g1 <- digamma(alpha) - digamma(alpha+beta) - mean(log(p))
  g2 <- digamma(beta) - digamma(alpha+beta) - mean(log(1-p))
  return(cbind(g1=g1,g2=g2))
}
