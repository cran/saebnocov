#' @title Matrix G in Newton Raphson method by J.N.K.Rao
#'
#' @param alpha An alpha estimate value on iterating process
#' @param beta A beta estimate value on iterating process
#' @param ni The number of sample in each area
#' @param yi The number of "success" value in each area
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
#' matrixRao(alpha = temp1$alpha_cap,
#'  beta = temp1$beta_cap, ni = temp$direst$ni,
#'   yi = temp$direst$yi)
#'
matrixRao <- function(alpha,beta,ni,yi){
  G <- matrix(0,ncol=2,nrow=2)
  G[1,1] <- sum(sapply(1:length(ni),
                       function(x) sum(1/(alpha+beta+0:ni[x]-1)^2)-sum(1/(alpha+0:yi[x]-1)^2)))
  G[1,2] <- G[2,1] <- sum(sapply(1:length(ni),function(x) sum(1/(alpha+beta+0:ni[x]-1)^2)))
  ni.yi <- ni-yi
  G[2,2] <- sum(sapply(1:length(ni),
                       function(x) sum(1/(alpha+beta+0:ni[x]-1)^2)-sum(1/(beta+0:ni.yi[x]-1)^2)))
  return(G)
}
