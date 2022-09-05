#' @title Weighted Sample Mean and Variance
#'
#' @param data the data must contain two or three columns : code, y, and weight data if exist.
#'
#' @return This function returns a list with following objects :
#' \item{direst}{an information about direct estimatior in each area}
#' \item{pcap}{pcap (the weighted sample mean), vardir (the weighted sample variance),yt (the total number of the "success" category from each area), and nt (the total number of sample from each area)}
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
#' pcapdir(dataEB[,-c(3)])
#'
#' ##load dataset with weight value
#' data(dataEB)
#' pcapdir(dataEB)
#'
pcapdir <- function(data){
  data.dir <- direct.est(data)
  nt <- sum(data.dir$ni)
  yt <- sum(data.dir$yi)
  if(ncol(data)==3){ ## if weight exists
    prw <- prop.table(data$weight)
    wsum <- aggregate(prw~data$code,FUN=sum)[,2]
    pcap <- weighted.mean(data.dir$p,wsum)
    pcap <- c(pcap=pcap,vardir=sum(wsum*(data.dir$p-pcap)^2),yt=yt,nt=nt)
  }
  else if(ncol(data)==2){
    pcap <- sum(data.dir$ni/nt*data.dir$p)
    pcap <- c(pcap=pcap,vardir=sum(data.dir$ni/nt*(data.dir$p-pcap)^2),yt=yt,nt=nt)
  }
  return(list(direst=data.dir, pcap=pcap))
}
