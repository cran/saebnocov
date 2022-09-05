#' @title Direct Estimates
#'
#' @param data the data must contain two or three columns : code, y, and weight data if exist.
#'
#' @return This function returns a data frame with following objects :
#' \item{p}{a proportion value of direct estimator in each area}
#' \item{yi}{the number of "success" category in each area (code 1)}
#' \item{ni}{the number of sample in each area}
#' \item{sumw}{the number of weight value in each area (if weight data exists)}
#' \item{prw}{the square of proportion table of weight data in each area (if weight data exists)}
#' \item{vardir}{variation value of direct estimator in each area}
#' \item{rrmse}{Relative Root Mean Square value of direct estimator in each area}
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
#' direct.est(dataEB[,-c(3)])
#'
#' ##load dataset with weight value
#' data(dataEB)
#' direct.est(dataEB)
#'

utils::globalVariables(c('code','ni','p','prw','vardir','weight'))

direct.est <- function(data){
  if(ncol(data)==3){
    data <- data%>%group_by(code)%>%mutate(wprop=prop.table(weight))
    ct = crosstab(data$code,data$y,data$weight,prop.r=T,plot=F)
    estdir <- data.frame(p=ct$prop.row[,2],yi=aggregate(y~code,data,sum)[,2],ni=aggregate(y~code,data,length)[,2])
    estdir <- estdir%>%mutate(sumw=aggregate(weight~code,data,sum)[,2],prw=aggregate(wprop^2~code,data,sum)[,2],
                              vardir=p*(1-p)*prw,rrmse=sqrt(vardir)/p*100)
  }
  else if(ncol(data)==2) {
    ct <- crosstab(data$code,data$y,prop.r=T,plot=F)
    estdir <- data.frame(p=ct$prop.row[,2],yi=aggregate(y~code,data,sum)[,2],ni=aggregate(y~code,data,length)[,2])
    estdir <- estdir%>%mutate(vardir=p*(1-p)/ni,rrmse=sqrt(vardir)/p*100)
  }
  return(estdir)
}
