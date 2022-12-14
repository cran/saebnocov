% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcapdir.R
\name{pcapdir}
\alias{pcapdir}
\title{Weighted Sample Mean and Variance}
\usage{
pcapdir(data)
}
\arguments{
\item{data}{the data must contain two or three columns : code, y, and weight data if exist.}
}
\value{
This function returns a list with following objects :
\item{direst}{an information about direct estimatior in each area}
\item{pcap}{pcap (the weighted sample mean), vardir (the weighted sample variance),yt (the total number of the "success" category from each area), and nt (the total number of sample from each area)}
}
\description{
Weighted Sample Mean and Variance
}
\examples{
## load dataset with no weight value
data(dataEB)
pcapdir(dataEB[,-c(3)])

##load dataset with weight value
data(dataEB)
pcapdir(dataEB)

}
