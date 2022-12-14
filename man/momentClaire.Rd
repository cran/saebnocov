% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/momentClaire.R
\name{momentClaire}
\alias{momentClaire}
\title{Estimates alpha and beta parameter with Moment method by Claire E.B.O.}
\usage{
momentClaire(data.dir, pcap)
}
\arguments{
\item{data.dir}{Direct estimates of the data from function pcapdir}

\item{pcap}{weighted sample mean and variance from function pcapdir}
}
\value{
This function returns a data frame with following objects :
\item{alpha_cap}{an alpha estimator by Moment method of Claire E.B.O.}
\item{beta_cap}{a beta estimator by Moment method of Claire E.B.O.}
}
\description{
Estimates alpha and beta parameter with Moment method by Claire E.B.O.
}
\examples{
## load dataset with no weight value
data(dataEB)
temp = pcapdir(dataEB[,-c(3)])
momentClaire(data.dir = temp$direst, pcap = temp$pcap)

##load dataset with weight value
data(dataEB)
temp = pcapdir(dataEB[,-c(3)])
momentClaire(data.dir = temp$direst, pcap = temp$pcap)

}
