\name{summary.ivbma}
\alias{summary.ivbma}
\title{Generic summary function for object of class ivbma}
\description{
Generic summary function for object of class \code{ivbma}
}
\usage{
\method{summary}{ivbma}(object,nms.U=NULL,nms.V=NULL,...)
}
\arguments{
\item{object}{\code{ivbma} object returned by \code{\link{ivbma}}}
\item{nms.U}{names for the parameter \eqn{\rho = (\delta,\tau)^T}}
\item{nms.V}{names for the parameter \eqn{\lambda = (\beta,\gamma)^T}}
\item{...}{additional parameters}
}
\seealso{
\code{\link{ivbma}}
}
\references{Anna Karl and Alex Lenkoski (2012). "Instrumental Variable Bayesian
Model Averaging via Conditional Bayes Factors" http://arxiv.org/abs/1202.5846}
