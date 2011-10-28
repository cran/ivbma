\name{ivbma.sample.theta}
\alias{ivbma.sample.theta}
\title{Updates the IVBMA parameters}
\description{
Runs one step of the Gibbs Sampler in IVBMA
}
\usage{
ivbma.sample.theta(theta,D,full)
}
\arguments{
\item{theta}{a list containing all the parameters in the IVBMA model}
\item{D}{dataset}
\item{full}{Indicator of whether model averaging is performed (FALSE)}
}

\value{
  Returns an updated theta object
}

\author{
  Alex Lenkoski (alex.lenkoski@uni-heidelberg.de)
  }
\seealso{
  \code{\link{ivbma}}
}
\references{Anna Karl and Alex Lenkoski (2012). "Instrumental Variable Bayesian
Model Averaging via Conditional Bayes Factors" http://arxiv.org/abs/1202.5846}
