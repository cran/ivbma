\name{ivbma.cv.study}
\alias{ivbma.cv.study}
\title{Conducts a leave-one-out cross validation study using IVBMA}
\description{
This function allows one to conduct a leave-one-out cross validation
study using IVBMA. It takes an appropriately constructed object and
proceeds to drop each observation, fit IVBMA using the remaining
observations, forms a posterior predictive distribution of the dropped
observation and scores the predictive distribution along a number of metrics.
}
\usage{
ivbma.cv.study(d, ...) 
}
\arguments{
\item{d}{Object containing data, d$Y, d$X, d$W, d$Z must all be defined}
\item{...}{additional parameters to be passed to ivbma.  In particular,
  you'll usually want to set s.}
}
\value{
  This returns an n by 4 matrix.  Row i of the matrix gives the squared
  error (SE), absolute error (AE), predictive variance (VAR) and continous ranked
  probability score (CRPS) of the posterior predictive distribution
  leaving observation i out and subsequently using i as the verifying
  observation.
  }
\seealso{
\code{\link{ivbma}}
}
\references{Anna Karl and Alex Lenkoski (2012). "Instrumental Variable Bayesian
Model Averaging via Conditional Bayes Factors"
http://arxiv.org/abs/1202.5846}

\examples{
set.seed(1)
data(growth)
s <- 2e1 ##To replicate KL, set this to 2e5
a.bma <- ivbma.cv.study(growth, s = s, print.every = s)
a.full <- ivbma.cv.study(growth, s = s, print.every = s, full = TRUE)
print(rbind(colMeans(a.bma), colMeans(a.full)))
}
