\name{ivbma}
\alias{ivbma}
\title{Instrumental Variable Bayesian Model Averaging via Conditional
  Bayes Factors}
\description{
This function estimates an Instrumental Variable (IV) system while
incorporating model uncertainty and performing model averaging using an
MC3-within-Gibbs Sampler.
}
\usage{
ivbma(Y, X, Z, W, s = 1000, b = round(s/10), full = FALSE, 
    odens = min(c(5000, s - b)), print.every = round(s/10), run.diagnostics = FALSE)
}
\arguments{
\item{Y}{n x 1 matrix. Response variable}
\item{X}{n x r matrix. Endogenous variables}
\item{W}{n x p matrix. Further explanatory variables.  You are
  responsible for including an intercept.}
\item{Z}{n x q matrix. Instrumental variables}
\item{s}{integer. Number of iterations}
\item{b}{integer. Number of iterations to discard as burn-in.}
\item{full}{If full is TRUE then model selection is not performed}
\item{odens}{Output density.  How many samples from the posterior should
  be returned?  Note that posterior expectations are taken over every
  sample after burn-in}
\item{print.every}{After how many iterations should the progress be printed?}
\item{run.diagnostics}{If TRUE, this will compute experimental
  diagnostics to assess the validity of the instruments in use.  Note
  that this adds a non-negligible amount of computing time.}
}
\details{
The function estimates the parameters based on the model
\deqn{Y = [X W] * \rho + \epsilon}
\deqn{X = [Z W] * \lambda + \eta} with
\deqn{(\epsilon_i, \eta_i)^T \sim N_2 ( 0,\Sigma)}
and its extension to multiple endogenous variables.
If \code{full} is set to \code{FALSE} model uncertainty is included
using conditional Bayes factors.
}
\value{
  \item{rho}{An  odens x (r + p) matrix with sampled values for the
    outcome stage.  Endogenous variables come first.} 
  \item{rho.bar}{Posterior expectation of the outcome stage taken over
    all iterations}
  \item{lambda}{A (p + q) x r x odens array with sampled values for
    the parameters of the first stage regressions.  Instruments come
    first.}
  \item{lambda.bar}{Posterior expectation of each first stage taken over
    all iterations}
  \item{Sigma}{odens sampled realizations of Sigma}
  \item{Sigma.bar}{Posterior expectation of Sigma taken over all
    iterations}
  \item{M}{Sampled first stage models}
  \item{M.bar}{Posterior first stage inclusion probabilities}
  \item{L}{Sampled second stage models}
  \item{L.bar}{Posterior second stage inclusion probabilities}

  If run.diagnostics was set to TRUE then you also receive

  \item{Sargan}{Model averaged Sargan p-values.  Lower values indicate
    lack of instrument validity}
  \item{Bayesian.Sargan}{An _Experimental_ Bayesian Sargan diagnostic based on
    Conditional Bayes Factors.  Same direction as above}

}

\author{
  Alex Lenkoski (alex.lenkoski@uni-heidelberg.de)
  }
\examples{
set.seed(1)
data(growth)
attach(growth)
## To replicate KL, set s to 1e5
a <- ivbma(Y, X, Z, W, s = 1e2)
summary(a, nms.U = c(names(Z), names(W)),nms.V = c(names(X), names(W)))
detach(growth)

set.seed(1)
data(margarine)
attach(margarine)
## To replicate KL, set s to 2.5e5
a <- ivbma(Y, X, Z, W, s=1e2)
summary(a, nms.U = c(names(Z), names(W)),nms.V = c(names(X), names(W)))
detach(margarine)

}

\seealso{
  \code{\link{summary.ivbma}}
  \code{\link{ivbma.cv.study}}
}
\references{Anna Karl and Alex Lenkoski (2012). "Instrumental Variable Bayesian
Model Averaging via Conditional Bayes Factors" http://arxiv.org/abs/1202.5846}
