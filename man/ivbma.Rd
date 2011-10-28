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
## Runs one example of the simulation study in Karl and Lenkoski (2012)
n <- 120
q <- 10
p <- 15
Sigma.true <- matrix(c(1,0.4,0.4,1),ncol=2)
epseta.true <- matrix(rnorm(n * 2),n,2) \%*\% chol(Sigma.true)
eps.true <- epseta.true[,1]
eta.true <- epseta.true[,2]
beta.true  <- 1.5
gamma.true <- c(2,0,0,1.4,0,0,0,2.7,1.25,0,0,0,3.3,0,0)
delta.true <- c(0,0,4.1,0,0,0,1.2,3,0,0.9)
tau.true <- c(0,2.5,0,0,0,0,0,0,1.7,0,0,0,0.8,0,0)
Z <- matrix(rnorm(n * q),ncol=q)
W <- matrix(rnorm(n * p),ncol=p)
X <- Z \%*\% delta.true + W \%*\% tau.true + eta.true
Y <- X \%*\% beta.true + W \%*\% gamma.true + eps.true
a <- ivbma(Y,X,W,Z)
summary(a)
}
\seealso{
\code{\link{summary.ivbma}}
}
\references{Anna Karl and Alex Lenkoski (2012). "Instrumental Variable Bayesian
Model Averaging via Conditional Bayes Factors" http://arxiv.org/abs/1202.5846}
