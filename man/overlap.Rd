\name{overlap}
\alias{overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Overlap}
\description{
Computes misclassification probabilities and pairwise overlaps for finite mixture models with Gaussian components. Overlap is defined as sum of two misclassification probabilities
}
\usage{
overlap(Pi, Mu, S, eps = 1e-06, acc = 1e-06, lim = 1e06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pi}{ vector of mixing proprtions (length K) }
  \item{Mu}{ matrix consisting of components' mean vectors (K x p) }
  \item{S}{ set of components' covariance matrices (p x p x K) }
  \item{eps}{ error bound for overlap computation }
  \item{acc}{ error bound for integration (Davies, 1980) }
  \item{lim}{ maximum number of integration terms (Davies, 1980) }
}

\value{
  \item{OmegaMap}{matrix of misclassification probabilities (K x K); OmegaMap[i,j] is the probability that X coming from the i-th component is classified to the j-th component}
  \item{BarOmega}{value of average overap}
  \item{MaxOmega}{value of maximum overap}
  \item{rcMax}{row and column numbers for the pair of components producing maximum overlap 'MaxOmega'}
}
\references{
Maitra, R. and Melnykov, V. (200?) "Simulating data to study performance of finite mixture modeling and clustering algorithms", The Journal of Computational and Graphical Statistics.

Davies, R. (1980) "The distribution of a linear combination of chi-square random variables", Applied Statistics, 29, 323-333.
}
\author{ Melnykov, V., Chen, W.-C., Maitra, R.}

\seealso{ MixSim, pdplot, simdataset }
\examples{
data(iris)

p <- dim(iris)[2] - 1
n <- dim(iris)[1]
K <- 3

id <- as.numeric(iris[,5])
Pi <- NULL
Mu <- NULL
S <- array(rep(0, p * p * K), c(p, p, K))

# estimate mixture parameters
for (k in 1:K){
	Pi <- c(Pi, sum(id == k) / n)
	Mu <- rbind(Mu, apply(iris[id == k,-5], 2, mean))
	S[,,k] <- var(iris[id == k,-5])
}

overlap(Pi = Pi, Mu = Mu, S = S)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }

