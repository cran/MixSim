\name{simdataset}
\alias{simdataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dataset Simulation }
\description{
  Simulates a datasets of sample size n given parameters of finite mixture model with Gaussian components
}
\usage{
simdataset(n, Pi, Mu, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ sample size }
  \item{Pi}{ vector of mixing proprtions (length K) }
  \item{Mu}{ matrix consisting of components' mean vectors (K x p) }
  \item{S}{ set of components' covariance matrices (p x p x K) }
}
\details{
  Numbers of observations in components are assigned as a realization from multinomial distribution with probabilities given by mixing proportions
}
\value{
  \item{x }{ simulated dataset (n x p)}
  \item{id }{ classification vector (length n)}
}
\references{
Maitra, R. and Melnykov, V. (200?) "Simulating data to study performance of finite mixture modeling and clustering algorithms", The Journal of Computational and Graphical Statistics.

Davies, R. (1980) "The distribution of a linear combination of chi-square random variables", Applied Statistics, 29, 323-333.
}
\author{ Melnykov, V., Chen, W.-C., Maitra, R.}

\seealso{ MixSim, overlap, pdplot }
\examples{
K <- 4
repeat{
   Q <- MixSim(BarOmega = 0.01, MaxOmega = 0.05, K = 4, p = 2)
   if (Q$fail == 0) break
}
A <- simdataset(n = 1000, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
colors <- c("red", "green", "blue", "brown")
plot(A$x, xlab = "x1", ylab = "x2", type = "n")
for (k in 1:K){
   points(A$x[A$id == k, ], col = colors[k], pch = 19, cex = 0.4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ datagen }% __ONLY ONE__ keyword per line
