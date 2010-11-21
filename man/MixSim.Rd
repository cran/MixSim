\name{MixSim}
\alias{MixSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixture Simulation }
\description{
Generates a finite mixture model with Gaussian components for prespecified levels of maximum and/or average overlaps
}
\usage{
MixSim(BarOmega = NULL, MaxOmega = NULL, K, p, sph = FALSE, hom = FALSE,
       ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 100,
       eps = 1e-06, lim = 1e06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BarOmega}{ value of desired average overlap }
  \item{MaxOmega}{ value of desired maximum overlap }
  \item{K}{ number of components }
  \item{p}{ number of dimensions }
  \item{sph}{ covariance matrix structure (FALSE - nonspherical, TRUE - spherical) }
  \item{hom}{ heterogeneous or homogeneous clusters (FALSE - heterogeneous, TRUE - homogeneous) }
  \item{ecc}{ maximum eccentricity }
  \item{PiLow}{ value of the smallest mixing proportion (if 'PiLow' is not reachable with respect to K, equal proportions are taken; PiLow = 1.0 implies equal proportions by default) }
  \item{int}{ mean vectors are simulated uniformly on a hypercube with sides specified by int = (lower.bound, upper.bound) }
  \item{resN}{ maximum number of mixture resimulations }
  \item{eps}{ error bound for overlap computation }
  \item{lim}{ maximum number of integration terms (Davies, 1980) }
}
\details{
  If 'BarOmega' is not specified, the function generates a mixture solely based on 'MaxOmega'; if 'MaxOmega' is not specified, the function generates a mixture solely based on 'BarOmega'.
}
\value{
  \item{Pi }{vector of mixing proprtions}
  \item{Mu }{matrix consisting of components' mean vectors (K x p)}
  \item{S }{set of components' covariance matrices (p x p x K)}
  \item{OmegaMap}{matrix of misclassification probabilities (K x K); OmegaMap[i,j] is the probability that X coming from the i-th component is classified to the j-th component}
  \item{BarOmega }{value of average overap}
  \item{MaxOmega }{value of maximum overap}
  \item{rcMax }{row and column numbers for the pair of components producing maximum overlap 'MaxOmega'}
  \item{fail }{flag value; 0 represents successful mixture generation, 1 represents failure}
}
\references{
Maitra, R. and Melnykov, V. (2010) "Simulating data to study performance of finite mixture modeling and clustering algorithms", The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Davies, R. (1980) "The distribution of a linear combination of chi-square random variables", Applied Statistics, 29, 323-333.
}
\author{ Melnykov, V., Chen, W.-C., Maitra, R.}

\seealso{ overlap, pdplot, simdataset}
\examples{
\dontrun{
set.seed(1234)

# controls average and maximum overlaps
MixSim(BarOmega = 0.05, MaxOmega = 0.15, K = 4, p = 5) 

# controls average overlap
MixSim(BarOmega = 0.05, K = 4, p = 5, hom = TRUE)

# controls maximum overlap
MixSim(MaxOmega = 0.15, K = 4, p = 5, sph = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ datagen }% __ONLY ONE__ keyword per line
