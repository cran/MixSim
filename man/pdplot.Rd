\name{pdplot}
\alias{pdplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallel Distribution Plot }
\description{
Constructs a parallel distribution plot for Gaussian finite mixture models
}
\usage{
pdplot(Pi, Mu, S, file = NULL, Nx = 5, Ny = 5, MaxInt = 1, marg = c(2,1,1,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pi}{ vector of mixing proprtions }
  \item{Mu}{ matrix consisting of components' mean vectors (K x p) }
  \item{S}{ set of components' covariance matrices (p x p x K) }
  \item{file}{ name of .pdf-file }
  \item{Nx}{ number of color levels for smoothing along the x-axis}
  \item{Ny}{ number of color levels for smoothing along the y-axis }
  \item{MaxInt}{ maximum color intensity }
  \item{marg}{ plot margins }
}
\details{
If 'file' is specified, produced plot will be saved as a .pdf-file
}

\references{
Maitra, R. and Melnykov, V. (2010) "Simulating data to study performance of finite mixture modeling and clustering algorithms", The Journal of Computational and Graphical Statistics, 2:19, 354-376.

}
\author{ Melnykov, V., Chen, W.-C., Maitra, R.}

\seealso{ MixSim, overlap, simdataset }
\examples{
\dontrun{
data(iris)

K <- 3
p <- dim(iris)[2] - 1
n <- dim(iris)[1]
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

pdplot(Pi = Pi, Mu = Mu, S = S)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ hplot }% __ONLY ONE__ keyword per line
