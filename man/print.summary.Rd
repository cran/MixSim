\name{print.object}
\alias{print.MixSim}
\alias{summary.MixSim}
\title{Functions for Printing or Summarizing Objects}
\description{
  A \code{MixSim} class is declared, and these are functions
  to print and summary objects.
}
\usage{
\method{print}{MixSim}(x, ...)
\method{summary}{MixSim}(object, ...)
}
\arguments{
  \item{x}{an object with the 'MixSim' class attributes.}
  \item{object}{an object with the 'MixSim' class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing and debugging.

  For other functions, they only show summaries of objects. Use
  \code{names} or \code{str} to explore the details.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', The Journal of Computational and Graphical Statistics, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, (accepted).
}
\author{ Volodymyr Melnykov, Wei-Chen Chen, and Ranjan Maitra. }

\seealso{
  \code{\link{MixSim}}.
}
\examples{
\dontrun{
# Functions applied by directly type the names of objects.
}
}
\keyword{programming}
