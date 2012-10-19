\name{QUnif}
\alias{QUnif}
\alias{sHalton}
%
\concept{Quasi Monte Carlo}
\concept{low discrepancy sequence}
\concept{space filling}
%
\title{Quasi Randum Numbers via Halton Sequences}
\description{
  These functions provide quasi random numbers or \emph{space filling} or
  \emph{low discrepancy} sequences in the \eqn{p}-dimensional unit cube.
}
\usage{
sHalton(n.max, n.min = 1, base = 2, leap = 1)
 QUnif (n, min = 0, max = 1, n.min = 1, p, leap = 1)
}
\arguments{
  \item{n.max}{maximal (sequence) number.}
  \item{n.min}{minimal sequence number.}
  \item{n}{number of \eqn{p}-dimensional points generated in
    \code{QUnif}.  By default, \code{n.min = 1, leap = 1} and
    the maximal sequence number is \code{n.max = n.min + (n-1)*leap}.}
  \item{base}{integer \eqn{\ge 2}{>= 2}: The base with respect to which
    the Halton sequence is built.}
  \item{min, max}{lower and upper limits of the univariate intervals.
    Must be of length 1 or \code{p}.}
  \item{p}{dimensionality of space (the unit cube) in which points are
    generated.}
  \item{leap}{integer indicating (if \eqn{> 1}) if the series should be
    leaped, i.e., only every \code{leap}th entry should be taken.}
}
\value{
  \code{sHalton(n,m)} returns a numeric vector of length \code{n-m+1} of
  values in \eqn{[0,1]}.

  \code{QUnif(n, min, max, n.min, p=p)} generates \code{n-n.min+1}
  p-dimensional points in \eqn{[min,max]^p} returning a numeric matrix
  with p columns.
}
\note{
  For \code{leap} Kocis and Whiten recommend values of
  \eqn{L=31,61,149,409}, and particularly the \eqn{L=409} for dimensions
  up to 400.
}
\references{
  James Gentle (1998)
  \emph{Random Number Generation and Monte Carlo Simulation}; sec.\ 6.3.
  Springer.

  Kocis, L. and Whiten, W.J. (1997)
  Computationsl Investigations of Low-Discrepancy Sequences.
  \emph{ACM Transactions of Mathematical Software} \bold{23}, 2, 266--294.
}
\author{Martin Maechler}
\examples{
32*sHalton(20, base=2)
QUnif(n=10,p=2,leap=409)
}
\keyword{math}
\keyword{multivariate}
\keyword{datagen}
