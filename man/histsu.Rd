\name{histsu}
\alias{histsu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Histogram with exact number of bins specified by the user }
\description{
The generic function histsu computes a histogram of the given data values.
}
\usage{
histsu(x, breaks = "Sturges", freq = NULL, probability = !freq, 
include.lowest = TRUE, right = TRUE, density = NULL, angle = 45, col = NULL, 
border = NULL, main = paste("Histogram of", xname), xlim = range(breaks), 
ylim = NULL, xlab = xname, ylab, axes = TRUE, plot = TRUE, labels = FALSE, 
nclass = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of values for which the histogram is desired. }
  \item{breaks}{ Either: 1) A vector giving the breakpoints between histogram 
  cells, OR 2) A single number giving the number of cells for the
  histogram, OR 3) A character string naming an algorithm to compute the number 
  of cells (see Details), OR 4) A function to compute the number of cells. }
  \item{freq}{ logical; if \code{TRUE}, the histogram graphic is a 
  representation of frequencies, the \code{counts} component of the result; 
  if \code{FALSE}, probability densities, component 'density', are plotted (so 
  that the histogram has a total area of one). 
  Defaults to \code{TRUE} iff 'breaks' are equidistant (and 'probability' is not 
  specified). }
  \item{probability}{ A logical value, \code{TRUE} means it is not a frequency 
  graph. }
  \item{include.lowest}{ If \code{TRUE}, an \code{x[i]} equal to the 'breaks' 
  value will be included in the first (or last, for \code{right=FALSE}) bar.  
  This will be ignored (with a warning) unless 'breaks' is a vector. }
  \item{right}{ If \code{TRUE}, the histograms cells are right-closed 
  (left open) intervals. }
  \item{density}{ The density of shading lines, in lines per inch. The default 
  value of \code{NULL} means that no shading lines are drawn.
  Non-positive values of 'density' also inhibit the drawing of shading lines. }
  \item{angle}{ The slope of shading lines, given as an angle in degrees 
  (counter-clockwise). }
  \item{col}{ A colour to be used to fill the bars. The default of \code{NULL} 
  yields unfilled bars. }
  \item{border}{ The color of the border around the bars.  The default is to use 
  the standard foreground color. }
  \item{main}{ Title of the graph. }
  \item{xlim}{ A two valued vector specifying the lower and upper limits of the 
  x axis. }
  \item{ylim}{ A two valued vector specifying the lower and upper limits of the 
  y axis. }
  \item{xlab}{ X axis labels. }
  \item{ylab}{ Y axis labels. }
  \item{axes}{ Logical value, if \code{TRUE}, axis will be drawn. }
  \item{plot}{ Logical value, if \code{TRUE}, plot will be drawn. }
  \item{labels}{ Logical or character.  Additionally draw labels on top of bars, 
  if not \code{FALSE}; see \code{\link{plot.histogram}}. }
  \item{nclass}{ Number of bins of the histogram. }
  \item{\dots}{ Other graphical parameters, see par for details. }
}
\details{
See \code{\link{hist}} help file. This function forces the number of class of 
histogram to that as specified by the user.
}
\value{
An object of class \code{"histogram"} which is a list with components:
  \item{breaks}{ The n+1 cell boundaries (=\code{breaks} if that was a vector). 
  }
  \item{counts}{ N integers; for each cell, the number of \code{x[]} inside. }
  \item{density}{ Values as estimated density values. If \code{all(diff(breaks) 
  == 1)}, they are the relative frequencies counts/n. } 
  \item{intensities}{ Same as density, deprecated. }
  \item{mids}{ The n cell midpoints. }
  \item{xname}{ A character string with the actual \code{x} argument name. }
  \item{equidist}{ Logical, indicating if the distances between \code{breaks} 
  are all the same.}
}
\references{ Venables, W. N. and Ripley. B. D. (2002) Modern Applied Statistics 
with S.  Springer.
 }
\author{ R development team with modifications by Steve Su }
\note{ Please see \code{hist} help file.
}
\seealso{ \code{\link{hist}} }
\examples{
# See hist for extended example:
junk<-rgamma(1000,5)
# Forcing the number of bins to be 10:
histsu(junk,nclass=10)
}
\keyword{hplot} 
