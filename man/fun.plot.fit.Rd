\name{fun.plot.fit}
\alias{fun.plot.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the univariate generalised lambda distribution fits on the data 
set. }
\description{ This function is designed for univariate generalised lambda 
distribution fits only.
}
\usage{
fun.plot.fit(fit.obj, data, nclass = 50, xlab = "", name = "", param.vec,
ylab="Density", main="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{ Fitted object from \code{\link{fun.data.fit.ml}}, 
  \code{\link{fun.data.fit.hs}}, \code{\link{fun.data.fit.hs.nw}}, 
  \code{\link{fun.RPRS.ml}}, \code{\link{fun.RMFMKL.ml}}, 
  \code{\link{fun.RPRS.hs}}, \code{\link{fun.RMFMKL.hs}}, 
  \code{\link{fun.RPRS.hs.nw}}, \code{\link{fun.RMFMKL.hs.nw}} }
  \item{data}{ Dataset to be plotted. }
  \item{nclass}{ Number of class of histogram, the default is 50. }
  \item{xlab}{ Label on the x axis. }
  \item{name}{ Naming the type of distribution fits. }
  \item{param.vec}{ A vector describing the type of generalised lambda 
  distribution used in the \code{fit.obj}. }
  \item{ylab}{ Label on the y axis. }
  \item{main}{ Title of the graph. }
}

\value{
A graphical output showing the data and the resulting distributional fits.
}

\author{ Steve Su }
\note{ If the distribution fits over fits the peak of the distribution, it can 
be difficult to see the actual data set.
}
\seealso{ \code{\link{fun.plot.fit.bm}}, \code{\link{fun.data.fit.ml}}, 
\code{\link{fun.data.fit.hs}}, \code{\link{fun.data.fit.hs.nw}}, 
\code{\link{fun.RPRS.ml}}, \code{\link{fun.RMFMKL.ml}}, 
\code{\link{fun.RPRS.hs}}, \code{\link{fun.RMFMKL.hs}}, 
\code{\link{fun.RPRS.hs.nw}}, \code{\link{fun.RMFMKL.hs.nw}} }
\examples{
\donttest{
# Generate Normally distribute random numbers as dataset
 junk<-rnorm(1000,3,2)

# Fit the data set using fun.data.fit.ml. 
# Also, fun.data.fit.hs or fun.data.fit.hs.nw can be used.
 obj.fit<-fun.data.fit.ml(junk)

# Plot the resulting fits
 fun.plot.fit(obj.fit,junk,xlab="x",name=".ML",param.vec=c("rs","fmkl","fmkl"))

# This function also works for singular fits such as those by fun.RPRS.ml,
# fun.RMFMKL.ml, fun.RPRS.hs, fun.RMFMKL.hs, fun.RPRS.hs.nw, fun.RMFMKL.hs.nw
 junk<-rnorm(1000,3,2)
 obj.fit<-fun.RPRS.ml(junk)
 fun.plot.fit(obj.fit,junk,xlab="x",name="RPRS.ML",param.vec=c("rs"))
}}
\keyword{hplot}