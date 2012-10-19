\name{fun.plot.many.gld}
\alias{fun.plot.many.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting many univariate generalised lambda distributions on one page. }
\description{
  This is a variant of the \code{\link{fun.plot.fit}} function.
}
\usage{
fun.plot.many.gld(fit.obj, data, xlab="", ylab="Density", main="", legd="",
param.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{ A matrix of generalised lambda distibutions parameters 
  from \code{\link{fun.data.fit.ml}}, \code{\link{fun.data.fit.hs}}, 
  \code{\link{fun.data.fit.hs.nw}}, \code{\link{fun.RPRS.ml}}, 
  \code{\link{fun.RMFMKL.ml}}, \code{\link{fun.RPRS.hs}}, 
  \code{\link{fun.RMFMKL.hs}}, \code{\link{fun.RPRS.hs.nw}}, 
  \code{\link{fun.RMFMKL.hs.nw}} functions. Or a matrix of generalised lambda
  distribution parameters. }
  \item{data}{ Dataset to be plotted or two values showing the ranges of value 
  to be compared. }
  \item{xlab}{ X-axis labels. }
  \item{ylab}{ Y-axis labels. }
  \item{main}{ Title for the plot. }
  \item{legd}{ Legend for the plot. }
  \item{param.vec}{ A vector showing the types of generalised lambda 
  distributions. This can be "rs" or "fmkl", only needed if you want to put your 
  own parameters for generalised lambda distributions which are not generated 
  from a fitting algorithm in this package. }
}
\value{
 A graph showing the different distributions on the same page.
}
\author{ Steve Su }
\note{ The data part of the function is not plotted, to see the dataset use the 
\code{\link{fun.plot.fit}} function.
}
\seealso{ \code{\link{fun.plot.fit}}, \code{\link{fun.plot.fit.bm}} }
\examples{

## Fit the dataset
# junk<-rnorm(1000,3,2)
# result.hs<-fun.data.fit.hs(junk,rs.default = "Y", fmkl.default = "Y", 
# rs.leap=3, fmkl.leap=3,rs.init = c(-1.5, 1.5), fmkl.init = c(-0.25, 1.5),
# no.c.rs=50,no.c.fmkl=50)

# par(mfrow=c(2,2))

## Plot the entire data range
# fun.plot.many.gld(result.hs,junk,"x","density","",
# legd=c("RPRS.hs", "RMFMKL.hs"))

## Plot and compare parts of the distributions
# fun.plot.many.gld(result.hs,c(1,2),"x","density","",legd=c("RPRS.hs", 
#"RMFMKL.hs"))
# fun.plot.many.gld(result.hs,c(0.1,0,2),"x","density","",legd=c("RPRS.hs", 
#"RMFMKL.hs"))
# fun.plot.many.gld(result.hs,c(3,4),"x","density","",legd=c("RPRS.hs", 
#"RMFMKL.hs"))

}

\keyword{hplot}