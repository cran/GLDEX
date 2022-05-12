\name{qqplot.gld}
\alias{qqplot.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Do a quantile plot on the univariate distribution fits. }
\description{
  This plots the theoretical and actual data quantiles to allow the user to 
  examine the adequacy of a single gld distribution fit.
}
\usage{
qqplot.gld(data, fit, param, len = 10000, name = "", 
corner = "topleft",type="",range=c(0,1),xlab="",main="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data fitted. }
  \item{fit}{ Parameters of distribution fit. }
  \item{param}{ Can be either \code{"rs"} or \code{"fmkl"}. }
  \item{len}{ Precision of the quantile calculatons. Default is 10000. This 
  means 10000 points are taken from 0 to 1.}
  \item{name}{ Name of the data set, added to the title of plot if \code{main}
  is missing. }
  \item{corner}{ Can be \code{"bottomright"}, \code{"bottom"}, 
  \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"},
   \code{"topright"}, \code{"right"}, \code{"center"} as in 
   \code{\link{legend}}. }
  \item{type}{ This can be "" or "str.qqplot", the first produces the raw 
  quantiles and the second plot them on a straight line. Default is "".}
  \item{range}{ This is the range for which the quantiles are to be plotted. 
  Default is \code{c(0,1)}. }
  \item{xlab}{ x axis label, if left blank, then default is "Data". }
  \item{main}{ Title of the plot, if left blank, a default title will 
  be added. }
}

\value{
 A plot is given.
}
\author{ Steve Su }

\seealso{ \code{\link{qqplot.gld.bi}} }
\examples{
\donttest{
 set.seed(1000)

 junk<-rweibull(300,3,2)

# Fit the function using fun.data.fit.ml
 obj.fit1.ml<-fun.data.fit.ml(junk)

# Do a quantile plot on the raw quantiles
 qqplot.gld(junk,obj.fit1.ml[,1],param="rs",name="RS ML fit")

# Or a qq plot to examine deviation from straight line
 qqplot.gld(junk,obj.fit1.ml[,1],param="rs",name="RS ML fit",type="str.qqplot")
}}
\keyword{ hplot }
