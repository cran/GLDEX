\name{qqplot.gld.bi}
\alias{qqplot.gld.bi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Do a quantile plot on the bimodal distribution fits. }
\description{
  This plots the theoretical and actual data quantiles to allow the user to 
  examine the adequacy of two gld distributions mixture fit.
}
\usage{
qqplot.gld.bi(data, fit, param1, param2, len = 10000, name = "", 
corner = "topleft",type="",range=c(0,1),xlab="",main="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{ Data fitted. }
  \item{fit}{ Parameters of distribution fit. }
  \item{param1}{ Can be either \code{"rs"} or \code{"fmkl"}. }
  \item{param2}{ Can be either \code{"rs"} or \code{"fmkl"}. }
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
  \item{xlab}{ x axis label, if left blank, then default is "Data" }
  \item{main}{ Title of the plot, if left blank, a default title will 
  be added. }
}
\value{
 A plot is given.
}
\author{ Steve Su }
\seealso{ \code{\link{qqplot.gld}} }
\examples{
# set.seed(1000)

# junk<-rweibull(300,3,2)

## Fitting mixture of generalised lambda distributions on the data set using 
## both the maximum likelihood and partition maximum likelihood and plot the 
## resulting fits
# junk<-fun.auto.bimodal.ml(faithful[,1],per.of.mix=0.1,clustering.m=clara,
# init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1.5),init2=c(-0.25,1.5),
# leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# name="Maximum likelihood using",xlab="faithful1",param.vec=c("rs","fmkl"))

## Do a quantile plot on the raw quantiles
# qqplot.gld.bi(faithful[,1],junk$par,param1="rs",param2="fmkl",
# name="RS FMKL ML fit")

## Or a qq plot to examine deviation from straight line
# qqplot.gld.bi(faithful[,1],junk$par,param1="rs",param2="fmkl",
# name="RS FMKL ML fit",type="str.qqplot")
}


\keyword{ hplot }
