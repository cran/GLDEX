\name{fun.plot.fit.bm}
\alias{fun.plot.fit.bm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting mixture of two generalised lambda distributions on the 
data set. }
\description{
This function is designed for mixture of two generalised lambda distributions
only.
}
\usage{
fun.plot.fit.bm(fit.obj, data, nclass = 50, xlab = "", name = "", main="", 
param.vec, ylab="Density")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{ Fitted object from \code{\link{fun.auto.bimodal.ml}}, 
  \code{\link{fun.auto.bimodal.pml}} }
  \item{data}{ Dataset to be plotted. }
  \item{nclass}{ Number of class of histogram, the default is 50. }
  \item{xlab}{ Label on the x axis. }
  \item{name}{ Legend, usually used to identify type of GLD used if \code{main}
  is provided. If \code{main} is not provided, then this is used in the title. }
  \item{main}{ Title of the graph. }
  \item{param.vec}{ A vector describing the type of generalised lambda 
  distribution used in the \code{fit.obj}. } 
  \item{ylab}{ Label on the y axis. }
  }
\value{
 A graphical output showing the data and the resulting distributional fits.
}

\author{ Steve Su }
\note{ If the distribution fits over fits the peak of the distribution, it can 
be difficult to see the actual data set.
}
\seealso{ \code{\link{fun.auto.bimodal.ml}}, \code{\link{fun.auto.bimodal.pml}}, 
\code{\link{fun.plot.fit}} }
\examples{

# par(mfrow=c(2,1))

## Fitting mixture of generalised lambda distributions on the data set using 
## both the maximum likelihood and partition maximum likelihood and plot 
## the resulting fits

# junk<-fun.auto.bimodal.ml(faithful[,1],per.of.mix=0.1,clustering.m=clara,
# init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),
# leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# name="Maximum likelihood using",xlab="faithful1",param.vec=c("rs","fmkl"))

# junk<-fun.auto.bimodal.pml(faithful[,1],clustering.m=clara,init1.sel="rprs",
# init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# name="Partition maximum likelihood using",xlab="faithful1",
# param.vec=c("rs","fmkl"))

# junk<-fun.auto.bimodal.ml(faithful[,1],per.of.mix=0.1,clustering.m=clara,
# init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),
# leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# main="Mixture distribution fit",
# name="RS and FMKL GLD",xlab="faithful1",param.vec=c("rs","fmkl"))

}
\keyword{hplot}