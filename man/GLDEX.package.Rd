\name{GLDEX-package}
\alias{GLDEX-package}
\alias{GLDEX}
\docType{package}
\title{
This package fits RS and FMKL generalised lambda distributions using 
various methods. It also provides functions for fitting bimodal 
distributions using mixtures of generalised lambda distributions.
}
\description{
The fitting algorithms considered in this package have two major objectives. One 
is to provide a smoothing device to fit distributions to data using the weight 
and unweighted discretised approach based on the bin width of the histogram. The 
other is to provide a definitive fit to the data set using the maximum 
likelihood estimation.

Copyright Information: To ensure the stability of this package, this package 
ports other functions from other open sourced packages directly so that any 
changes in other packages will not cause this package to malfunction.

All functions obtained from other sources have been acknowledged by the author 
in the authorship or the description sections of the help files and they are 
freely available online for all to use. Please contact the author 
for any copyright issues.

Specifically the following functions have been modified from R:

hist.su, ks.gof, pretty.su 

The following functions are taken from other open source packages in R:

runif.pseudo, rnorm.pseudo, runif.halton, rnorm.halton, runif.sobol, rnorm.sobol 
by Diethelm Wuertz distributed under GPL.

digitsBase, QUnif and sHalton written by Martin Maechler distributed under GPL.

dgl, pgl, rgl, qgl, starship.adpativegrid, starship.obj and starship written by 
Robert King and some functions modified by Steve Su distributed under GPL.

Lmoments and t1lmoments written by Juha Karvanen distributed under GPL.

}
\details{
\tabular{ll}{
Package: \tab GLDEX\cr
Type: \tab Package\cr
Version: \tab 2.0.0.0\cr
Date: \tab 2014-03-25\cr
License: \tab GPL (>= 3.0)\cr
}
This package allows a direct fitting method onto the data set using 
\code{\link{fun.RMFMKL.ml}}, \code{\link{fun.RMFMKL.ml.m}}, \code{\link{fun.RMFMKL.hs}}, 
\code{\link{fun.RMFMKL.hs.nw}}, \code{\link{fun.RPRS.ml}}, \code{\link{fun.RPRS.ml.m}},
\code{\link{fun.RPRS.hs}}, \code{\link{fun.RPRS.hs.nw}},
\code{\link{fun.RMFMKL.qs}}, \code{\link{fun.RPRS.qs}}, 
\code{\link{fun.RMFMKL.mm}}, \code{\link{fun.RPRS.mm}}, 
\code{\link{fun.RMFMKL.lm}}, \code{\link{fun.RPRS.lm}} and in the case of 
bimodal data set: \code{\link{fun.auto.bimodal.qs}},
\code{\link{fun.auto.bimodal.ml}}, 
\code{\link{fun.auto.bimodal.pml}} functions.
The resulting fits can be graphically gauged by \code{\link{fun.plot.fit}} or 
\code{\link{fun.plot.fit.bm}} (for bimodal data), or examined by numbers
using the Kolmogorov-Smirnoff resample tests (\code{\link{fun.diag.ks.g}}) and 
\code{\link{fun.diag.ks.g.bimodal}}). For unimodal data fits, it is also
possible to compare the mean, variance, skewness and kurtosis of the fitted 
distribution with the data set using \code{\link{fun.comp.moments.ml}} and 
\code{\link{fun.comp.moments.ml.2}} functions. Similarly, for bimodal data fits,
this is done via \code{\link{fun.theo.bi.mv.gld}} and \code{\link{fun.moments.r}}. 
Additionally, L moments for single generalised lambda distribution can be 
obtained using \code{\link{fun.lm.theo.gld}}. For graphical display of 
goodness of fit, quantile plots can be used, these can be done using
\code{\link{qqplot.gld}} and \code{\link{qqplot.gld.bi}} for univariate 
and bimodal generalised lambda distribution fits respectively.
}
\author{
Steve Su <allegro.su@gmail.com>
}
\references{

Asquith, W. (2007), L-moments and TL-moments of the generalized lambda 
distribution, Computational Statistics and Data Analysis 51(9): 4484-4496. 

Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), A study of the 
generalized tukey lambda family, Communications in Statistics - Theory and 
Methods *17*, 3547-3567. 

Gilchrist, Warren G. (2000), Statistical Modelling with Quantile Functions, 
Chapman & Hall 

Karian, Z.A., Dudewicz, E.J., and McDonald, P. (1996), The extended generalized 
lambda distribution system for fitting distributions to data: history, 
completion of theory, tables, applications, the ``Final Word'' on Moment fits, 
Communications in Statistics - Simulation and Computation *25*, 611-642. 

Karian, Zaven A. and Dudewicz, Edward J. (2000), Fitting statistical 
distributions: the Generalized Lambda Distribution and Generalized Bootstrap 
methods, Chapman & Hall 

Karvanen, J. and A. Nuutinen (2008), 
Characterizing the generalized lambda distribution by L-moments, 
Computational Statistics and Data Analysis 52(4): 1971-1983. 

King, R.A.R. & MacGillivray, H. L. (1999), A starship method for fitting the 
generalised lambda distributions, Australian and New Zealand Journal of 
Statistics, 41, 353-374 

Ramberg, J. S. & Schmeiser, B. W. (1974), An approximate method for generating 
asymmetric random variables, Communications of the ACM *17*, 78-82. 

Su, S. (2005). A Discretized Approach to Flexibly Fit Generalized Lambda 
Distributions to Data. Journal of Modern Applied Statistical Methods (November): 
408-424.

Su, S. (2007). Nmerical Maximum Log Likelihood Estimation for Generalized Lambda 
Distributions.  Computational Statistics and Data Analysis: *51*, 8, 3983-3998.

Su, S. (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. Journal of 
Statistical Software: *21* 9.

Su, S. (2009). Confidence Intervals for Quantiles Using Generalized Lambda 
Distributions. Computational Statistics and Data Analysis *53*, 9, 3324-3333.

Su, S. (2010). Chapter 14: Fitting GLDs and Mixture of GLDs to Data using 
Quantile Matching Method. Handbook of Distribution Fitting Methods with R. 
(Karian and Dudewicz) 557-583.

Su, S. (2010). Chapter 15: Fitting GLD to data using GLDEX 1.0.4 in R. Handbook 
of Distribution Fitting Methods with R. (Karian and Dudewicz) 585-608.

Su, S. (In Press) "Flexible Parametric Quantile Regression Model" Statistics &
Computing  

}
\seealso{
\pkg{GLDreg} package in R for GLD regression models.
}
\examples{

#### Univariate distributions example:

set.seed(1000)

junk<-rweibull(300,3,2)

### A faster ML estimtion 

# junk.fit1<-fun.RMFMKL.ml.m(junk)
# junk.fit2<-fun.RPRS.ml.m(junk)

# qqplot.gld(junk.fit1,data=junk,param="fmkl")
# qqplot.gld(junk.fit2,data=junk,param="rs")

### Using discretised approach, with 50 classes

## Using discretised method with weights
# obj.fit1.hs<-fun.data.fit.hs(junk)

## Plot the resulting fit. The fun.plot.fit function also works for singular fits 
## such as those by fun.plot.fit(obj.fit1.ml,junk,nclass=50,
## param=c("rs","fmkl","fmkl"),xlab="x")

# fun.plot.fit(obj.fit1.hs,junk,nclass=50,param=c("rs","fmkl"),xlab="x")

## Compare the mean, variance, skewness and kurtosis of the fitted distribution 
## with actual data

# fun.theo.mv.gld(obj.fit1.hs[1,1],obj.fit1.hs[2,1],obj.fit1.hs[3,1],
# obj.fit1.hs[4,1],param="rs")
# fun.theo.mv.gld(obj.fit1.hs[1,2],obj.fit1.hs[2,2],obj.fit1.hs[3,2],
# obj.fit1.hs[4,2],param="fmkl")
# fun.moments.r(junk)

## Conduct resample KS tests
# fun.diag.ks.g(obj.fit1.hs[,1],junk,param="rs")
# fun.diag.ks.g(obj.fit1.hs[,2],junk,param="fmkl")

### Try another fit, say 15 classes

# obj.fit2.hs<-fun.data.fit.hs(junk,rs.default="N",fmkl.default="N",no.c.rs = 15, 
# no.c.fmkl = 15)

# fun.plot.fit(obj.fit2.hs,junk,nclass=50,param=c("rs","fmkl"),xlab="x")

# fun.theo.mv.gld(obj.fit2.hs[1,1],obj.fit2.hs[2,1],obj.fit2.hs[3,1],
# obj.fit2.hs[4,1],param="rs")
# fun.theo.mv.gld(obj.fit2.hs[1,2],obj.fit2.hs[2,2],obj.fit2.hs[3,2],
# obj.fit2.hs[4,2],param="fmkl")
# fun.moments.r(junk)

# fun.diag.ks.g(obj.fit2.hs[,1],junk,param="rs")
# fun.diag.ks.g(obj.fit2.hs[,2],junk,param="fmkl")

### Uses the maximum likelihood estimation method

## Fit the function using fun.data.fit.ml
# obj.fit1.ml<-fun.data.fit.ml(junk)

## Plot the resulting fit
# fun.plot.fit(obj.fit1.ml,junk,nclass=50,param=c("rs","fmkl","fmkl"),xlab="x",
# name=".ML")

## Compare the mean, variance, skewness and kurtosis of the fitted distribution 
## with actual data
# fun.comp.moments.ml(obj.fit1.ml,junk)

## Do a quantile plot

# qqplot.gld(junk,obj.fit1.ml[,1],param="rs",name="RS ML fit")

## Run a KS resample test on the resulting fit

# fun.diag2(obj.fit1.ml,junk,1000)

## It is possible to use say fun.data.fit.ml(junk,rs.leap=409,fmkl.leap=409,
## FUN="QUnif") to find solution under a different set of low discrepancy number 
## generators.

#### Bimodal distributions example:

## Fitting mixture of generalised lambda distributions on the data set using both 
## the maximum likelihood and partition maximum likelihood and plot the resulting 
## fits

# par(mfrow=c(2,1))

# junk<-fun.auto.bimodal.ml(faithful[,1],per.of.mix=0.01,clustering.m=clara,
# init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1.5),init2=c(-0.25,1.5),
# leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# name="Maximum likelihood using",xlab="faithful1",param.vec=c("rs","fmkl"))

## Do a quantile plot

# qqplot.gld.bi(faithful[,1],junk$par,param1="rs",param2="fmkl",
# name="RS FMKL ML fit",range=c(0.001,0.999))

# junk<-fun.auto.bimodal.pml(faithful[,1],clustering.m=clara,init1.sel="rprs",
# init2.sel="rmfmkl",init1=c(-1.5,1.5),init2=c(-0.25,1.5),leap1=3,leap2=3)
# fun.plot.fit.bm(nclass=50,fit.obj=junk,data=faithful[,1],
# name="Partition maximum likelihood using",xlab="faithful1",
# param.vec=c("rs","fmkl"))

## Fit the faithful[,1] data from the dataset library using sobol sequence 
## generator for the first distribution fit and halton sequence for the second 
## distribution fit.

# fit1<-fun.auto.bimodal.ml(faithful[,1],init1.sel="rmfmkl",init2.sel="rmfmkl",
# init1=c(-0.25,1.5),init2=c(-0.25,1.5),leap1=3,leap2=3,fun1="runif.sobol",
# fun2="runif.halton")

## Run diagnostic KS tests

# fun.diag.ks.g.bimodal(fit1$par[1:4],fit1$par[5:8],prop1=fit1$par[9],
# data=faithful[,1],param1="fmkl",param2="fmkl")

## Find the theoretical moments of the fit

# fun.theo.bi.mv.gld(fit1$par[1],fit1$par[2],fit1$par[3],fit1$par[4],"fmkl",
# fit1$par[5],fit1$par[6],fit1$par[7],fit1$par[8],"fmkl",fit1$par[9])

## Compare this with the empirical moments from the data set.

# fun.moments.r(faithful[,1])

## Do a quantile plot

# qqplot.gld.bi(faithful[,1],fit1$par,param1="fmkl",param2="fmkl",
# name="FMKL FMKL ML fit")

## Quantile matching method

## Fitting faithful data from the dataset library, with the clara clustering 
## regime. The first distribution is RS and the second distribution is fmkl. 
## The percentage of data mix is 1\%.

## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.qs(junk)

# fun.auto.bimodal.qs(faithful[,1],per.of.mix=0.01,clustering.m=clara,
# init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),
# leap1=3,leap2=3)

# L Moment matching

## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.lm(junk)

# Moment matching method

## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.mm(junk)

# Example on fitting mixture of normal distributions

#data1<-c(rnorm(1500,-1,2/3),rnorm(1500,1,2/3))

#junk<-fun.auto.bimodal.ml(data1,per.of.mix=0.01,clustering.m=data1>0,
#init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1.5),init2=c(-0.25,1.5),
#leap1=3,leap2=3)

#fun.plot.fit.bm(nclass=50,fit.obj=junk,data=data1,
#name="Maximum likelihood using",xlab="faithful1",param.vec=c("rs","fmkl"))

#qqplot.gld.bi(data1,junk$par,param1="rs",param2="fmkl",
#name="RS FMKL ML fit",range=c(0.001,0.999))

# Generate random observations from FMKL generalised lambda distributions with 
# parameters (1,2,3,4) and (4,3,2,1) with 50\% of data from each distribution.
fun.simu.bimodal(c(1,2,3,4),c(4,3,2,1),prop1=0.5,param1="fmkl",param2="fmkl")
}
\keyword{distribution} 
\keyword{smooth}