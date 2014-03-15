\name{GLD functions}
\alias{dgl}
\alias{pgl}
\alias{qdgl}
\alias{qgl}
\alias{rgl}
\title{ The Generalised Lambda Distribution Family }
\description{
     Density, quantile density, distribution function, quantile 
     function and random generation for the generalised lambda
     distribution  (also known as the asymmetric lambda, or Tukey
     lambda).  Works for  both the \code{"fmkl"} and \code{"rs"}
     parameterisations. These functions originate from the \pkg{gld} 
     library by Robert King and they are modified in this package
     to allow greater functionality and adaptability to new fitting
     methods. It does not give an error message for invalid distributions
     but will return NAs instead. To allow comparability with the pkg(gld) 
     package, this package uses the same notation and description as those 
     written by Robert King.
}
\usage{
     dgl(x, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
       param = "fmkl", inverse.eps = 1e-08, 
       max.iterations = 500)
     pgl(q, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
       param = "fmkl",  inverse.eps = 1e-08, 
       max.iterations = 500)
     qgl(p, lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL,
       param = "fmkl")
     rgl(n, lambda1=0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL,
       param = "fmkl")
}

\arguments{
  \item{x}{ Vector of actual values for \code{dgl} }
  \item{p}{ Vector of probabilities for \code{qgl} or \code{qdgl} } 
  \item{q}{ Vector of quantiles for \code{pgl} } 
  \item{n}{ Number of observations to be generated for \code{rgl} }
  \item{lambda1}{ This can be either a single numeric value or a vector. If it 
  is a vector, it must be of length 4 for parameterisations  \code{"fmkl"} or 
  \code{"rs"} and of length 5 for parameterisation \code{"fm5"} and the other 
  'lambda' arguments must be left as \code{NULL}. The numbering of the lambda 
  parameters for the \code{"fmkl"} parameterisation is different to that used by 
  Freimer, Mudholkar, Kollia and Lin (1988).
 }
  \item{lambda2}{ Scale parameter }
  \item{lambda3}{ First shape parameter }
  \item{lambda4}{ Second shape parameter }
  \item{param}{ \code{"fmkl"} or \code{"fkml"} uses Freimer, Kollia, Mudholkar,
          and Lin (1988) and it is the default setting. \code{"rs"} uses 
          Ramberg and Schmeiser (1974) }
  \item{inverse.eps}{ Accuracy of calculation for the numerical determination of
           F(x), defaults to 1e-8 }
  \item{max.iterations}{ Maximum number of iterations in the numerical 
          determination of F(x), defaults to 500 }
}
\details{
     The generalised lambda distribution, also known as the asymmetric
     lambda,  or Tukey lambda distribution, is a distribution with a
     wide range of shapes.   The distribution is defined by its
     quantile function, the inverse of the distribution function.  The
     'gld' package implements three parameterisations of the
     distribution.   The default parameterisation (the FMKL) is that
     due to Freimer  Mudholkar, Kollia and Lin (1988) (see references
     below), with a quantile  function:
     
     \deqn{F^{-1}(u)= \lambda_1 + { { \frac{u^{\lambda_3}-1}{\lambda_3} -    
     \frac{(1-u)^{\lambda_4}-1}{\lambda_4} } \over \lambda_2 }}{F inverse (u) = 
     lambda1 + ( (u^lambda3 -1)/lambda3 - ((1-u)^lambda4-1) 
      /lambda4 ) / lambda 2}
      for \eqn{\lambda_2 > 0}{lambda2 >0}.

     A second parameterisation, the RS, chosen by setting \code{param="rs"} is 
     that due to Ramberg and Schmeiser (1974), with the quantile function:
     \deqn{F^{-1}(u)= \lambda_1 + \frac{u^{\lambda_3} - (1-u)^{\lambda_4}}   
     {\lambda_2} }{F inverse (u) = lambda1 + ( u^lambda3 - (1-u)^lambda4 ) / 
     lambda 2}

     This parameterisation has a complex series of rules determining
     which values  of the parameters produce valid statistical
     distributions. See \code{gl.check.lambda} for details.}

\value{
  
     \item{dgl}{ gives the density (based on the quantile density and a 
     numerical solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}, }

     \item{qdgl}{ gives the quantile density, }

     \item{pgl}{ gives the distribution function (based on a numerical
     solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}, }

     \item{qgl}{ gives the quantile function, and }

     \item{rgl}{ generates random observations. }
}
\references{
     Freimer, M., Kollia, G., Mudholkar, G. S. & Lin, C. T. (1988),  A
     study of the generalized tukey lambda family, Communications  in
     Statistics - Theory and Methods *17*, 3547-3567.

     Gilchrist, Warren G. (2000), Statistical Modelling with Quantile
     Functions, Chapman & Hall

     Karian, Z.A., Dudewicz, E.J., and McDonald, P. (1996), The
     extended  generalized lambda distribution system for fitting
     distributions to data:  history, completion of theory, tables,
     applications, the ``Final Word''  on Moment fits, Communications
     in Statistics - Simulation and Computation  *25*, 611-642.

     Karian,  Zaven A. and Dudewicz, Edward J. (2000), Fitting
     statistical distributions: the Generalized Lambda Distribution and
     Generalized Bootstrap methods, Chapman & Hall

     Ramberg, J. S. & Schmeiser, B. W. (1974), An approximate method
     for generating asymmetric random variables, Communications of the
     ACM *17*,  78-82. }
\examples{
 qgl(seq(0,1,0.02),0,1,0.123,-4.3)
 pgl(seq(-2,2,0.2),0,1,-.1,-.2,param="fmkl",inverse.eps=.Machine$double.eps)
}
\keyword{distribution}