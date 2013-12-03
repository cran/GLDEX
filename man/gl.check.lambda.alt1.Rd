\name{gl.check.lambda.alt1}
\alias{gl.check.lambda.alt1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Checks whether the parameters provided constitute a valid generalised 
lambda distribution. }
\description{
A replacement to the gl.check.lambda function in \pkg{gld} package.
}
\usage{
gl.check.lambda.alt1(l1, l2 = NULL, l3 = NULL, l4 = NULL, 
param = "fmkl", vect = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l1}{ Lambda 1. }
  \item{l2}{ Lambda 2. }
  \item{l3}{ Lambda 3. }
  \item{l4}{ Lambda 4. }
  \item{param}{ \code{"rs"} or \code{"fmkl"} generalised lambda distribution. }
  \item{vect}{ A logical, set this to \code{TRUE} if the parameters are given in 
  the vector form (it turns off checking of the format of 'lambdas' and the 
  other lambda arguments
 }
}
\details{
This is a modified \code{gl.check.lambda} function in replace of \pkg{gld} 
library's \code{gl.check.lambda} function to allow for 5 parameters FMKL 
distributions and vector input of parameter values into this function.
}
\value{
A logical value, \code{TRUE} or \code{FALSE}. \code{TRUE} indicates the 
parameters given is a valid probability distribution.
}
\references{  Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), A
     study of the generalized tukey lambda family, Communications in
     Statistics - Theory and Methods *17*, 3547-3567.

     Karian, Z.E., Dudewicz, E.J., and McDonald, P. (1996), The
     extended  generalized lambda distribution system for fitting
     distributions to data: history, completion of theory, tables,
     applications, the ``Final Word''  on Moment fits, Communications
     in Statistics - Simulation and Computation  *25*, 611-642.

     Ramberg, J. S. & Schmeiser, B. W. (1974), An approximate method
     for generating asymmetric random variables, Communications of the
     ACM *17*, 78-82.
}
\author{ Steve Su }
\seealso{ \code{\link{gl.check.lambda.alt}} }
\examples{
gl.check.lambda.alt1(c(0,1,.23,4.5),param="fmkl",vect=TRUE) 
## TRUE, Using vector input of parameter values.
gl.check.lambda.alt1(0,-1,.23,4.5,param="fmkl") ## FALSE 
gl.check.lambda.alt1(0,1,0.5,-0.5,param="rs") ## FALSE
}
\keyword{distribution} 