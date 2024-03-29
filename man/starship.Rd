\name{starship}
\alias{starship}
\title{Carry out the ``starship'' estimation method for the generalised
lambda distribution} 
\description{
Calculates estimates for the FMKL parameterisation of the generalised lambda
distribution on the basis of data, using the starship method.
The starship method is built on the fact that the 
generalised lambda distribution 
is a transformation of the uniform distribution.  This method finds the 
parameters that transform the data closest to the uniform distribution.
This function uses a grid-based search to find a suitable starting point (using
\code{\link{starship.adaptivegrid}}) then uses \code{\link{optim}} to find
the parameters that do this.
}
\usage{
starship(data, optim.method = "Nelder-Mead", initgrid = NULL, param="FMKL",
optim.control=NULL)
}
\arguments{
  \item{data}{Data to be fitted, as a vector}
  \item{optim.method}{Optimisation method for \code{\link{optim}} to use, 
defaults to Nelder-Mead}
  \item{initgrid}{Grid of values of \eqn{\lambda_3}{lambda 3} and
\eqn{\lambda_4}{lambda 4}
to try, in \code{\link{starship.adaptivegrid}}.  This should be a list with
elements, 
  \code{lcvect}, a vector of values for \eqn{\lambda_3}{lambda 3},
  \code{ldvect}, a vector of values for \eqn{\lambda_4}{lambda 4} and
  \code{levect}, a vector of values for \eqn{\lambda_5}{lambda 5} 
  (\code{levect} is only required if \code{param} is \code{fm5}).
  
If it is left as NULL, the default grid depends on the parameterisation.
For \code{fmkl}, both lcvect and ldvect default to:
\tabular{ccccccccccc}{-1.5\tab -1\tab -.5\tab -.1\tab 0\tab .1\tab .2\tab .4\tab 
.8\tab 1\tab 1.5}
(\code{levect} is NULL).

For \code{rs}, both lcvect and ldvect default to:
\tabular{cccccc}{.1\tab .2\tab .4\tab .8\tab 1\tab 1.5}
(\code{levect} is NULL).

For \code{fm5}, both lcvect and ldvect default to:
\tabular{ccccccccccc}{-1.5\tab -1\tab -.5\tab -.1\tab 0\tab .1\tab .2\tab .4\tab 
.8\tab 1\tab 1.5}
and \code{levect} defaults to:
\tabular{ccccc}{-0.5\tab 0.25\tab 0\tab 0.25\tab 0.5}
	}
\item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}
  \code{fm5} uses the 5 parameter version of the FMKL parameterisation 
  (paper to appear)}
\item{optim.control}{List of options for the optimisation step.  See
\code{\link{optim}} for details.  If left as NULL, the parscale 
control is set to scale \eqn{\lambda_1}{lambda 1}
and \eqn{\lambda_2}{lambda 2} by the absolute value of their starting points.
}
}
\details{
The starship method is described in King and MacGillivray, 1999 (see
references). It is built on the fact that the 
generalised lambda distribution 
is a transformation of the uniform distribution.  Thus the inverse of this
transformation is the distribution function for the gld.  The starship method
applies different values of the parameters of the distribution to the
distribution function, calculates the depths \emph{q} corresponding to the data
and chooses the parameters that make the depths closest to a uniform
distribution.

The closeness to the uniform is assessed by calculating the Anderson-Darling
goodness-of-fit test on the transformed data against the uniform, for a
sample of size \code{length(data)}.

This is implemented in 2 stages in this function.  First a grid search is
carried out, over a small number of possible parameter values
(see \code{\link{starship.adaptivegrid}} for details).  Then the minimum from
this search is given as a starting point for an optimisation of the 
Anderson-Darling value using optim, with method given by \code{optim.method}

See references for details on 
parameterisations.
}
\value{
Returns a list, with 
  \item{lambda}{A vector of length 4, giving 
the estimated parameters, in order,
  \eqn{\lambda_1}{lambda 1} - location parameter
  \eqn{\lambda_2}{lambda 2} - scale parameter
  \eqn{\lambda_3}{lambda 3} - first shape parameter
  \eqn{\lambda_4}{lambda 4} - second shape parameter }
  \item{grid.results}{output from the grid search - see
\code{\link{starship.adaptivegrid}} for details}
  \item{optim}{output from the optim search -
\code{\link{optim}} for details}
}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.
 
King, R.A.R. & MacGillivray, H. L. (1999), \emph{A starship method for
fitting the generalised \eqn{\lambda}{lambda} distributions},
Australian and New Zealand Journal of
Statistics \bold{41}, 353--374

Owen, D. B. (1988), \emph{The starship}, Communications in Statistics - 
Computation and Simulation \bold{17}, 315--323.


}
\author{Robert King, Darren Wraith}
\seealso{ 
\code{\link{starship.adaptivegrid}}, 
\code{\link{starship.obj}} 
}
\examples{
\donttest{
 data <- rgl(100,0,1,.2,.2)
 starship(data,optim.method="Nelder-Mead",initgrid=list(lcvect=(0:4)/10,
 ldvect=(0:4)/10))
}}
\keyword{distribution}
