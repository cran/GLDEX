\name{fun.simu.bimodal}
\alias{fun.simu.bimodal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a mixture of two generalised lambda distributions. }
\description{
This function allows the user to simulate observations from a mixture of two 
generalised lambda distributions. It 
can be very useful for sensitivity analysis.
}
\usage{
fun.simu.bimodal(result1, result2, prop1, prop2, len = 1000, 
no.test = 1000, param1, param2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result1}{ A vector comprising four values for the first generalised 
  lambda distribution. }
  \item{result2}{ A vector comprising four values for the second generalised 
  lambda distribution. }
  \item{prop1}{ Proportion of the first generalised lambda distribution }
  \item{prop2}{ 1-prop1, this can be left unspecified. }
  \item{len}{ Length of object for each simulation run. }
  \item{no.test}{ Number of simulation run. }
  \item{param1}{ This can be \code{"rs"} or \code{"fmkl"}, specifying the type 
  of the first generalised lambda distribution. }
  \item{param2}{ This can be \code{"rs"} or \code{"fmkl"}, specifying the type 
  of the second generalised lambda distribution. }
}
\details{
The length of object in \code{len} means how many observations should
be generated in each simulation run, with the number of simulation runs governed 
by \code{no.test}.
}
\value{
A list with length equal to the number of simulation runs. Each subset of the 
list has random observations equal to the the number specified in
\code{len}.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.theo.bi.mv.gld}}, \code{\link{fun.moments.bimodal}}, 
\code{\link{fun.rawmoments}} }
\examples{
# Generate random observations from FMKL generalised lambda distributions with 
# parameters (1,2,3,4) and (4,3,2,1) with 50\% of data from each distribution.
junk<-fun.simu.bimodal(c(1,2,3,4),c(4,3,2,1),prop1=0.5,param1="fmkl",
param2="fmkl")

# Calculate the maximum number from each simulation run
sapply(junk,max)

# Calculate the median from each simulation run
sapply(junk,median)
}
\keyword{datagen}