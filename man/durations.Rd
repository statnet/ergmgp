\name{durations}
\alias{durations}

\title{
Obtain Edge Spell Durations from an ERGM Generating Process Trajectory
}

\description{
Given an input trajectory (in networkDynamic form, or network form with additional attributes), return the set of all edge durations (along with censoring information, if desired).
}

\usage{
durations(net, censor = c("obs", "omit"), return.censoring = TRUE)
}


\arguments{
  \item{net}{
a \code{network} or \code{networkDynamic} object containing the trajectory information.
}
  \item{censor}{
how should censoring be handled?  (Currently, only returning observed spell lengths and omitting censored spells are supported.)
}
  \item{return.censoring}{
logical; return censoring information?
}
}

\details{
This function extracts information on edge spells (periods of time in which edges are present) from the input network, and returns the spell durations (optionally, together with censoring information).  The durations should not be assumed to be in any particular order; this function is generally invoked to examine duration distributions.

If \code{net} is a \code{networkDynamic} object, both spell and censoring information are extracted from its edge activities.  If \code{net} is a \code{network} object, then its \code{"LastChangeTime"} network attribute is used to obtain spell information.  (Both can be obtained from \code{\link{simEGP}} with appropriate settings.)  Currently, \code{network} objects with \code{"EventHistory"} attributes are not supported - use the \code{networkDynamic} output type to examine complete event histories.  For the \code{network} case, the observation period is assumed to span the interval from 0 to \code{net\%n\%"Time"} (so be sure that temporal offsets were not used if employing that data type).

Spells may be left-censored, right-censored, or both.  \code{censor=="obs"} results in lengths being reported as-is (subject to truncation to the observation period), and \code{censor=="omit"} results in censored spells being omitted.  Censoring indicators are also included when \code{return.censoring==TRUE}.  Note that if \code{"LastChangeTime"} information is being used, all spells are censored (we see only the onset times for edges that were present when the simulation was terminated), so the \code{"omit"} option will return a zero-length vector.  

When using \code{durations} to estimate equilibrium duration distributions, it is important to bear in mind that EGP trajectories stopped by event count are not terminated at a random time, and hence will provide biased estimates.  Consider using \code{\link{EGPRateEst}} to calibrate a reasonable simulation time, and sampling with a temporal stopping rule.
}

\value{
A vector of spell durations (order not guaranteed), or a matrix containing said durations and censoring indicators (0=uncensored, 1=right-censored, 2=left-censored, and 3=interval censored).
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}

\seealso{
\code{\link{simEGP}}
}

\examples{
\donttest{  #Examples are a bit slow, so not automatically run

#Generate a simple CD-CSTERGM trajectory; equilibrium mean outdegree
#is 2, dissolution rate is 1/3
set.seed(1331)
n <- 25
net <- simulate(network.initialize(n)~edges, coef=log(2/(n-3)))
traj <- simEGP(net~edges, coef=list(formation=log(2/(n-3)*1/3),
    dissolution=log(1/3)), time=5000, process="CDCSTERGM", 
    return.changetime=TRUE, verbose=FALSE)
network.edgecount(traj)/(n-1)             #Mean degree apx 2
dur <- durations(traj)                    #Get durations
head(dur)                                 #All are right-censored
mean(dur[,1])                             #Apx 3 (despite censoring)

#Repeat, but now using a networkDynamic object
set.seed(1331)
net <- simulate(network.initialize(n)~edges, coef=log(2/(n-3)))
traj <- simEGP(net~edges, coef=list(formation=log(2/(n-3)*1/3),
    dissolution=log(1/3)), time=500, process="CDCSTERGM", 
    return.networkDynamic=TRUE, verbose=FALSE)
slice <- traj \%t\% 499                     #Take a slice near the end
network.edgecount(slice)/(n-1)            #Mean degree apx 2
dur <- durations(traj)                    #Get durations
head(dur)                                 #More of a mix
mean(dur[,1])                             #Apx 3
hist(dur[,1], xlab="Time", main="Duration Distribution")  #Visualize

}
}

 \keyword{ manip }
 \keyword{ graphs }
 \keyword{ survival }

