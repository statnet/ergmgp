\name{EGPRateEst}
\alias{EGPRateEst}

\title{
Estimate Event Rates for an ERGM Generating Process
}

\description{
Given an EGP, estimate either the expected time required for a specified number of transitions to occur, or the expected number of transitions within a specified time period.
}

\usage{
EGPRateEst(formula, coef, process = c("LERGM", "CRSAOM", "CI", "DS",
    "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), time.target = NULL,
    event.target = NULL, reps = 25, cores = 1, rate.factor = 1, 
    verbose = FALSE, ...)
}

\arguments{
  \item{formula}{
an ERGM formula for the EGP (or a list with \code{formation} and \code{dissolution} formulas, for \code{CSTERGM} processes).  The left-hand side is used as the initial state.
}
  \item{coef}{
a vector of EGP coefficients, or a list of vectors with elements \code{formation} and \code{dissolution} for \code{CSTERGM} and variants.
}
  \item{process}{
the ERGM generating process to use.
}
  \item{time.target}{
if specified, the length of the time period for which trajectories should be simulated (in which case the estimand is the number of events within this period). 
}
  \item{event.target}{
if specified, the number of transition events over which trajectories should be simulated (in which case the estimand is the time required for this number of events to be realized). 

}
  \item{reps}{
number of replicate trajectories to use.
}
  \item{cores}{
number of cores to use for simultaneous simulation of trajectories.
}
  \item{rate.factor}{
rate or pacing factor (sets the time scale).
}
  \item{verbose}{
logical; show progress information?
}
  \item{\dots}{
additional arguments to \code{\link{simEGP}}.
}
}

\details{
This function can be used to estimate the expected amount of time needed for a specific number of transitions to be realized (in which case \code{event.target} should be supplied) or the expected number of transition events occurring within a specified time period (in which case \code{time.target} should be supplied).  Either way, one of \code{time.target} and \code{event.target} must be given.  The function works by simulating \code{reps} trajectories (using \code{simEGP}) for the specified time/number of events, and returning the mean outcome (along with some other associated statistics).

A typical use case for this function is to calibrate the simulation time needed to obtain a reasonable number of transitions from some starting point (e.g., to ensure burn-in).  Simply simulating a fixed number of transition events will result in a biased system state; however, one can avoid this problem by using this function to determine the average duration needed for the desired number of events to be realized, and then using this duration as a stopping rule for subsequent simulations.  Alternately, another use is to estimate the rate at which events accrue, e.g. to estimate compute time or memory requirements for a longer simulation study.  Some processes are particularly prone to entering regimes in which they produce very large numbers of events per unit phenomenological time, and it can be useful to identify this issue before committing resources to simulating a long trajectory.

Note that, at present, all trajectories have the same starting point (the network on the left-hand side of the input formula).  They are hence coupled by the initial condition (despite being otherwise independent).  When equilibrium estimates from short sequences are desired, it may be wise to call this function more than once with different starting networks and integrate the results.
}

\value{
A vector containing the mean outcome (time or event count), its standard error, the standard deviation of the outcome, and the number of replicates used.
}

\references{
Butts, Carter T.  (2023).  \dQuote{Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.} \emph{Journal of Mathematical Sociology}.  \doi{10.1080/0022250X.2023.2180001}
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}

\seealso{
\code{\link{ergmgp}} for information on EGPs, \code{\link[ergm]{ergm}} for information on ERGM specifications, \code{\link{simEGP}}
}

\examples{
#Simulate a small network with triadic dependence
n <- 25
set.seed(1331)
co <- c(log(2.5/(n-3.5)), -0.75)
net <- simulate(network.initialize(n, directed=FALSE) ~ edges + esp(0),
    coef = co)

#Estimate the time needed for 500 events in a LERGM
etime <- EGPRateEst(net ~ edges + esp(0), coef = co, process = "LERGM",
    event.target = 500)
etime

#Estimate the mean number of events in the above time
eevents <- EGPRateEst(net ~ edges + esp(0), coef = co, process = "LERGM",
    time.target = etime[1])
eevents   #Expectation should be close to 500

}

 \keyword{ graphs }
 \keyword{ models }


