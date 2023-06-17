\name{EGPHazard}
\alias{EGPHazard}

\title{
Calculate Transition Hazards for an ERGM Generating Process
}

\description{
Given an EGP and an initial state, calculate the transition rates to one or more neighboring states.
}

\usage{
EGPHazard(form, coef, toggles = NULL, rate.factor = 1, process = c("LERGM",
    "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"))
}

\arguments{
  \item{form}{
an ERGM formula for the EGP (or a list with \code{formation} and \code{dissolution} formulas, for \code{CSTERGM} processes).  The left-hand side is used as the current state when computing transition rates.
}
  \item{coef}{
a vector of EGP coefficients, or a list of vectors with elements \code{formation} and \code{dissolution} for \code{CSTERGM} and variants.
}
  \item{toggles}{
edge variables to evaluate; passing \code{"all"} or \code{NULL} leads to all edge variables being evaluated, \code{"edges"} evaluates only dissolution events, \code{"nulls"} evaluates only formation events, and passing a two-column matrix of IDs (tail, head order) evaluates the selected dyads.
}
  \item{rate.factor}{
rate or pacing factor (sets the time scale).
}
  \item{process}{
the ERGM generating process to use.
}
}

\details{
An ERGM generating process (EGP) is a continuous time graph process with an equilibrium distribution having a known ERGM form.  See \code{\link{ergmgp}} for an overview of EGPs, including the specifications supported here.

\code{EGPHazard} calculates the log transition rates (i.e., hazards) from an initial or current state (specified by the left-hand side of the input formula) to one or more target states.  These states are specified by the edge variables whose states would change (often called \dQuote{toggles} in ERGM nomenclature).  By default, all possible transitions are evaluated; this can also be obtained by setting \code{toggles=="all"}.  Dissolution rates for all current edges can be obtained by setting \code{toggles=="edges"}, and formation rates for all current nulls can be obtained by setting \code{toggles=="nulls"}.  Otherwise, the \code{toggles} argument expects a two-column matrix of tail and head vertex IDs indicating the edge variables to be evaluated.  Note that only instantaneous rates from the origin state are computed; toggles are not cumulative.

EGP specifications are as per \code{\link{simEGP}}.  Transition rates for all currently implemented EGPs follow the specifications of Butts (2023), with the trivial addition of a pacing constant for all families (which simply sets the timescale).
}

\value{
a matrix containing the toggles, indicators for whether each event would have been a formation event, and the log event hazards (one row per toggle).
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

#Compute all rates under a LERGM
lr <- EGPHazard(net ~ edges + esp(0), coef = co, process = "LERGM")
head(lr)  #Sender, receiver, formation (1=yes), log rate

#Use a toggle matrix to obtain the same outcome
lrt <- EGPHazard(net ~ edges + esp(0), coef = co, toggles = lr[,1:2],
    process = "LERGM")
all(lrt == lr)   #TRUE

#Examine edge dissolution rates
ldissr <- EGPHazard(net ~ edges + esp(0), coef = co, toggles = "edges", 
    process = "LERGM")
a <- function(z){(z-min(z))/diff(range(z))}
plot(net, edge.col = hsv(a(ldissr[,4])*0.6))  #Blue=fast, red=slow

}

 \keyword{ graphs }
 \keyword{ models }

