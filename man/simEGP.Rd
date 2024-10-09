\name{simEGP}
\alias{simEGP}
\alias{simEGPTraj}

\title{
Simulate Trajectories from an ERGM Generating Process
}

\description{
Given an \code{\link[ergm]{ergm}} formula, simulate trajectories from a continuous time graph process having the specified ERGM as a limiting distribution.  A number of different processes are supported, and termination may be specified either by phenomenological time or event counts.
}

\usage{
simEGP(form, coef, events = 1, time = NULL, rate.factor = 1, 
    time.offset = 0, event.offset = 0, process = c("LERGM", "CRSAOM", 
    "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), 
    use.logtime = FALSE, return.changetime = FALSE, 
    changetime.offset = NULL, return.history = FALSE, 
    return.networkDynamic = FALSE, verbose = TRUE,  trace.interval = 100,
    ...)

simEGPTraj(form, coef, events = 1, time = NULL, checkpoints = 1, 
    rate.factor = 1, trajectories = 1, mc.cores = 1, 
    log.sampling = FALSE, process = c("LERGM", "CRSAOM", "CI", "DS",
    "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), use.logtime = FALSE,
    return.changetime = FALSE, return.history = FALSE, verbose = TRUE,
    trace.interval = 100, statsonly = FALSE, monitor = NULL) 
}

\arguments{
  \item{form}{
an \code{\link[ergm]{ergm}} formula defining terms for the EGP; the left-hand side must be a network object, whose properties are used to determine the state space.  For the \code{CSTERGM} process, a list containing two such formulas must be used, with named elements \code{formation} (for the formation model) and \code{dissolution} (for the dissolution model).
}
  \item{coef}{
vector of coefficients for the EGP; for the CSTERGMs, this should be a containing named elements \code{formation} and \code{dissolution}, each of which should be the coefficient vector for its respective model.
}
  \item{events}{
optionally, the number of simulated events to draw (if \code{time==NULL}); if \code{time} is specified, this is ignored.
}
  \item{time}{
optionally, the temporal length of the simulation; if not supplied, \code{events} is used instead to determine when to stop.
}
  \item{rate.factor}{
a multiplicative factor scaling the time evolution of the system; higher values correspond to faster dynamics.
}
  \item{time.offset}{
optionally, an initial \dQuote{clock} offset for the start of a trajectory; this allows time 0 (the start of the simulation interval) to be set to an arbitrary time point.  This is only used for book-keeping (e.g., when a trajectory is run as multiple segments), and does not affect e.g. the meaning of the \code{time} argument (which is always interpreted as units after the start time).
}
  \item{event.offset}{
optionally, an initial offset to the step or event count for the start of a trajectory (e.g., for trajectories being run in segments).  As with \code{time.offset}, this only affects book-keeping, and has no other effect.
}
  \item{process}{
the ERGM generating process to use (described below).
}
  \item{use.logtime}{
logical; internally, use logarithmic timescale?  This can potentially protect against overflow or underflow when rates are extreme, but can reduce precision and adds some overhead.
}
  \item{return.changetime}{
logical; should we return a matrix with the last update times for each edge variable as a network attribute?
}
  \item{changetime.offset}{
optionally, an \code{n x n} matrix of last change times (for trajectories being resumed in process).
}
  \item{return.history}{
logical; return the entire event history as a network attribute?
}
  \item{return.networkDynamic}{
logical; retain the entire event history and return as a \code{\link[networkDynamic]{networkDynamic}} object?
}
  \item{verbose}{
logical; provide trace messages regarding simulation progress?
}
  \item{trace.interval}{
for verbose output, the interval at which messages should be printed (in events).
}
  \item{checkpoints}{
number of checkpoints at which the trajectory should be sampled (in addition to the initial state).
}
  \item{trajectories}{
number of independent trajectories to simulate (all start from the seed network, but evolve independently).
}
  \item{mc.cores}{
number of cores to use when simulating trajectories.
}
  \item{log.sampling}{
logical; should time points to sample be logarithmically spaced?
}
  \item{statsonly}{
logical; should only network statistics be retained (and not the graphs themselves)?
}
  \item{monitor}{
optionally, an \code{\link[ergm]{ergm}} formula with additional statistics to track.
}
  \item{\dots}{
additional arguments (currently unused).
}
}

\details{
An ERGM generating process (EGP) is a continuous time graph process with an equilibrium distribution having a known ERGM form.  See \code{\link{ergmgp}} for an overview of EGPs, including the specifications supported here.

\code{simEGP} generates a single trajectory from an EGP, with the EGP being specified via its graph potential (as a \code{\link[ergm]{ergm}} formula or pair thereof and associated coefficients) and its initial state being given by the left-hand side of the input formula.  The trajectory length can be specified either in terms of the number of transitions to be simulated (\code{events}) or the length of the trajectory in phenomenological time (\code{time}); only the latter leads to the specified ERGM equilibrium (since event times are not \dQuote{random} times, stopping after a fixed number of events biases the final state).  If desired for bookkeeping purposes, an offset can be added to the simulation clock (which otherwise starts at 0), event count (likewise), and most recent change times (also likewise).  By default, the return value is a \code{\link[network]{network}} object containing the final graph state, with network attributes giving the final time (\code{"Time"}), event count (\code{"Events"}), ERGM potential (\code{"Potential"}).  A square matrix containing the time of the most recent transition experienced by each edge variable can be returned as a network attribute (\code{"LastChangeTime"}) if \code{return.changetime} is selected.  By default, the entire event history is not stored (as it can become extremely large).  However, if \code{return.history==TRUE}, a matrix containing the event history is saved and returned as a network attribute (\code{"EventHistory"}). Alternately, setting \code{return.networkDynamic=TRUE} will lead to the event history being stored and the entire trajectory being returned as a \code{\link[networkDynamic]{networkDynamic}} object, with edge activities set based on the observed transitions.  This format may be easier to use for visualization, or to query the state of the network at an arbitrary point in the trajectory.  The \code{\link{durations}} function can be used to extract edge durations from such objects, as well.

For models with extreme transition rates, the option \code{use.logtime} may be useful for avoiding overflow or underflow; this only affects internal calculation, and not reported event times.  Note that logscale calculations add some overhead, and may be less precise in some cases than the default, so this option is not suggested unless specifically needed.  

\code{simEGPTraj} is a wrapper for \code{simEGP}, which adds additional capabilities for simulation of multiple trajectories and/or sampling of longer trajectories.  Each returned trajectory contains the initial state, as well as \code{simEGP} output from \code{checkpoints} points along the trajectory (including the end).  The default behavior (\code{checkpoints==1}) returns the initial and final states.  Checkpoints are evenly spaced (with termination criteria indicated as per \code{simEGP}) by default, or logarithmically spaced if \code{log.sampling==TRUE}.  Multiple independent trajectories can be simulated by setting \code{trajectories>1}; these can be run in parallel by setting \code{mc.cores>1}.  If desired, the model statistics can be returned without the graph state by choosing \code{statsonly=TRUE}, and a one-sided \code{monitor} formula can likewise be used to calculate additional statistics if desired (with similar functionality to the \code{ergm} \code{simulate} method).  Otherwise, \code{network.list} objects are returned containing the states in the respective trajectories.

Simulation itself follows the discrete event approach described in Butts (2023).  Transition hazards are computed for all edge variables (making the scaling no better then \eqn{O(N^2)} for each update, and are used to draw both the next event and the event time.  Because the cost of computing each transition is unrelated to waiting time, this algorithm can be quite efficient at simulating long time periods when events are sparse (unlike, e.g. a discrete-time algorithm that updates in every period).  By turns, however, trajectories can become quite expensive (per unit phenomenological time) when event rates are high.  This issue is especially pronounced for the \code{CRSAOM} and \code{DS} processes, which can both generate very high transition rates in some cases.  Unless otherwise specified, event histories are not stored, so storage costs are by default unrelated to trajectory length.  Care should be taken when storing event histories, as they can become quite large when transition rates are high. 

To obtain equilibrium graph distributions from an EGP, it is generally (much) more efficient to use the \code{simulate} functions in the \code{\link[ergm]{ergm}} package: they employ MCMC algorithms that are unconstrained by the need to follow realistic trajectories, and that are optimized for rapid mixing.  (In particular, note that many systems can become \emph{kinetically trapped}, spending very long periods in metastable states that are far from equilibrium.  This can be a real-world phenomenon, but is not always desirable from a computational point of view.  Functions such as \code{simEGP} are intended to faithfully reproduce such dynamics, while MCMC algorithms are intended to avoid them.)  Comparison of late-phase draws from a \code{simERGMPot} trajectory with equilibrium ERGM draws can be used to evaluate convergence to equilibrium behavior (where desired); alternately, \code{simEGP} can be seeded with ERGM draws to follow trajectories from equilibrated states.  Consult the \code{ergm} package documentation for details.
}

\value{
  For \code{simEGP}, a network object containing the final graph state, with network attributes \code{Time}, \code{Events}, and \code{Potential} listing the time, event count, and ERGM potential at the end of the simulation interval.  See above for additional attributes that may be added if history retention is activated.  If \code{return.networkDynamic==TRUE}, then the return value is instead a \code{\link[networkDynamic]{networkDynamic}} object containing the event history as edge activity data; be aware that an edge will exist in this object if any corresponding edge is ever active, so the raw graph state should not be used to access the final system state.  Instead, use the \code{\link[networkDynamic]{network.extract}} method to query the network state at the desired time point.
  
  For \code{simERGMPotTraj}, a list containing the simulated trajectories.  These are either \emph{network.list} objects containing the networks at each checkpoint (with time, step, and potential attribute as described above), or else matrices of trace statistics (if \code{statsonly==TRUE}).  Note that the statistics are in any event included as an attribute to each network list, so the effect of \code{statsonly==TRUE} is simply not to retain the graph states.
}

\references{
Butts, Carter T.  (2023).  \dQuote{Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.} \emph{Journal of Mathematical Sociology}.  \doi{10.1080/0022250X.2023.2180001}
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}

\note{
Using \code{steps} to control trajectory termination will lead to biased samples (sometimes severely so); this is because transitions are not random times. If your goal is to obtain equilibrium draws (or draws en route thereto), use \code{time} to set the stopping point.  See \code{\link{EGPRateEst}} for a simple tool for calibrating simulation times.
}

\seealso{
\code{\link{ergmgp}} for information on EGPs, \code{\link[ergm]{ergm}} for information on ERGM specifications, \code{\link{EGPHazard}}, \code{\link{EGPRateEst}}, \code{\link[networkDynamic]{networkDynamic}}
}

\examples{
\donttest{
#Small example of 2-ribbon generation
n<-100
set.seed(1331)
net<-network.initialize(n,directed=FALSE)
sim<-simEGP(net~edges+kstar(2)+nsp(1:2), 
    coef=c(109-log(n),-25,-1.25,3.25), time=100, process="LERGM",
    verbose = TRUE)
plot(sim) #Return value is a single network

#Generate a trajectory showing the process at multiple stages
set.seed(1331)
sim<-simEGPTraj(net~edges+kstar(2)+nsp(1:2), 
    coef=c(109-log(n),-25,-1.25,3.25), time=100, checkpoints = 5, 
    trajectories = 2, mc.cores = 1, log.sampling = TRUE, 
    process = "LERGM", verbose = TRUE) 
length(sim)==2        #One entry per simulated trajectory
op<-par(mfrow=c(2,3))
for(i in 1:6)         #Show the first trajectory
  plot(sim[[1]][[i]],main=paste("Time",round(sim[[1]][[i]]\%n\%"Time",2)))
summary(sim[[2]]~edges+kstar(2))  #Show selected stats from the second
attributes(sim[[1]])  #Show precomputed statistics
par(op)

#A simple example with statsonly
set.seed(1331)
sim<-simEGPTraj(net~edges+esp(0), coef = c(log(2)-log(n), -1), time = 200,
    checkpoints = 25, process = "LERGM", statsonly = TRUE, 
    monitor = ~triangle)
sim                   #Note the monitor stat
op<-par(mfrow=c(1,1))
plot(sim[,"Time"], sim[,"edges"], type = "l")  #Time by edge count
lines(sim[,"Time"], sim[,"esp0"], col = 2)     #Add ESP(0)s
par(op)
}
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ graphs }% __ONLY ONE__ keyword per line

