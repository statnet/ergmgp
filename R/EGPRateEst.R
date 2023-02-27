######################################################################
#
# EGPRateEst.R
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/16/23
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/ergmgp package
#
# This file contains routines related to rate estimation from EGPs.
#
#######################################################################

#Given a specified EGP, estimate the expected number of events within a given
#period or the expected time to reach a given event count.
#
#Arguments:
#  formula - an EGP formula, per simEGPTraj
#  coef - EGP coefficients, per simEGPTraj
#  process - the EGP to use, per simEGPTraj
#  time.target - a simulation period to use as a target
#  event.target - an event count to use as a target
#  reps - replicate trajectories to simulate
#  cores - cores to use during simulation
#  rate.factor - rate factor for the EGP
#  verbose - verbosity argument for trajectory simulation
#  ... - additional arguments to simEGPTraj
#
#Return value:
#  A vector containing the mean outcome (time or event count), its standard error,
#  the standard deviation of the outcome, and the number of replicates used.
#
#Note: this is really just a wrapper for simEGPTraj.  We either use it to simulate
#trajectories for a fixed period (time.target) or for a fixed number of events
#(event.target), and then examine the respective distribution of event counts or
#simulation times.  This is extremely useful for calibrating a model that one
#wants to observe in equilibrium, since choosing a large event target will often
#give reasonable mixing, but terminating simulations using event counts yields
#biased states.  One can e.g. use EGPRateEst to determine how long a trajectory
#must be to, on average, generate a particular number of events, and then use 
#this as the simulation time for subsequent simulation studies.
#
EGPRateEst<-function(formula, coef, process=c("LERGM", "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), time.target=NULL, event.target=NULL, reps=25, cores=1, rate.factor=1, verbose=FALSE, ...){
  #Ensure that we were given a target
  if(is.null(time.target)+is.null(event.target)!=1)
    stop("Exactly one of time.target and event.target must be specified.\n")
  #Simulate trajectories with the chosen target, and return
  if(is.null(time.target)){  #Estimate expected time as a function of event count
    sim<-simEGPTraj(formula, coef=coef, process=process, events=event.target, trajectories = reps, mc.cores = cores, rate.factor=rate.factor, verbose=verbose, ...)
    time<-sapply(sim,function(z){z[[2]]%n%"Time"})  #Extract times
    #Return the summary statistics
    c(mean.time=mean(time), mean.time.se=sd(time)/sqrt(reps), sd.time=sd(time), replications=reps)
  }else{                     #Estimate expected event count as a function of time
    sim<-simEGPTraj(formula, coef=coef, process=process, time=time.target, trajectories = reps, mc.cores = cores, rate.factor=rate.factor, verbose=verbose, ...)
    events<-sapply(sim,function(z){z[[2]]%n%"Events"})  #Extract event counts
    #Return the summary statistics
    c(mean.events=mean(events), mean.event.se=sd(events)/sqrt(reps), sd.events=sd(events), replications=reps)
  }
}

