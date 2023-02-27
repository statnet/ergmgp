######################################################################
#
# simEGP.R
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/9/23
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/ergmgp package
#
# This file contains various routines related to simulation of EGP
# trajectories.
#
#######################################################################


#Internal function for initializing EGP simulation calls
#
#Arguments:
#  form - ERGM formula, or a list with elements "formation" and "dissolution"
#    containing such formulas (for CSTERGM)
#  coef - a coefficient vector, or a list with elements "formation" and "dissolution"
#    containing such formulas (for the CSTERGM variants)
#  process - the type of EGP to simulate
#
#Return value:
#  A list with elements needed by the EGP simulation functions.
#
EGP_init<-function(form, coef=NULL, process=c("LERGM", "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM")){
  proc<-match.arg(process)
  procnum<-switch(process,  #Process number, for C code
    LERGM=0,
    CRSAOM=1,
    CI=2,
    DS=3,
    CDCSTERGM=4,
    CFCSTERGM=5,
    CSTERGM=6,
    CTERGM=7
  )
  if(proc=="CSTERGM"){           #Have to deal with formation/dissolution potentials
    nw<-ergm.getnetwork(form[[1]])                 #Get the graph
    f<-c(as.character(form$formation)[length(as.character(form$formation))],as.character(form$dissolution)[length(as.character(form$dissolution))])  #Get the RHSs
    fcomb<-paste(f[1],"+",f[2])             #Combined formula (because reasons)
    model<-ergm_model(as.formula(paste("nw~",fcomb)),nw)        #Combined model
    state<-ergm_state(nw,model=model)
    state<-update(state,stats=summary(state))
    np<-nparam(model)                              #Total parameters (combined)
  }else{                         #One potential will rule them all (mostly)
    nw<-ergm.getnetwork(form)                      #Get the graph
    model<-ergm_model(form,nw)                     #Set up the model state
    state<-ergm_state(nw,model=model)
    state<-update(state,stats=summary(state))
    np<-nparam(model)                              #Get the number of parameters
    f<-as.character(form)[3]                       #Save the formula RHS
    fcomb<-f
  }
  if(!is.null(coef)){
    if(proc%in%c("LERGM","CRSAOM","CI","DS","CTERGM")){
      if(length(coef)!=np){
        stop("Uh-oh - you only provided ",length(coef)," coefficients, but there are ",np," model parameters.  Stopping.\n")
      }
      pot<-c(as.numeric(summary(as.formula(paste("nw",f,sep="~")))%*%coef),0) #Second element not used
      cooffset<-0                               #Coef offset
    }else if(proc=="CSTERGM"){
      if(length(coef[[1]])+length(coef[[2]])!=np){
        stop("Uh-oh - you only provided ",length(coef[[1]]) + length(coef[[2]])," coefficients, but there are ",np," model parameters.  Stopping.\n")
      }
      pot<-c(as.numeric(summary(as.formula(paste("nw",f[1],sep="~")))%*%coef$formation),  as.numeric(summary(as.formula(paste("nw",f[2],sep="~")))%*%coef$dissolution))
      cooffset<-length(coef$formation)          #Formation coef offset
      coef<-c(coef$formation,coef$dissolution)  #Combine
    }else if(proc=="CDCSTERGM"){
      if(length(coef$formation)!=np){
        stop("Uh-oh - you only provided ",length(coef$formation)," formation coefficients, but there are ",np," model parameters.  Stopping.\n")
      }
      pot<-c(as.numeric(summary(as.formula(paste("nw",f,sep="~")))%*%coef$formation),  as.numeric(summary(as.formula("nw~edges"))%*%coef$dissolution))
      cooffset<-1                               #Coef offset
      coef<-c(coef$dissolution,coef$formation)  #Constant goes first
    }else if(proc=="CFCSTERGM"){
      if(length(coef$dissolution)!=np){
        stop("Uh-oh - you only provided ",length(coef$dissolution)," dissolution coefficients, but there are ",np," model parameters.  Stopping.\n")
      }
      pot<-c(as.numeric(summary(as.formula("nw~edges"))%*%coef$formation), as.numeric(summary(as.formula(paste("nw",f,sep="~")))%*%coef$dissolution))
      cooffset<-1                               #Coef offset
      coef<-c(coef$formation,coef$dissolution)  #Constant foes first
    }
  }else{
    cooffset<-NULL
    pot<-NULL
  }
  #Return the processed inputs
  list(proc=proc,procnum=procnum,nw=nw,np=np,model=model,state=state,coef=coef,cooffset=cooffset,pot=pot, f=f,fcomb=fcomb)
}

#Simulate a trajectory from a continuous time ERGM generating process (EGP)
#
#Arguments:
#  form - an ergm formula, or a vector of two formulas (for CSTERGM)
#  coef - ergm coefficients vector; for the CSTERGMs, this should be a list of
#    two vectors for formation and dissolution
#  events - optionally, the number of simulated events to draw (if time=NULL); if
#    time is specified, this is ignored
#  time - optionally, the temporal length of the simulation; if not supplied, events
#    is used instead to determine when to stop
#  rate.factor - rate or pacing factor (sets the time scale)
#  time.offset - an offset to be added to the simulation time (only affects the
#    time listed in the output object)
#  event.offset - an offset to be added to the simulation event count (only affects
#    cumulative event count listed in the output object)
#  process - the ERGM generating process to use
#  use.logtime - logical; internally, use logarithmic timescale?  This can
#    potentially protect against overflow or underflow when rates are extreme,
#    but adds some overhead
#  return.changetime - logical; should we return a matrix with the last update
#    times for each edge variable?
#  changetime.offset - optionally, an n x n matrix of last change times (for
#    trajectories being resumed in process)
#  return.history - logical; return the entire event history?
#  return.networkDynamic - logical; retain history and return in networkDynamic
#    form?
#  verbose - logical; report on progress during the simulation?
#  trace.interval - for verbose output, interval at which event updates should be
#    given (in events)
#  ... - additional arguments (not currently used)
#
#
#Return value: a network object with the state of the system after the specified
#  number of events or time period (depending on whether the time argument was used).
#  The final time (plus offset, if any) is listed in the "Time" attribute of the 
#  network object, the cumulative number of events (plus any offset) is listed in the
#  "Events" attribute, and the ERGM potential (theta*t(Y)) is listed in the
#  "Potential" attribute.  If indicated, last change times are given in the 
#  "LastChangeTime" attribute, and event histories in the "EventHistory" attribute.
#
simEGP<-function(form, coef, events=1, time=NULL, rate.factor=1, time.offset=0, event.offset=0, process=c("LERGM", "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), use.logtime=FALSE, return.changetime=FALSE, changetime.offset=NULL, return.history=FALSE, return.networkDynamic=FALSE, verbose=TRUE, trace.interval=100, ...){
  #Set things up
  ini<-EGP_init(form=form, coef=coef, process=process) #Pass to EGP_init for setup
  proc<-ini$proc
  procnum<-ini$procnum
  nw<-ini$nw
  model<-ini$model
  state<-ini$state
  coef<-ini$coef
  pot<-ini$pot
  np<-ini$np
  cooffset<-ini$cooffset
  #Determine the key simulation parameters
  if(is.null(time)||(!is.finite(time))){
    tmax<- -1
    evmax<-events
  }else{
    tmax<-time
    evmax<- -1
  }
  lratefact<-log(rate.factor)
  #If returning as a networkDynamic object, activate history tracking
  if(return.networkDynamic)
    return.history<-TRUE
  #Run the simulation, using the specified ERGM generating process
  if(verbose)
    verbose<-max(1,trace.interval)
  z<-.Call("simEGP_R", procnum, state, coef, cooffset, tmax, evmax, lratefact, pot, as.integer(use.logtime), return.changetime, changetime.offset, return.history, as.integer(verbose), PACKAGE="ergmgp")
  #Construct the response
  if(return.networkDynamic){  #Fold the event history into a networkDynamic object
    #Create a base object (ignore warnings about not having data - we add it later)
    suppressWarnings(nw<-networkDynamic(base.net=nw, net.obs.period=list(mode="continuous", time.increment=NA, time.unit="unknown", observations=list(c(0, ifelse(tmax>0,tmax,0))))))
    activate.edges(nw,onset=-Inf,terminus=Inf,e=1:length(nw$mel))  #Ensure edges initially active
    #Add edges and events
    alladd<-z$evhist[z$evhist[,4]>0,2:3,drop=FALSE]
    if(NROW(alladd)>0){
      alladd<-unique(alladd)  #All edges formed
      sel<-!is.adjacent(nw,alladd[,1],alladd[,2])
      if(any(sel)){
        nw<-add.edges(nw,tail=alladd[sel,1],head=alladd[sel,2])
      }
    }
    if(z$evcount>0){
      for(i in 1:NROW(z$evhist)){
        eid<-get.edgeIDs(nw, v=z$evhist[i,2],alter=z$evhist[i,3])
        if(z$evhist[i,4]>0){  #Onset
          activate.edges(nw, onset=z$evhist[i,1], terminus=Inf, e=eid)
        }else{                #Terminus
          deactivate.edges(nw, onset=z$evhist[i,1], terminus=Inf, e=eid)
        }
      }
    }
  }else{                      #Use a standard network object
    nw<-z$state$nw
    if(NROW(z$state$el)>0)
      nw<-add.edges(nw,as.numeric(z$state$el[[1]]),as.numeric(z$state$el[[2]]))
    if(return.history){
      colnames(z$evhist)<-c("Time","Snd","Rec","Onset")
      nw%n%"EventHistory"<-z$evhist
    }
  }
  #Store some useful things
  if(proc%in%c("CDCSTERGM","CFCSTERGM","CSTERGM"))
    nw%n%"Potential"<-z$potential
  else
    nw%n%"Potential"<-z$potential[1]
  nw%n%"Time"<-z$simtime+time.offset
  nw%n%"Events"<-z$evcount+event.offset
  if(return.changetime)
    nw%n%"LastChangeTime"<-z$lasttog
  nw
} 


#A convenient wrapper for simEGP, that simulates multiple trajectories (in
#parallel, if desired) using temporal or eventwise checkpointing, and returns the
#collated results.
#
#Arguments:
#  form - an ergm formula, or a vector of two formulas (for CSTERGM)
#  coef - ergm coefficients vector; for the CSTERGMs, this should be a list of
#    two vectors for formation and dissolution
#  events - total number of events in each trajectory, if time is not specified
#  time - total time for each trajectory (overrides events)
#  checkpoints - number of checkpoints for which values should be returned
#  rate.factor - rate or pacing factor (sets the time scale)
#  trajectories - number of independent trajectories to simulate (all start from
#    the seed network in form, but evolve independently)
#  mc.cores - the number of cores to use when computing trajectories
#  log.sampling - logical; should time/event checkpoints be determined on log scale?
#    (This is helpful for tracking systems that slow down exponentially fast.)
#  process - the ERGM generating process to use
#  use.logtime - logical; internally, use logarithmic timescale?  This can
#    potentially protect against overflow or underflow when rates are extreme,
#    but adds some overhead
#  return.changetime - logical; should we return a matrix with the last update
#    times for each edge variable?
#  return.history - logical; return the entire event history?
#  verbose - logical; report on progress during the simulation?  (This can slow
#    things down rather a lot)
#  trace.interval - for verbose output, interval at which event updates should be
#    given (in events)
#  statsonly - logical; return only summary statistics?  (Saves on memory, but
#    will not otherwise improve performance.)
#  monitor - optionally, a RHS ergm formula with additional stats to monitor.  If
#    set to null, only the base formula stats are calculated
#
#Return value:
#  If statsonly=TRUE, a list containing matrices with the summary statistics (as well
#  as the time, event count, and energy at each checkpoint.  Otherwise, a list of
#  network.list objects is returned, with statistics information saved as a "stats"
#  attribute.  If only one trajectory is computed, a single matrix or network.list
#  is returned (instead of a list thereof).
# 
simEGPTraj<-function(form, coef, events=1, time=NULL, checkpoints=1, rate.factor=1, trajectories=1, mc.cores=1, log.sampling=FALSE, process=c("LERGM", "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), use.logtime=FALSE, return.changetime=FALSE, return.history=FALSE, verbose=TRUE, trace.interval=100, statsonly=FALSE, monitor=NULL){
  #Set things up
  ini<-EGP_init(form=form, coef=coef, process=process) #Pass to EGP_init for setup
  proc<-ini$proc
  bn<-ini$nw
  ipot<-ini$pot
  f<-ini$f
  fcomb<-ini$fcomb
  #Compute the trajectories
  traj<-mclapply(1:trajectories,function(z,bn,f,coef,mon,ipot){
    #Set things up for the simulation
    if(!is.null(mon)){                            #Set up the monitor formula
      fm<-as.character(mon)
      fm<-paste(fcomb,fm[length(fm)],sep="+")
     }else{
      fm<-fcomb
    }
    if(!statsonly){  #We will be saving the actual draws, so set things up accordingly
      sim<-vector(mode="list",length=checkpoints+1)
    }
    simstats<-summary(as.formula(paste("bn",fm,sep="~"))) #Get baseline stats
    #Get initial potential
    pot<-ipot
    if(proc%in%c("CDCSTERGM","CDCSTERGM","CSTERGM")){
      pnam<-c("Potential.Form","Potential.Diss")
    }else{
      pnam<-c("Potential")
    }
    snam<-names(simstats)                                 #Get stat names
    simstats<-c(0,0,pot[1:(1+(proc%in%c("CDCSTERGM","CDCSTERGM","CSTERGM")))],simstats)                         #Append time, etc.
    names(simstats)<-c("Time","Events",pnam,snam)         #Set names
    net<-bn                                               #Create the working copy
    net%n%"Time"<-0                                       #Initialize network props
    net%n%"Events"<-0
    net%n%"Potential"<-pot
    if(!statsonly)
      sim[[1]]<-net
    if(!is.null(time)){  #Set up the time or event sampling points
      if(log.sampling){
        timeinc<-as.list(diff(c(0,exp((1:checkpoints)/checkpoints*log1p(time))-1)))
      }else{
        timeinc<-as.list(rep(time/checkpoints,checkpoints))
      }
      evinc<-rep(1,checkpoints)
    }else{
      timeinc<-replicate(checkpoints,NULL)
      if(log.sampling){
        evinc<-pmax(round(diff(c(0,exp(seq(from=0,to=log(events),length=checkpoints))))), 1)
      }else{
        evinc<-rep(max(round(events/checkpoints),1),checkpoints)
      }
    }
    #Walk through the checkpoints, simulating as we go
    for(i in 1:checkpoints){
      #Get the next network state
      if(proc=="CSTERGM")
        nform<-list(formation=as.formula(paste("net",f[1],sep="~")), dissolution=as.formula(paste("net",f[2],sep="~")))
      else
        nform<-as.formula(paste("net",f,sep="~"))
      net<-simEGP(nform, coef=coef, events=evinc[i], time=timeinc[[i]], rate.factor=rate.factor, time.offset=net%n%"Time", event.offset=net%n%"Events", changetime.offset=net%n%"LastChangeTime", process=process, use.logtime=use.logtime, return.changetime=return.changetime, return.history=return.history, verbose=verbose, trace.interval=trace.interval)
      #If desired, save the network
      if(!statsonly)
        sim[[i+1]]<-net
      #Save the network statistics
      simstats<-rbind(simstats,c(net%n%"Time",net%n%"Events",net%n%"Potential", summary(as.formula(paste("net",fm,sep="~")))))
    }
    rownames(simstats)<-0:checkpoints
    #Return the result for this trajectory
    if(statsonly){
      simstats
    }else{
      class(sim)<-"network.list"
      attr(sim,"stats")<-simstats
      sim
    }
  },bn=bn,f=f,coef=coef,mon=monitor,ipot=ipot,mc.cores=mc.cores)
  #Return the result
  if(trajectories==1)
    traj[[1]]
  else
    traj
}

