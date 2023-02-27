######################################################################
#
# EGPHazards.R
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/16/23
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/ergmgp package
#
# This file contains various routines related to calculation of event
# hazards for EGPs.
#
#######################################################################


#Simulate a trajectory from a continuous time ERGM generating process (EGP)
#
#Arguments:
#  form - an ergm formula, or a vector of two formulas (for CSTERGM)
#  coef - ergm coefficients vector; for the CSTERGMs, this should be a list of
#    two vectors for formation and dissolution
#  toggles - edge variables to evaluate; "all" or NULL leads to all edge variables
#    being evaluated, "edges" evaluates only dissolution events, "nulls" evaluates
#    only formation events, and passing a two-column matrix of IDs evaluates the 
#    selected dyads
#  rate.factor - rate or pacing factor (sets the time scale)
#  process - the ERGM generating process to use
#
#
#Return value: a matrix containing the toggles, indicators for whether each event
#  would have been a formation event, and the log event hazards (one row per toggle)
#
EGPHazard<-function(form, coef, toggles=NULL, rate.factor=1, process=c("LERGM", "CRSAOM", "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM")){
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
  lratefact<-log(rate.factor)
  n<-network.size(nw)
  #Determine the key simulation parameters
  if(is.null(toggles)||(is.character(toggles)&&(toggles=="all"))){ #By default, get all toggles
    el<-cbind(as.vector(row(diag(n))),as.vector(col(diag(n))))
    if(is.directed(nw))
      el<-el[el[,1]!=el[,2],,drop=FALSE]
    else
      el<-el[el[,1]<el[,2],,drop=FALSE]
  }else if(is.character(toggles)&&(toggles=="edges")){            #Toggle only edges
    el<-as.edgelist(nw)
  }else if(is.character(toggles)&&(toggles=="nulls")){            #Toggle only nulls
    g<-1-as.sociomatrix(nw)
    diag(g)<-0
    el<-cbind(row(diag(n))[g==1],col(diag(n))[g==1])
    if(is.directed(nw))
      el<-el[el[,1]!=el[,2],,drop=FALSE]
    else
      el<-el[el[,1]<el[,2],,drop=FALSE]
  }else{                                 #Toggle what we were asked to toggle (two-col matrix)
    if(length(dim(toggles))==0)
      toggles<-as.matrix(toggles,ncol=2)
    toggles<-as.matrix(toggles)
    if((!is.matrix(toggles))||(NCOL(toggles)!=2))
      stop("toggles must be given either as a directive, or a two-column matrix of vertex IDs.\n")
    if(any(is.na(toggles))||any(toggles<1)||any(toggles>n))
      stop("Illegal toggle specified; toggles must be valid vertex IDs.\n")
    if(!is.directed(nw)){
      toggles<-cbind(pmin(toggles[,1],toggles[,2]),pmax(toggles[,1],toggles[,2]))
    }
    el<-unique(toggles)
  }
  #Obtain the log hazards
  z<-.Call("EGPHazard_R", procnum, state, coef, cooffset, lratefact, pot, el[,1], el[,2], PACKAGE="ergmgp")
  #Assemble the results
  out<-cbind(el,1-is.adjacent(nw,el[,1],el[,2]),z)
  colnames(out)<-c("Snd","Rec","IsForm","logHaz")
  out
} 



