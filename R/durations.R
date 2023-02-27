######################################################################
#
# durations.R
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/26/23
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/ergmgp package
#
# This file contains routines related to edge duration estimation from 
# EGPs.
#
#######################################################################

#Given an input trajectory (in networkDynamic form, or network form with
#additional attributes), return the set of all edge durations (correcting
#for censoring, if desired)
#
#Arguments:
#  net - a networkDynamic object containing the trajectory, or a network
#        object with a LastChangeTime attribute
#  censor - censoring correction procedure to use
#  return.censoring - return a censoring indicator?
#Return value:
#  A vector of durations (one per spell), or a matrix with durations and
#  a censoring indicator (0=uncensored, 1=right-censored, 2=left-censored,
#  3=interval censored).
#
#Currently, the two ways of handling censoring are to return the observed
#lengths ("obs") or to omit the censored spells ("omit").  The "omit"
#method is not useful if LastChangeTime data is being used (all spells
#are censored), and in the latter case it is assumed that the observation
#period goes from 0 to the "Time" network attribute.  If desired, a 
#censoring indicator is also returned, facilitating other analyses.
#
#(Early versions of this function supported a primitive censoring
#correction method, but it did not work at all well and was often worse
#than doing nothing at all.  Such a method may be added later, if one
#can be found that gives reasonable and robust performance without
#making too many assumptions.)
#
durations<-function(net, censor=c("obs","omit"), return.censoring=TRUE){
  #Function to process spell matrices, returning durations
  procSpells<-function(z,obs.int,censor){
    #First, figure out which spells are censored
    iscen<-(z[,2]>obs.int[2])
    iscen[(z[,1]<obs.int[1])]<-2+iscen[(z[,1]<obs.int[1])]
    #Cut to fit
    z[,1]<-pmax(obs.int[1],z[,1])
    z[,2]<-pmin(obs.int[2],z[,2])
    #Get durations
    dur<-z[,2]-z[,1]
    #If needed, deal with censoring
    if(censor=="omit"){   #Drop censored spells
      dur<-dur[!iscen]
      iscen<-iscen[!iscen]  #I know, I know - done for simpler processing later
    }
    #Return the results
    cbind(dur,iscen)
  }
  #Check the input
  if(inherits(net,"networkDynamic")){ #Complete spell info, as edge activity spells
    #Obtain all of the durations
    censor<-match.arg(censor)                            #Censoring type
    obs.int<-(net%n%"net.obs.period")$observations[[1]]  #Simulation interval
    dur<-sapply(get.edge.activity(net),function(z){
      procSpells(z,obs.int=obs.int,censor=censor)
    })
    dur<-cbind(unlist(sapply(dur,"[",,1)),unlist(sapply(dur,"[",,2)))
    if(return.censoring){   #Return both durations and censoring info?
      colnames(dur)<-c("Duration","Censoring")
      dur
    }else                   #Else, just the durations
      dur[,1]
  }else if("LastChangeTime"%in%list.network.attributes(net)){ #Last change time for obs edges
    g<-as.sociomatrix(net)
    ct<-net%n%"LastChangeTime"  #Change time
    if(is.directed(net)){
      ct<-ct[g>0]
    }else{
      ct<-ct[upper.tri(ct)][g[upper.tri(g)]>0]
    }
    #All cases are censored...
    iscen<-rep(1,length(ct))                             #Right by default
    iscen[ct<=0]<-3                                      #Could also be interval
    censor<-match.arg(censor)                            #How to deal with censoring
    endtime<-net%n%"Time"
    ct<-pmax(0,ct)
    if(censor=="omit")
      return(numeric(0))
    dur<-endtime-ct
    if(return.censoring){                               #Including censoring info?
      dur<-cbind(dur,iscen)
      colnames(dur)<-c("Duration","Censoring")
    }
    dur
  }
}

