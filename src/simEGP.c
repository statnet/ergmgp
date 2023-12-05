/*
######################################################################
#
# simEGP.c
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/4/23
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# This file contains functions for simulating realizations from ERGM
# generating processes.
#
######################################################################
*/

#include "egp.h"
#include "ergm_model.h"
#include "ergm_state.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_util.h"
#include <R_ext/Utils.h>
#include "utils.h"

/*
  Simulate a single trajectory from an ERGM generating process (EGP); several
  different types of EGPs are supported.  This function simulates a trajectory 
  of specified length and returns the final state (along with, optionally, history
  information); length may be specified by the number of events drawn (time will 
  be random) or by the time taken (event count will be random).  It is important
  to note that event updates are *not* random times, and hence trajectories
  terminated by event count will not be in equilibrium.
  
  Arguments:
    segp - an ERGM generating process type code, saying which process should be used
           for the simulation
    mstate - an input ergm state, from ergm_state(), containing the model, graph state,
             and anything else needed
    scoef - coefficient vector for the model; for constant rate CSTERGMs, the first
            element should be the constant rate (formation or dissolution).  For
            general CSTERGMs, this should be a vector of formation coefficients
            followed by the dissolution coefficients
    scooffset - coefficient offsets for CSTERGMs.  For constant rate CSTERGMs, this
                should be 1, and the first element should be whichever rate is
                constant; for general CSTERGMs, this should be the number of formation
                parameters
    stmax - maximum simulation time; if constraining by events, this should be negative
            (meaning unlimited or infinite)
    sevmax - maximum number of events to simulate; if constraining by time, this should be
               negative (meaning unlimited or infinite)
    slcrate - log "collision rate" or intrinsic rate coefficient; this is the maximum 
              transition rate as the log potential difference goes to infinity for some
              processes, and for others simply sets the timescale of dynamics
    spot - "negative effective energy" (aka ERGM potential) at onset; this allows
              the potential to be tracked
    slogtime - logical; track time on log scale?  This is purely internal; using it
               can potentially avoid overflow/underflow errors for models with extreme
               rates, but also increases computational burden somewhat
    slasttog - logical; or every edge variable, track the time of its last toggle
               (i.e., from off to on, or on to off)?  If TRUE, an n x n matrix is 
               included in the return value conaining the last toggle time
    slttoff - n x n matrix of initial last toggle times, in case we are resuming a
               simulation already in progress
    skeephist - logical; should we keep and return the entire event history?  (Not for
                the faint of heart or those short of RAM.)
    sverbose - if positive, we will return progress messages; events are shown every
               verbose events (so e.g., verbose=100 is a good idea, and verbose=1 is not)
  
  Return value:
    A named list with elements 
      state: the final state, in ergm_state() form
      evcount: the final event or step count (number of observed transitions)
      simtime: the final simulation time (which may not be the time of the last event,
               depending on whether events or total time were left random)
      potential: final ERGM potential
      lasttog: if used, an n x n matrix containing the last toggle time for each edge
               variable
      evhist: if used, the complete event history (as an event by 4 matrix, with cols
               sender, receiver, time, and formation indicator (1=onset, 0=terminus)
*/
SEXP simEGP_R(SEXP segp, SEXP mstate, SEXP scoef, SEXP scooffset, SEXP stmax, SEXP sevmax, SEXP slcrate, SEXP spot, SEXP slogtime, SEXP slasttog, SEXP slttoff, SEXP skeephist, SEXP sverbose){
  int pc=0,verbose,i,j,egp,cooffset,isedge,logtime,lasttog,keephist,selisedge=0;
  Vertex t,h,selth[2]={0,0};
  double *coef,tmax,evmax,lcrate,ecount,simtime,pot[2],lHneigh,*ltt,*dp;
  double tm,relpot[2],ltotrat=R_NegInf,selval=R_NegInf,selpot[2],selmax=R_NegInf;
  double *pehist,*pdp;
  element *ehist;
  SEXP outl,stats,sltt,sehist;
  R_xlen_t nev,xi;
  selpot[0]=selpot[1]=R_NegInf;
  ehist=NULL;
  sltt=NULL;
  ltt=NULL;
  
  /*Ensure that inputs are processed correctly*/
  PROTECT(segp=coerceVector(segp,INTSXP));           /*EGP type*/
  egp=INTEGER(segp)[0];
  UNPROTECT(1);
  PROTECT(scoef=coerceVector(scoef,REALSXP)); pc++;  /*Model coefficients*/
  coef=REAL(scoef);
  PROTECT(stmax=coerceVector(stmax,REALSXP));        /*Max simulation time*/
  tmax=REAL(stmax)[0];
  UNPROTECT(1);
  PROTECT(sevmax=coerceVector(sevmax,REALSXP));      /*Max events*/
  evmax=REAL(sevmax)[0];
  UNPROTECT(1);
  PROTECT(slcrate=coerceVector(slcrate,REALSXP));    /*Log collision rate*/
  lcrate=REAL(slcrate)[0];
  UNPROTECT(1);
  PROTECT(spot=coerceVector(spot,REALSXP));          /*ERGM potential*/
  pot[0]=REAL(spot)[0];                                /*Either total or formation pot*/
  pot[1]=REAL(spot)[1];                                /*Either unused or dissolution pot*/
  UNPROTECT(1);
  PROTECT(sverbose=coerceVector(sverbose,INTSXP));   /*Verbosity*/
  verbose=INTEGER(sverbose)[0];
  UNPROTECT(1);
  PROTECT(slogtime=coerceVector(slogtime,INTSXP));   /*Use log time*/
  logtime=INTEGER(slogtime)[0];
  UNPROTECT(1);
  PROTECT(slasttog=coerceVector(slasttog,INTSXP));   /*Track last toggles*/
  lasttog=INTEGER(slasttog)[0];
  UNPROTECT(1);
  PROTECT(skeephist=coerceVector(skeephist,INTSXP)); /*Keep event history*/
  keephist=INTEGER(skeephist)[0];
  UNPROTECT(1);
  PROTECT(scooffset=coerceVector(scooffset,INTSXP));  /*CSTERGM form coef offset*/
  cooffset=INTEGER(scooffset)[0];
  UNPROTECT(1);
  
  /*Initialize other key things*/
  GetRNGstate();
  ErgmState *s=ErgmStateInit(mstate, ERGM_STATE_NO_INIT_PROP);  /*Set up the ERGM state*/
  Model *m = s->m;              /*Extract the model from the state*/
  Network *nwp=s->nwp;           /*Extract the network from the state*/
  int directed_flag = nwp->directed_flag;  /*Is this directed?*/
  Vertex n_nodes = nwp->nnodes;
  if(tmax<0.0)                 /*Deal w/R's prohibition on passing Inf values*/
    tmax=R_PosInf;
  if(evmax<0.0)
    evmax=R_PosInf;
  if((egp==EGP_CDCSTERGM)||(egp==EGP_CFCSTERGM)) /*Set coef offset for constant rate CSTERGMs*/
    cooffset=1;
  lHneigh=log(n_nodes*(n_nodes-1.0)/(2.0-directed_flag));  /*Log Hamming neighbors*/
  if(lasttog){     /*If storing last toggle times, initialize the n x n timing matrix*/
    if(slttoff==R_NilValue){  /*Start from scratch*/
      PROTECT(sltt=allocMatrix(REALSXP,n_nodes,n_nodes)); pc++;
      ltt=REAL(sltt);
      for(i=0;i<n_nodes;i++)
        for(j=0;j<n_nodes;j++)
          ltt[i+j*n_nodes]=0.0;
    }else{                      /*Use the one we were given*/
      PROTECT(sltt=coerceVector(slttoff,REALSXP)); pc++;
      if(LENGTH(sltt)!=n_nodes*n_nodes){
        error("We were given a last toggle time matrix of total length %ld, but it should have been %ld.  Stopping.\n", (long)LENGTH(sltt), (long)(n_nodes*n_nodes));
      }
      ltt=REAL(sltt);
    }
  }
  //Rprintf("FWIW, model claims to have directed flag %d, and to have %d stats\n",directed_flag, m->n_stats);

  /*Memory to store one row of changescores*/
  SEXP snewrow,scurcs;
  PROTECT(snewrow=allocVector(REALSXP,m->n_stats)); pc++;
  PROTECT(scurcs=allocVector(REALSXP,m->n_stats)); pc++;
  double *newRow = REAL(snewrow);
  double *curcs = REAL(scurcs);
  memset(newRow,0.0,m->n_stats*sizeof(double));
  memset(curcs,0.0,m->n_stats*sizeof(double));
  
  /*Memory to store a matrix of log rates*/
  SEXP slrates;
  PROTECT(slrates=allocVector(REALSXP,n_nodes*n_nodes)); pc++;
  double *lrates = REAL(slrates);
  memset(lrates,0.0,n_nodes*n_nodes*sizeof(double));

  /*Memory for ergm stat state (needs to be protected)*/
  //Rprintf("\tGoing to allocVector for %d stats\n",m->n_stats);
  PROTECT(stats=allocVector(REALSXP,m->n_stats)); pc++;
  //Rprintf("\tGoing to memcpy allocated memory (%d bytes) from s->stats to stats\n",m->n_stats*sizeof(double));
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));
  //Rprintf("\tOK, that happened without obvious incident\n");
  
  /*Run the event simulation until time expires or we run out of events*/
  if(logtime){
    tm=R_NegInf;   /*Store time on log scale, to deal with extreme rate cases*/
    if(tmax<R_PosInf)
      tmax=log(tmax);
    simtime=R_NegInf;
  }else{
    tm=0.0;        /*...or just use standard scale, which is faster*/         
    simtime=0.0;
  }
  ecount=0.0;
  if(verbose){
    Rprintf("Initializing simulation: max events=%.0f, max time=%0f, initial pot=(%f,%f)\n", evmax, tmax, pot[0],pot[1]);
  }
  while((ecount<evmax)&&(tm<tmax)){
    if((verbose)&&(((int)(ecount))%verbose==0)){
      Rprintf("event=%.0f, t=%f, pot=(%f,%f)\n",ecount,tm,pot[0],pot[1]);
    }
    R_CheckUserInterrupt();
    /*Flip every edge, and compute the potential difference*/
    ltotrat=R_NegInf;
    selmax=R_NegInf;
    for(t=1; t<=n_nodes; t++){  /*Note: tail/head values are 1-indexed*/
      for(h=(directed_flag ? 1 : t+1); h<=n_nodes; h++)
        if(t!=h){
          isedge=IS_OUTEDGE(t,h);  /*Is the (t,h) edge present?*/
          /*Compute the potential difference for the (t,h) toggle*/
          relpot[0]=relpot[1]=0.0;
          memset(newRow,0.0,m->n_stats*sizeof(double));
          EXEC_THROUGH_TERMS_INTO(m, newRow, {
            if(mtp->c_func){
              ZERO_ALL_CHANGESTATS();
              (*(mtp->c_func))(t, h, mtp, nwp, IS_OUTEDGE(t,h));
            }else if(mtp->d_func){
              (*(mtp->d_func))(1, &t, &h, mtp, nwp);
            }
            addonto(dstats, mtp->dstats, N_CHANGE_STATS);
          });
          //Rprintf("\tdstats: ");
          //for(j=0;j<m->n_stats;j++)
          //  Rprintf("%f ",curcs[j]);
          //Rprintf("\n");
          /*Compute relative potential (or potentials, for CSTERGMs)*/
          if(egp==EGP_CDCSTERGM){        /*Constant dissolution CSTERGMs*/
            for(j=0;j<m->n_stats;j++)
              relpot[0]+=coef[j+1]*newRow[j];
            relpot[1]=coef[0];
          }else if(egp==EGP_CFCSTERGM){  /*Constant formation CSTERGMs*/
            for(j=0;j<m->n_stats;j++)
              relpot[1]+=coef[j+1]*newRow[j];
            relpot[0]=coef[0];
          }else if(egp==EGP_CSTERGM){     /*For CSTERGMs, consider form and diss*/
            for(j=0;j<cooffset;j++)
              relpot[0]+=coef[j]*newRow[j];
            for(j=cooffset;j<m->n_stats;j++)
              relpot[1]+=coef[j]*newRow[j];
          }else{                          /*Standard calculation for most models*/
            for(j=0;j<m->n_stats;j++)    
              relpot[0]+=coef[j]*newRow[j];
          }
          /*Store the rate information*/
          switch(egp){
            /*Longitudinal ERGM*/
            /*log rate = log(A) - log1p(-potdiff)*/
            case EGP_LERGM:
              lrates[t-1+(h-1)*n_nodes] = lcrate - logspace_add(0.0,-relpot[0]);
              break;
            /*Competing rate SAOM*/
            /*log rate = log(A) + targpot*/
            case EGP_CRSAOM:
              lrates[t-1+(h-1)*n_nodes] = lcrate + relpot[0] + pot[0];
              break;
            /*Change inhibition process*/
            /*log rate = log(A) + min(0, potdiff)*/
            case EGP_CI:
              lrates[t-1+(h-1)*n_nodes] = lcrate + MIN(0.0,relpot[0]);
              break;
            /*Differential stability process*/
            /*log rate = log(A) - log |H| - curpot*/
            case EGP_DS:
              lrates[t-1+(h-1)*n_nodes] = lcrate - lHneigh - pot[0];
              break;
            /*Constant dissolution rate continuum STERGM*/
            /*log rate = log(A) + (1-G_th) potdiff + G_th coef_d*/
            case EGP_CDCSTERGM:
              lrates[t-1+(h-1)*n_nodes] = lcrate + (isedge ? relpot[1] : relpot[0]);
              break;
            /*Constant formation rate continuum STERGM*/
            /*log rate = log(A) + (1-G_th) coef_f + G_th potdiff*/
            case EGP_CFCSTERGM:
              lrates[t-1+(h-1)*n_nodes] = lcrate + (isedge ? relpot[1] : relpot[0]);
              break;
            /*Continuum STERGM*/
            /*log rate = log(A) + G_th potdiff_d + (1-G_th) potdiff_f*/
            case EGP_CSTERGM:
              lrates[t-1+(h-1)*n_nodes] = lcrate + (isedge ? relpot[1]: relpot[0]);
              break;
            /*Continuum TERGM*/
            /*log rate = log(A) + potdiff*/
            case EGP_CTERGM:
              lrates[t-1+(h-1)*n_nodes] = lcrate + relpot[0];
              break;
          }
          ltotrat=logspace_add(ltotrat,lrates[t-1+(h-1)*n_nodes]);
          //Rprintf("(%ld,%ld) relpot[0]=%f, lrate=%f, ltotrat=%f\n",t,h,relpot[0],lrates[t-1+(h-1)*n_nodes],ltotrat);
          /*Determine the selection value for the Gumbel trick: we select events by drawing a Gumbel(0) deviate
            for each transition, and then finding arg max_ij log(r_ij) + G(0)_ij.  Because we can do this as we
            go, we can do it in one step.  We still have to pay two log operations to generate the Gumbel draw,
            alas, but things could be worse.*/
          selval=lrates[t-1+(h-1)*n_nodes]-log(-log(runif(0.0,1.0)));
          if(selval>selmax){  /*This is our current winner!*/
            selmax=selval;
            selisedge=isedge;
            selth[0]=t;
            selth[1]=h;
            selpot[0]=relpot[0];
            selpot[1]=relpot[1];
          }
        }
    }
    /*Draw the next event time*/
    if(logtime)
      tm=logspace_add(tm, log(-log(runif(0.0,1.0)))-ltotrat); /*Log version*/
    else
      tm+=rexp(exp(-ltotrat)); /*(Non log version) Note that they use scale rather than rate here...*/
    //Rprintf("\tNext event happens at %f\n",tm);
    /*If the next event would occur before the end of the period, make the transition*/
    if(tm<tmax){
      ecount++;                   /*Increment the event count*/
      simtime=tm;                 /*Update the time of the last event*/
      pot[0] += selpot[0];        /*Formation (or general) potential*/
      pot[1] += selpot[1];        /*Dissolution potential*/
      TOGGLE(selth[0],selth[1]);  /*Toggle the winning edge!*/
      if(lasttog){                /*Store the update time, if desired*/
        ltt[selth[0]-1+(selth[1]-1)*n_nodes]=(logtime ? exp(tm) : tm);
      }
      if(keephist){
        dp=(double *)R_alloc(3,sizeof(double));
        dp[0]=(double)selth[0];
        dp[1]=(double)selth[1];
        dp[2]=1.0-(double)selisedge;
        ehist=enqueue(ehist,(logtime ? exp(tm) : tm), (void *)dp);
      }
    }
  }
  if(tm>tmax){  /*If we truncated at time tmax, update simulation time accordingly*/
    simtime=tmax;
  }
  if(logtime)  /*If we used logtime, transform back*/
    simtime=exp(simtime);
    
  /*Save our stuffs*/
  //Rprintf("Survived, trying to save output.\n");
  const char *outputnams[] = {"state","evcount","simtime","potential","lasttog","evhist",""};
  PROTECT(outl=mkNamed(VECSXP,outputnams)); pc++;
  //Rprintf("\tCreated output vector, now writing ERGM save state\n");
  s->stats=REAL(stats);
  SET_VECTOR_ELT(outl,0,ErgmStateRSave(s));
  //Rprintf("\tSaving other sundries\n");
  SEXP secount,ssimtime;
  PROTECT(secount=allocVector(REALSXP,1)); pc++;
  PROTECT(ssimtime=allocVector(REALSXP,1)); pc++;
  REAL(secount)[0]=ecount;
  REAL(ssimtime)[0]=simtime;
  REAL(spot)[0]=pot[0];
  REAL(spot)[1]=pot[1];
  SET_VECTOR_ELT(outl,1,secount);
  SET_VECTOR_ELT(outl,2,ssimtime);
  SET_VECTOR_ELT(outl,3,spot);
  if(lasttog){                         /*Save the toggle matrix*/
    SET_VECTOR_ELT(outl,4,sltt);
  }
  if(keephist){                        /*Save the event history matrix*/
    nev=(R_xlen_t)ecount;
    PROTECT(sehist=allocMatrix(REALSXP,nev,4)); pc++;
    pehist=REAL(sehist);
    for(xi=nev-1;xi>=0;xi--){
      pehist[xi]=ehist->val;                          /*Event time*/
      pdp=(double *)(ehist->dp);
      pehist[xi+nev]=pdp[0];                          /*Tail*/
      pehist[xi+2*nev]=pdp[1];                        /*Head*/
      pehist[xi+3*nev]=pdp[2];                        /*Onset/terminus indicator*/
      ehist=ehist->next;
    }
    SET_VECTOR_ELT(outl,5,sehist);
  }
  /*Clean up and return*/
  //Rprintf("Starting cleanup\n");
  ErgmStateDestroy(s);
  PutRNGstate();
  UNPROTECT(pc);
  return outl;
}


