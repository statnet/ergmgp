/*
######################################################################
#
# EGPHazards.c
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/10/23
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# This file contains functions for calculating event hazards from 
# EGPs.
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

/*
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
    slcrate - log "collision rate" or intrinsic rate coefficient; this is the maximum 
              transition rate as the log potential difference goes to infinity for some
              processes, and for others simply sets the timescale of dynamics
    spot - "negative effective energy" (aka ERGM potential) at onset; this is needed for
           some processes
    stoggt - a vector of vertex IDs indicating the tails of edge variables to be toggled
    stoggh - a vector of vertex IDs indicating the heads of edge variables to be toggled

  Return value: a vector of log event hazards for the selected toggles.

*/
SEXP EGPHazard_R(SEXP segp, SEXP mstate, SEXP scoef, SEXP scooffset, SEXP slcrate, SEXP spot, SEXP stoggt, SEXP stoggh){
  int pc=0,i,j,egp,cooffset,isedge,*toggt,*toggh,ntog;
  Vertex t,h;
  double *coef,lcrate,pot[2],relpot[2],lHneigh;
  double *lrates;
  SEXP stats,slrates;
  
  /*Ensure that inputs are processed correctly*/
  PROTECT(segp=coerceVector(segp,INTSXP));            /*EGP type*/
  egp=INTEGER(segp)[0];
  UNPROTECT(1);
  PROTECT(scoef=coerceVector(scoef,REALSXP)); pc++;   /*Model coefficients*/
  coef=REAL(scoef);
  PROTECT(slcrate=coerceVector(slcrate,REALSXP));     /*Log collision rate*/
  lcrate=REAL(slcrate)[0];
  UNPROTECT(1);
  PROTECT(scooffset=coerceVector(scooffset,INTSXP));  /*CSTERGM form coef offset*/
  cooffset=INTEGER(scooffset)[0];
  UNPROTECT(1);
  PROTECT(stoggt=coerceVector(stoggt,INTSXP)); pc++;  /*Tails to be toggled*/
  toggt=INTEGER(stoggt);
  PROTECT(stoggh=coerceVector(stoggh,INTSXP)); pc++;  /*Heads to be toggled*/
  toggh=INTEGER(stoggh);
  PROTECT(spot=coerceVector(spot,REALSXP));          /*ERGM potential*/
  pot[0]=REAL(spot)[0];                                /*Either total or formation pot*/
  pot[1]=REAL(spot)[1];                                /*Either unused or dissolution pot*/
  UNPROTECT(1);
  ntog=LENGTH(stoggh);
  
  /*Initialize other key things*/
  GetRNGstate();
  ErgmState *s=ErgmStateInit(mstate, ERGM_STATE_NO_INIT_PROP);  /*Set up the ERGM state*/
  Model *m = s->m;               /*Extract the model from the state*/
  Network *nwp=s->nwp;           /*Extract the network from the state*/
  int directed_flag = nwp->directed_flag;  /*Is this directed?*/
  Vertex n_nodes = nwp->nnodes;
  if((egp==EGP_CDCSTERGM)||(egp==EGP_CFCSTERGM)) /*Set coef offset for constant rate CSTERGMs*/
    cooffset=1;
  lHneigh=log(n_nodes*(n_nodes-1.0)/(2.0-directed_flag));      /*Log Hamming neighbors*/
  PROTECT(slrates=allocVector(REALSXP,ntog)); pc++;            /*Log event hazards*/
  lrates=REAL(slrates);

  /*Memory to store one row of changescores*/
  SEXP snewrow,scurcs;
  PROTECT(snewrow=allocVector(REALSXP,m->n_stats)); pc++;
  PROTECT(scurcs=allocVector(REALSXP,m->n_stats)); pc++;
  double *newRow = REAL(snewrow);
  double *curcs = REAL(scurcs);
  memset(newRow,0.0,m->n_stats*sizeof(double));
  memset(curcs,0.0,m->n_stats*sizeof(double));
  
  /*Memory for ergm stat state (needs to be protected)*/
  PROTECT(stats=allocVector(REALSXP,m->n_stats)); pc++;
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));
  
  /*Run through the list of toggles*/
  for(i=0;i<ntog;i++){
    R_CheckUserInterrupt();
    /*Flip selected edges, and compute the potential difference*/
    t=(Vertex)(toggt[i]);       /*Reminder: tail/head values are 1-indexed*/
    h=(Vertex)(toggh[i]);
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
        lrates[i] = lcrate - logspace_add(0.0,-relpot[0]);
        break;
      /*Competing rate SAOM*/
      /*log rate = log(A) + targpot*/
      case EGP_CRSAOM:
        lrates[i] = lcrate + relpot[0] + pot[0];
        break;
      /*Change inhibition process*/
      /*log rate = log(A) + min(0, potdiff)*/
      case EGP_CI:
        lrates[i] = lcrate + MIN(0.0,relpot[0]);
        break;
      /*Differential stability process*/
      /*log rate = log(A) - log |H| - curpot*/
      case EGP_DS:
        lrates[i] = lcrate - lHneigh - pot[0];
        break;
      /*Constant dissolution rate continuum STERGM*/
      /*log rate = log(A) + (1-G_th) potdiff + G_th coef_d*/
      case EGP_CDCSTERGM:
        lrates[i] = lcrate + (isedge ? relpot[1] : relpot[0]);
        break;
      /*Constant formation rate continuum STERGM*/
      /*log rate = log(A) + (1-G_th) coef_f + G_th potdiff*/
      case EGP_CFCSTERGM:
        lrates[i] = lcrate + (isedge ? relpot[1] : relpot[0]);
        break;
      /*Continuum STERGM*/
      /*log rate = log(A) + G_th potdiff_d + (1-G_th) potdiff_f*/
      case EGP_CSTERGM:
        lrates[i] = lcrate + (isedge ? relpot[1]: relpot[0]);
        break;
      /*Continuum TERGM*/
      /*log rate = log(A) + potdiff*/
      case EGP_CTERGM:
        lrates[i] = lcrate + relpot[0];
        break;
    }
  }
    
  /*Clean up and return*/
  ErgmStateDestroy(s);
  PutRNGstate();
  UNPROTECT(pc);
  return slrates;
}


