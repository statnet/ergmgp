/*
######################################################################
#
# egp.h
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 1/31/23
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# This file contains general macros and definitions for ERGM generating
# process code.
#
######################################################################
*/
#ifndef EGP_H
#define EGP_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#define EGP_LERGM       0  /*Longitudinal ERGM*/
#define EGP_CRSAOM      1  /*Competing rate SAOM*/
#define EGP_CI          2  /*Change inhibition process*/
#define EGP_DS          3  /*Differential stability process*/
#define EGP_CDCSTERGM   4  /*Constant dissolution continuum STERGM*/
#define EGP_CFCSTERGM   5  /*Constant formation continuum STERGM*/
#define EGP_CSTERGM     6  /*Continuum STERGM*/
#define EGP_CTERGM      7  /*Continuum TERGM*/

#endif

