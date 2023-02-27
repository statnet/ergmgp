/*
######################################################################
#
# utils.h
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/8/23
# Licensed under the GNU General Public License version 3 or later.
# Portions taken from the sna library by Carter T. Butts,
#  (self-licensed under GPL).
#
# Part of the R/ergmgp package
#
# This file contains headers for  various utility routines to be 
# called by other ergmgp functions.
#
######################################################################
*/

#ifndef UTILS_H
#define UTILS_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

#ifndef MAX
  #define MAX(a,b) ((a)>(b) ? (a) : (b))
#endif
#ifndef MIN
  #define MIN(a,b) ((a)<(b) ? (a) : (b))
#endif
#ifndef IISNA
  #define IISNA(a) (a == NA_INTEGER)
#endif

/*The element datatype; contains a double value, an abstract pointer, and a
pointer to the next element. The purpose of the element is to serve in generic
stacks and queues, which are needed for a number of purposes (including the
BFS algorithms used here).*/

typedef struct elementtype{
   double val;
   void *dp;
   struct elementtype *next;
} element;


/*Element datatypes for skip lists.  These are fairly similar to "elements," except that they add an extra "depth" parameter to keep track of the length of the pointer vector contained in next (as opposed to the single pointer for a standard element).  The list itself will ultimately look something like this:

 head -> el1 -> el2 -> el3 -> el4 -> el5 -> el6 -> el7
  |       |             |      |      |             |
 ph1 --> p11 --------> p31 -> p41 -> p51 --------> p71
  |       |             |             |
 ph2 --> p12 --------> p32 --------> p52
  |                     |             |
 ph3 ----------------> p33 --------> p53

Note that the initial pointer to the next element is intended to be the first element of next, and is carried by all list members.  Additional pointers may also be carried, this being (of course) random.  In practice, skip lists always require a head pointer with depth equal to the maximum list depth; it may be helpful to maintain the list length in its val entry (and this is done here).

Important implementation note: as used here, depth=length(next vector)-1 (i.e., the mandatory first element is not included).  This allows el->next[depth] to be the outermost pointer, with depth==0 denoting the minimal case.  (Otherwise, we'd have an extra subtraction operation carried through all of our lookup routines, for no particularly good reason.)  Just bear in mind that the length of the next vector is depth+1, in case this is important for some application or other.
*/

typedef struct slelementtype{
  double val;
  void *dp;
  struct slelementtype **next;
  int depth;
} slelement;


/*INTERNAL ROUTINES---------------------------------------------------------*/

/*STACK/QUEUE/LIST ROUTINES*/

int isInSList(slelement *head, double val);

slelement *slistDelete(slelement *head, double val);

slelement *slistInsert(slelement *head, double val, void *dp);

void slistPrint(slelement *head);

slelement *slistSearch(slelement *head, double val);

int isInList(element *head, double val);

element *listInsert(element *head, double val, void *dp);

element pop(element *head);

element *push(element *head, double val, void *dp);

element *pushCalloc(element *head, double val, void *dp);

element *clearstack(element *head);

long int stacklen(element *head);

char isinstack(element *head,double val);

element stackdel(element *head,double val);

element dequeue(element *head);

element *enqueue(element *head, double val, void *dp);

element *clearqueue(element *head);

long int queuelen(element *head);

char isinqueue(element *head,double val);

element queuedel(element *head,double val);


#endif
