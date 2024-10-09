/*
######################################################################
#
# utils.c
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 10/8/24
# Licensed under the GNU General Public License version 3 or later.
# Portions taken from the sna library by Carter T. Butts,
#  (self-licensed under GPL).
#
# Part of the R/ergmgp package
#
# This file contains various utility routines to be called by other
# ergmgp functions.
#
######################################################################
*/

#include "utils.h"


/*QUEUE/STACK/LIST ROUTINES-------------------------------------------------*/


/*Skip Lists (Ascending)*/

int isInSList(slelement *head, double val)
/*Is val in the skip list pointed to by head?*/
{
  slelement *ep;
  int i;
  
  /*Return immediately if no list*/  
  if(head==NULL)
    return 0;
  /*Otherwise, go looking for trouble*/
  ep=head;
  for(i=head->depth;i>=0;i--)
    for(;(ep->next[i]!=NULL)&&(ep->next[i]->val<val);ep=ep->next[i]);
  /*We have reached the end of the line; check to see where we are*/
  if((ep->next[0]==NULL)||(ep->next[0]->val>val))
    return 0;
  return 1;    
}


slelement *slistDelete(slelement *head, double val)
/*Remove the first element matching val from the list, if present (otherwise, leave things intact).  A pointer to the deleted element is returned (or NULL if not found); note that we do not have to return a pointer to the list, since the head will remain unchanged regardless.  (This means that empty lists are not NULLed out, for what that's worth.)*/
{
  slelement *ep,**epp,**tochange,*rp;
  int i,olddepth;

//  Rprintf("\tTrying to delete item with val %.1f\n",val);

  /*Return immediately if no list*/  
  if(head==NULL)
    return NULL;

  /*Try to delete val*/
  tochange=(slelement **)R_alloc(head->depth+1,sizeof(slelement *));
  ep=head;
  for(i=head->depth;i>=0;i--){  /*Filter down*/
    for(;(ep->next[i]!=NULL)&&(ep->next[i]->val<val);ep=ep->next[i]);
    tochange[i]=ep;    /*Record the last element at each level*/
  }

  /*We have reached the end of the line; is there a value here to delete?*/
  if((ep->next[0]==NULL)||(ep->next[0]->val>val))
    return NULL;
//  Rprintf("\t\tStopped search at %.1f\n",ep->next[0]->val);
  rp=ep->next[0];
  /*Apparently, ep->next[0] should be scheduled for demolition*/
  for(i=0;i<=head->depth;i++){  /*Update pointers from search trace*/
    if(tochange[i]->next[i]!=rp)  /*Nothing deeper goes here*/
      break;
    tochange[i]->next[i]=rp->next[i];  /*Always depth-safe*/
  }

  /*Update the maximum list depth and the list length*/
//  Rprintf("\t\tNew length %.0f, old depth %d\n",head->val-1,head->depth);
  head->val--;
  olddepth=head->depth;
  for(;(head->depth>0)&&(head->next[head->depth]==NULL);head->depth--);
//  Rprintf("\t\t\tNew depth %d\n",head->depth);
  if(head->depth!=olddepth){
    epp=head->next;
    head->next=(slelement **)R_alloc(head->depth+1,sizeof(slelement *));
    for(i=0;i<=head->depth;i++)
      head->next[i]=epp[i];
  }
  
  /*Return the item pointer*/
//  Rprintf("\t\tAbout to return item %.1f\n",rp->val);
  return rp;
}


slelement *slistInsert(slelement *head, double val, void *dp)
/*Add the indicated item to a list, returning a pointer to the updated head (which might have been changed, if called with NULL).  Note that the return value will simply be head unless except in the case noted above.*/
{
  slelement *ep,*new,**tochange,**epp;
  int i;
  
  /*Create the new element*/
  new=(slelement *)R_alloc(1,sizeof(slelement));
  new->depth=(int)rgeom(0.5);
  new->next=(slelement **)R_alloc(new->depth+1,sizeof(slelement *));
  new->val=val;
  new->dp=dp;
  
  /*Add it to the list*/
  if(head==NULL){  /*If no list, create from whole cloth....*/
    head=(slelement *)R_alloc(1,sizeof(slelement));
    head->val=1.0;
    head->dp=NULL;
    head->depth=new->depth;
    head->next=(slelement **)R_alloc(head->depth+1,sizeof(slelement *));
    for(i=0;i<=head->depth;i++){   /*Head -> new, new -> NULL*/
      head->next[i]=new;
      new->next[i]=NULL;
    }
  }else{           /*Otherwise, insert in place*/
    head->val++;    /*Increment the list length indicator*/
    tochange=(slelement **)R_alloc(MAX(new->depth,head->depth)+1, sizeof(slelement *));
    ep=head;
    for(i=head->depth;i>=0;i--){
      for(;(ep->next[i]!=NULL)&&(ep->next[i]->val<val);ep=ep->next[i]);
      tochange[i]=ep;    /*Record the last element at each level*/
    }
    if(new->depth>head->depth){  /*For new levels, head is last element*/
      epp=head->next;
      head->next=(slelement **)R_alloc(new->depth+1,sizeof(slelement *));
      for(i=0;i<=head->depth;i++)
        head->next[i]=epp[i];
      for(i=head->depth+1;i<=new->depth;i++){
        tochange[i]=head;
        head->next[i]=NULL;
      }
      head->depth=new->depth;
    }
    for(i=0;i<=new->depth;i++){  /*Adjust pointers for leftward elements*/
      new->next[i]=tochange[i]->next[i];
      tochange[i]->next[i]=new;
    }
  }
  
  /*Return the possibly updated head pointer*/
  return head;
}


void slistPrint(slelement *head)
/*Troubleshooting utility to print the contents of an slist.*/
{
  slelement *ep,*ep2;
  int count=0,i,j;
  
  Rprintf("SkipList Printout:\n");
  if(head==NULL)
    Rprintf("\tEmpty list.\n");
  else{
    for(ep=head;ep!=NULL;ep=ep->next[0]){
      Rprintf("  %d: [%.1f] ",count++,ep->val);
      for(i=0;i<=ep->depth;i++){
        for(j=0,ep2=head;(ep2!=NULL)&&(ep->next[i]!=ep2);ep2=ep2->next[0])
          j++;
        Rprintf("--%03d",j);
      }
      Rprintf("\n");
    }
  }
  Rprintf("--------------------\n");
}


slelement *slistSearch(slelement *head, double val)
/*Return a pointer to the first element matching val, or else NULL.*/
{
  slelement *ep;
  int i;

  /*Return immediately if no list*/  
  if(head==NULL)
    return NULL;
  /*Otherwise, go looking for trouble*/
  ep=head;
  for(i=head->depth;i>=0;i--)
    for(;(ep->next[i]!=NULL)&&(ep->next[i]->val<val);ep=ep->next[i]);
  /*We have reached the end of the line; check to see where we are*/
  if((ep->next[0]==NULL)||(ep->next[0]->val>val))
    return NULL;
  return ep->next[0];    
}


/*Linked Lists (Ascending)*/

int isInList(element *head, double val)
/*Is val in the sorted list pointed to by head?*/
{
  element *ep;
  
  for(ep=head;(ep!=NULL)&&(ep->val<val);ep=ep->next);
  if(ep==NULL)
    return 0;
  if(ep->val==val)
    return 1;
  return 0;
}


element *listInsert(element *head, double val, void *dp)
/*Add a new element to a sorted list, returning a pointer to the updated
list.*/
{
  element *elem,*ep;

  /*Initialize the element*/
  elem=(element *)R_alloc(1,sizeof(struct elementtype));
  elem->val=val;
  elem->dp=dp;
  elem->next=NULL;


  if(head==NULL){  /*If this is the only element, make it the head*/
    return elem;
  }else if(head->val>val){  /*If this is first, make it the head*/
    elem->next=head;
    return elem;
  }else{          /*Otherwise, traverse until we get to the right spot*/
    for(ep=head;(ep->next!=NULL)&&(ep->next->val<val);ep=ep->next);
    if(ep->next==NULL){   /*We ran out of list, apparently*/
      ep->next=elem;
      return head;
    }else{                /*We need to add elem after ep*/
      elem->next=ep->next;
      ep->next=elem;
      return head;
    }
  }
}


/*Stacks*/

element pop(element *head)
/*Pop an element from the stack pointed to by head*/
{
element rval;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=head;
}else{
   if(head->next==NULL){
      rval.val=head->val;
      rval.dp=head->dp;
      head=NULL;
      rval.next=NULL;
   }else{
      rval.next=head->next;
      rval.val=head->val;
      rval.dp=head->dp;
      head=rval.next;
   }
}

return rval;
}


element *push(element *head, double val, void *dp)
/*Adds element with value val to the stack, returning the head 
pointer.*/
{
element *newnode;

/*Create the new entry*/
newnode=(element *)R_alloc(1,sizeof(struct elementtype));

newnode->val=val;   /*Set the element value*/
newnode->dp=dp;

/*Set the next pointer equal to the current first entry (if any)*/
newnode->next=head;

/*Place the new node at the head of the stack*/
head=newnode;

return head;
}


element *pushCalloc(element *head, double val, void *dp)
/*Adds element with value val to the stack, returning the head 
pointer. This function uses R_alloc for memory allocation, so
memory allocated is automatically freed by R at the end of the .C
call.*/
{
element *newnode;

/*Create the new entry*/
newnode=(element *)R_alloc(1,sizeof(struct elementtype));

newnode->val=val;   /*Set the element value*/
newnode->dp=dp;

/*Set the next pointer equal to the current first entry (if any)*/
newnode->next=head;

/*Place the new node at the head of the stack*/
head=newnode;

return head;
}


long int stacklen(element *head)
/*Returns the length of the stack pointed to by head*/
{
element *p;
int count=0;

for(p=head;p!=NULL;p=p->next)
   count++;

return count;
}


char isinstack(element *head,double val)
/*Returns a 1 if val is in the stack pointed to by head, otherwise 0*/
{
element *p;

for(p=head;p!=NULL;p=p->next)
   if(p->val==val)
      return 1;

return 0;
}


element stackdel(element *head,double val)
/*! Find the element val in the stack pointed to by head and delete it,
returning the deleted element.*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=NULL;
}else if(head->val==val){
   rval.val=head->val;
   rval.dp=head->dp;
   rval.next=head->next;
   head=rval.next;
}else{
   for(p=head;(p->next!=NULL)&&(p->next->val!=val);p=p->next);
   if(p->next==NULL){
      rval.val=-1.0;
      rval.dp=NULL;
      rval.next=NULL;
   }else{
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      rval.next=p->next->next;
      p->next=rval.next;
   }
}
      
return rval;
}


/*Queues*/

element dequeue(element *head)
/*Dequeue an element from the queue pointed to by head*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=head;
}else{
   if(head->next==NULL){
      rval.val=head->val;
      rval.dp=head->dp; 
      head=NULL;
      rval.next=NULL;
   }else{
      for(p=head;p->next->next!=NULL;p=p->next);
      rval.next=NULL;
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      p->next=NULL;
   }
}
      
return rval;
}


element *enqueue(element *head, double val, void *dp)
/*Adds element with value val to the queue, returning the head 
pointer.*/
{
element *newnode;

/*Create the new entry*/
newnode=(element *)R_alloc(1,sizeof(struct elementtype));

newnode->val=val;   /*Set the element value*/
newnode->dp=dp;

/*Set the next pointer equal to the current first entry (if any)*/
newnode->next=head;

/*Place the new node at the head of the queue*/
head=newnode;

return head;
}

 
long int queuelen(element *head)
/*Returns the length of the queue pointed to by head*/
{
element *p; 
int count=0;

for(p=head;p!=NULL;p=p->next)
   count++;

return count;
}


char isinqueue(element *head,double val)
/*Returns a 1 if val is in the queue pointed to by head, otherwise 0*/
{
element *p;

for(p=head;p!=NULL;p=p->next)
   if(p->val==val)
      return 1;

return 0;
}


element queuedel(element *head,double val)
/*Find the element val in the queue pointed to by head and delete it,
returning the deleted element.*/
{
element rval,*p;

if(head==NULL){
   rval.val=-1.0;
   rval.dp=NULL;
   rval.next=NULL;
}else if(head->val==val){
   rval.val=head->val;
   rval.dp=head->dp;
   rval.next=head->next;
   head=rval.next;
}else{
   for(p=head;(p->next!=NULL)&&(p->next->val!=val);p=p->next);
   if(p->next==NULL){
      rval.val=-1.0;
      rval.dp=NULL;
      rval.next=NULL;
   }else{
      rval.val=p->next->val;
      rval.dp=p->next->dp;
      rval.next=p->next->next;
      p->next=rval.next;
   }
}
      
return rval;
}

