// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.

//   A. L. Delcher
//
//     File:  ~delcher/TIGR/delcher.h
//  Version:  1.01  31 Jul 97
//
//  Common generic routines.
//


#ifndef  __DELCHER_H_INCLUDED
#define  __DELCHER_H_INCLUDED


#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <errno.h>


#define  TRUE  1
#define  FALSE  0
#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif


FILE *  File_Open  (const char *, const char *);
template <class DT>
  DT  Max  (DT, DT);
template <class DT>
  DT  Min  (DT, DT);
void *  Safe_malloc  (size_t);
void *  Safe_realloc  (void *, size_t);
char *  strdup  (char * &, const char *);
template <class DT>
void  Swap  (DT &, DT &);


FILE *  File_Open  (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }



template <class DT>
DT  Max  (DT A, DT B)

/* Return the larger of  A  and  B . */

  {
   if  (A > B)
       return  A;
     else
       return  B;
  }



template <class DT>
DT  Min  (DT A, DT B)

/* Return the smaller of  A  and  B . */

  {
   if  (A < B)
       return  A;
     else
       return  B;
  }



void *  Safe_malloc  (size_t Len)

/* Allocate and return a pointer to  Len  bytes of memory.
*  Exit if fail. */

  {
   void  * P;

   P = malloc (Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  malloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_realloc  (void * Q, size_t Len)

/* Reallocate memory for  Q  to  Len  bytes and return a
*  pointer to the new memory.  Exit if fail. */

  {
   void  * P;

   P = realloc (Q, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  realloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



template <class DT>
void  Swap  (DT & A, DT & B)

/* Swap the values in  A  and  B . */

  {
   DT  Save;

   Save = A;
   A = B;
   B = Save;

   return;
  }



char *  strdup  (char * & S1, const char * S2)

/* Allocate memory in  S1  for a copy of string  S2  and copy
*  it.  Return a pointer to  S1 . */

  {
   S1 = (char *) Safe_malloc (1 + strlen (S2));
   strcpy (S1, S2);

   return  S1;
  }



#endif
