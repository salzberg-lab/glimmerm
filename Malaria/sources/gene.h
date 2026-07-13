// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.
//
//   A. L. Delcher
//
//     File:  ~delcher/TIGR/gene.h
//  Version:  1.01  31 Jul 97
//
//  Common routines for DNA sequences.
//


#ifndef  __GENE_H_INCLUDED
#define  __GENE_H_INCLUDED


const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const int  MAX_LINE = 300;


char  Complement  (char);
char  Filter  (char);
int  Read_String  (FILE *, char * &, long int &, char [], int);
int  Read_Multi_String  (FILE *, char * &, long int &, char [], int,int);



char  Complement  (char Ch)

/* Returns the DNA complement of  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      case  'r' :          // a or g
        return  'y';
      case  'y' :          // c or t
        return  'r';
      case  's' :          // c or g
        return  's';
      case  'w' :          // a or t
        return  'w';
      case  'm' :          // a or c
        return  'k';
      case  'k' :          // g or t
        return  'm';
      case  'b' :          // c, g or t
        return  'v';
      case  'd' :          // a, g or t
        return  'h';
      case  'h' :          // a, c or t
        return  'd';
      case  'v' :          // a, c or g
        return  'b';
      default :            // anything
        return  'n';
     }
  }



char  Filter  (char Ch)

//  Return a single  a, c, g or t  for  Ch .  Choice is to minimize likelihood
//  of a stop codon on the primary strand.

  {
   switch  (tolower (Ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  Ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        return  'c';
    }
  }



int  Read_String  (FILE * fp, char * & T, long int & Size, char Name [],
                   int Partial)

/* Read next string from  fp  (assuming FASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read.  If  Partial  is
*  true, then the first line must have 2 integers indicating positions
*  in the string and only those positions will be put into  T .  If
*  Partial  is false, the entire string is put into  T .  Sets  Name
*  to the first string after the starting '>' character. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct, Found = FALSE;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch == EOF)
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }

int  Read_Multi_String  (FILE * fp, char * & T, long int & Size, char Name [],
                   int Partial, int no)

/* Read next string #no from  fp  (assuming MULTIFASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read.  If  Partial  is
*  true, then the first line must have 2 integers indicating positions
*  in the string and only those positions will be put into  T .  If
*  Partial  is false, the entire string is put into  T .  Sets  Name
*  to the first string after the starting '>' character. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct, Found = FALSE;
   int i;

   i=0;

   while(i<no) {
     while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
       ;

     if  (Ch == EOF)
       return  FALSE;
     if (Ch == '>') i++;
   }

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }


#endif
