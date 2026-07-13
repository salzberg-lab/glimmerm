//Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.

//  A. L. Delcher
//
//  Written:  12 May 97
//  Revised:  19 May 97   Convert Delta & Lambda to arrays, combine output files,
//                        and condense so that Lambda isn't used.
//  Revised:  31 Jul 97   Make outputting the simple model file optional
//                        to save space.  Rename to  build-imm
//
//     File:  ~delcher/TIGR/build-imm.c
//  Version:  1.01  31 Jul 97
//
//  Create and output the lambda and delta files for the nonuniform model
//  based on the data in the input file.  Uses the simple, statistical
//  significance training model.


#include  "delcher.h"


const int  MODEL_LEN = 9;             //  n  in the paper
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const int  MAX_STRING_LEN = 6000;
const double  MINIMUM_DELTA_VALUE = 0.1;
const int  WINDOW_SIZE = 48;
const int  DEFAULT_OUTPUT_SIMPLE_MODEL = FALSE;


#include  "context.h"


void  Find_Stop_Codons  (char [], int, int []);
void  Process_Options  (int, char * []);
int  Read_String  (FILE *, char * &, long int &, char []);


int  Output_Simple_Model = DEFAULT_OUTPUT_SIMPLE_MODEL;



main  (int argc, char * argv [])

  {
   FILE  * Infile;
   String_Array  * Delta [6], * Lambda [6];
   char  * Line, File_Name [MAX_NAME_LEN], Name [MAX_NAME_LEN];
   double  Incr_Val;
   long int  i, j, Input_Size, Num_Strings, T;
   
   if  (argc < 2)
       {
        fprintf (stderr, "USAGE:  %s <string-file>\n", argv [0]);
        exit (-1);
       }

   Process_Options (argc, argv);

   Infile = File_Open (argv [1], "r");

   for  (i = 0;  i < 6;  i ++)
     {
      Delta [i] = new String_Array (MODEL_LEN, ALPHABET_SIZE);
      Delta [i] -> Set (0.0);
     }

   Num_Strings = 0;
   Incr_Val = 1.0;
   Line = (char *) Safe_malloc (INIT_SIZE);
   Line [0] = ' ';
   Input_Size = INIT_SIZE;

   while  (Read_String (Infile, Line, Input_Size, Name))
     {
      T = strlen (Line + 1);
      assert (T < Input_Size - 1);

{
 int  Stop [7];

 Find_Stop_Codons  (Line, T, Stop);
 if  (Stop [0])
     printf ("*** String %ld has stop codon\n", Num_Strings + 1);
}  

//  j%3  is the model determined by the frame.  Add  Incr_Val  to that
//  model at the entry corresponding to the string ending at position  j .

      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          Delta [j % 3] -> Incr (Line + j - i, i + 1, Incr_Val);

      Reverse_Complement (Line, T);

//  Add similarly for the reverse models

      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          Delta [3 + j % 3] -> Incr (Line + j - i, i + 1, Incr_Val);

      Num_Strings ++;
     }

//  Do each of the six frame models
   for  (i = 0;  i < 6;  i ++)
     {
      Lambda [i] = new String_Array (MODEL_LEN - 1, ALPHABET_SIZE);

//  Set IMM interpolation parameters based on counts in  Delta
      Lambda [i] -> Set_Lambda (* Delta [i]);

//  Convert counts obtained above into transition probabilities in model
      Delta [i] -> Normalize ();
     }

   if  (Output_Simple_Model)
       {
        strcpy (File_Name, argv [1]);
        strcat (File_Name, ".deltas.simple");
        Delta [0] -> Write (File_Name);
        for  (i = 1;  i < 6;  i ++)
          Delta [i] -> Append (File_Name);
       }

//  Convert IMM into equivalent simple model
   for  (i = 0;  i < 6;  i ++)
     Delta [i] -> Condense (* Lambda [i]);

   strcpy (File_Name, argv [1]);
   strcat (File_Name, ".deltas.context");
   Delta [0] -> Write (File_Name);
   for  (i = 1;  i < 6;  i ++)
     Delta [i] -> Append (File_Name);

   printf ("Number of strings = %ld\n", Num_Strings);
   return  0;
  }



void  Find_Stop_Codons  (char X [], int T, int Stop [])

//  Set  Stop [0 .. 6]  TRUE  or  FALSE   according to whether
//  X [1 .. T] has a stop codon in the corresponding reading frame.
//  Stop [6]  is always set  FALSE .

  {
   static int  Next [10] [4] =
     {{ 0,  2,  0,  1},     //  0  a, g
      { 3,  4,  5,  6},     //  1  t
      { 0,  2,  0,  6},     //  2  c
      { 7,  2,  7,  1},     //  3  ta
      { 8,  2,  0,  6},     //  4  tc
      { 7,  2,  0,  1},     //  5  tg
      { 9,  4,  5,  6},     //  6  ct, tt
      { 0,  2,  0,  1},     //  7  taa, tag, tga    Forward stop
      { 0,  2,  0,  1},     //  8  tca              Reverse stop
      { 7,  2,  7,  1}};    //  9  cta, tta         Reverse stop
   int  i, State;

   for  (i = 0;  i < 7;  i ++)
     Stop [i] = 0;

   State = 0;
   for  (i = 1;  i <= T;  i ++)
     {
      switch  (tolower (X [i]))
        {
         case  'a' :
           State = Next [State] [0];
           break;
         case  'c' :
           State = Next [State] [1];
           break;
         case  'g' :
           State = Next [State] [2];
           break;
         case  't' :
           State = Next [State] [3];
           break;
         default :
           fprintf (stderr, "Unexpected character %c\n", X [i]);
           State = 0;
        }
      if  (State == 7)
          Stop [i % 3] = TRUE;
      else if  (State > 7)
          Stop [3 + i % 3] = TRUE;
     }

   return;
  }



void  Process_Options  (int argc, char * argv [])

//  Process command-line options and set corresponding global switches
//  and parameters.
//
//    -s     Output simple probability model

  {
   int  i;

   for  (i = 2;  i < argc;  i ++)
     {
      switch  (argv [i] [0])
        {
         case  '-' :
           switch  (argv [i] [1])
             {
              case  's' :       // output simple probability model
                Output_Simple_Model = TRUE;
                break;
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
           break;
         default :
           fprintf (stderr, "Unrecognized option %s\n", argv [i]);
        }
     }

   return;
  }



int  Read_String  (FILE * fp, char * & T, long int & Size, char Name [])

/* Read next string from  fp  (assuming one string per line with each string
*  preceded by a name tag) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF). */

  {
   long int  Len;
   int  Ch;

   if  (fscanf (fp, "%s", Name) == EOF)
       return  FALSE;

   while  ((Ch = fgetc (fp)) != EOF && (Ch == ' ' || Ch == '\t'))
     ;
   if  (Ch == EOF)
       return  FALSE;

   ungetc (Ch, fp);
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '\n')
     {
      if  (isspace (Ch))
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
           break;
         case  's' :
           Ch = 'c';
           break;
         case  'w' :
           Ch = 'a';
           break;
         case  'r' :
           Ch = 'a';
           break;
         case  'y' :
           Ch = 'c';
           break;
         case  'm' :
           Ch = 'a';
           break;
         case  'k' :
           Ch = 'g';
           break;
         case  'b' :
           Ch = 'c';
           break;
         case  'd' :
           Ch = 'a';
           break;
         case  'h' :
           Ch = 'a';
           break;
         case  'v' :
           Ch = 'a';
           break;
         case  'n' :
           Ch = 'a';
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'a';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';

   return  TRUE;
  }



