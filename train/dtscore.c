// Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,Maryland, U.S.A.  All rights reserved.


#include  <stdio.h>
#include <stdlib.h>

const int ALPHABET_SIZE = 4;
const int MODEL_LEN = 9;
const int MAX_STRING_LEN = 500000;
const int SIMPLE_MODEL_LEN = 6;
const double  MAX_LOG_DIFF = -20.0;
const double epsilon=0.001;

#include "strarray.h"
#include "context_prob.h"
#include "delcher.h"
#include "gene.h"
#include "oc1_prob.h"

char TRAIN_DIR[500]="";
String_Array  * Delta [6];
double  Ch_Ct [4];


int String_Score(char X [], int T);
void  Indep_Eval  (char X [], int T, double P [], double & Prob_X);
void  Find_Stop_Codons  (char X [], int T, int Stop []);
double Perc(char *Data,int start,int end);
double comp_entropy(char *Data,int len,int start,int incr);


// ********************* DECISION TREES DECLARATIONS *****************

#define NO_OF_VAR 2 	/* number of variables for the decision trees */
#define NO_OF_TREES 10     /* number of decision trees  

/* this file contains freq's of all in-frame hexamers in all
   training sequences. */
#define IF_6MER_TRAIN        "exon.hexfreq"

/* this file contains the freq's of all hexamers in all training
   sequences, including exons, introns, and intergenic DNA */
#define TRAIN_6MERS          "all.hexfreq"

int no_of_trees = NO_OF_TREES;
int no_of_dimensions = NO_OF_VAR;
int no_of_categories = 2;
struct tree_node **exonroots = NULL;
struct tree_node **startroots = NULL;
struct tree_node **intronroots = NULL;
struct tree_node **lastroots = NULL;
struct tree_node **snglroots = NULL;

double *HexamerVec;    /* 4096 vector of in-frame hexamer frequencies */
int get_if_hexamer(char *,char *,int,double *);
double ifhexamer(char *Data,int start, int stop);
void classify(POINT *,struct tree_node **);
extern struct tree_node *read_tree(char *);
void decision_tree();

main  (int argc, char * argv [])
{
  FILE  * fp;
  char line[MAX_STRING_LEN],exon[MAX_STRING_LEN];
  int i,len,score,seqno,ind;
  char seq[MAX_STRING_LEN];
  char Name[50];
  POINT point;
  double exscore;
  double score1,score2,hexamerval;
  char label[2];

  // usage -------------------------------------------------------------
  
  if (argc != 4 ) { 
    fprintf (stderr, "USAGE:  %s <dtorfs-file> <deltas_file> <train-dir> \n", argv [0]);
    exit (-1);
  }
  strcpy(TRAIN_DIR,argv[3]);

  // read delta's file -------------------------------------------------

  strcpy (line, TRAIN_DIR);
  strcat (line, argv[2]);
   
  fp = File_Open (line, "r");
  for  (i = 0;  i < 6;  i ++) {                  // Create six frame models
     
    Delta [i] = new String_Array (MODEL_LEN,ALPHABET_SIZE);
    
    // Read the model that  build-imm  created
    Delta [i] -> Read (fp);      
    
    Delta [i] -> Convert_To_Logs ();
  }
  fclose (fp);

 // ---------------------------- INIT_DECISION_TREES ------------------

  decision_tree();  


  // read exons file ----------------------------------------------

  fp = fopen(argv[1], "r");
  if  (fp == NULL) {
    fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
    exit (0);
  }
  
  while(fgets(line,MAX_STRING_LEN,fp)!=NULL) {
    sscanf(line,"%s %s %d %s %lf %lf",label,Name,&seqno,seq,&score1,&score2);
    len=strlen(seq);

    Ch_Ct[0]=Ch_Ct[1]=Ch_Ct[2]=Ch_Ct[3]=0;
    for(i=0;i<len;i++) {
      switch(seq[i]) {
      case  'a' :
      case  't' :
	Ch_Ct [0] += 1.0;
	break;
      case  'c' :
      case  'g' :
	Ch_Ct [1] += 1.0;
	break;
      }
    }
    Ch_Ct [2] = Ch_Ct [1];               // Counts are for *both* strands
    Ch_Ct [3] = Ch_Ct [0];
    for  (i = 0;  i < 4;  i ++)
      // Convert to log of proportion of at vs. gc
      Ch_Ct [i] = log (Ch_Ct [i] / (2.0 * len));  

    ind=1;

    switch (label[0]) {
    case 'S':
      hexamerval=ifhexamer(seq,0,len-1-3);
      point.dimension = (float *) malloc (sizeof(float) * 2);
      no_of_dimensions = 2;
      point.dimension[0] = (double) len;
      point.dimension[1] = hexamerval;
      classify(&point, snglroots);
      break;
    case 'L':
      if(len<=9) hexamerval=0;
      hexamerval=ifhexamer(seq,0,len-1-3);
      point.dimension = (float *) malloc (sizeof(float) * 3);
      no_of_dimensions = 3;
      point.dimension[0] = score1;
      point.dimension[1] = (double) len;
      point.dimension[2] = hexamerval;
      classify(&point, lastroots);
      if(score1==-99) ind=0;
      break;
    case 'F':
      if(len<=6) hexamerval=0;
      hexamerval=ifhexamer(seq,0,len-1);
      point.dimension = (float *) malloc (sizeof(float) * 3);
      no_of_dimensions = 3;
      point.dimension[0] = score1;
      point.dimension[1] = (double) len;
      point.dimension[2] =  hexamerval;
      classify(&point, startroots);
      if(score1==-99) ind=0;
      break;
    case 'E':
      if(len<=6) hexamerval=0;
      hexamerval=ifhexamer(seq,0,len-1);
      point.dimension = (float *) malloc (sizeof(float) * 4);
      no_of_dimensions = 4;
      point.dimension[0] = score1;
      point.dimension[1] = score2;
      point.dimension[2] = (double) len;
      point.dimension[3] = hexamerval;
      classify(&point, exonroots);
      if(score1==-99||score2==-99) ind=0;
      break;
    case 'I':
      if(len<=6) hexamerval=0;
      hexamerval=ifhexamer(seq,0,len-1);
      point.dimension = (float *) malloc (sizeof(float) * 4);
      no_of_dimensions = 4;
      point.dimension[0] = score1;
      point.dimension[1] = score2;
      point.dimension[2] = (double) len;
      point.dimension[3] = hexamerval;
      classify(&point, intronroots);
      if(score1==-99||score2==-99) ind=0;
      break;
    }

    exscore = point.prob[0];
    free(point.dimension); 
    
    strcpy(exon,"a");
    strcat(exon,seq);
    if(label[0]=='L' || label[0]=='S') { seq[len-3]='\0';len-=3;}
    score=String_Score(exon,len);
    if(label[0]=='L' || label[0]=='S') len+=3;
    printf("%s %s %d %d %f %d\n",label,Name,len,score,exscore,ind);
    

  }

}

double comp_entropy(char *Data,int len,int start,int incr)
{
  double freq[4];
  int i;
  double entropy,sum;

  freq[0]=freq[1]=freq[2]=freq[3]=0;

  i=start;
  while(i<len) {
    switch (Data[i]) {
    case 'a':
    case 'A': freq[0]++;break;
    case 'c':
    case 'C': freq[1]++;break;
    case 'g':
    case 'G': freq[2]++;break;
    case 't':
    case 'T': freq[3]++;break;
    }
    i+=incr;
  }

  freq[0]+=epsilon;
  freq[1]+=epsilon;
  freq[2]+=epsilon;
  freq[3]+=epsilon;
  
  sum=freq[0]+freq[1]+freq[2]+freq[3];
  freq[0]=freq[0]/sum;
  freq[1]=freq[1]/sum;
  freq[2]=freq[2]/sum;
  freq[3]=freq[3]/sum;

  entropy=-freq[0]*log(freq[0])-freq[1]*log(freq[1])-freq[2]*log(freq[2])-freq[3]*log(freq[3]);

  return(entropy);
}

void decision_tree()
{

  FILE *treenames;              /* file hold all trees' names */
  int i;
  char *s;
  char treefile[LINESIZE];
  char filename[550];
  char allfilename[550];

  struct tree_node
    **exontrees = NULL;       /* trees for classifing internal exons */
  struct tree_node
    **introntrees = NULL;       /* trees for classifing introns */
  struct tree_node
    **starttrees = NULL;       /* trees for classifying initial exons */
  struct tree_node
    **lasttrees = NULL;       /* trees for classifying final exons */
  struct tree_node
    **sngltrees = NULL;       /* trees for classifying single exons */

  double *if_6mervec;     /* in-frame hexamer freqs */

  /****************************************/
  /* read in the decision trees           */
  /****************************************/
  
  /* read in sets of trees for: 
     1. exons
     2. introns
     3. initial exons
     4. final exons  */
  
  /* allocate space for exonroots */
  exontrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);
  introntrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);
  starttrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);
  lasttrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);
  sngltrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);

  /* the decision tree names must be in a file called "treenames" */
  sprintf(filename,"%streenames",TRAIN_DIR);
  if((treenames = fopen(filename, "r")) == NULL) {
    fprintf (stderr, "Decision tree names file can not be opened.");
    exit(-1);
  }
  for (i=0; i< (no_of_trees * 5); i++) {
    if ( fgets(treefile, LINESIZE, treenames) == NULL) {
      fprintf (stderr, "Not enough trees.");
      exit(-1);
    }
    s = strchr(treefile, '\n');
    *s = '\0';

    sprintf(filename,"%s%s",TRAIN_DIR,treefile);
    /* first 10 trees are for initial exons -- starttrees */
    if (i < 10) {
      if ((starttrees[i] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);	
      }
    }
    /* next 10 trees are for internal exons -- exontrees */
    else if (i < 20) {
      if ((exontrees[i-10] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);	
      }
    }
    /* next 10 trees are for final exons -- lasttrees */
    else if (i < 30) {
      if ((lasttrees[i-20] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);
      }
    }
    /* next 10 trees are for introns */
    else if (i < 40) {
      if ((introntrees[i-30] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);
      }
    }
/* next 10 trees are for single exons */
    else if (i < 50) {
      if ((sngltrees[i-40] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);
      }
    }
  }
  fclose(treenames);

  exonroots = exontrees;
  intronroots = introntrees;
  startroots = starttrees;
  lastroots = lasttrees;
  snglroots = sngltrees;


  /****************************************/
  /* read in the hexamer frequency vector */
  /****************************************/

  if_6mervec = (double *) malloc (sizeof(double) * 4096); 
  //get_if_hexamer(IF_6MER_TRAIN, TRAIN_6MERS, 4096, if_6mervec);
  sprintf(filename,"%s%s",TRAIN_DIR,IF_6MER_TRAIN);
  sprintf(allfilename,"%s%s",TRAIN_DIR,TRAIN_6MERS);
  get_if_hexamer(filename, allfilename, 4096, if_6mervec);
  HexamerVec    = if_6mervec;
}


int String_Score(char X [], int T) 
{
  double  Max, Min, Sum, S [7];
  int  i, Has_Stop [7];
  int Score[7];

  Fast_Evaluate (X, T, MODEL_LEN, * Delta [0], * Delta [1], * Delta [2],S [0]);
  Fast_Evaluate (X, T, MODEL_LEN, * Delta [2], * Delta [0], * Delta [1],S [1]);
  Fast_Evaluate (X, T, MODEL_LEN, * Delta [1], * Delta [2], * Delta [0],S [2]);
  Fast_Evaluate (X, T, MODEL_LEN, * Delta [3], * Delta [4], * Delta [5],S [3]);
  Fast_Evaluate (X, T, MODEL_LEN, * Delta [5], * Delta [3], * Delta [4],S [4]);
  Fast_Evaluate (X, T, MODEL_LEN, * Delta [4], * Delta [5], * Delta [3],S [5]);

  //  if  (Use_Independent)
  Indep_Eval (X, T, Ch_Ct, S [6]);
  //else    S [6] = MIN_LOG_PROB_FACTOR * T;

  Find_Stop_Codons (X, T, Has_Stop);


  Max = - DBL_MAX;
  Min = DBL_MAX;
  for  (i = 0;  i < 7;  i ++)
    {
      if  (! Has_Stop [i])
	{
	  if  (S [i] > Max)
	    Max = S [i];
	  if  (S [i] < Min)
	    Min = S [i];
	}
    }

  assert (Max != - DBL_MAX && Min != DBL_MAX);

  if  (Min < Max + MAX_LOG_DIFF)
    Min = Max + MAX_LOG_DIFF;
  
  for  (i = 0;  i < 7;  i ++)
    if  (Has_Stop [i])
      S [i] = Min + MAX_LOG_DIFF;
    else if  (S [i] < Min)
      S [i] = Min;
  
  Sum = 0.0;
  for  (i = 0;  i < 7;  i ++)
    {
      S [i] = exp (S [i] - Min);
      Sum += S [i];
    }
  

  for  (i = 0;  i < 7;  i ++)
    {
      S [i] /= Sum;
      
      if  (Has_Stop [i]) {
	Score [i] = -1;
	
      }
      else
	{
	  Score [i] = int (100.0 * S [i]);
	  if  (Score [i] > 99)
	    Score [i] = 99;
	}
    }

   return(Score[0]);
}

void  Find_Stop_Codons  (char X [], int T, int Stop [])
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


void  Indep_Eval  (char X [], int T, double P [], double & Prob_X)
{
  int  i;
  
  Prob_X = 0.0;
  
  for  (i = 1;  i <= T;  i ++)
    switch  (X [i])
      {
      case  'a' :
	Prob_X += P [0];
	break;
      case  'c' :
	Prob_X += P [1];
	break;
      case  'g' :
	Prob_X += P [2];
	break;
      case  't' :
	Prob_X += P [3];
	break;
      }
  
  Prob_X = Max (Prob_X, MIN_LOG_PROB_FACTOR * T);
  
  return;
}

double Perc(char *Data,int start, int end)
{
  int i,gc;
  double ret;

  gc=0;

  for(i=start;i<=end;i++) {
    switch(Data[i]) {
    case 'C':
    case 'c':
    case 'G':
    case 'g':gc++;break;
    }
  }
  
  ret=gc*100/(end-start+1);

  return(ret);

}

int get_if_hexamer(char *file1,char *file2,int nvals,double *vector)
{ 
   FILE *fp1;
   FILE *fp2;
   int index;
   double hexfreq;

   /* Open files.  The first file stores the in-frame hexamer
      frequencies.  The second file stores the hexamer frequencies
      for all frames of exons, introns, and intergenic DNA */
   if (((fp1=fopen(file1,"r")) == NULL)  ||
       ((fp2=fopen(file2,"r")) == NULL))
     {
       printf ("vector file -> %s\n", file1);
       fprintf(stderr,"Could not open one of %s or %s\n",file1,file2);
       exit(-1);
     }
   /* the file format is hexamer and frequency, e.g., 
      AGATCA  0.000245.  We just want to keep the numbers; in
      particular the log ratio between the frequencies. */
   for ( index=0; index<nvals; ++index )
     {
       fscanf(fp1,"%*s %lf",&(vector[index]));
       fscanf(fp2,"%*s %lf",&hexfreq);
       if (vector[index] != 0) /* avoid taking log of 0 */
	 vector[index] = log(vector[index]/hexfreq);
       else
	 vector[index] = -15.0;
     }
   fclose(fp1);
   fclose(fp2);
   return(0);
}


/* here is the ifhexamer statistic computation.  Sum the log ratios
   (stored in HexamerVec) of the hexamers that actually occur in the
   current window.  Consider all 3 frames and sum separately for each
   frame.  Return the max.  NOTE: THIS SHOULD NOT NECESSARILY BE THE
   MAX -- PERHAPS RETURNING THE FRAME THAT IS KNOWN TO BE THE CORRECT
   ONE WILL IMPROVE PERFORMANCE.  TRY THIS LATER. */

double ifhexamer(char *Data, int start, int stop)
{ 
   int index,weight,val;                /* location of hexamer in HexamerVec */
   int pos1, phase,i;
   static double *hex;            /* compute in 3 frames, return max */
   static int firsttime=1; 
   double retval;

   /* Initialize */
   if ( firsttime )
   {
      firsttime=0;
      hex = (double *) malloc( sizeof(double) * 3 );
   }

   if(start<stop) {
     for (phase=0; phase<3; phase++) {
       hex[phase] = 0;
       /* loop through the window in-phase, 3 at a time */
       for ( pos1=phase+start; pos1+5<stop+1; pos1 += 3 )
	 { 
	   index = 0;
	   weight = 1;
	   /* compute index into hexamerVec based on current hexamer */
	   for ( i=5; i>=0; --i )
	     {
	       /* pos1 is the current position in in the window */
	       switch(Data[pos1+i]) {
	       case 'a': val=0; break;
	       case 'c': val=1; break;
	       case 'g': val=2; break;
	       case 't': val=3; break;
	       default: val=1;
	       }
	       
	       index += weight*val;
	       weight *= 4;
	     }
	   /* we store "acgt" as "0123" and compute a unique number
	      based on this string, which is the index into HexamerVec.
	      HexamerVec stores log ratios, so we can add them. */
	   hex[phase] += HexamerVec[index];
	 }
     }
   }
   else {
     for (phase=0; phase<3; phase++) {
       hex[phase] = 0;
       /* loop through the window in-phase, 3 at a time */
       for ( pos1=start-phase; pos1-5>stop-1; pos1 -= 3 )
	 { 
	   index = 0;
	   weight = 1;
	   /* compute index into hexamerVec based on current hexamer */
	   for ( i=5; i>=0; --i )
	     {
	       /* pos1 is the current position in in the window */
	       switch(Data[pos1-i]) {
	       case 'a': val=3; break;
	       case 'c': val=2; break;
	       case 'g': val=1; break;
	       case 't': val=0; break;

		 /*case 'a': val=0; break;
	       case 'c': val=1; break;
	       case 'g': val=2; break;
	       case 't': val=3; break;*/

	       default: val=1;
	       }
	       
	       index += weight*val;
	       weight *= 4;
	     }
	   /* we store "acgt" as "0123" and compute a unique number
	      based on this string, which is the index into HexamerVec.
	      HexamerVec stores log ratios, so we can add them. */
	   hex[phase] += HexamerVec[index];
	 }
     }
   }

   retval = hex[0];   /* return value, max of 3 phases */
   if (hex[1] > retval)
     retval = hex[1];
   if (hex[2] > retval)
     retval = hex[2];
   return(retval);
}


void classify (POINT *point,struct tree_node **roots)
{
  int j, t;
  struct tree_node *cur_node;
  double sum;
  double probtmp;

  /* initialize Exon and Intron prob */
  point->prob[0]=0;
  point->prob[1]=0;
  /* get length and skip if <= 2 */
  for (t=1; t<=no_of_trees; t++) {
    //printf("Tree %d\n",t);
    cur_node = roots[t-1];
    while (cur_node != NULL) {
      sum = cur_node->coefficients[no_of_dimensions];

      for (j=1;j<=no_of_dimensions;j++)
	sum += cur_node->coefficients[j-1] * point->dimension[j-1];

      if (sum < 0) {
	if (cur_node->left != NULL) 
	  cur_node = cur_node->left;
	else {
	  /* New for prob. classification, added by Xin Chen */
	  probtmp =
	    (double)cur_node->left_count[(cur_node->left_cat)-1]/
	    cur_node->left_total;
	  
	  if(cur_node->left_cat == 1) {
	    point->prob[0] += probtmp;
	    point->prob[1] += 1-probtmp;
	  }
	  else {
	    point->prob[0] += 1-probtmp;
	    point->prob[1] += probtmp;
	  }
	  
	  /****/
	  /*printf("left cat= %d  left_count[%d]= %d no_points= %d prob=%f\n",
		 cur_node->left_cat,
		 cur_node->left_cat,
		 cur_node->left_count[cur_node->left_cat],
		 cur_node->left_total,point->prob[0]);*/
	  
	  break;
	}
      }
      else {
	if (cur_node->right != NULL) 
	  cur_node = cur_node->right;
	else {
	  /* New for prob. classification, added by Xin Chen */
	  probtmp =
	    (double)cur_node->right_count[(cur_node->right_cat)-1]/
	    cur_node->right_total;
	  
	  if(cur_node->right_cat == 1) {
	    point->prob[0] += probtmp;
	    point->prob[1] += 1-probtmp;
	  }
	  else {
	    point->prob[0] += 1-probtmp;
	    point->prob[1] += probtmp;
	  }
	  
	  /****/
	  /*printf("right cat= %d  right_count[%d]= %d no_points= %d prob=%f\n",
		 cur_node->right_cat,
		 cur_node->right_cat,
		 cur_node->right_count[cur_node->right_cat],
		 cur_node->right_total,point->prob[0]);*/
	    
	  break;
	}
      }
    }
  }

  point->prob[0] /= (double) no_of_trees;
  point->prob[1] /= (double) no_of_trees;
  if(point->prob[0] >= point->prob[1]){
    point->category = 1;
  }
  else{
    point->category = 2;
  }

}

