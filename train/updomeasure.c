/* Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
 * Maryland, U.S.A.  All rights reserved.
 * 
 * sitesk.c compute a score for splice sites based on karlin's paper 
*/

#include  <stdio.h>
#include <string.h>
#include "delcher.h"
#include "gene.h"

const int  MODEL_LEN = 12;

char UPSTREAM_MODEL[MAX_LINE];
char DOWNSTREAM_MODEL[MAX_LINE];
char CODING_MODEL[MAX_LINE];

const double  MIN_LOG_PROB_FACTOR = -6.0;
#include "icm.h"

double *C[4];

void  Read_Models  ();
void compute_codings(int start,int end, double *cod);
void  Indep_Eval  (char X [], int T, double P [], double & Prob_X);

void main ( int argc, char * argv [])
{ 
  char *sequence;
  int seqlen;
  int scorelen;
  int start;
  int i,j;
  double *IU,*ID;
  int stop[3];
  int model;
  double score,logscore,logdiff;
  double cod[3];
  int max1=1;
  int measure;
  double Ch_All[4];

  if  (argc < 12) {
    fprintf (stderr, "USAGE:  %s <sequence> <start> <len> <upstream_model> <downstream_model> <coding_model> <which_measure> <freq_a> <freq_c> <freq_g> <freq_t>\n",   /* measure is 1 for upstream and -1 for downstream */
	     argv [0]);
    exit (EXIT_FAILURE);
  } 

  strcpy(UPSTREAM_MODEL,argv[4]);
  strcpy(DOWNSTREAM_MODEL,argv[5]);
  strcpy(CODING_MODEL,argv[6]);

  Read_Models();

  seqlen=strlen(argv[1]);
  sequence = (char *) malloc((seqlen+2)*sizeof(char));
  if (sequence == NULL) {
    fprintf(stderr,"Memory allocation for input sequence failure.\n"); 
    abort();
  }
  sequence[0]='a';
  strcpy(sequence+1,argv[1]);

  start=atoi(argv[2]);
  scorelen=atoi(argv[3]);

  measure=atoi(argv[7]);
  for(i=8;i<12;i++) 
    Ch_All[i-8]=log(atof(argv[i]));
  
  // upstream measure   
  if(measure>0) {
    IU=(double *) malloc((seqlen+2)*sizeof(double));
    if (IU == NULL) {
      fprintf(stderr,"Memory allocation for upstream failure.\n");
      abort();
    }
  }

  // downstream measure   
  if(measure<0) {
    ID=(double *) malloc((seqlen+2)*sizeof(double));
    if (ID == NULL) {
      fprintf(stderr,"Memory allocation for downstream failure.\n");
      abort();
    }
  }

  // coding measure
  for(i=0;i<4;i++) {
    C[i]=(double *) malloc((seqlen+2)*sizeof(double));
    if (C[i] == NULL) {
      fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
      abort();
    }
  }

  if(measure>0) IU[start]=0;
  if(measure<0) ID[start]=0;
  C[3][start]=0;

  stop[0]=0;
  stop[1]=1;
  stop[2]=2;

  for(i=1;i<=seqlen;i++) {

    if(i<MODEL_LEN) {
      if(measure>0) IU[i]=IU[i-1]+log(get_prob_of_window2(i-1,UMODEL,sequence+1));
      if(measure<0) ID[i]=ID[i-1]+log(get_prob_of_window2(i-1,DMODEL,sequence+1));
    }
    else {
      if(measure>0) IU[i]=IU[i-1]+log(get_prob_of_window1(i-MODEL_LEN+1,UMODEL,sequence));
      if(measure<0) ID[i]=ID[i-1]+log(get_prob_of_window1(i-MODEL_LEN+1,DMODEL,sequence));
    }

    C[0][i]=0;
    C[1][i]=0;
    C[2][i]=0;

    if(i>=3) { 
      if((sequence[i-2]=='t' && sequence[i-1]=='g' && sequence[i]=='a')||
	 (sequence[i-2]=='t' && sequence[i-1]=='a' && sequence[i]=='a')||
	 (sequence[i-2]=='t' && sequence[i-1]=='a' && sequence[i]=='g')) {
	stop[i%3]=i;
      }
    }
     
    Indep_Eval(sequence+i-1,1,Ch_All,score);
    C[3][i]=C[3][i-1]+score;

    for(j=0;j<3;j++) {
      if(i>stop[j]) {
	
	if(i>=stop[j]+MODEL_LEN) {
	  model = (int)(i-MODEL_LEN-stop[j])%3;
	  score=get_prob_of_window1(i-MODEL_LEN+1,MODEL[model],sequence);
	}
	else { 
	  model = ((int)(i-stop[j]-1)%3 +1)%3;
	  score=get_prob_of_window2(i-stop[j]-1,MODEL[model],sequence+stop[j]+1);
	}

	assert(score!=0);
	C[j][i]=log(score);
      }
    }

  }

  if(measure>0) logscore = IU[start+scorelen]-IU[start];
  if(measure<0) logscore = ID[start+scorelen]-ID[start];
  logdiff=logscore-C[3][start+scorelen]+C[3][start];

  compute_codings(start+1,start+scorelen,cod);

  if(logscore<cod[0] || logscore<cod[1] || logscore<cod[2]) { max1=0;}

  printf("%f %d",logdiff,max1);
}
  
    
void  Read_Models  ()

//  Read in the probability models of the non/coding regions.

{
  char File_Name[MAX_LINE];

  // upstream
  FILE  * fp;
  fp = File_Open (UPSTREAM_MODEL, "r");   // maybe rb ?
  UMODEL=Read_NonCoding_Model (fp);
  fclose (fp);

  // downstream
  fp = File_Open (DOWNSTREAM_MODEL, "r");   // maybe rb ?
  DMODEL=Read_NonCoding_Model (fp);
  fclose (fp);

  // coding
  fp = File_Open (CODING_MODEL, "r");   // maybe rb ?
  Read_Scoring_Model (fp);
  fclose (fp);
  
  return;
}


void compute_codings(int start,int end, double *cod)
{
  int i;
  int j;

  for(j=0;j<3;j++) { cod[j]=0;}

  if(start<end) {
    for(i=start;i<=end;i++) {
      for(j=0;j<3;j++) {
	if(C[j][i]==0) cod[j]-=10;
	else cod[j]+=C[j][i];
      }
    }
  }
  else {
    for(i=start;i>=end;i--) {
      for(j=0;j<3;j++) {
	if(C[4+j][i]==0) cod[j]-=10;
	else cod[j]+=C[4+j][i];
      }
    }
  }
}
  

void  Indep_Eval  (char X [], int T, double P [], double & Prob_X)

//  Set  Prob_X  to the log of the probability of generating DNA string
//  X [1 .. T]  using the independent logs of probabilities of single
//  characters in  P [] .

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
