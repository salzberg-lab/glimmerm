//Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
// Maryland, U.S.A.  All rights reserved.

//read param from config_file
/*
*   malgene.cpp was designed by Mihaela PERTEA starting from glimmer.cpp
*   This program intends to find open reading frames in the file named
*  on the command line and scores them using the delta information
*  in the files whose names are prefixed by the second command-line
*  parameter.
*  It was designed to find also introns 
*/


#include  "delcher.h"
#include  "gene.h"
#include "oc1.h"

// ********************** DEFINE SECTION **************************

// Define GT_AG only to print the GT AG of the sequence file
// Define MESSAGE to make a very verbose version of the program!!
//#define GT_AG
//#define MESSAGE
//#define GENESCORE
//#define ATG
//#define STOP

#define INTERGENE 50; // minimum length of intergenic region


// ********************* CONSTANT DECLARATIONS **********************

const int  DEFAULT_MIN_GENE_LEN = 60;

const int  DEFAULT_THRESHOLD_SCORE = 90;
const int  DEFAULT_USE_INDEPENDENT = TRUE;
const int  ORF_SIZE_INCR = 1000;


const int  MODEL_LEN = 12;
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const double  MAX_LOG_DIFF = -20.0;

#include  "context.h"
#include "icm.h"

const unsigned int  OK = 0x0;

// ########################### WEIGHT PARAMETERS ########################

const double WEIGHT_FEx=0.95;
const double WEIGHT_Ex=1;
const double WEIGHT_LEx=1;
const double WEIGHT_SEx=1; 
const double WEIGHT_In=0.9; 
const double WEIGHT_Interg=0.08;

// ******************** STRUCTURE DEFINITIONS ***********************

struct Gene_Info
{
  char line[MAX_LINE];
  Gene_Info * link;

};

struct Splice_Site
{
  long int poz; // pozition of the site
  char type; // 4=for.start; 1=for.gt; 2=for.ag; 3=for.stop; 8=rev.stop; 7=rev.gt; 6=rev.ag; 5 =rev.start
  double score[3];              // best score in all frames
  struct Splice_Site * interg;  // next gene (in order) from node
  struct Splice_Site * next;    // next node in List
  struct Splice_Site * link[3]; // link to best score in all frames
  long int end; // for for.stop indicates position of for.start, and for
                // for rev.start indicates position of rev.stop 
  unsigned char bad_link[3];
};

struct intron
{
  long int gt,ag;
  struct intron *leg;
};

struct Node {
    long int val;
    int mark;
};


struct myint
{
  int yn; // 1 for GT, 2 for AG, -2 for CT, -1 for AC, 3 for ATG, -3 for CAT
  double score;
};



// ******************** FUNCTION DECLARATIONS ***********************


void  Find_Stop_Codons  (char [], int, int []);
void  Indep_Eval  (char [], int, double [], double &);
int  If_Stop  (char *);
void  Permute  (int [], int);
void  Process_Options  (int, char * []);
void  Score_String  (char *, long int, long int, long int, double [ALPHABET_SIZE], int [7]);
double String_Score_S0(char X [], int T);
double String_Score_S1(char X [], int T);
double String_Score_SMax(char X [], int T);
void  Simple_Score  (char [], int, int, double [ALPHABET_SIZE], int [7]);
void  Transfer  (char * , long int, int Len,char * , long int );
int Score_Gene_F(int Frame,  long int start, long int i,intron *Cap);
int Score_Gene_R(int Frame,  long int start, long int i,intron *Cap);
int StopCod(long int i);
void printreverse(intron *Cap);
void printgene(long int Hi,long int Lo,int Frame,intron *Cap);
void printsites(intron *Cap);
int  Is_Acceptor  (const int *, double *, double, int);
int  Is_Donor  (const int *, double *, double, int);
int  Is_Atg  (const int *, double *,double);
int  Is_Stop  (const int *, double *,double);
void DP_Left(long int stop);
char inverse(char);
int in_frame(long int stop, Splice_Site * span,long int * laststop);
double enhance_score(double score,double accthr,double donthr);
void  Read_Probability_Model(char *);
void  Read_NC_Model  ();

void compute_vect(char *seq,double *r,int frame);
double scalarprod(double *r1,double *r2);

// ********************* VARIABLE DECLARATIONS **********************
double DT_FEx_Min;
double DT_Ex_Min;
double DT_LEx_Min;
double DT_SEx_Min;
double DT_In_Min;

double DT_Down=0.6;
double DT_Up=0.8;
int NONCOD_LEN=99;

// default values
long int MAX_INTRON_LENGTH = 4000;
long int MIN_INTRON_LENGTH = 5;
long int MIN_EXON_LENGTH = 20;
long int MAX_EXON_LENGTH = 3000;
long int MAX_GENE_LENGTH = 7000;


double  Ch_Ct [ALPHABET_SIZE] = {0.0};   // Count of number of occurrences of acgt's in entire genome
double Ch_All[ALPHABET_SIZE]={0.0};
char  * Data, *AllData;
long int  Data_Len, AllData_Len;
long int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
char  * Orf_Buffer;
long int  Orf_Buffer_Len;
int  Threshold_Score = DEFAULT_THRESHOLD_SCORE;
int  Use_Independent = DEFAULT_USE_INDEPENDENT;
myint *IF;
myint *Stop;
int intergene;
long int last,no;
char TRAIN_DIR[500]="";

int Use_Start_Site = 1;
int Use_Stop_Site = 1;
int is_stop=0;
int is_start=0;
double Acc_Thr = -99;
double Don_Thr = -99;
double Acc_Thr_Last = -99;
double Don_Thr_First = -99;
double Atg_Thr = -99;
double Stop_Thr = -99;
int Use_Filter = 1;
int Acc_Win_Len;
int Don_Win_Len;

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
struct tree_node **uproots = NULL;
struct tree_node **downroots = NULL;

double *HexamerVec;    /* 4096 vector of in-frame hexamer frequencies */
int get_if_hexamer(char *,char *,int,double *);
double ifhexamer(long int start, long int stop);
void classify(POINT *,struct tree_node **);
extern struct tree_node *read_tree(char *);
void decision_tree();

int offset;
int ngene;

// ************************** IMM scoring
double *C[8];
double *IU,*ID;
double *RU,*RD;
//double *IM;
double SimEx(long int start, long int end, int frame);
void compute_codings(long int start,long int end, double *cod);
int basetoint(char c);

// ************************** Main PROGRAM ***************************

main  (int argc, char * argv [])
{
  FILE  * fp;
  char  File_Name [MAX_LINE], Name [MAX_LINE];
  char atgstr[MAX_LINE];
  int  ret;
  long int  i, j, k,  Input_Size;
  int B[200];
  double score;
  int istacc, istdon;

  // ---------------------------- OPTIONS ------------------------------

  if  (argc < 2) {
    fprintf (stderr,
	     "USAGE:  %s <genome-file> [options] \n",
	     argv [0]);
    exit (-1);
  }


  Process_Options (argc, argv);   // Set global variables to reflect status of
                                  // command-line options.


  // ----------------------------------
  //         read config_file
  // ----------------------------------

  strcpy(File_Name,TRAIN_DIR);
  strcat(File_Name,"config_file");
  
  fp = fopen (File_Name, "r");
  if  (fp == NULL)
    {
      fprintf (stderr, "ERROR:  Unable to open file %s\n", File_Name);
      exit (0);
    }

  // read min intron
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    errno = 0;
    i=strtol(Name,NULL,10);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad min intron length in config_file: %s\n", Name);
    else
      // set different intron length if wanted
      ;
  }
  
  // read max intron
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    errno = 0;
    i=strtol(Name,NULL,10);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad max intron length in config_file: %s\n", Name);
    else
      if(i>MAX_INTRON_LENGTH) MAX_INTRON_LENGTH=i+50;
  }

  // read min exon
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    errno = 0;
    i=strtol(Name,NULL,10);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad min exon length in config_file: %s\n", Name);
    else
      if(i<MIN_EXON_LENGTH) MIN_EXON_LENGTH =i;
      else // set different min exon length if wanted
	;
      if(i<3) MIN_EXON_LENGTH =3;
  }

  // read max exon
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    errno = 0;
    i=strtol(Name,NULL,10);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad max exon length in config_file: %s\n", Name);
    else
      if(i>MAX_EXON_LENGTH) MAX_EXON_LENGTH=i+100;
  }
  
  // read max gene
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    errno = 0;
    i=strtol(Name,NULL,10);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad max exon length in config_file: %s\n", Name);
    else
      if(i>MAX_GENE_LENGTH) MAX_GENE_LENGTH=i+1000;
  }

  // read if trees were used
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    istacc=atoi(Name);
    assert(istacc == 0 || istacc ==1);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    istdon=atoi(Name);
    assert(istdon == 0 || istacc ==1);
  }

  // read splice site scores

  if(fgets(Name,MAX_LINE,fp)==NULL) { // acceptor - no filter
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(Acc_Thr==-99 && !Use_Filter) {
      errno = 0;
      score = strtod (Name, NULL);
      if  (errno == ERANGE)
	fprintf (stderr, "ERROR:  Bad acceptor threshold score in config_file: %s\n", Name);
      else
	Acc_Thr = score;
    }
  }

  if(fgets(Name,MAX_LINE,fp)==NULL) { // donor - no filter
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(Don_Thr==-99 && !Use_Filter) {
      errno = 0;
      score = strtod (Name, NULL);
      if  (errno == ERANGE)
	fprintf (stderr, "ERROR:  Bad donor threshold score in config_file: %s\n", Name);
      else
	Don_Thr = score;
    }
  }

  // read filter windows length
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    Acc_Win_Len=atoi(Name);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    Don_Win_Len=atoi(Name);
  }
  

  if(fgets(Name,MAX_LINE,fp)==NULL) { // acceptor - with filter
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(Acc_Thr==-99 && Use_Filter) {
      errno = 0;
      score = strtod (Name, NULL);
      if  (errno == ERANGE)
	fprintf (stderr, "ERROR:  Bad acceptor threshold score in config file: %s\n", Name);
      else
	Acc_Thr = score;
    }
  }

  if(fgets(Name,MAX_LINE,fp)==NULL) { // donor - with filter
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(Don_Thr==-99 && Use_Filter) {
      errno = 0;
      score = strtod (Name, NULL);
      if  (errno == ERANGE)
	fprintf (stderr, "ERROR:  Bad donor threshold score in config_file: %s\n", Name);
      else
	Don_Thr = score;
    }
  }

  // read start site score if necessary
  if(fgets(Name,MAX_LINE,fp)==NULL) { 
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(atoi(Name)) {
      is_start=1;
      if(fgets(Name,MAX_LINE,fp)==NULL) { 
	fprintf(stderr,"config_file is not in corresponding format.\n");
	abort();
      }
      else {
	errno = 0;
	score = strtod (Name, NULL);
	if  (errno == ERANGE)
	  fprintf (stderr, "ERROR:  Bad start site threshold score in config_file: %s\n", Name);
	else
	  Atg_Thr = score;
      }
    }
    else {
      Use_Start_Site=0;
    }
  }

  
  //readind min DT
 
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    DT_FEx_Min=atof(Name);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    DT_Ex_Min=atof(Name);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
   DT_LEx_Min =atof(Name);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
    DT_SEx_Min=atof(Name);
  }
  if(fgets(Name,MAX_LINE,fp)==NULL) {
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else { 
   DT_In_Min =atof(Name);
  }

  // read first and last splice site scores

  if(fgets(Name,MAX_LINE,fp)==NULL) { // short last exon acceptor score
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    errno = 0;
    score = strtod (Name, NULL);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad acceptor threshold score in config_file: %s\n", Name);
    else
      Acc_Thr_Last = score;
  }
  
  if(fgets(Name,MAX_LINE,fp)==NULL) { // short first exon donor score
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    errno = 0;
    score = strtod (Name, NULL);
    if  (errno == ERANGE)
      fprintf (stderr, "ERROR:  Bad acceptor threshold score in config_file: %s\n", Name);
    else
      Don_Thr_First = score;
  }

  // read stop site score if necessary
  if(fgets(Name,MAX_LINE,fp)==NULL) { 
    fprintf(stderr,"config_file is not in corresponding format.\n");
    abort();
  }
  else {
    if(atoi(Name)) {
      is_stop=1;
      if(fgets(Name,MAX_LINE,fp)==NULL) { 
	fprintf(stderr,"config_file is not in corresponding format.\n");
	abort();
      }
      else {
	errno = 0;
	score = strtod (Name, NULL);
	if  (errno == ERANGE)
	  fprintf (stderr, "ERROR:  Bad start site threshold score in config_file: %s\n", Name);
	else
	  Stop_Thr = score;
      }
    }
    else {
      Use_Stop_Site=0;
    }
  }

  fclose(fp);

  //-----------------------------
  //      end config_file
  //-----------------------------
  

  intergene=INTERGENE;

  // ---------------------------- INIT_DECISION_TREES ------------------

  decision_tree();  

  // --------------------------- READ DATA ------------------------------

  fp = File_Open (argv [1], "r");

  AllData = (char *) Safe_malloc (INIT_SIZE);
  Input_Size = INIT_SIZE;
  
  // Read entire genome into  Data [1 .. ]
  // it seams that if there is not enough space, the space is increased
  Read_String (fp, AllData, Input_Size, Name, FALSE);   
  
  fclose (fp);
  
  AllData_Len = strlen (AllData + 1);

  printf("GlimmerM (Version 3.0)\n");    
  printf("Sequence name: %s\n",Name); 

  printf("Sequence length: %ld bp\n",AllData_Len);

  for  (i = 1;  i <= AllData_Len;  i ++) {

    // Converts all characters to  acgt
    AllData [i] = Filter (tolower (AllData [i]));
    
    switch  (AllData [i]) {
      
    case  'a' :
      Ch_Ct [0] += 1.0;
      Ch_All [0] += 1.0;
      break;
    case  't' :
      Ch_Ct [0] += 1.0;
      Ch_All [3] += 1.0;
      break;
    case  'c' :
      Ch_Ct [1] += 1.0;
      Ch_All [1] += 1.0;
      break;
    case  'g' :
      Ch_Ct [1] += 1.0;
      Ch_All [2] += 1.0;
      break;
    }
  }

 // --------------------- SOME DATA STATISTICS -------------------
  
  Ch_Ct [2] = Ch_Ct [1];               // Counts are for *both* strands
  Ch_Ct [3] = Ch_Ct [0];
  for  (i = 0;  i < 4;  i ++) {
    // Convert to log of proportion of at vs. gc
    Ch_Ct [i] = log (Ch_Ct [i] / (2.0 * AllData_Len));  
    Ch_All[i] = log(Ch_All[i]/AllData_Len);
  }

  // the second file is used in order to find the name of deltas.context
  strcpy (File_Name, TRAIN_DIR);
  strcat (File_Name, "icm.model");
   
  Read_Probability_Model (File_Name);

  // read the noncoding models

  Read_NC_Model();

   // Create a buffer to hold orfs for further processing

 Orf_Buffer_Len = ORF_SIZE_INCR;
 Orf_Buffer = (char *) Safe_malloc (Orf_Buffer_Len);
 Orf_Buffer [0] = ' ';

 //  Prepare first codon which will be characters 1, 2 and 3
 //  for wraparound, since genome is circular.

 ngene=0;

 printf("\nPredicted genes/exons\n\n");
 printf("Gene Exon Strand  Exon            Exon Range      Exon\n");
 printf("   #    #         Type                           Length\n\n");

 
 if((float)AllData_Len/500000>1) Data_Len=500000;
 else Data_Len=AllData_Len%500000;

 Data=(char *) malloc((Data_Len+2)*sizeof(char));
 if (Data == NULL) {
   fprintf(stderr,"Memory allocation for Data failure.\n"); 
   abort();
 }

 IF=(myint *) malloc((Data_Len+2)*sizeof(myint));
 if (IF == NULL) {
   fprintf(stderr,"Memory allocation for intron failure.\n"); 
   abort();
 }

 Stop=(myint *) malloc((Data_Len+2)*sizeof(myint));
 if (Stop == NULL) {
   fprintf(stderr,"Memory allocation for intron failure.\n"); 
   abort();
 }
 
 for(i=0;i<8;i++) {
   C[i]=(double *) malloc((Data_Len+2)*sizeof(double));
   if (C[i] == NULL) {
     fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
     abort();
   }
 }
   
 // intron measure   
 /*IM=(double *) malloc((Data_Len+2)*sizeof(double));
   if (IM == NULL) {
   fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
   abort();
   }*/
 

 // upstream measure   
 IU=(double *) malloc((Data_Len+2)*sizeof(double));
 if (IU == NULL) {
   fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
   abort();
 }


 // downstream measure   
 ID=(double *) malloc((Data_Len+2)*sizeof(double));
 if (ID == NULL) {
   fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
   abort();
 }


 // upstream measure   
 RU=(double *) malloc((Data_Len+2)*sizeof(double));
 if (RU == NULL) {
   fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
   abort();
 }

 // downstream measure   
 RD=(double *) malloc((Data_Len+2)*sizeof(double));
 if (RD == NULL) {
   fprintf(stderr,"Memory allocation for coding %d failure.\n",i);
   abort();
 }
     



 for(offset=0;offset<=AllData_Len/500000;offset++) {

   if(offset<AllData_Len/500000) Data_Len=500000;
   else Data_Len=AllData_Len%500000;

   Data[0]=' ';
   strncpy(Data+1,AllData+offset*500000+1,Data_Len);
   Data[Data_Len+1]='\0';

   // ----------------------- Coding measure -----------------------------

   
   long int stop[3];
   int model;

   stop[0]=0;
   stop[1]=1;
   stop[2]=2;

   //IM[0]=0;
   IU[0]=0;
   ID[0]=0;
   C[3][0]=0;

   for(i=1;i<Data_Len+1;i++) {


     if(i<MODEL_LEN) {

       //IM[i]=get_prob_of_window2(i-1,IMODEL,Data+1);
       IU[i]=IU[i-1]+log(get_prob_of_window2(i-1,UMODEL,Data+1));
       ID[i]=ID[i-1]+log(get_prob_of_window2(i-1,DMODEL,Data+1));
     }
     else {
       //IM[i]=get_prob_of_window1(i-MODEL_LEN+1,IMODEL,Data);
       IU[i]=IU[i-1]+log(get_prob_of_window1(i-MODEL_LEN+1,UMODEL,Data));
       ID[i]=ID[i-1]+log(get_prob_of_window1(i-MODEL_LEN+1,DMODEL,Data));

       }

     C[0][i]=0;
     C[1][i]=0;
     C[2][i]=0;

     if(i>=3) { 
       if((Data[i-2]=='t' && Data[i-1]=='g' && Data[i]=='a')||
	  (Data[i-2]=='t' && Data[i-1]=='a' && Data[i]=='a')||
	  (Data[i-2]=='t' && Data[i-1]=='a' && Data[i]=='g')) {
	 stop[i%3]=i;
       }
     }
     
     
     Indep_Eval(Data+i-1,1,Ch_All,score);
     C[3][i]=C[3][i-1]+score;

     for(j=0;j<3;j++) {
       if(i>stop[j]) {

	 if(i>=stop[j]+MODEL_LEN) {
	   model = (int)(i-MODEL_LEN-stop[j])%3;
	   score=get_prob_of_window1(i-MODEL_LEN+1,MODEL[model],Data);
	 }
	 else { 
	   model = ((int)(i-stop[j]-1)%3 +1)%3;
	   score=get_prob_of_window2(i-stop[j]-1,MODEL[model],Data+stop[j]+1);
	 }

	 assert(score!=0);
	 C[j][i]=log(score);
       }
     }
     
     //     printf("%d  %f %f  %f %f  %f %f  %f %f %f\n",i,C[0][i],exp(C[3][i]),C[1][i],exp(C[3][i]),C[2][i],exp(C[3][i]),IU[i],IM[i],ID[i]);exit(0);
   }
   
   // reverse complement

   char *Copy;

   Copy=(char *) malloc((Data_Len+2)*sizeof(char));
   if (Copy == NULL) {
     fprintf(stderr,"Memory allocation for copy failure.\n");
     abort();
   }
   
   strcpy(Copy,Data);
   Reverse_Complement(Copy,Data_Len);
     
   stop[0]=0;
   stop[1]=1;
   stop[2]=2;
   
   RU[Data_Len+1]=0;
   RD[Data_Len+1]=0;
   C[7][Data_Len+1]=0;

   for(i=1;i<Data_Len+1;i++) {

     if(i<MODEL_LEN) {
       
       //IM[i]=get_prob_of_window2(i-1,IMODEL,Data+1);
       RU[Data_Len-i+1]=RU[Data_Len-i+2]+log(get_prob_of_window2(i-1,UMODEL,Copy+1));
       RD[Data_Len-i+1]=RD[Data_Len-i+2]+log(get_prob_of_window2(i-1,DMODEL,Copy+1));
     }
     else {
       //IM[i]=get_prob_of_window1(i-MODEL_LEN+1,IMODEL,Data);
       RU[Data_Len-i+1]=RU[Data_Len-i+2]+log(get_prob_of_window1(i-MODEL_LEN+1,UMODEL,Copy));
       RD[Data_Len-i+1]=RD[Data_Len-i+2]+log(get_prob_of_window1(i-MODEL_LEN+1,DMODEL,Copy));
       
     }


     C[4][Data_Len-i+1]=0;
     C[5][Data_Len-i+1]=0;
     C[6][Data_Len-i+1]=0;
     
     if(i>=3) { 
       if((Copy[i-2]=='t' && Copy[i-1]=='g' && Copy[i]=='a')||
	  (Copy[i-2]=='t' && Copy[i-1]=='a' && Copy[i]=='a')||
	  (Copy[i-2]=='t' && Copy[i-1]=='a' && Copy[i]=='g')) {
	 stop[i%3]=i;
       }
     }

     Indep_Eval(Copy+i-1,1,Ch_All,score);
     C[7][Data_Len-i+1]=C[7][Data_Len-i+2]+score;
     
     for(j=0;j<3;j++) {
       if(i>stop[j]) {
	 
	 if(i>=stop[j]+MODEL_LEN) {
	   model = ((int)(i-MODEL_LEN-stop[j])%3 )%3;
	   score=get_prob_of_window1(i-MODEL_LEN+1,MODEL[model],Copy);
	 }
	 else {
	   model = ((int)(i-stop[j]-1)%3 +1)%3;
	   score=get_prob_of_window2(i-stop[j]-1,MODEL[model],Copy+stop[j]+1);
	 }
	 assert(score!=0);
	 C[4+j][Data_Len-i+1]=log(score);
       }
     }
     //     printf("%d  %f %f  %f %f  %f %f  %f %f %f\n",i,C[0][i],exp(C[3][i]),C[1][i],exp(C[3][i]),C[2][i],exp(C[3][i]),IU[i],IM[i],ID[i]);exit(0);     
   }
   
   free(Copy);


   // ---------------------- SPLICE SITES --------------------------------
   
   // forward direction

   for(i=0;i<Data_Len+1;i++)    { IF[i].yn=0; Stop[i].yn=0; }
  
   // first look fot gt and ag and start sites
   for(i=80;i<=Data_Len-82;i++){

     if((Data[i]=='t' && Data[i+1]=='a' && Data[i+2]=='a') ||
	(Data[i]=='t' && Data[i+1]=='g' && Data[i+2]=='a') ||
	(Data[i]=='t' && Data[i+1]=='a' && Data[i+2]=='g') ) { // Deal w/ stop sites
       k=0;
       for(j=i-4;j<i+15;j++) {
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=0;break;
	 case 'C':
	 case 'c': B[k]=1;break;
	 case 'G':
	 case 'g': B[k]=2;break;
	 case 'T':
	 case 't': B[k]=3;break;
	 default: B[k]=1;
	 }
	 k++;
       }
#ifdef GT_AG	
       ret=Is_Stop(B,&score,Stop_Thr);
       printf("Found forw stop at %ld w/ score %f\n",i,score);
#endif

       if((!Use_Stop_Site) || Is_Stop(B,&score,Stop_Thr)) {
	 // add acceptor to list;
	  
	 Stop[i].yn=1;
	 if(!Use_Stop_Site) Stop[i].score=1;
	 else Stop[i].score=score;
	 
       }

     }
     
     if(Data[i]=='a' && Data[i+1]=='t' && Data[i+2]=='g') { // Deal w/ start sites
       k=0;
       for(j=i-12;j<i+7;j++) {
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=0;break;
	 case 'C':
	 case 'c': B[k]=1;break;
	 case 'G':
	 case 'g': B[k]=2;break;
	 case 'T':
	 case 't': B[k]=3;break;
	 default: B[k]=1;
	 }
	 k++;
       }
       
#ifdef GT_AG	
       ret=Is_Atg(B,&score,Atg_Thr);
       printf("Found start at %ld w/ score %f\n",i,score);
#endif

#ifdef ATG

       if(i<Data_Len-201) {
	 //if(i>200 && i<Data_Len-201) {
	 strncpy(atgstr,Data+i-1,201);
	 
	 atgstr[201]='\0';
      
	 score=String_Score_S0(atgstr,200);
	 
	 //score=String_Score_S1(atgstr,200);

	 //strncpy(atgstr,Data+i-201,201);
	 //atgstr[201]='\0';
	 //score=score-String_Score_SMax(atgstr,200);

	 if(score>0) {
	   IF[i].yn=3;
	   IF[i].score=-10;
	 }
       }
#endif
       
       if((!Use_Start_Site) || Is_Atg(B,&score,Atg_Thr)) {
	 // add acceptor to list;
	  
#ifdef GT_AG
	 //printf("Found forward start at %ld w/ score %f\n",i,score);
#endif

	 IF[i].yn=3;
	 if(!Use_Start_Site) IF[i].score=1;
	 else IF[i].score=score;
	 
       }
     }

     if(Data[i]=='a' && Data[i+1]=='g') { // Deal with acceptors
       k=0;
       for(j=i-80;j<i+82;j++){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=0;break;
	 case 'C':
	 case 'c': B[k]=1;break;
	 case 'G':
	 case 'g': B[k]=2;break;
	 case 'T':
	 case 't': B[k]=3;break;
	 default: B[k]=1;
	 }
	 k++;
       }

#ifdef GT_AG	
       ret=Is_Acceptor(B,&score, Acc_Thr, istacc);
       printf("Found ag at %ld w/ score %f\n",i,score);
#endif

       if(Is_Acceptor(B,&score, Acc_Thr, istacc)) {
	 // add acceptor to list;
	  
#ifdef GT_AG
	 //printf("Found forward ag at %ld w/ score %f\n",i,score);
#endif

	 IF[i].yn=2;
	 IF[i].score=score;
       }
     }

     if(Data[i]=='g' && Data[i+1]=='t') { // Deal with donors
       k=0;
       for(j=i-80;j<i+82;j++){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=0;break;
	 case 'C':
	 case 'c': B[k]=1;break;
	 case 'G':
	 case 'g': B[k]=2;break;
	 case 'T':
	 case 't': B[k]=3;break;
	 default: B[k]=1;
	 }
	 k++;
       }
       
#ifdef GT_AG
       ret=Is_Donor(B,&score, Don_Thr, istdon);
       printf("Found gt at %ld w/ score %f\n",i,score);
#endif
       
       if(Is_Donor(B,&score, Don_Thr, istdon)) {
	 // add donor to list;

#ifdef GT_AG
	 //printf("Found forward gt at %ld w/ score %f\n",i,score);
#endif
	
	 IF[i].yn=1;
	 IF[i].score=score;
       }
     }
   }


   // reversed direction

   // first look for gt and ag
   for(i=Data_Len-81;i>=82;i--){


     if((Data[i]=='a' && Data[i-1]=='t' && Data[i-2]=='c') ||
	(Data[i]=='a' && Data[i-1]=='t' && Data[i-2]=='t') ||
	(Data[i]=='a' && Data[i-1]=='c' && Data[i-2]=='t')) { // Deal w/ stop sites
       k=0; 
       for(j=i+4;j>i-15;j--){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=3;break;
	 case 'C':
	 case 'c': B[k]=2;break;
	 case 'G':
	 case 'g': B[k]=1;break;
	 case 'T':
	 case 't': B[k]=0;break;
	 default: B[k]=2;
	 }
	 k++;
       }

#ifdef GT_AG	
       ret=Is_Stop(B,&score,Stop_Thr);
       printf("Found rev stop at %ld w/ score %f\n",i,score);
#endif

       if((!Use_Stop_Site) || Is_Stop(B,&score,Stop_Thr)) {
	 Stop[i].yn=-1;
	 if(!Use_Stop_Site) Stop[i].score=1;
	 else Stop[i].score=score;
	 }
     }


     if(Data[i]=='t' && Data[i-1]=='a' && Data[i-2]=='c') { // Deal w/ start sites
       k=0;
       for(j=i+12;j>i-7;j--){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=3;break;
	 case 'C':
	 case 'c': B[k]=2;break;
	 case 'G':
	 case 'g': B[k]=1;break;
	 case 'T':
	 case 't': B[k]=0;break;
	 default: B[k]=2;
	 }
	 k++;
       }

#ifdef GT_AG
       ret=Is_Atg(B,&score,Atg_Thr);
       printf("Found start at %ld w/ score %f\n",i,score);
#endif
       
#ifdef ATG

       if(i>200) {
	 //if(i>200 && i<Data_Len-201) {

	 strncpy(atgstr,Data+i-199,201);
	 atgstr[201]='\0';
	 
	 Reverse_Complement(atgstr,200);
	
	 score=String_Score_S0(atgstr,200);
	 //score=String_Score_S1(atgstr,200);

	 //strncpy(atgstr,Data+i-1,201);
	 //atgstr[201]='\0';
	 //score-=String_Score_SMax(atgstr,200);

	 if(score>0) {
	   IF[i].yn=-3;
	   IF[i].score=-10;
	 }
       }
#endif

       if((!Use_Start_Site) || Is_Atg(B,&score,Atg_Thr)) {
	 // add acceptor to list;
	
#ifdef GT_AG
	 //printf("Found reverse atg at %ld w/ score %f\n",i,score);
#endif
	 
	 IF[i].yn=-3;
	 if(!Use_Start_Site) IF[i].score=1;
	 else IF[i].score=score;
	 }
     }
  

     if(Data[i]=='t' && Data[i-1]=='c') { // Deal with acceptors
       k=0;
       for(j=i+80;j>i-82;j--){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=3;break;
	 case 'C':
	 case 'c': B[k]=2;break;
	 case 'G':
	 case 'g': B[k]=1;break;
	 case 'T':
	 case 't': B[k]=0;break;
	 default: B[k]=2;
	 }
	 k++;
       }

#ifdef GT_AG
       ret=Is_Acceptor(B,&score, Acc_Thr,istacc);
       printf("Found ag at %ld w/ score %f\n",i,score);
#endif

       if(Is_Acceptor(B,&score, Acc_Thr,istacc)) {
	 // add acceptor to list;
	
#ifdef GT_AG
	 //printf("Found reverse ag at %ld w/ score %f\n",i,score);
#endif

	 IF[i].yn=-2;
	 IF[i].score=score;
       }
     }
  

     if(Data[i]=='c' && Data[i-1]=='a') { // Deal with donors
       k=0;
       for(j=i+80;j>i-82;j--){
	 switch (Data[j]){
	 case 'A':
	 case 'a': B[k]=3;break;
	 case 'C':
	 case 'c': B[k]=2;break;
	 case 'G':
	 case 'g': B[k]=1;break;
	 case 'T':
	 case 't': B[k]=0;break;
	 default: B[k]=1;
	 }
	 k++;
       }

#ifdef GT_AG
       ret=Is_Donor(B,&score, Don_Thr,istdon);
       printf("Found gt at %ld w/ score %f\n",i,score);
#endif

       if(Is_Donor(B,&score, Don_Thr,istdon)) {
	 // add donor to list;

#ifdef GT_AG
	 //printf("Found reverse gt at %ld w/ score %f\n",i,score);
#endif

	 IF[i].yn=-1;
	 IF[i].score=score;
       }
     }
   }

   // keep the best score within filter lengths if the Use_Filter is active

   if(Use_Filter) {
     // filter acceptor sites
     for(i=1+Acc_Win_Len/2; i<Data_Len-Acc_Win_Len/2;i++) 
       if(IF[i].yn==2 || IF[i].yn==-2) 
	 for(j=i-Acc_Win_Len/2;j<i+Acc_Win_Len/2;j++) 
	   if(i!=j) {
	     if(IF[i].yn==IF[j].yn) {
	       if(IF[i].score<IF[j].score) { IF[i].yn=0; break; }
	       if(IF[j].score<IF[i].score) {IF[j].yn=0; }
	     }
	     if(IF[i].yn==-IF[j].yn) {
	       if(IF[i].score<IF[j].score) { IF[i].yn=0;  break; }
	       if(IF[j].score<IF[i].score) {IF[j].yn=0; }
	     }
	   }
     // filter donor sites
     for(i=1+Don_Win_Len/2; i<Data_Len-Don_Win_Len/2;i++) 
       if(IF[i].yn==1 || IF[i].yn==-1) 
	 for(j=i-Don_Win_Len/2;j<i+Don_Win_Len/2;j++) 
	   if(i!=j) {
	     if(IF[i].yn==IF[j].yn) {
	       if(IF[i].score<IF[j].score) { IF[i].yn=0; break; }
	       if(IF[j].score<IF[i].score) {IF[j].yn=0; }
	     }
	     if(IF[i].yn==-IF[j].yn) {
	       if(IF[i].score<IF[j].score) { IF[i].yn=0;  break; }
	       if(IF[j].score<IF[i].score) {IF[j].yn=0; }
	     }
	   }
   }

  /* for(i=31;i<Data_Len-30;i++) {
     //    if(IF[i].yn) to apply the process to all type of scores
     if(IF[i].yn!=0 && IF[i].yn!=3 && IF[i].yn!=-3) // only for splice sites
     for(j=i-30;j<i+30;j++) 
     if(i!=j) {
     if(IF[i].yn==IF[j].yn) {
     if(IF[i].score<IF[j].score) { IF[i].yn=0; break; }
     if(IF[j].score<IF[i].score) {IF[j].yn=0; }
     }
     if(IF[i].yn==-IF[j].yn) {
     if(IF[i].score<IF[j].score) { IF[i].yn=0;  break; }
     if(IF[j].score<IF[i].score) {IF[j].yn=0; }
     }
     }
     }
  */
	       
  //print found splice sites
  
#ifdef GT_AG 

   printf("Forward:\n"); 
   for(i=1;i<Data_Len;i++){
     if(IF[i].yn==3) printf("atg: %d %f\n",i,IF[i].score);
     if(IF[i].yn==1) printf("gt: %d %f\n",i,IF[i].score);
     if(IF[i].yn==2)  { 
       printf("ag: %d %f\n",i,IF[i].score);
      /*if(i>20 && i<Data_Len-24) {
     	for(j=i+2;j<i+MAX_EXON_LENGTH;j++) {
	  if(j<Data_Len-22 && IF[j].yn==1) {

	    for(k=0;k<20;k++) {
	      switch(Data[k-18+i]) {
	      case 'A':
	      case 'a': B[k]=0;break;
	      case 'C':
	      case 'c': B[k]=1;break;
	      case 'G':
	      case 'g': B[k]=2;break;
	      case 'T':
	      case 't': B[k]=3;break;
	      default: B[k]=1;
	      }
	    }

	    for(k=20;k<40;k++){
	      switch(Data[k-20+j]) {
	      case 'A':
	      case 'a': B[k]=0;break;
	      case 'C':
	      case 'c': B[k]=1;break;
	      case 'G':
	      case 'g': B[k]=2;break;
	      case 'T':
	      case 't': B[k]=3;break;
	      default: B[k]=1;
	      }
	      }
	   
	    //	    Is_AD(B,&score);

	    //	    printf("ag-gt: %d %d %f\n",i,j,score);
	  }
	}
	}*/
     }
   }

   printf("Reverse:\n");
   for(i=Data_Len;i>1;i--){
     if(IF[i].yn==-1) printf("gt: %d %f\n",i,IF[i].score);
     if(IF[i].yn==-2) printf("ag: %d %f\n",i,IF[i].score);
     if(IF[i].yn==-3) printf("atg: %d %f\n",i,IF[i].score);
   }
   exit(1);

#endif


 
 

   // ------------------------ MAIN CALLS ------------------------------

   DP_Left(Data_Len);

 }
   
 free(IF);
 for(i=0;i<8;i++) free(C[i]);
 free(ID);
 free(IU);
 free(RU);
 free(RD);
 free(Data);

 exit(0);


}

void printsites(intron *Cap)
{
  double poz[100];
  int l,i;
  
  l=-1;

  while(Cap!=NULL) {
    poz[++l]=IF[Cap->ag].score;
    poz[++l]=IF[Cap->gt].score;
    Cap=Cap->leg;
  }

  printf("          start ");

  for(i=l;i>=0;i--)
    printf("%f ",poz[i]);
 
  printf("stop");
}

void printreverse(intron *Cap)
{
  long int poz[100];
  int l,i;
  
  l=-1;

  while(Cap!=NULL) {
    poz[++l]=(Cap->ag)-2;
    poz[++l]=(Cap->gt)+1;
    Cap=Cap->leg;
  }

  for(i=l;i>=0;i--) {
    printf("%8ld ",poz[i]);fflush(stdout);
    last=poz[i];

  }
      
}

void DP_Left(long int stop)
{
  long int i,length,go_back,exonend[1000];
  Splice_Site * List, * site, *span, * max_link[3], *lastinterg;
  int temp_frame;
  intron *Cap,*node,*prev_node;
  double max_score[3],temp_score;
  POINT point;
  double hexamerval;
  long int laststopf[3],laststopr[3];
  int j,k,nex,exlength;
  char line[MAX_LINE];
  Gene_Info *geneCap,*geneNode;
  //  long int myalloc;
  double big_score;
  unsigned char Blink[3], bad_link;
  int stopvect;
#ifdef STOP
  char gene[100000],exon[100000],vect1[300],vect2[300];
  double score;
#endif

  //  fprintf(stderr,"In DP_Left!\n");fflush(stderr);fflush(stdout);

  //while(1);

  for(j=0;j<3;j++) {
    laststopf[j]=0;
    laststopr[j]=0;
  }

  // create new node in list of splice site List

  List=NULL;
  lastinterg=NULL;

  // done creation of node

  for(i=1;i<=stop;i++) {

    //fprintf(stderr,"i=%ld\n",i);fflush(stdout);

    // ** stop in forward
    if(i>=3 && StopCod(i)) {
      if (Stop[i-2].yn==1) { // found a real stop

	//      printf("1\n");fflush(stdout);
	
#ifdef MESSAGE
	printf("DP: found stop at %ld %c%c%c\n",i,Data[i-2],Data[i-1],Data[i]);fflush(stdout);
#endif

	max_score[0]=-HUGE_VAL;
	
	span=List;
	while(span != NULL && i-span->poz-1 < MAX_EXON_LENGTH) {
	  if(span->type == 4) {  // span is a start codon
	    length=i-span->poz + 1;
	    

	    if(length >= Min_Gene_Len && length %3 == 0 && in_frame(i-2,span,laststopf)) {
	      hexamerval=ifhexamer(span->poz,i-3);
	      point.dimension = (float *) malloc (sizeof(float) * 2);
	      no_of_dimensions = 2;
	      point.dimension[0] = length;
	      point.dimension[1] = hexamerval;
	      classify(&point, snglroots);
	      temp_score = point.prob[0];
	      free(point.dimension); 
	      
	      if(temp_score>DT_SEx_Min) {

		//	      	      printf("Exon %ld %ld: %f Score: %f Length: %d\n",span->poz,i,temp_score,temp_score*length+span->score[0],length);
		
#ifdef MESSAGE
		printf("  exon %ld %ld : score %f total_score %f\n",span->poz,i,temp_score,log(temp_score*WEIGHT_SEx)*length+span->score[0]);
#endif

		temp_score=log(temp_score*WEIGHT_SEx)*length+span->score[0]; // ??????? sau = temp_score*length*5;
		//temp_score=temp_score+span->score[0]; 
		
		if(temp_score > max_score[0]) {
		  max_score[0]=temp_score;
		  max_link[0]=span;
		}
	      }
	    }
	  }
	
	  if(span->type == 2) {
	    length=i-span->poz-1;
	    temp_frame=(3-length%3)%3;
	    
	    if(span->score[temp_frame]>-HUGE_VAL && length >2 &&in_frame(i-2,span,laststopf) ) {

	      if(length<=6) hexamerval=0;
	      else hexamerval=ifhexamer(span->poz+2,i);
	      
	      point.dimension = (float *) malloc (sizeof(float) * 3);
	      no_of_dimensions = 3;
	      point.dimension[0] = IF[span->poz].score;
	      point.dimension[1] = length;
	      point.dimension[2] = hexamerval;
	      classify(&point, lastroots);
	      bad_link=0;
	      temp_score = enhance_score(point.prob[0],IF[span->poz].score,-100);
	      //	    printf("temp_score=%f span->poz=%ld span->bad_link=%d\n",temp_score,span->poz,span->bad_link);
	      if(temp_score<0.2) { // try different ranges here
		bad_link=1;
		if(span->bad_link[temp_frame]>=2) { bad_link=2; }
	      }
	      free(point.dimension);
	      //	    printf("bad_link=%d\n",bad_link);
	      
	      if(length<=15 && IF[span->poz].score <Acc_Thr_Last) temp_score=0;

	      if(bad_link<2 && temp_score>DT_LEx_Min) {	    

		//	      	      printf("Exon %ld %ld: %f Score: %f Length: %d\n",span->poz+2,i,temp_score,temp_score*length+span->score[temp_frame],length);

#ifdef MESSAGE
		printf("  exon %ld %ld : score %f total_score %f\n",span->poz+2,i,temp_score,log(temp_score*WEIGHT_LEx)*length+span->score[temp_frame]);
#endif
	  
		temp_score=log(temp_score*WEIGHT_LEx)*length+span->score[temp_frame];
		//temp_score=temp_score+span->score[temp_frame];

		if(temp_score > max_score[0]) {
		  max_score[0]=temp_score;
		  max_link[0]=span;
		}
	      }
	    }
	  }
	  span=span->next;

       
	}

	length=0;go_back=i;
#ifdef STOP
	strcpy(gene,"");
#endif
	if(max_score[0]!= -HUGE_VAL) {

	  Cap=NULL;
	  span=max_link[0];
	  
	  //	printf("span->poz=%ld span->type=%d\n",span->poz,span->type);
	  big_score=0;nex=0;
	  while(span->type !=4) {
	    node=(intron *)malloc(sizeof(intron));  
	    //myalloc+=sizeof(intron);
	    //printf("Alloc+=%d: %ld\n",sizeof(intron),myalloc);
	    
	    if (node == NULL) {
	      fprintf(stderr,"Memory allocation for node failure.\n"); 
	      abort();
	    }

	    // here look for an exon with score >= 0.5
	    if(go_back==i) {
	      if(go_back-span->poz-1<=6) hexamerval=0;
	      else hexamerval=ifhexamer(span->poz+2,go_back);
	      
	      point.dimension = (float *) malloc (sizeof(float) * 3);
	      no_of_dimensions = 3;
	      point.dimension[0] = IF[span->poz].score;
	      point.dimension[1] = go_back-span->poz-1;
	      point.dimension[2] = hexamerval;
	      classify(&point, lastroots);
	      big_score += point.prob[0];
	      free(point.dimension);
	      nex++;
	    }
	    else {
	      hexamerval=ifhexamer(span->poz+2,go_back);
	      point.dimension = (float *) malloc (sizeof(float) * 4);
	      no_of_dimensions = 4;
	      point.dimension[0] = IF[span->poz].score;
	      point.dimension[1] = IF[go_back+1].score;
	      point.dimension[2] = go_back-span->poz-1;
	      point.dimension[3] = hexamerval;
	      classify(&point, exonroots);
	      big_score += point.prob[0];
	      free(point.dimension);
	      nex++;
	    }

	    length+=go_back-span->poz-1;
	    temp_frame=(3-length%3)%3;

#ifdef STOP
	    strncpy(exon,Data+span->poz+2,go_back-span->poz-1);
	    exon[go_back-span->poz-1]='\0';
	    strcat(exon,gene);
	    strcpy(gene,"");
	    strcat(gene,exon);
	    
# endif

	    node->ag=span->poz;
	    span=span->link[temp_frame];
	    node->gt=span->poz;
	    node->leg=Cap;
	    go_back=span->poz-1;
	    Cap=node;
	    span=span->link[temp_frame];

	    //	  printf("span->poz=%ld span->type=%d\n",span->poz,span->type);

	  }

	  stopvect=1;
#ifdef STOP

	  strncpy(exon,Data+span->poz,go_back-span->poz+1);
	  exon[go_back-span->poz+1]='\0';
	  strcat(exon,gene);
	  strcpy(gene,"");
	  strcat(gene,exon);
	  length=strlen(gene);

	  if(length>=Min_Gene_Len) {

	    strcpy(vect1,gene+length-Min_Gene_Len-1);
	    score=String_Score_S1(vect1,Min_Gene_Len-3);

	    strncpy(vect2,Data+i,Min_Gene_Len+1);
	    vect2[Min_Gene_Len+1]='\0';
	    length=strlen(vect2);
	    if(length==Min_Gene_Len) {
	      score-=String_Score_SMax(vect2,Min_Gene_Len);
	      if(score<0) stopvect=0;
	    }
	  }


# endif

	  if(go_back==i) { big_score=1; } // if only one exon consider it
	  else {
	    length=go_back-span->poz+1;
	    if(length<=6) hexamerval=0;
	    else hexamerval=ifhexamer(span->poz,go_back);
	    point.dimension = (float *) malloc (sizeof(float) * 3);
	    no_of_dimensions = 3;
	    point.dimension[0] = IF[go_back+1].score;
	    point.dimension[1] = length;
	    point.dimension[2] = hexamerval;
	    classify(&point, startroots);
	    big_score += point.prob[0];
	    nex++;
	    free(point.dimension);
	    big_score = big_score/nex;
	  }

#ifdef MESSAGE
	  printf("Try score_gene_f: link at %ld big_score=%f nex=%d\n",max_link[0]->poz,big_score,nex);
#endif	


	  if(stopvect && big_score>0.4 && Score_Gene_F((span->poz)%3,span->poz+2,i,Cap)) {
#ifdef MESSAGE
	    printf("linked to ag/start at %ld\n",max_link[0]->poz); fflush(stdout);
#endif


	    site=(Splice_Site *)malloc(sizeof(Splice_Site)); 
	    //myalloc+=sizeof(Splice_Site);
	    //printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	    if (site == NULL) {
	      fprintf(stderr,"Memory allocation for site failure.\n"); 
	      abort();
	    }

	    site->poz=i;
	    site->type=3; 
	    site->end=span->poz-2;
	    for(j=1;j<3;j++) { 
	      site->score[j]=-HUGE_VAL;
	      site->bad_link[j]=0;
	    }  
	    site->score[0]=max_score[0];

	    // *** here compute the downstream value
	    if(i>Data_Len-NONCOD_LEN || !is_stop || !is_start) { 
	      site->score[1]=WEIGHT_Interg; 
	      // or maybe	site->score[1]=DT_Down; ? 
	    }
	    else {
	      double logscore=ID[i+NONCOD_LEN]-ID[i];
	      double logdiff=logscore-C[3][i+NONCOD_LEN]+C[3][i];
	      double cod[3];
	      compute_codings(i+1,i+NONCOD_LEN,cod);
	      float max1=1;
	      if(logscore<cod[0] || logscore<cod[1] || logscore<cod[2]) { max1=0;}
	      point.dimension = (float *) malloc (sizeof(float) * 3);
	      no_of_dimensions = 3;
	      point.dimension[0] = Stop[i-2].score;
	      point.dimension[1] = logdiff;
	      point.dimension[2] = max1;
	      classify(&point, downroots);
	      site->score[1]=point.prob[0];
	      free(point.dimension); 
	    }
	    
	    point.dimension = (float *) malloc (sizeof(float) * 3);
	    no_of_dimensions = 3;
	    point.dimension[0] = IF[go_back+1].score;
	    point.dimension[1] = length;
	    point.dimension[2] = hexamerval;
	    classify(&point, startroots);
	    big_score += point.prob[0];
	    free(point.dimension); 
	    
	    site->interg=lastinterg;
	    lastinterg=site;
	    site->next=List;
	    for(j=0;j<3;j++) site->link[j]=NULL;
	    site->link[0]=max_link[0];
	    List=site;
	    
#ifdef GENESCORE
	    printf("Gene_score: %f\n",max_score[0]); fflush(stdout);
#endif

	}
	

	  // free Cap here
	  prev_node=Cap;
	  if(prev_node !=NULL) node=prev_node->leg;
	  while(prev_node !=NULL) {
	    free(prev_node);
	    //myalloc-=sizeof(intron);
	    //printf("Alloc-=%d: %ld\n",sizeof(intron),myalloc);
	    prev_node=node;
	    if(node !=NULL) node=prev_node->leg;
	  }
	}

	//###
      }
      j=i%3;
      laststopf[j]=i;
      //###
    }


      
    // ** gt in forward

    if(IF[i].yn==1) { // found a gt

#ifdef MESSAGE
      printf("DP: ##gt at %ld %c%c\n",i,Data[i],Data[i+1]);fflush(stdout);
#endif


      span=List;
      max_score[0]=-HUGE_VAL;
      max_score[1]=-HUGE_VAL;
      max_score[2]=-HUGE_VAL;

      while(span != NULL && i-span->poz < MAX_EXON_LENGTH) {

	if(span->type == 4) {
	  length=i-span->poz; 
	  temp_frame=length%3;
	  if(in_frame(i-temp_frame,span,laststopf)) {
	    if(length>=5) {

	      if(length<=6) hexamerval=0;
	      else hexamerval=ifhexamer(span->poz,i-1);
	      
	      point.dimension = (float *) malloc (sizeof(float) * 3);
	      no_of_dimensions = 3;
	      point.dimension[0] = IF[i].score;
	      point.dimension[1] = length;
	      point.dimension[2] = hexamerval;
	      classify(&point, startroots);
	      bad_link=0;
	      temp_score = enhance_score(point.prob[0],-100,IF[i].score);
	      if(temp_score<0.2) { // try different ranges here
		bad_link=1;
	      }
	      free(point.dimension);
	      
	      if(length<=15 && IF[i].score<Don_Thr_First ) temp_score=0;
	      if(temp_score>DT_FEx_Min) {
		//if(temp_score>0.02) {

		//				printf("Exon %ld %ld: %f Score: %f Length: %d\n",span->poz,i-1,temp_score,temp_score*length+span->score[0],length);

#ifdef MESSAGE
		printf("  exon %ld %ld : score %f total_score %f\n",span->poz,i-1,temp_score,log(temp_score*WEIGHT_FEx)*length+span->score[0]);
#endif
		
		temp_score=log(temp_score*WEIGHT_FEx)*length+span->score[0];
		//temp_score=temp_score+span->score[0];

		if(temp_score > max_score[temp_frame]) {
		  max_score[temp_frame]=temp_score;
		  max_link[temp_frame]=span;
		  Blink[temp_frame]=bad_link;
		}
	      }
	    }
	  }
	}

	if(span->type == 2 ) {

	  length=i-1-span->poz-1;

	  for(temp_frame=0;temp_frame<3;temp_frame++) {

	    if(length>=MIN_EXON_LENGTH) {
	      j=(3-(length-temp_frame)%3)%3;
	      if(span->score[j]>-HUGE_VAL && in_frame(i-temp_frame,span,laststopf)) {

		if(length<=6) hexamerval=0;
		else hexamerval=ifhexamer(span->poz+2,i);
		point.dimension = (float *) malloc (sizeof(float) * 4);
		no_of_dimensions = 4;
		point.dimension[0] = IF[span->poz].score;
		point.dimension[1] = IF[i].score;
		point.dimension[2] = length;
		point.dimension[3] = hexamerval;
		classify(&point, exonroots);
		bad_link=0;
		temp_score = enhance_score(point.prob[0],IF[span->poz].score,IF[i].score);
		if(temp_score<0.2) { // try different ranges here
		  bad_link=1;
		  if(span->bad_link[j]>=2) { bad_link=2; }
		}

		free(point.dimension);

		if(bad_link<2 && temp_score>DT_Ex_Min) {
		  //if(temp_score>0.02) {

		  //		  		  printf("Exon %ld %ld: %f Score: %f Length: %d\n",span->poz+2,i-1,temp_score,temp_score*length+span->score[j],length);

#ifdef MESSAGE
		  printf("  exon %ld %ld : score %f total_score %f\n",span->poz+2,i-1,temp_score,log(temp_score*WEIGHT_Ex)*length+span->score[j]);
#endif
		  
		  temp_score=log(temp_score*WEIGHT_Ex)*length+span->score[j];
		  //temp_score=temp_score+span->score[j];

		  if(temp_score > max_score[temp_frame]) {
		    max_score[temp_frame]=temp_score;
		    max_link[temp_frame]=span;
		    if(bad_link) Blink[temp_frame]=bad_link+span->bad_link[j]; 
		    else Blink[temp_frame]=0;
		  }
		}
	      }
	    }
	  }
	}
	span=span->next;
      }

      if(max_score[0]!= -HUGE_VAL || max_score[1] != -HUGE_VAL || max_score[2] != -HUGE_VAL) {

#ifdef MESSAGE
	for(j=0;j<3;j++) 
	  if(max_score[j]!= -HUGE_VAL) printf("frame %d: linked to ag/start at %ld\n",j,max_link[j]->poz); 
	fflush(stdout);
#endif

	site=(Splice_Site *)malloc(sizeof(Splice_Site)); 
	//myalloc+=sizeof(Splice_Site);
	//printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	if (site == NULL) {
	  fprintf(stderr,"Memory allocation for site failure.\n"); 
	  abort();
	}
	
	
	site->poz=i;
	site->type=1; 
	for(j=0;j<3;j++) {
	  site->score[j]=max_score[j];
	  site->bad_link[j]=Blink[j];  
	}
	site->interg=lastinterg;
	site->next=List;
	for(j=0;j<3;j++) 
	  if(max_score[j]!=-HUGE_VAL) site->link[j]=max_link[j];
	  else site->link[j]=NULL;
	List=site;
      }
      
    }


    // ** ag in forward

    if(IF[i].yn==2) { // found an ag

      //      printf("3\n");fflush(stdout);

#ifdef MESSAGE
      printf("DP: ##found ag at %ld %c%c\n",i,Data[i],Data[i+1]);fflush(stdout);
#endif
      
      max_score[0]=-HUGE_VAL;
      max_score[1]=-HUGE_VAL;
      max_score[2]=-HUGE_VAL;
      span=List;
      while(span!=NULL && i-span->poz < MAX_INTRON_LENGTH) {
	if(span->type==1) 
	  for(temp_frame=0;temp_frame<3;temp_frame++) 
	    if(span->score[temp_frame]>-HUGE_VAL) {
	      length=i+1-span->poz+1;
	      if(length>=5) {
		if(length<=6) hexamerval=0;
		else hexamerval=ifhexamer(span->poz,i+1);

		point.dimension = (float *) malloc (sizeof(float) * 4);
		no_of_dimensions = 4;
		point.dimension[0] = IF[span->poz].score;
		point.dimension[1] = IF[i].score;
		point.dimension[2] = length;
		point.dimension[3] = hexamerval;
		classify(&point, intronroots);
		bad_link=0;
		temp_score = enhance_score(point.prob[0],IF[i].score,IF[span->poz].score);
		if(temp_score<0.2) { // try different ranges here
		  bad_link=1;
		}
		free(point.dimension); 
	    
		if(temp_score >DT_In_Min  ) {

		  //		  printf("Intron %ld %ld: %f Length=%d Hex=%f\n",span->poz,i+1,temp_score,length,hexamerval);


		  //if(length>150) length=150;	
		  

#ifdef MESSAGE
		  printf("  intron %ld %ld : score %f total_score %f\n",span->poz,i,temp_score,log(temp_score*WEIGHT_In)*length+span->score[temp_frame]);
#endif
	      
		  temp_score=log(temp_score*WEIGHT_In)*length+span->score[temp_frame];
		  //temp_score=temp_score*0.9+span->score[temp_frame];

		  if(temp_score > max_score[temp_frame]) {
		    max_score[temp_frame]=temp_score;
		    max_link[temp_frame]=span;
		    if(bad_link) Blink[temp_frame]=bad_link+span->bad_link[temp_frame];
		    else Blink[temp_frame]=0;
		  }
		}
	      }
	    }
      	span=span->next;
      }


      if(max_score[0]!= -HUGE_VAL || max_score[1] != -HUGE_VAL || max_score[2] != -HUGE_VAL) {

#ifdef MESSAGE
	for(j=0;j<3;j++) 
	  if(max_score[j]!= -HUGE_VAL) printf("linked to gt at %ld\n",max_link[j]->poz); 
	fflush(stdout);
#endif

	site=(Splice_Site *)malloc(sizeof(Splice_Site));  
	//myalloc+=sizeof(Splice_Site);
	//printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);

	if (site == NULL) {
	  fprintf(stderr,"Memory allocation for site failure.\n"); 
	  abort();
	}

	site->poz=i;
	site->type=2; 
	site->interg=lastinterg;
	for(j=0;j<3;j++) { 
	  site->score[j]=max_score[j];
	  site->bad_link[j]=Blink[j];
	}
	site->next=List;
	for(j=0;j<3;j++) 
	  if(max_score[j]!=-HUGE_VAL) site->link[j]=max_link[j];
	  else site->link[j]=NULL;
	List=site;
      }

    }

    
    // ** start in forward
     
    if(IF[i-2].yn==3) {

      //      printf("4\n");fflush(stdout);

#ifdef MESSAGE
      printf("DP: Add start at %ld %c%c%c\n",i,Data[i-2],Data[i-1],Data[i]); fflush(stdout);
#endif

      span=List;

      // here compute the upstream value;
      double upstreamVal=WEIGHT_Interg;
      if( is_start && is_stop && i-2-NONCOD_LEN>=1) { 
	double logscore=IU[i-3]-IU[i-2-NONCOD_LEN-1];
	double logdiff=logscore-C[3][i-3]+C[3][i-2-NONCOD_LEN-1];
	double cod[3];
	compute_codings(i-2-NONCOD_LEN,i-3,cod);
	float max1=1;
	if(logscore<cod[0] || logscore<cod[1] || logscore<cod[2]) { max1=0;}
	point.dimension = (float *) malloc (sizeof(float) * 3);
	no_of_dimensions = 3;
	point.dimension[0] = IF[i-2].score;
	point.dimension[1] = logdiff;
	point.dimension[2] = max1;
	classify(&point, uproots);
	upstreamVal=point.prob[0];
	free(point.dimension); 
#ifdef MESSAGE
	printf(" atgscore=%f logscore=%f logdiff=%f max1=%f\n",IF[i-2].score,logscore,logdiff,max1);
#endif 
      }

      max_score[0]=log(0.5*(upstreamVal+WEIGHT_Interg)/2)*i;
      max_link[0]=NULL;
      go_back=0;
      while(span!=NULL && span->poz-go_back>0) {

	//	printf("span->poz=%ld, span->type=%d, span->score[0]=%f\n",span->poz,(int)span->type,span->score[0]);

	if(i-span->poz >= intergene)  // prima modificare
	  if(span->type==3 || span->type==5) {
	    if(span->end>go_back) go_back=span->end;
	    temp_score=log(0.5*(upstreamVal+span->score[1])/2)*(i-2-span->poz+1)+span->score[0];
	    if(temp_score > max_score[0]) {
	      max_score[0]=temp_score;
	      max_link[0]=span;
	    }
	  }
	span=span->interg;
      }

#ifdef MESSAGE
            printf(" score=%f upstreamVal=%f\n",max_score[0],upstreamVal);
#endif MESSAGE

      site=(Splice_Site *)malloc(sizeof(Splice_Site)); 
      //myalloc+=sizeof(Splice_Site);
      //printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
      if (site == NULL) {
	fprintf(stderr,"Memory allocation for site failure.\n"); 
	abort();
      }
      site->poz=i-2;
      site->type=4; 
      for(j=0;j<3;j++) {
	site->score[j]=-HUGE_VAL; 
	site->bad_link[j]=0;
      }
      site->score[0]=max_score[0];// a doua modificare: se reseteaza scorul la o gena noua
      site->interg=lastinterg;
      site->next=List;
      site->link[0]=max_link[0];

      List=site;

#ifdef MESSAGE
      //      printf(" linked with start/stop at position %ld\n",max_link[0]->poz);
#endif
    }

    // ** stop in reversed

    if(StopCod(-i) && i>=3) { // found a stop
      if(Stop[i].yn==-1) {
	//      printf("5\n");fflush(stdout);
	
#ifdef MESSAGE
	printf("DPr: Add stop at %ld %c%c%c\n",i-2,Data[i-2],Data[i-1],Data[i]);fflush(stdout);
#endif

	span=List;

	// here compute the downstream value;
	double downstreamVal=WEIGHT_Interg;
	if(i-2-NONCOD_LEN>=1 && is_stop && is_start) { 
	  double logscore=RD[i-2-NONCOD_LEN]-RD[i-2];
	  double logdiff=logscore+C[7][i-2]-C[7][i-2-NONCOD_LEN];
	  double cod[3];
	  compute_codings(i-3,i-2-NONCOD_LEN,cod);
	  float max1=1;
	  if(logscore<cod[0] || logscore<cod[1] || logscore<cod[2]) { max1=0;}
	  point.dimension = (float *) malloc (sizeof(float) * 3);
	  no_of_dimensions = 3;
	  point.dimension[0] = Stop[i].score;
	  point.dimension[1] = logdiff;
	  point.dimension[2] = max1;
	  classify(&point, downroots);
	  downstreamVal=point.prob[0];
	  free(point.dimension); 
	}


	max_score[0]=log(0.5*(WEIGHT_Interg+downstreamVal)/2)*i;
	max_link[0]=NULL;
	go_back=0;
	while(span!=NULL && span->poz-go_back>0) {
	  if(i-span->poz >= intergene) // prima modificare
	    

	    //	printf("span->poz=%ld, span->type=%d, span->score[0]=%f go_back=%ld\n",span->poz,(int)span->type,span->score[0],go_back);

	    if(span->type==3 || span->type==5) {
	      if(span->end>go_back) go_back=span->end;
	      temp_score=log(0.5*(downstreamVal+span->score[1])/2)*(i-span->poz+1)+span->score[0]; // give a small probability to intergenic region
	      if(temp_score > max_score[0]) {
		max_score[0]=temp_score;
		max_link[0]=span;
	      }
	    }
	  span=span->interg;
	}

	site=(Splice_Site *)malloc(sizeof(Splice_Site));  
	//myalloc+=sizeof(Splice_Site);
	//printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	if (site == NULL) {
	  fprintf(stderr,"Memory allocation for site failure.\n"); 
	  abort();
	}

	site->poz=i;
	site->type=8; 
	for(j=0;j<3;j++) {
	  site->score[j]=-HUGE_VAL; 
	  site->bad_link[j]=0;
	}
	site->score[0]=max_score[0]; // a doua modificare: se reseteaza scorul la o gena noua
	site->next=List;
	site->interg=lastinterg;
	site->link[0]=max_link[0];
	List=site;
      }
      j=i%3;
      laststopr[j]=i; 
    }


    // ** gt in reversed

    if(IF[i].yn==-1) { // found a gt = ag's equivalent in forward

      //      printf("6\n");fflush(stdout);

#ifdef MESSAGE
      printf("DPr: ##gt at %ld %c%c\n",i,Data[i],Data[i-1]);fflush(stdout);
#endif


      span=List;
      while(span!=NULL && i-span->poz<3) span=span->next;
      max_score[0]=-HUGE_VAL;
      max_score[1]=-HUGE_VAL;
      max_score[2]=-HUGE_VAL;
      
      while(span!=NULL && i-span->poz < MAX_INTRON_LENGTH) {
	if(span->type==6) 
	  for(temp_frame=0;temp_frame<3;temp_frame++)
	    if(span->score[temp_frame]>-HUGE_VAL) {
	      length=i-span->poz+2;
	      if(length>=5) {
		if(length<=6) hexamerval=0;
		else hexamerval=ifhexamer(i,span->poz-1);

		point.dimension = (float *) malloc (sizeof(float) * 4);
		no_of_dimensions = 4;
		point.dimension[0] = IF[i].score;
		point.dimension[1] = IF[span->poz].score;
		point.dimension[2] = length;
		point.dimension[3] = hexamerval;
		classify(&point, intronroots);
		bad_link=0;
		temp_score = enhance_score(point.prob[0],IF[span->poz].score,IF[i].score);
		if(temp_score<0.2) { // try different ranges here
		  bad_link=1;
		}
		free(point.dimension);  
		
		if(temp_score > DT_In_Min) {
		  //		  printf("Intron %ld %ld: %f Length=%d Hex=%f\n",i,span->poz-1,temp_score,length,hexamerval);
		  //if(length>150) length=150;
			  
#ifdef MESSAGE
		  printf("  intron %ld %ld : score %f total_score %f\n",i,span->poz-1,temp_score,log(temp_score*WEIGHT_In)*length+span->score[temp_frame]);
#endif
	    
		  temp_score=log(temp_score*WEIGHT_In)*length+span->score[temp_frame];
		  //temp_score=temp_score*0.9+span->score[temp_frame];

		  if(temp_score > max_score[temp_frame]) {
		    max_score[temp_frame]=temp_score;
		    max_link[temp_frame]=span;
		    if(bad_link) Blink[temp_frame]=bad_link+span->bad_link[temp_frame];
		    else Blink[temp_frame]=0;
		  }
		}
	      }
	    }
      	span=span->next;
      }

      if(max_score[0]!= -HUGE_VAL || max_score[1] != -HUGE_VAL || max_score[2] != -HUGE_VAL) {

#ifdef MESSAGE
	for(j=0;j<3;j++)
	  if(max_score[j]!= -HUGE_VAL) printf("  linked to ag at %ld\n",max_link[j]->poz);
	fflush(stdout);
#endif

	site=(Splice_Site *)malloc(sizeof(Splice_Site)); 
	//myalloc+=sizeof(Splice_Site);
	//printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	if (site == NULL) {
	  fprintf(stderr,"Memory allocation for site failure.\n"); 
	  abort();
	}

	site->poz=i;
	site->type=7; 
	for(j=0;j<3;j++) { 
	  site->score[j]=max_score[j];
	  site->bad_link[j]=Blink[j];
	}
	site->interg=lastinterg;
	for(j=0;j<3;j++) 
	  if(max_score[j]!=-HUGE_VAL) site->link[j]=max_link[j];
	  else site->link[j]=NULL;
	site->next=List;
	List=site;
      }
    }

    // ** ag in reversed

    if(IF[i].yn==-2) { // found an ag = gt's equivalent in forward

     
#ifdef MESSAGE
	printf("DPr: ##found ag at %ld %c%c\n",i,Data[i],Data[i-1]);fflush(stdout);
#endif


      span=List;
      max_score[0]=-HUGE_VAL;
      max_score[1]=-HUGE_VAL;
      max_score[2]=-HUGE_VAL;
      

      while(span!=NULL && i-span->poz < MAX_EXON_LENGTH) {
	
	if(span->type==8) {
	  length=i-2-span->poz+3;
	  temp_frame=length%3;
	  if(in_frame(-(i-1-temp_frame),span,laststopr)) {
	    if(length>=3) {
	      if(length<=6) hexamerval=0;
	      else hexamerval=ifhexamer(i-2,span->poz-2);

	      point.dimension = (float *) malloc (sizeof(float) * 3);
	      no_of_dimensions = 3;
	      point.dimension[0] = IF[i].score;
	      point.dimension[1] = length;
	      point.dimension[2] = hexamerval;
	      classify(&point, lastroots);
	      bad_link=0;
	      temp_score = enhance_score(point.prob[0],IF[i].score,-100);
	      if(temp_score<0.2) { // try different ranges here
		bad_link=1;
	      }
	      free(point.dimension);
	      //free((float*)point.dimension +1);

	      if(length<=15 && IF[span->poz].score <Acc_Thr_Last) temp_score=0;

	      if(temp_score>DT_LEx_Min) {
#ifdef MESSAGE
		printf("  exon %ld %ld : score %f total_score %f\n",i-2,span->poz+1,temp_score,log(temp_score*WEIGHT_LEx)*length+span->score[0]);
#endif

		temp_score=log(temp_score*WEIGHT_LEx)*length+span->score[0];
		//temp_score=temp_score+span->score[0];

		if(temp_score > max_score[temp_frame]) {
		  max_score[temp_frame]=temp_score;
		  max_link[temp_frame]=span;
		  if(bad_link) Blink[temp_frame]=bad_link+span->bad_link[0];else Blink[temp_frame]=0;
		}
	      }
	    }
	  }
	}

	  
	if(span->type==7) 
	  for(temp_frame=0;temp_frame<3;temp_frame++) {
	    length=i-2-span->poz;
	    if(length>=MIN_EXON_LENGTH) {
	      j=(3-(length-temp_frame)%3)%3;
	      if(span->score[j]>-HUGE_VAL && in_frame(-(i-1-temp_frame),span,laststopr)) {

		if(length<=6) hexamerval=0;
		else hexamerval=ifhexamer(i-2,span->poz+1);

		point.dimension = (float *) malloc (sizeof(float) * 4);
		no_of_dimensions = 4;
		point.dimension[0] = IF[i].score;
		point.dimension[1] = IF[span->poz].score;
		point.dimension[2] = length;
		point.dimension[3] = hexamerval;
		classify(&point, exonroots);
		bad_link=0;
		temp_score = enhance_score(point.prob[0],IF[i].score,IF[span->poz].score);
		if(temp_score<0.2) { // try different ranges here
		  bad_link=1;
		  if(span->bad_link[j]>=2) { bad_link=2; }
		}
		free(point.dimension);

		if(bad_link<2 && temp_score>DT_Ex_Min) {
		//if(temp_score>0.02) {
#ifdef MESSAGE
		  printf("  exon %ld %ld : score %f total_score %f\n",i-2,span->poz+1,temp_score,log(temp_score*WEIGHT_Ex)*length+span->score[j]);
#endif
	      
		  temp_score=log(temp_score*WEIGHT_Ex)*length+span->score[j];
		  //temp_score=temp_score+span->score[j];

		  if(temp_score > max_score[temp_frame]) {
		    max_score[temp_frame]=temp_score;
		    max_link[temp_frame]=span;
		    if(bad_link) Blink[temp_frame]=bad_link+span->bad_link[j];else Blink[temp_frame]=0;
		  }
		}
	      }
	    }
	  }
	span=span->next;
      }
      
      if(max_score[0]!= -HUGE_VAL || max_score[1] != -HUGE_VAL || max_score[2] != -HUGE_VAL) {
	
#ifdef MESSAGE
	for(j=0;j<3;j++) 
	  if(max_score[j]!= -HUGE_VAL) printf("frame %d:  linked with gt/stop at %ld\n",j,max_link[j]->poz);
	fflush(stdout);
#endif


	  site=(Splice_Site *)malloc(sizeof(Splice_Site));  
	  //myalloc+=sizeof(Splice_Site);
	  //printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	  if (site == NULL) {
	    fprintf(stderr,"Memory allocation for site failure.\n"); 
	    abort();
	  }
	  site->poz=i;
	  site->type=6; 
	  for(j=0;j<3;j++) {
	    site->score[j]=max_score[j];
	    site->bad_link[temp_frame]=Blink[temp_frame];
	  }
	  site->interg=lastinterg;
	  for(j=0;j<3;j++) 
	    if(max_score[j]!=-HUGE_VAL) site->link[j]=max_link[j];
	    else site->link[j]=NULL;
	  site->next=List;
	  List=site;
      }
    }


    // ** start in reversed

    if(IF[i].yn==-3) {

      //      printf("8\n");fflush(stdout);

#ifdef MESSAGE
      printf("DPr: found start at %ld %c%c%c\n",i,Data[i-2],Data[i-1],Data[i]); fflush(stdout);
#endif

      

      span=List;
      max_score[0]=-HUGE_VAL;


      while(span!=NULL && i-span->poz < MAX_EXON_LENGTH) {

	if(span->type==8) { // span is a stop
	  length=i-span->poz+3;
	  if(length >= Min_Gene_Len && length %3 == 0 && in_frame(-(i+1),span,laststopr)) {
	    hexamerval=ifhexamer(i,span->poz+1);
	    point.dimension = (float *) malloc (sizeof(float) * 2);
	    no_of_dimensions = 2;
	    point.dimension[0] = length;
	    point.dimension[1] = hexamerval;
	    classify(&point, snglroots);
	    temp_score = point.prob[0];
	    free(point.dimension);

	    if(temp_score>DT_SEx_Min) {
#ifdef MESSAGE
	      printf("  exon %ld %ld : score %f total_score %f\n",span->poz,i,temp_score,log(temp_score*WEIGHT_SEx)*length+span->score[0]);
#endif

	      temp_score=log(temp_score*WEIGHT_SEx)*length+span->score[0];
	      //temp_score=temp_score+span->score[0];

	      if(temp_score > max_score[0]) {
		max_score[0]=temp_score;
		max_link[0]=span;
	      }
	    }
	  }
	}

	if(span->type==7) {
	  length=i-span->poz;
	  temp_frame=(3-length%3)%3;
	  if(span->score[temp_frame]>-HUGE_VAL && length >=5 && in_frame(-(i+1),span,laststopr)) {
	    if(length<=6) hexamerval=0;
	    else hexamerval=ifhexamer(i,span->poz+1);
		
	    point.dimension = (float *) malloc (sizeof(float) * 3);
	    no_of_dimensions = 3;
	    point.dimension[0] = IF[span->poz].score;
	    point.dimension[1] = length;
	    point.dimension[2] = hexamerval;
	    classify(&point, startroots);
	    bad_link=0;
	    temp_score = enhance_score(point.prob[0],-100,IF[span->poz].score);
	    if(temp_score<0.2) { // try different ranges here
	      bad_link=1;
	      if(span->bad_link[temp_frame]>=2) { bad_link=2; }
	    }
	    free(point.dimension);
	    if(length<=15 && IF[span->poz].score<Don_Thr_First ) temp_score=0;
	    if(bad_link<2 && temp_score>DT_FEx_Min) {
	    //if(temp_score>0.02) {
#ifdef MESSAGE
	      printf("  exon %ld %ld : score %f total_score %f\n",i,span->poz+1,temp_score,log(temp_score*WEIGHT_FEx)*length+span->score[temp_frame]);
#endif
		
	      temp_score=log(temp_score*WEIGHT_FEx)*length+span->score[temp_frame];
	      //temp_score=temp_score+span->score[temp_frame];
	      
	      if(temp_score>max_score[0]) {
		max_score[0]=temp_score;
		max_link[0]=span;
	      }
	    }
	  }
	}
	span=span->next;
      }

      length=0; go_back=i;
#ifdef STOP
      strcpy(gene,"");
#endif
      if(max_score[0]!=-HUGE_VAL) {


	Cap=NULL;
	span=max_link[0];
	big_score=0;nex=0;
	while(span->type !=8) {
	  node=(intron *)malloc(sizeof(intron));  
	  //myalloc+=sizeof(intron);
	  //printf("Alloc+=%d: %ld\n",sizeof(intron),myalloc);
	  if (node == NULL) {
	    fprintf(stderr,"Memory allocation for node failure.\n"); 
	    abort();
	  }

	  if(go_back==i) {
	    if(go_back-span->poz<=6) hexamerval=0;
	    else hexamerval=ifhexamer(go_back,span->poz);

	    point.dimension = (float *) malloc (sizeof(float) * 3);
	    no_of_dimensions = 3;
	    point.dimension[0] = IF[span->poz].score;
	    point.dimension[1] = go_back-span->poz;
	    point.dimension[2] = hexamerval;
	    classify(&point, startroots);
	    big_score += point.prob[0];
	    free(point.dimension);
	    nex++;
	  }
	  else {
	    hexamerval=ifhexamer(go_back,span->poz+1);
	    point.dimension = (float *) malloc (sizeof(float) * 4);
	    no_of_dimensions = 4;
	    point.dimension[0] = IF[go_back+2].score;
	    point.dimension[1] = IF[span->poz].score;
	    point.dimension[2] = go_back-span->poz;
	    point.dimension[3] = hexamerval;
	    classify(&point, exonroots);
	    big_score += point.prob[0];
	    nex++;
	    free(point.dimension);
	  }
	  

	  length+=go_back-span->poz;
	  temp_frame=(3-length%3)%3;

#ifdef STOP
	  strncpy(exon,Data+span->poz+1,go_back-span->poz);
	  exon[go_back-span->poz]='\0';
	  strcat(exon,gene);
	  strcpy(gene,"");
	  strcat(gene,exon);
# endif

	  node->gt=span->poz;
	  span=span->link[temp_frame];
	  node->ag=span->poz;
	  go_back=span->poz-2;
	  node->leg=Cap;
	  Cap=node;
	  //	  prev_node->leg=node;  
	  
	  span=span->link[temp_frame];
	  // prev_node=node;
	}
#ifdef MESSAGE
	//printf("Try score_gene_r(%d,%d,%d)\n",i%3,i,span->poz);
	printf("Try score_gene_r big_score=%f\n",big_score);
#endif

	stopvect=1;
#ifdef STOP

	strncpy(exon,Data+span->poz-2,go_back-span->poz+3);
	exon[go_back-span->poz+3]='\0';
	strcat(exon,gene);
	strcpy(gene,"a");
	strcat(gene,exon);
	length=strlen(gene);


	Reverse_Complement(gene,length-1);

	if(length>=Min_Gene_Len) {
	
	  strcpy(vect1,gene+length-Min_Gene_Len-1);
	  score=String_Score_S1(vect1,Min_Gene_Len-3);

	  if(span->poz-2-Min_Gene_Len>=1) {
	    strncpy(vect2,Data+span->poz-2-Min_Gene_Len-2,Min_Gene_Len+1);
	    vect2[Min_Gene_Len+1]='\0';
	    score-=String_Score_SMax(vect2,Min_Gene_Len);
	    if(score<0) stopvect=0;
	  }
	}


# endif


	if(go_back==i) big_score=1;
	else {
	  if(go_back-span->poz+1<=6) hexamerval=0;
	  else hexamerval=ifhexamer(go_back,span->poz);

	  point.dimension = (float *) malloc (sizeof(float) * 3);
	  no_of_dimensions = 3;
	  point.dimension[0] = IF[go_back+2].score;
	  point.dimension[1] = go_back-span->poz+1;
	  point.dimension[2] = hexamerval;
	  classify(&point, lastroots);
	  big_score += point.prob[0];
	  free(point.dimension);
	  nex++;
	  big_score=big_score/nex;
	}    

	if(stopvect && big_score>0.4 && Score_Gene_R(i%3,i,span->poz,Cap)) {
#ifdef MESSAGE
	  printf("  linked with gt/stop at %ld\n",max_link[0]->poz);fflush(stdout);
#endif

	  site=(Splice_Site *)malloc(sizeof(Splice_Site)); 
	  //myalloc+=sizeof(Splice_Site);
	  //printf("Alloc+=%d: %ld\n",sizeof(Splice_Site),myalloc);
	  if (site == NULL) {
	    fprintf(stderr,"Memory allocation for site failure.\n"); 
	    abort();
	  }

	  site->poz=i;
	  site->type=5; 
	  site->end=span->poz-2;
	  for(j=1;j<3;j++) {
	    site->score[j]=-HUGE_VAL;
	    site->bad_link[j]=0;
	  }
	  site->score[0]=max_score[0];

	  // *** here compute the upstream value
	  if(i>Data_Len-NONCOD_LEN || !is_stop || !is_start) { 
	    site->score[1]=WEIGHT_Interg; 
	    // or maybe	site->score[1]=DT_Up; ? 
	  }
	  else {
	    double logscore=RU[i]-RU[i+NONCOD_LEN];
	    double logdiff=logscore+C[7][i+NONCOD_LEN]-C[7][i];
	    double cod[3];
	    compute_codings(i+NONCOD_LEN,i+1,cod);
	    float max1=1;
	    if(logscore<cod[0] || logscore<cod[1] || logscore<cod[2]) { max1=0;}
	    point.dimension = (float *) malloc (sizeof(float) * 3);
	    no_of_dimensions = 3;
	    point.dimension[0] = IF[i].score;
	    point.dimension[1] = logdiff;
	    point.dimension[2] = max1;
	    classify(&point, uproots);
	    site->score[1]=point.prob[0];
	    free(point.dimension); 
	  }

	  site->interg=lastinterg;
	  lastinterg=site;
	  site->next=List;
	  for(j=0;j<3;j++) site->link[j]=NULL;
	  site->link[0]=max_link[0];
	  List=site;

#ifdef GENESCORE
	  printf("Gene_score: %f\n",max_score); fflush(stdout);
#endif
	  
	}
	

	// free Cap here
	prev_node=Cap;
	if(prev_node !=NULL) node=prev_node->leg;
	while(prev_node !=NULL) {
	  free(prev_node);
	  //myalloc-=sizeof(intron);
	  //printf("Alloc-=%d: %ld\n",sizeof(intron),myalloc);
	  prev_node=node;
	  if(node !=NULL) node=prev_node->leg;
	}

      }
#ifdef MESSAGE
      printf("\n"); fflush(stdout);
#endif

    }

  }

  // print the results
  
  
  max_score[0]=-HUGE_VAL;
  span=List;
  go_back=0;
  /*  if((tempfile=fopen("glimm.tempfile","w"))==NULL) {
      fprintf(stderr,"Could not create temporary file!\n");
      exit(-1);
      }
      fprintf(tempfile,"\n");*/
  
  geneCap=NULL;

  max_link[0]=NULL;

  while(span!=NULL && span->poz-go_back>0) {
    if(span->type==3 || span->type == 5) {
      if(span->end>go_back) go_back=span->end;
      if(span->score[0]+log(0.5*(WEIGHT_Interg+span->score[1])/2)*(Data_Len-span->poz+1)>max_score[0]) {
	max_score[0]=span->score[0]+log(0.5*(WEIGHT_Interg+span->score[1])/2)*(Data_Len-span->poz+1);
	max_link[0]=span;
      }
    }
    span=span->next;
  }

  span=max_link[0];

  k=ngene+1;


  // printf("Print results:\n");

  while(span!=NULL) {

    switch(span->type) {
    case 1: // donor forward
      exonend[nex++]=span->poz-1;
      go_back=span->poz-1;
      break;
    case 2: // acceptor forward
      exonend[nex++]=span->poz+2;
      length+=go_back-span->poz-1;
      break;
    case 3: // stop forward
      nex=0;
      exonend[nex++]=span->poz;
      length=0; go_back=span->poz;
      break;
    case 4: // start forward
      exonend[nex++]=span->poz;
      ngene++;
      //      fprintf(tempfile,"\n");
      geneNode=(Gene_Info *)malloc(sizeof(Gene_Info)); 
      strcpy(geneNode->line,"\n");
      geneNode->link=geneCap;
      geneCap=geneNode;
      if(nex==2) {
	exlength=(int)(exonend[0]-exonend[1]+1);
	//	fprintf(tempfile,"   1  +  Single   %10ld %10ld  %7d\n",exonend[1],exonend[0],exlength);
	sprintf(line,"   1  +  Single   %10ld %10ld  %7d\n",exonend[1]+500000*offset,exonend[0]+500000*offset,exlength); 
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info)); 
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
      }
      else {
	exlength=(int)(exonend[0]-exonend[1]+1);
	//fprintf(tempfile,"%4d  +  Terminal %10ld %10ld  %7d\n",(nex+1)/2,exonend[1],exonend[0],exlength);
	sprintf(line,"%4d  +  Terminal %10ld %10ld  %7d\n",(nex+1)/2,exonend[1]+500000*offset,exonend[0]+500000*offset,exlength); 
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info));
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
	for(j=2;j<nex-2;j+=2) {
	  exlength=(int)(exonend[j]-exonend[j+1]+1);
	  //  fprintf(tempfile,"%4d  +  Internal %10ld %10ld  %7d\n",(nex-j+1)/2,exonend[j+1],exonend[j],exlength);
	  sprintf(line,"%4d  +  Internal %10ld %10ld  %7d\n",(nex-j+1)/2,exonend[j+1]+500000*offset,exonend[j]+500000*offset,exlength); 
	  geneNode=(Gene_Info *)malloc(sizeof(Gene_Info)); 
	  strcpy(geneNode->line,line);
	  geneNode->link=geneCap;
	  geneCap=geneNode;
	}
	exlength=(int)(exonend[nex-2]-exonend[nex-1]+1);
	//fprintf(tempfile,"   1  +  Initial  %10ld %10ld  %7d\n",exonend[nex-1],exonend[nex-2],exlength);
	sprintf(line,"   1  +  Initial  %10ld %10ld  %7d\n",exonend[nex-1]+500000*offset,exonend[nex-2]+500000*offset,exlength); 
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info));  
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
      }
      length+=go_back-span->poz+1;
      break;
    case 5: // start reverse
      nex=0;
      exonend[nex++]=span->poz;
      go_back=span->poz;length=0;
      break;
    case 6: // acceptor reverse
      exonend[nex++]=span->poz-2;
      go_back=span->poz-2;
      break;
    case 7: // donor reverse
      exonend[nex++]=span->poz+1;
      length+=go_back-span->poz;
      break;
    case 8: // stop reverse
      exonend[nex++]=span->poz-2;
      ngene++;
      //fprintf(tempfile,"\n");
      geneNode=(Gene_Info *)malloc(sizeof(Gene_Info)); 
      strcpy(geneNode->line,"\n");
      geneNode->link=geneCap;
      geneCap=geneNode;
      if(nex==2) {
	exlength=(int)(exonend[0]-exonend[1]+1);
	//fprintf(tempfile,"   1  -  Single   %10ld %10ld  %7d\n",exonend[1],exonend[0],exlength);
	sprintf(line,"   1  -  Single   %10ld %10ld  %7d\n",exonend[1]+500000*offset,exonend[0]+500000*offset,exlength);
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info));
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
      }
      else {
	exlength=(int)(exonend[0]-exonend[1]+1);
	//fprintf(tempfile,"   1  -  Initial  %10ld %10ld  %7d\n",exonend[1],exonend[0],exlength);
	sprintf(line,"   1  -  Initial  %10ld %10ld  %7d\n",exonend[1]+500000*offset,exonend[0]+500000*offset,exlength); 
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info));
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
	for(j=2;j<nex-2;j+=2) {
	  exlength=(int)(exonend[j]-exonend[j+1]+1);
	  //fprintf(tempfile,"%4d  -  Internal %10ld %10ld  %7d\n",(j+2)/2,exonend[j+1],exonend[j],exlength);
	  sprintf(line,"%4d  -  Internal %10ld %10ld  %7d\n",(j+2)/2,exonend[j+1]+500000*offset,exonend[j]+500000*offset,exlength);
	  geneNode=(Gene_Info *)malloc(sizeof(Gene_Info)); 
	  strcpy(geneNode->line,line);
	  geneNode->link=geneCap;
	  geneCap=geneNode;
	}
	exlength=(int)(exonend[nex-2]-exonend[nex-1]+1);
	//fprintf(tempfile,"%4d  -  Terminal %10ld %10ld  %7d\n",(nex+1)/2,exonend[nex-1],exonend[nex-2],exlength);
	sprintf(line,"%4d  -  Terminal %10ld %10ld  %7d\n",(nex+1)/2,exonend[nex-1]+500000*offset,exonend[nex-2]+500000*offset,exlength);
	geneNode=(Gene_Info *)malloc(sizeof(Gene_Info));  
	strcpy(geneNode->line,line);
	geneNode->link=geneCap;
	geneCap=geneNode;
      }
      length+=go_back-span->poz+3;
      break;
    }

    temp_frame=(3-length%3)%3;
    span=span->link[temp_frame];
  }

  //  fclose(tempfile);
 
  j=k;
 
  while(geneCap!=NULL) {
    if(geneCap->line[0]=='\n') { 
      j++;
      printf("\n");
    }
    else printf("%4d %s",j,geneCap->line);
    geneNode=geneCap;
    geneCap=geneNode->link;

    free(geneNode);
  }
  
  // free List here
  site=List;
  if(site) span=site->next;
  while(site != NULL) {
    free(site);
    site=span;
    if(span!=NULL) span=site->next;
  }

}

int in_frame(long int stop, Splice_Site * span, long int *laststop) {
  long int j;
  int i,frame;
  char buf[3];
  
  if(stop>0) {
    if(span->type==4) { 
      j=span->poz;
      frame=0;
    }
    if(span->type==2) {
      j=span->poz+2;
      frame=(3-(stop-1-span->poz-1)%3)%3;
    }

    switch(frame) {
    case 1: buf[0]=Data[(span->link[frame])->poz-1];
      buf[1]=Data[j++];buf[2]=Data[j++]; 
      if(If_Stop(buf)) return(0);
      break;
    case 2: buf[0]=Data[(span->link[frame])->poz-2];buf[1]=Data[(span->link[frame])->poz-1];
      buf[2]=Data[j++]; 
      if(If_Stop(buf)) return(0);
      break;
    }
    
    //### see if there is stop in the reading frame
    j+=2;
    i=j%3;
    if(j<=laststop[i]) return(0);
  }
  else {
    stop=-stop;
    if(span->type==8) {
      j=span->poz+1;
      frame=0;
    }
    if(span->type==7) {
      j=span->poz+1;
      frame=(3-(stop-1-span->poz)%3)%3;
    }

    switch(frame) {
    case 1: buf[2]=inverse(Data[(span->link[frame])->poz-2]);
      buf[1]=inverse(Data[j++]);buf[0]=inverse(Data[j++]); 
      if(If_Stop(buf)) return(0);
      break;
    case 2: buf[2]=inverse(Data[(span->link[frame])->poz-3]);
      buf[1]=inverse(Data[(span->link[frame])->poz-2]);
      buf[0]=inverse(Data[j++]); 
      if(If_Stop(buf)) return(0);
      break;
    }
    j+=2;
    i=j%3;
    if(j<=laststop[i]) return(0);
  }

  return(1);
}

char inverse(char ch)
{
  switch(ch) {
  case 'a': return('t');
  case 'A': return('T');
  case 'c': return('g');
  case 'C': return('G');
  case 'g': return('c');
  case 'G': return('C');
  case 't': return('a');
  case 'T': return('A');
  default: return('c');
  }
}

 
int StopCod(long int i)
{
  if(i>0) {
    if(Data[i-2]=='t' && Data[i-1]=='a' && Data[i]=='a') return(1);
    if(Data[i-2]=='t' && Data[i-1]=='g' && Data[i]=='a') return(1);
    if(Data[i-2]=='t' && Data[i-1]=='a' && Data[i]=='g') return(1);
    return(0);
  }
  else {
    i=-i;
    if(Data[i]=='a' && Data[i-1]=='t' && Data[i-2]=='t') return(1);
    if(Data[i]=='a' && Data[i-1]=='c' && Data[i-2]=='t') return(1);
    if(Data[i]=='a' && Data[i-1]=='t' && Data[i-2]=='c') return(1);
    return(0);
  }
}

int Score_Gene_R(int Frame,long int start,long int i,intron *Cap)
{
  int Score[7];
  long int Gene_Len,Start,j,end,Copy_Len,Rev_Start;
  int In_Frame_Score;
  char *Copy;
  intron *temp;

  // here create the copy to eliminate the introns
  Copy=(char *) malloc((Data_Len+2)*sizeof(char));
  if (Copy == NULL) {
    fprintf(stderr,"Memory allocation for copy failure.\n");
    abort();
  }
 
  j=1;
  Copy[1]='\0';
  Copy_Len=Data_Len;

  temp=Cap;
  end=i-2;

  while(temp!=NULL) {

    // eliminate the portion between gt and ag;
    strncat(Copy+1,Data+j,temp->ag-j-1);
    j=temp->gt+1;
    end=temp->ag-2;
    temp=temp->leg;
  }

  strcat(Copy+1,Data+j);
  
  Copy_Len=strlen(Copy+1);

  Rev_Start=start;

  if  (Rev_Start != 0) {
	 
    Gene_Len = Rev_Start - i -(Data_Len-Copy_Len) ;
    Start=Rev_Start;

    Gene_Len = Start - i -(Data_Len-Copy_Len);

    if  (Gene_Len >= Min_Gene_Len && end<=Start) {
      
#ifdef GENESCORE

      printf("Score gene:\n");fflush(stdout);
      printgene(Start,i+1,-Frame-1,Cap);
#endif

      Score_String (Copy, Start-(Data_Len-Copy_Len),
		    i  + 1, Copy_Len,
		    Ch_Ct, Score);
      In_Frame_Score = Score [0];

      Permute (Score, - Frame - 1);


#ifdef GENESCORE

      printf("Score is %d\n",In_Frame_Score);fflush(stdout);

#endif
      
      if  (In_Frame_Score >= Threshold_Score) {
	free(Copy);
	return(1);
      }
    }
  }
  free(Copy);  
  return(0);
}



int Score_Gene_F(int Frame,long int start, long int i,intron *Cap) 
{
  int Score[7];
  long int Gene_Len,Start,j,end,Copy_Len,For_Start;
  int In_Frame_Score;
  char *Copy;
  intron *temp;

  // here create the copy to eliminate the introns
  Copy=(char *) malloc((Data_Len+2)*sizeof(char));
  if (Copy == NULL) {
    fprintf(stderr,"Memory allocation for copy failure.\n"); 
    abort();
  }

  j=1;
  Copy[1]='\0';
  end=i-2;

  temp=Cap;

  while(temp!=NULL) {
    // eliminate the portion between gt and ag;

    

    if(j==1) end=Cap->gt;

    strncat(Copy+1,Data+j,temp->gt-j);


    j=temp->ag+2;
    temp=temp->leg;

  }

  strcat(Copy+1,Data+j);
	   
  Copy_Len=strlen(Copy+1);


  // found start codon with no stop codons
  For_Start=start;


  if  (For_Start != 0){  // If have seen a start codon in this orf
    
    For_Start-=2;

    // see if the gene is a good one


	 
    Gene_Len = 1 + i - 3 - For_Start-(Data_Len-Copy_Len);
    
    // modify choose start in order to work for eukariots 
    Start=For_Start;

    // Allows use of a fancy start codon selection alg
    Gene_Len = 1 + i - Start-(Data_Len-Copy_Len);

    //printf("Gene_Len=%d\n",Gene_Len);

    if  (Gene_Len >= Min_Gene_Len && Start<end) {
      
      // Score this orf in all 6 frames plus the random model

      Score_String (Copy, Start, i - 3-(Data_Len-Copy_Len),
		    Copy_Len, Ch_Ct,
		    Score);

	
      In_Frame_Score = Score [0]; //Frame is relative to beginning of orf


      // my variant
      
      Permute (Score, Frame + 1);

      // end my variant

      
#ifdef GENESCORE

      printf("Score gene:\n");
      printgene(i-3,Start,Frame+1,Cap);
      printf("Score is %d\n",In_Frame_Score);fflush(stdout);

#endif 
      
      if  (In_Frame_Score >= Threshold_Score) {
	free(Copy);
	return(1);
      }

    }
  }
  free(Copy);
  return(0);

}


void printgene(long int Hi,long int Lo,int Frame,intron *Cap)
{
  intron *temp;
 
  if  (Frame > 0) {
    printf ("%8ld ", Lo);fflush(stdout);

    temp=Cap;
    while(temp !=NULL) {
      printf("%8ld %8ld ",(temp->gt)-1,(temp->ag)+2);fflush(stdout);
      temp=temp->leg; 
    }

    printf("%8ld",Hi+3);fflush(stdout);
  }
  else {
    printf("%8ld ",Hi);fflush(stdout);
    temp=Cap; no=0;last=Hi;
    if(temp!=NULL) printreverse(temp);fflush(stdout);
    printf ("%8ld", Lo-3);fflush(stdout);
    
  }
     
  putchar ('\n');
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


void  Permute  (int S [], int F)

/* Permute the values in  S  which are in frame  F  to make them
*  consistent with scores from frame  +1 . */

  {
   int  Save;

   switch  (F)
     {
      case  1 :
        Save = S [4];
        S [4] = S [5];
        S [5] = Save;
        return;
      case  2 :
        Save = S [0];
        S [0] = S [2];
        S [2] = S [1];
        S [1] = Save;
        Save = S [3];
        S [3] = S [5];
        S [5] = Save;
        return;
      case  3 :
        Save = S [0];
        S [0] = S [1];
        S [1] = S [2];
        S [2] = Save;
        Save = S [3];
        S [3] = S [4];
        S [4] = Save;
        return;
      case  -1 :
        Save = S [0];
        S [0] = S [3];
        S [3] = Save;
        Save = S [1];
        S [1] = S [5];
        S [5] = S [2];
        S [2] = S [4];
        S [4] = Save;
        return;
      case  -2 :
        Save = S [0];
        S [0] = S [4];
        S [4] = S [2];
        S [2] = S [5];
        S [5] = Save;
        Save = S [1];
        S [1] = S [3];
        S [3] = Save;
        return;
      case  -3 :
        Save = S [0];
        S [0] = S [5];
        S [5] = S [1];
        S [1] = S [4];
        S [4] = Save;
        Save = S [2];
        S [2] = S [3];
        S [3] = Save;
        return;
     }
   return;
  }


void  Process_Options  (int argc, char * argv [])

//  Process command-line options and set corresponding global switches
//  and parameters.
//
//    -d dir Set the training files' directory
//    -g n   Set minimum gene length to n
//    -r     Don't use independent probability score column
//    +r     Use independent probability score column
//    -t n   Set threshold score for calling as gene to n.  If the in-frame
//           score >= n, then the region is given a number and considered
//           a potential gene.
//    -f     Don't use filtering of the splice sites
//    +f     Use filtering of the splice sites
//    -s     Don't use start site model
//    +s     Use start site model
//    -b     Don't use stop site model
//    +b     Use stop site model
//    -5 th  Use threshold th for the acceptor sites
//    -3 th  Use threshold th for the donor sites

  {
   char  * P;
   long int  W;
   double D;
   int  i;

   for  (i = 2;  i < argc;  i ++)
     {
      switch  (argv [i] [0])
        {
         case  '-' :
           switch  (argv [i] [1])
             {
	     case 'd' :   //  the training files' directory
	       strcpy(TRAIN_DIR,argv[++i]);strcat(TRAIN_DIR,"/");
	       //printf("train dir = %s\n",TRAIN_DIR);fflush(stdout);
	       break;
	     case  'g' :       // minimum gene length
	       errno = 0;
	       if  (argv [i] [2] != '\0')
		 P = argv [i] + 2;
	       else
		 P = argv [++ i];
	       W = strtol (P, NULL, 10);
	       if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum gene length %s\n", P);
	       else
		 Min_Gene_Len = W;
	       assert (Min_Gene_Len > 0);
	       break;
	     case  'r' :       // don't use random/independent score column
	       Use_Independent = FALSE;
	       break;
	     case 'f' :
	       Use_Filter=0;
	       break;
	     case 's' :
	       Use_Start_Site=0;
	       break;
	     case  't' :       // threshold score for calling as gene
	       errno = 0;
	       if  (argv [i] [2] != '\0')
		 P = argv [i] + 2;
	       else
		 P = argv [++ i];
	       W = strtol (P, NULL, 10);
	       if  (errno == ERANGE)
		 fprintf (stderr, "ERROR:  Bad threshold score for gene %s\n", P);
	       else
		 Threshold_Score = W;
	       assert (Threshold_Score > 0 && Threshold_Score < 100);
	       break;
	     case '5':
	       errno = 0;
	       if  (argv [i] [2] != '\0')
		 P = argv [i] + 2;
	       else
		 P = argv [++ i];
	       D = strtod (P, NULL);
	       if  (errno == ERANGE)
		 fprintf (stderr, "ERROR:  Bad acceptor threshold score %s\n", P);
	       else
		 Acc_Thr = D;
	       break;
	     case '3':
	       errno = 0;
	       if  (argv [i] [2] != '\0')
		 P = argv [i] + 2;
	       else
		 P = argv [++ i];
	       D = strtod (P, NULL);
	       if  (errno == ERANGE)
		 fprintf (stderr, "ERROR:  Bad acceptor threshold score %s\n", P);
	       else
		 Don_Thr = D;
	       break;
	     default :
	       fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
           break;
	case  '+' :
	  switch  (argv [i] [1])
	    {
	    case  'r' :       // use random/independent score column
	      Use_Independent = TRUE;
	      break;
	    case 'f' :
	      Use_Filter = TRUE;
	      break;
	    case 's' :
	      Use_Start_Site = TRUE;
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



void  Score_String  (char * Copy, long int Start, long int Stop, long int Copy_Len,
                     double Ch_Ct [ALPHABET_SIZE],
                    int Score [7])

//  Sets  Score  array to the scores in each frame for the string
//  Copy [Start .. Stop]  using the Markov model in  Delta
//  and the simple independent log probabilities in  Ch_Ct .
//  Each score is an integer 0 .. 99 and the sum of all scores is 99
//  (modulo round-off errors).  The model is the simple Markov equivalent
//  of the IMM.
//  Copy_Len  is the last position in  Copy  to compute wraparounds.
//  Score [6] gets the scaled score for the random model that uses  Ch_Ct .

  {
   long int  Len;


   if  (Start < Stop ) {
	      
     Len = 1 + Stop - Start;
     
     if  (2 + Len > Orf_Buffer_Len)
       {
	 Orf_Buffer_Len = Max (2 + Len, Orf_Buffer_Len + ORF_SIZE_INCR);
	 Orf_Buffer = (char *) Safe_realloc (Orf_Buffer, Orf_Buffer_Len);
       }
     Transfer (Orf_Buffer + 1, Start, Len,Copy,Copy_Len);
   }
   else {
     
     Len = 1 + Start - Stop;
     
     if  (2 + Len > Orf_Buffer_Len)
       {
	 Orf_Buffer_Len = Max (2 + Len, Orf_Buffer_Len + ORF_SIZE_INCR);
	 Orf_Buffer = (char *) Safe_realloc (Orf_Buffer, Orf_Buffer_Len);
       }
     Transfer (Orf_Buffer + 1, Start, - Len,Copy,Copy_Len);
   }

   Simple_Score (Orf_Buffer, Len, MODEL_LEN, Ch_Ct, Score);

   return;
  }



void  Simple_Score  (char X [], int T, int Model_Len,
                     double Ch_Ct [ALPHABET_SIZE],
                     int Score [])

/* Set  Score  to the probabilites of string  X [1 .. T]  being
*  generated in each of the 3 forward and 3 reverse reading frames
*  using simple nonhomogeneous Markov models in  Delta []  with
*  model length equal to  Model_Len . */

{
  double  Max, Min, Sum, S [7], W [7];
  double  Weak_Max, Weak_Min, Weak_Sum;
  int  i, Has_Stop [7], Offset;
  int Weak_Score;

  Find_Stop_Codons (X, T, Has_Stop);

  Max = - DBL_MAX;
  Min = DBL_MAX;

  if  (! Has_Stop [0]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [0],
		   MODEL [1], MODEL [2], S [0]);
    if  (S [0] > Max)
      Max = S [0];
    if  (S [0] < Min)
      Min = S [0];
  }
  if  (! Has_Stop [1]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [2],
		   MODEL [0], MODEL [1], S [1]);
    if  (S [1] > Max)
      Max = S [1];
    if  (S [1] < Min)
      Min = S [1];
  }
  if  (! Has_Stop [2]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [1],
		   MODEL [2], MODEL [0], S [2]);
    if  (S [2] > Max)
      Max = S [2];
    if  (S [2] < Min)
      Min = S [2];
  }

  Offset = T % 3;
  Reverse_Complement (X, T);
  if  (! Has_Stop [3 + Offset]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [0],
		   MODEL [1], MODEL [2], S [3 + Offset]);
    if  (S [3 + Offset] > Max)
      Max = S [3 + Offset];
    if  (S [3 + Offset] < Min)
      Min = S [3 + Offset];
  }
  Offset = (Offset + 1) % 3;
  if  (! Has_Stop [3 + Offset]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [1],
		   MODEL [2], MODEL [0], S [3 + Offset]);
    if  (S [3 + Offset] > Max)
      Max = S [3 + Offset];
    if  (S [3 + Offset] < Min)
      Min = S [3 + Offset];
  }
  Offset = (Offset + 1) % 3;
  if  (! Has_Stop [3 + Offset]) {
    Fast_Evaluate (X + 1, T, Model_Len, MODEL [2],
		   MODEL [0], MODEL [1], S [3 + Offset]);
    if  (S [3 + Offset] > Max)
      Max = S [3 + Offset];
    if  (S [3 + Offset] < Min)
      Min = S [3 + Offset];
  }

  Weak_Max = Max;
  Weak_Min = Min;
  
  Has_Stop [6] = ! Use_Independent;
  Has_Stop [6] = 0;

  if  (Use_Independent)  {
       
    Indep_Eval (X, T, Ch_Ct, S [6]);
    if  (S [6] > Max)
      Max = S [6];
    if  (S [6] < Min)
      Min = S [6];
  }

  assert (Max != - DBL_MAX && Min != DBL_MAX);

  if  (Min < Max + MAX_LOG_DIFF)
    Min = Max + MAX_LOG_DIFF;
  if  (Weak_Min < Weak_Max + MAX_LOG_DIFF)
    Weak_Min = Weak_Max + MAX_LOG_DIFF;

  Sum = 0.0;
  for  (i = 0;  i < 7;  i ++)
    if  (Has_Stop [i])
      W [i] = -1.0;
    else if  (S [i] >= Min)
      {
	W [i] = exp (S [i] - Min);
	Sum += W [i];
      }
    else
      W [i] = 0.0;
  assert (Sum > 0.0);

  for  (i = 0;  i < 7;  i ++)
    if  (Has_Stop [i])
      Score [i] = -1;
    else
      {
	Score [i] = int (100.0 * (W [i] / Sum));
	if  (Score [i] > 99)
	  Score [i] = 99;
      }

  if  (! Use_Independent)
  Weak_Score = Score [0];
  else
  {
    Weak_Sum = 0.0;
    for  (i = 0;  i < 6;  i ++)
      if  (! Has_Stop [i] && S [i] >= Weak_Min)
	{
	  W [i] = exp (S [i] - Weak_Min);
	  Weak_Sum += W [i];
	}
    assert (Weak_Sum > 0.0);
    Weak_Score = int (100.0 * (W [0] / Weak_Sum));
    if  (Weak_Score > 99)
      Weak_Score = 99;
  }

  if  (Score [6] < 0)
       Score [6] = 0;

   return;
}



void  Transfer  (char * S, long int Start, int Len,char * Copy, long int Copy_Len)

/* Transfer  |Len|  characters from  Copy [Start ...]  to
*  S  and add null terminator.  If  Len > 0 go in forward direction;
*  otherwise, go in reverse direction and use complements.
*  Allow for wraparound at  Copy_Len . */

  {
   long int  i, j;

   if  (Len > 0)
       {
        for  (i = 0;  i < Len;  i ++)
          {
           j = Start + i;
           if  (j > Copy_Len)
               j -= Copy_Len;
           else if  (j < 1)
               j += Copy_Len;
           S [i] = Filter (tolower (Copy [j]));
          }
        S [i] = '\0';
       }
     else
       {
        for  (i = 0;  i < - Len;  i ++)
          {
           j = Start - i;
           if  (j > Copy_Len)
               j -= Copy_Len;
           else if  (j < 1)
               j += Copy_Len;
           S [i] = Filter (tolower (Complement (Copy [j])));
          }
        S [i] = '\0';
       }

   return;
  }

float Perc(long int start,long int end, long int *atp)
{
  long int i,at;
  float ret;

  at=0;

  for(i=start;i<=end;i++) {
    switch(Data[i]) {
    case 'A':
    case 'a':
    case 'T':
    case 't':at++;break;
    }
  }
  
  ret=at*100/(end-start+1);

  (*atp)=at;

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

double ifhexamer(long int start, long int stop)
{ 
   int index,weight,val;                /* location of hexamer in HexamerVec */
   long int pos1, phase,i;
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
  struct tree_node
    **uptrees = NULL;       /* trees for classifying upstream regions */
  struct tree_node
    **downtrees = NULL;       /* trees for classifying downstream exons */

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
  uptrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);
  downtrees = (struct tree_node **)
    malloc (sizeof(struct tree_node *) * no_of_trees);

  /* the decision tree names must be in a file called "treenames" */
  sprintf(filename,"%streenames",TRAIN_DIR);
  if((treenames = fopen(filename, "r")) == NULL) {
    fprintf (stderr, "Decision tree names file can not be opened.");
    exit(-1);
  }
  for (i=0; i< (no_of_trees * 7); i++) {
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
    /* next 10 trees are for upstream region */
    else if (i < 60) {
      if ((uptrees[i-50] = read_tree(filename)) == NULL) {
	fprintf(stderr,"Mktree: Cannot read %s.\n",treefile);
	exit(-1);
      }
    }
    /* next 10 trees are for downstream region */
    else if (i < 70) {
      if ((downtrees[i-60] = read_tree(filename)) == NULL) {
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
  uproots = uptrees;
  downroots = downtrees;

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
	  /*printf("left cat= %d  left_count[%d]= %d no_points= %d\n",
	     cur_node->left_cat,
	     cur_node->left_cat,
	     cur_node->left_count[cur_node->left_cat],
	     cur_node->left_total);*/
	  
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
	  /*printf("right cat= %d  right_count[%d]= %d no_points= %d\n",
	       cur_node->right_cat,
	       cur_node->right_cat,
	       cur_node->right_count[cur_node->right_cat],
	       cur_node->right_total);*/
	    
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



int  If_Stop  (char * S)

/* Return  TRUE  iff  S  is a stop codon. */

  {
   return  (strncmp (S, "taa", 3) == 0
              || strncmp (S, "tag", 3) == 0
              || strncmp (S, "tga", 3) == 0);
  }


double enhance_score(double score,double accthr,double donthr)
{
  if(accthr>=30) score+=0.2;

  if(donthr>=30) score+=0.2;


  return(score);
}

void  Read_Probability_Model  (char * Param)

//  Read in the probability model indicated by  Param .

{
  FILE  * fp;
  
  fp = File_Open (Param, "r");   // maybe rb ?
  
  Read_Scoring_Model (fp);
  
  fclose (fp);
  
  return;
}

void  Read_NC_Model  ()

//  Read in the probability models of the non-coding regions.

{
  char File_Name[MAX_LINE];

  // upstream

  strcpy(File_Name,TRAIN_DIR);
  strcat(File_Name,"upstream.model");

  FILE  * fp;
  
  fp = File_Open (File_Name, "r");   // maybe rb ?
  
  UMODEL=Read_NonCoding_Model (fp);
  
  fclose (fp);

  // downstream

  strcpy(File_Name,TRAIN_DIR);
  strcat(File_Name,"downstream.model");

  fp = File_Open (File_Name, "r");   // maybe rb ?
  
  DMODEL=Read_NonCoding_Model (fp);
  
  fclose (fp);

  // introns

  /*strcpy(File_Name,TRAIN_DIR);
    strcat(File_Name,"introns.model");
    
    fp = File_Open (File_Name, "r");   // maybe rb ?
  
    IMODEL=Read_NonCoding_Model (fp);
  
  fclose (fp);*/
  
  return;
}

double String_Score_S0(char X [], int T) 
{
  double S[2];

  Fast_Evaluate (X + 1, T, MODEL_LEN, MODEL [0],MODEL [1], MODEL [2], S [0]);

  Indep_Eval (X, T, Ch_Ct, S [1]);

  return(S[0]-S[1]);
}

double String_Score_S1(char X [], int T) 
{
  double S;

  Fast_Evaluate (X + 1, T, MODEL_LEN, MODEL [0],MODEL [1], MODEL [2], S);

  return(S);
}

double String_Score_SMax(char X [], int T) 
{
  double S,SMax;

  SMax=- DBL_MAX;

  Fast_Evaluate (X + 1, T, MODEL_LEN, MODEL [0],MODEL [1], MODEL [2], S);
  if(S>=SMax) SMax=S;
  Fast_Evaluate (X + 1, T, MODEL_LEN, MODEL [2],MODEL [0], MODEL [1], S);
  if(S>=SMax) SMax=S;
  Fast_Evaluate (X + 1, T, MODEL_LEN, MODEL [1],MODEL [2], MODEL [0], S);
  if(S>=SMax) SMax=S;

  return(SMax);
}

double scalarprod(double *r1,double *r2)
{
  double val=0,sum1=0,sum2=0,sum=0;
  int i;

  for(i=0;i<64;i++) {
    sum+=r1[i]*r2[i];
    sum1+=r1[i]*r1[i];
    sum2+=r2[i]*r2[i];
  }

  val=sum/sqrt(sum1*sum2);

  return(val);
}

void compute_vect(char *seq,double *r,int frame)
{
  int i,j;
  int len;
  int val;

  for(i=0;i<64;i++) {
    r[i]=0;
  }

  len=strlen(seq);

  i=0;

  while(i+3<=len) {
    val=0;
    for(j=i+2;j>=i;j--) 
      switch (seq[j]) {
      case 'a': 
      case 'A': break;
      case 'c':
      case 'C': val+=1*(int)pow(4,2-j+i);break;
      case 'g':
      case 'G': val+=2*(int)pow(4,2-j+i);break;
      case 't':
      case 'T': val+=3*(int)pow(4,2-j+i);break;
      }
    r[val]++;
    if(frame) i+=1;
    else i+=3;
  }

  for(i=0;i<64;i++) r[64]+=r[i]*r[i];
  r[64]=sqrt(r[64]);
}
    
    
double SimEx(long int start, long int end, int frame) 
{
  long int i;
  double sim=0;
  
  if(start<end) {
    for(i=start;i<=end;i++) {
      if(C[frame][i]==0) return(-HUGE_VAL);
      sim+=C[frame][i];
    }
  }
  else {
    for(i=start;i>=end;i--) {
      if(C[4+frame][i]==0) return(-HUGE_VAL);
      sim+=C[4+frame][i];
    }
  }
   
  return(sim);
}

void compute_codings(long int start,long int end, double *cod)
{
  long int i;
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
  
int basetoint(char c)
{
  switch(c) {
  case 'a':
  case 'A': return(0);
  case 'c':
  case 'C': return(1);
  case 'g':
  case 'G': return(2);
  case 't':
  case 'T': return(3);
  }
  return(-1);
}
