/*Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
Maryland, U.S.A.  All rights reserved.*/

/* 
*   glimmerm.cpp was designed by Mihaela PERTEA starting from glimmer.cpp
*   This program intends to find open reading frames in the file named
*  on the command line and scores them using the delta information
*  in the files whose names are prefixed by the second command-line
*  parameter.
*  It was designed to find also introns 
*/


#include  "delcher.h"
#include  "gene.h"

// Define GT_AG only to print the GT AG of the sequence file
//#define GT_AG
//#define MESSAGE


//const int  DEFAULT_MIN_GENE_LEN = 90;
const int  DEFAULT_MIN_GENE_LEN = 175;

const double  DEFAULT_MIN_OLAP_PERCENT = 0.10;
const int  DEFAULT_THRESHOLD_SCORE = 99;
const int  DEFAULT_MIN_OLAP = 30;
const int  DEFAULT_CHOOSE_FIRST_START_CODON = TRUE;
const int  DEFAULT_USE_INDEPENDENT = TRUE;
const int  PATTERN_LEN = 8;
const double  ATG_THRESHOLD = 2.0;
const double  BIG_NEGATIVE = -1000.0;
const int  MAX_FREE_LEN = 3;
const int  ORF_SIZE_INCR = 1000;
const int  UPSTREAM_LEN = 14;
const int  UPSTREAM_OFFSET = 3;

const unsigned  ATG_MASK = 0x184;
const unsigned  CAA_MASK = 0x211;
const unsigned  CAC_MASK = 0x212;
const unsigned  CAG_MASK = 0x214;
const unsigned  CAT_MASK = 0x218;
const unsigned  CAY_MASK = 0x21a;
const unsigned  CTA_MASK = 0x281;
const unsigned  CTG_MASK = 0x284;
const unsigned  GTG_MASK = 0x484;
const unsigned  RTG_MASK = 0x584;
const unsigned  TAA_MASK = 0x811;
const unsigned  TAG_MASK = 0x814;
const unsigned  TAR_MASK = 0x815;
const unsigned  TCA_MASK = 0x821;
const unsigned  TGA_MASK = 0x841;
const unsigned  TRA_MASK = 0x851;
const unsigned  TTA_MASK = 0x881;
const unsigned  TTG_MASK = 0x884;
const unsigned  TYA_MASK = 0x8a1;
const unsigned  YTA_MASK = 0xa81;
const unsigned  SHIFT_MASK = 0xFF;

const int  MODEL_LEN = 9;
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const int  WINDOW_SIZE = 48;
const double  MAX_LOG_DIFF = -20.0;


#include  "context.h"


const unsigned int  OK = 0x0;
const unsigned int  SCORES_WORSE = 0x1;
const unsigned int  SHORTER = 0x2;
const unsigned int  SHADOWED_BY = 0x4;
const unsigned int  SHADOWS_ANOTHER = 0x8;
const unsigned int  REJECT_MASK = 0x3;
const char  SCORES_WORSE_CHAR = 'B';
const char  SHORTER_CHAR = 'S';
const char  SHADOWED_BY_CHAR = 'W';
const char  REJECT_CHAR = 'R';

char TRAIN_DIR[500]="";

struct  Problem_Node
  {
   char  Problem_Code;
   long int  From, Olap, Delay;
   int  Score;
   Problem_Node  * Next;
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

struct  Gene_Ref
  {
    long int  Lo, Hi, Max_Hi, Min_Lo;
    int  Frame;
    unsigned int  Status;
    char  Reject_Code;
    //Problem_Node  * Problem_List;
    intron * Intron_List;
    int Score[7];
    long int newno;
  };

struct  ED_Struct
  {
   float  Free_i, Free_j, Both_Free, Match;
   int  Free_i_Len, Free_j_Len, Both_i_Len, Both_j_Len;
  }  ED_Score [1 + PATTERN_LEN] [1 + UPSTREAM_LEN];

struct myint
{
  int yn;
  double score;
};



#define MAX_INTRON_LENGTH  850 // found one of 801!!!
#define MIN_INTRON_LENGTH  50
#define MIN_FIRST_EXON_LENGTH 3
#define MIN_EXON_LENGTH 10
#define MIN_LASTEXON_LENGTH 13
#define MAX_EXON_LENGTH 9500
#define MAX_GENE_LENGTH 20000
#define  ALLOVER 100;   // the minimum length to reject the shortest gene


void  Add_Problem (char, Gene_Ref &, long int, long int, int, long int);
double  Bulge_Cost  (int);
unsigned  Ch_Mask  (char);
void  Check_Previous  (Gene_Ref * &, long int *, long int, long int, int,intron *Cap, int [7]);
int  Choose_Score  (int [7], int);
long int  Choose_Start  (long int, long int);
int  Cmp  (const void *, const void *);
double  Doublet_Score  (char, char, char, char);
double  Edit_Distance  (const char *, const char *);
void  Find_Stop_Codons  (char [], int, int []);
void  Indep_Eval  (char [], int, double [], double &);
int  Is_Forward_Start  (unsigned);
int  Is_Forward_Stop  (unsigned);
int  Is_Reverse_Start  (unsigned);
int  Is_Reverse_Stop  (unsigned);
int  Is_Start  (char *);
int  Is_Stop  (char *);
double  Loop_Cost  (int, int);
int  Match  (char, char);
void  Permute  (int [], int);
void  Process_Options  (int, char * []);
void  Score_String  (char *, long int, long int, long int,
                    String_Array * [6], double [ALPHABET_SIZE], int [7]);
void  Simple_Score  (char [], int, int, String_Array * [6],
                     double [ALPHABET_SIZE], int [7]);
void  Transfer  (char * , long int, int Len,char * , long int );
void Score_Gene_F(int Frame, long int For_Prev, 
		  Gene_Ref * &Previous,long int *Id_Num,
		  long int i,long int Len,intron *Cap);
void Score_Gene_R(int Frame, long int Rev_Prev,
		  Gene_Ref * &Previous,long int *ID_Num,
		  long int i,long int Len,intron *Cap);
int StopCod(long int i);
void Intron_Left(long int poz, int Frame,
		 long int *For_Prev,Gene_Ref *&Previous,
		 long int *ID_Num,long int iL, long int Len,intron *Cap);
void Intron_Right(long int poz, int Frame,
		 long int *Rev_Prev,Gene_Ref *&Previous,
		 long int *ID_Num,long int iL, long int Len,intron *Cap);
/* void printreverse(intron *Cap,long int *last,float percent[50], long int *no);*/
void printreverse(intron *Cap);

int Compare_length(intron *Cap1,long int lo1,long int hi1,int f1,intron *Cap2,long int lo2,long int hi2,int f2);
void replace(Gene_Ref * &,long int,long int,long int,int,intron *,long int, int[7]);
int comp_record(const void *a, const void *b);
int Simple_Overlap(long int lo1,long int hi1,long int lo2, long int hi2);
long int Overlap(long int lo1,long int hi1,intron *Cap1,int Frame1,long int lo2, long int hi2,intron *Cap2,int Frame2);
void printoverlap(Gene_Ref Prev,long int Hi,long int Lo,int Frame,intron *Cap);
int CmpGenes(long int,long int,intron *,int ,long int ,long int ,intron *,int,double *,double *);
double  One_Score  (char [],int,int,String_Array * [6],double [ALPHABET_SIZE]);
int Delayed(Gene_Ref * & Previous,long int i,long int *Hi,long int *Lo,intron *Cap,int Frame);
void printgene(long int Hi,long int Lo,int Frame,intron *Cap);
int has_right_exon_f(intron *Cap);
int has_right_exon_r(intron *Cap);
void printsites(intron *Cap);
float Perc(long int start,long int end,long int *);
int Perc_OK(int,long int,long int,intron *);

double  Ch_Ct [ALPHABET_SIZE] = {0.0};   // Count of number of occurrences of
                                         // acgt's in entire genome
int  Choose_First_Start_Codon = DEFAULT_CHOOSE_FIRST_START_CODON;
String_Array  * Context_Delta [6];
char  * Data;
long int  Data_Len;
long int  Min_Gene_Len = DEFAULT_MIN_GENE_LEN;
long int Min_Last_Exon=MIN_LASTEXON_LENGTH;
int  Min_Olap = DEFAULT_MIN_OLAP;
double  Min_Olap_Percent = DEFAULT_MIN_OLAP_PERCENT;
char  * Orf_Buffer;
long int  Orf_Buffer_Len;
int  Threshold_Score = DEFAULT_THRESHOLD_SCORE;
int  Use_Independent = DEFAULT_USE_INDEPENDENT;
myint *IR,*IF;
int allover;

float percent[50];
long int last,no;
int Use_Filter=1;
int Acc_Win_Len=60;
int Don_Win_Len=30;

int  Is_Acceptor  (const int *, double *);
int  Is_Donor  (const int *, double *);
int Is_AD(const int *, double *);
int  Is_Atg  (const int * , double *);

Gene_Ref  * Previous;

main  (int argc, char * argv [])
{
  FILE  * fp;
  //Gene_Ref  * Previous, * * Ptr;
  long int * Ptr;
  Problem_Node  * Pred, * Prob;
  char  File_Name [MAX_LINE], Name [MAX_LINE];
  int  Frame,ret;
  unsigned  Codon;
  long int  For_Prev [3] = {1,1,1};
  long int  Rev_Prev [3] = {LONG_MAX, LONG_MAX, LONG_MAX};
  long int  For_Start [3] = {0};
  long int  Rev_Start [3] = {0};
  long int  ID_Num = 0;
  long int  Delay_Len, Virtual_Start, Virtual_End;
  long int  Lo, Hi;
  long int  i, j, k, l,Ct, Input_Size, Len,t;
  int B[200];
  double score,S1,S2;
  intron *Cap;


  if  (argc < 3) {
    fprintf (stderr,
	     "USAGE:  %s <genome-file> <delta-prefix> [options] \n",
	     argv [0]);
    exit (-1);
  }

  Process_Options (argc, argv);   // Set global variables to reflect status of
                                  // command-line options.
  allover=ALLOVER;


  fp = File_Open (argv [1], "r");

  Data = (char *) Safe_malloc (INIT_SIZE);
  Input_Size = INIT_SIZE;
  
  // Read entire genome into  Data [1 .. ]
  // it seams that if there is not enough space, the space is increased
  Read_String (fp, Data, Input_Size, Name, FALSE);   
  
  fclose (fp);

  Data_Len = strlen (Data + 1);

  for  (i = 1;  i <= Data_Len;  i ++) {

    // Converts all characters to  acgt
    Data [i] = Filter (tolower (Data [i]));
    
    switch  (Data [i]) {
      
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

  // find introns

  // forward direction

  IF=(myint *) malloc((Data_Len+2)*sizeof(myint));
  if (IF == NULL) {
    fprintf(stderr,"Memory allocation for intron failure.\n"); 
    abort();
  }

  for(i=0;i<Data_Len+1;i++)    IF[i].yn=0;
  
  // first look fot gt and ag
  for(i=80;i<=Data_Len-82;i++){
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
      ret=Is_Acceptor(B,&score);
      printf("Found ag at %ld w/ score %f\n",i,score);
#endif GT_AG

      if(Is_Acceptor(B,&score)) {
	// add acceptor to list;
	
#ifdef GT_AG
	printf("Found forward ag at %ld w/ score %f\n",i,score);
#endif GT_AG
	
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
      ret=Is_Donor(B,&score);
      printf("Found gt at %ld w/ score %f\n",i,score);
#endif GT_AG
	
      if(Is_Donor(B,&score)) {
	// add donor to list;
	
#ifdef GT_AG
	printf("Found forward gt at %ld w/ score %f\n",i,score);
#endif GT_AG
	
	IF[i].yn=1;
	IF[i].score=score;
      }
    }
    
  }


  // reversed direction

  IR=(myint *) malloc((Data_Len+2)*sizeof(myint));
  if (IR == NULL) {
    fprintf(stderr,"Memory allocation for intron failure.\n"); 
    abort();
  }

  for(i=0;i<Data_Len+1;i++) IR[i].yn=0;

  // first look fot gt and ag
  for(i=Data_Len-81;i>=82;i--){
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
      ret=Is_Acceptor(B,&score);
      printf("Found ag at %ld w/ score %f\n",i,score);
#endif GT_AG

      if(Is_Acceptor(B,&score)) {
	// add acceptor to list;
	
#ifdef GT_AG
	printf("Found reverse ag at %ld w/ score %f\n",i,score);
#endif GT_AG

	IR[i].yn=2;
	IR[i].score=score;
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
      ret=Is_Donor(B,&score);
      printf("Found gt at %ld w/ score %f\n",i,score);
#endif GT_AG

      if(Is_Donor(B,&score)) {
	// add donor to list;
	
#ifdef GT_AG
	printf("Found reverse gt at %ld w/ score %f\n",i,score);
#endif GT_AG
	
	IR[i].yn=1;
	IR[i].score=score;
      }
    }
  }

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

  //print found splice sites
  
#ifdef GT_AG 

  printf("Forward:\n");
  for(i=1;i<Data_Len;i++){
    if(IF[i].yn==1) printf("gt: %d\n",i);
    if(IF[i].yn==2) printf("ag: %d\n",i);
  }
     
  printf("Reverse:\n");
  for(i=Data_Len;i>1;i--){
    if(IR[i].yn==1) printf("gt: %d\n",i);
    if(IR[i].yn==2) printf("ag: %d\n",i);
  }
  exit(1);

#endif GT_AG


  
  Ch_Ct [2] = Ch_Ct [1];               // Counts are for *both* strands
  Ch_Ct [3] = Ch_Ct [0];
  for  (i = 0;  i < 4;  i ++)
    // Convert to log of proportion of at vs. gc
    Ch_Ct [i] = log (Ch_Ct [i] / (2.0 * Data_Len));  
      

  // the second file is used in order to find the name of deltas.context
  strcpy (File_Name, argv [2]);
  strcat (File_Name, ".deltas.context");
  
  fp = File_Open (File_Name, "r");
  for  (i = 0;  i < 6;  i ++) {                  // Create six frame models
     
    Context_Delta [i] = new String_Array (MODEL_LEN, ALPHABET_SIZE);
    
    // Read the model that  build-imm  created
    Context_Delta [i] -> Read (fp);      
    
    Context_Delta [i] -> Convert_To_Logs ();
  }
  fclose (fp);


   // Create a buffer to hold orfs for further processing

 Orf_Buffer_Len = ORF_SIZE_INCR;
 Orf_Buffer = (char *) Safe_malloc (Orf_Buffer_Len);
 Orf_Buffer [0] = ' ';

 // Previous is an array of structures representing potential
 // genes that have been found so far.

 Previous = (Gene_Ref *) Safe_malloc (sizeof (Gene_Ref));

 printf ("Minimum gene length = %d\n", Min_Gene_Len);
 printf ("Minimum overlap length = %d\n", Min_Olap);
 printf ("Minimum overlap percent = %.1f%%\n", 100.0 * Min_Olap_Percent);
 printf ("Threshold score = %d\n", Threshold_Score);
 printf ("Use independent scores = %s\n", Use_Independent ? "True" : "False");
 printf ("Use first start codon = %s\n",
	 Choose_First_Start_Codon ? "True" : "False");
 putchar ('\n');
 

 //  Prepare first codon which will be characters 1, 2 and 3
 //  for wraparound, since genome is circular.

 Codon = Ch_Mask (Data [1]) << 4 | Ch_Mask (Data [2]);

 Frame = 0;


 //  Main loop for each base in genome.

 for  (i = 3;  Data [i] != '\0';  i ++) {
     
   Codon = (Codon & SHIFT_MASK) << 4;
   Codon |= Ch_Mask (Data [i]);
   Frame = (Frame + 1) % 3;

   // if a stop codon here
   if  (Is_Forward_Stop (Codon)) {

#ifdef MESSAGE
     printf("Found stop at %ld\n",i);
#endif MESSAGE
    
     // printf("Frame=%d For_prev=%ld\n",Frame,For_Prev[Frame]);


     Len = i - For_Prev [Frame] - 3;     // Length of forward orf ending here


     // printf("Len=%ld\n",Len);
    
     if  (Len >= Min_Last_Exon) {

       // Here the intron model should be introduce, like following:
       Cap=NULL;
       Score_Gene_F(Frame,For_Prev[Frame],Previous,&ID_Num,i,Len,Cap);

       // find all genes with introns
       Intron_Left(i,Frame,For_Prev,Previous,&ID_Num,i,Len,Cap);

     }


     For_Prev [Frame] = i;    // Remember this stop codon location
     For_Start [Frame] = 0;   // Reset position of first start codon to empty
   }

   // If is first start codon in this frame
   // Modify Is_Forward_Start for the mask of the start codon
   if  (Is_Forward_Start (Codon) && For_Start [Frame] == 0)    
     For_Start [Frame] = i - 2;                              // remember it

   // why is i-2 ?

   // Comment for the moment the reverse fragment. Maybe it's better to 
   // use only direct strand, then generate reverse and look for genes 
   // there too?

   // Here is a comment

   if  (Is_Reverse_Stop (Codon)) {

     Len = i - Rev_Prev [Frame] - 3;
     if  (Len >= Min_Last_Exon) {
       Cap=NULL;

       
       Score_Gene_R(Frame,Rev_Prev[Frame],Previous,&ID_Num,Rev_Prev[Frame],Len,Cap);


       Intron_Right(Rev_Prev[Frame],Frame,Rev_Prev,Previous,&ID_Num,Rev_Prev[Frame],Len,Cap);

     }

     Rev_Prev [Frame] = i;
     Rev_Start [Frame] = 0;
   }

   if  (Is_Reverse_Start (Codon))    // Remember all reverse starts since want the rightmost
     Rev_Start [Frame] = i;  


	// Here ends the comment
 }
	


 // Why is this assert ?
 assert (Data_Len == i - 1);

 // what is .Min_Lo ?
 Previous [ID_Num] . Min_Lo = Previous [ID_Num] . Lo;
 for  (i = ID_Num - 1;  i > 0;  i --)
   Previous [i] . Min_Lo = Min (Previous [i] . Lo, Previous [i + 1] . Min_Lo);
 
   
 // End of processing  
 //printf ("End = %ld\n", Data_Len);fflush(stdout);

 // post process the data
 for(i=1;i<=ID_Num-1;i++) 
   if(Previous[i].Reject_Code!='R') {

     /* see if the gene is ok by percent */
     // if(Perc_OK(Previous[i].Frame,Previous[i].Lo,Previous[i].Hi,Previous[i].Intron_List)) {

     for(j=i+1;j<=ID_Num;j++) {

       if(Previous[j].Reject_Code!='R') {
	 if(Simple_Overlap(Previous[i].Lo,Previous[i].Hi,Previous[j].Lo,Previous[j].Hi)){

	   if(Overlap(Previous[i].Lo,Previous[i].Hi,Previous[i].Intron_List,Previous[i].Frame,Previous[j].Lo,Previous[j].Hi,Previous[j].Intron_List,Previous[j].Frame))
	      {
		ret=Compare_length(Previous[j].Intron_List,Previous[j].Lo,Previous[j].Hi,Previous[j].Frame,Previous[i].Intron_List,Previous[i].Lo,Previous[i].Hi,Previous[i].Frame);

		if(ret==0) { Previous[i].Reject_Code='R'; Previous[i].Status=j;} 
		if(ret==1) { Previous[j].Reject_Code='R'; Previous[j].Status=i;}
		if(ret==3) Previous[i].Status=j;
		if(ret==2) Previous[j].Status=i;

	      }
	 }
       }
     }
     
   }


 


  /*  for(i=1;i<=ID_Num-1;i++) 
      if(Previous[i].Reject_Code!='R') {
      
      for(j=i+1;j<=ID_Num;j++) {
      if(Previous[j].Reject_Code!='R') {
      if(Simple_Overlap(Previous[i].Lo,Previous[i].Hi,Previous[j].Lo,Previous[j].Hi)){
      
      if(Compare_length(Previous[j].Intron_List,Previous[j].Lo,Previous[j].Hi,Previous[j].Frame,Previous[i].Intron_List,Previous[i].Lo,Previous[i].Hi,Previous[i].Frame))  
      if( (Previous[j].Reject_Code!='o' ||Previous[i].Reject_Code!='O') &&
      (Previous[i].Reject_Code!='o' ||Previous[j].Reject_Code!='O'))
      Previous[j].Reject_Code='R';
      
      }
      }
      }
      }*/
    

 /* print info comparisons about data */

 /* for(i=1;i<=ID_Num;i++)
   if(Previous[i].Reject_Code!='R') {
     
     for(j=i+1;j<=ID_Num;j++) {
       if(Previous[j].Reject_Code!='R') {
	 if(Simple_Overlap(Previous[i].Lo,Previous[i].Hi,Previous[j].Lo,Previous[j].Hi)){
	   
	   if(CmpGenes(Previous[i].Lo,Previous[i].Hi,Previous[i].Intron_List,Previous[i].Frame,Previous[j].Lo,Previous[j].Hi,Previous[j].Intron_List,Previous[j].Frame,&S1,&S2)) {
	     printf("Gene %ld is better then gene %ld : %f-%f\n",Previous[j].newno,Previous[i].newno,S2,S1);
	   }
	   else 
	     printf("Gene %ld is better then gene %ld : %f-%f\n",Previous[i].newno,Previous[j].newno,S1,S2);
	 }
       }
     }
   }*/
  
 

 Ptr = (long int *) Safe_malloc (ID_Num * sizeof (long int));

 for  (i = 1;  i <= ID_Num;  i ++)
   Ptr [i - 1] = i;

 qsort (Ptr, ID_Num, sizeof (long int), Cmp);   // sort by decreasing length


 k=0;

 for(i=0;i<ID_Num;i++) {
   if(Previous[Ptr[i]].Reject_Code!='R') {
     k++;
     Previous[Ptr[i]].newno=k;
   }
   else Previous[Ptr[i]].newno=0;
 }



 
 /* for  (i = 0;  i < ID_Num;  i ++)
    if  (Ptr [i] -> Problem_List != NULL) {
    
    Ct = 0;
    Prob = Ptr [i] -> Problem_List;
    do {
    
    Prob = Prob -> Next;
    Ct ++;
    }  while  (Prob != Ptr [i] -> Problem_List);
    

    Ptr [i] -> Reject_Code = ' ';
    Prob = Ptr [i] -> Problem_List;
    for  (j = 0;  j < Ct;  j ++) {
       
    Pred = Prob;
    Prob = Prob -> Next;
    if  (Previous [Prob -> From] . Reject_Code == 'R') {
	 
    if  (Pred == Prob) {
    
    free (Prob);
    Prob = NULL;
    }
    else {
    
    Pred -> Next = Prob -> Next;
    free (Prob);
    Prob = Pred;
    }
    }
    else if  (Prob -> Problem_Code == REJECT_CHAR) {
    
    if  (Prob -> Delay > 0 && Ptr [i] -> Reject_Code != 'R') {
    
    Ptr [i] -> Reject_Code = 'D';
    }
    else
    Ptr [i] -> Reject_Code = 'R';
    }
    }
    Ptr [i] -> Problem_List = Prob;
    }*/

 printf("\n\n\nGene scores:\n\n");fflush(stdout);
 printf("GeneNo   F1 F2 F3 F4 F5 F6 IndScore     Splite sites scores\n");fflush(stdout);
 printf("_____________________________________________________________________________________________\n");fflush(stdout);
 k=0;
 for(i=0;i<ID_Num;i++)
   if  (Previous[Ptr[i]]. Reject_Code != 'R') {

     k++;
     printf(" %4ld    ",k);fflush(stdout);
     for(j=0;j<7;j++)
       if(Previous[Ptr[i]].Score[j]<0) {
	 printf(" - ");fflush(stdout);
       }
       else printf("%2ld ",Previous[Ptr[i]].Score[j]);	fflush(stdout);

     if(Previous[Ptr[i]].Frame >0) {
	Cap=Previous[Ptr[i]].Intron_List;
	ret=0;
	if(Cap!=NULL) { printf("          start ");fflush(stdout); ret=1;}
	while(Cap!=NULL) {
	   printf("%f %f ",IF[Cap->gt].score,IF[Cap->ag].score);fflush(stdout);
	   Cap=Cap->leg;
	}
	if(ret) {printf("stop");fflush(stdout);}
     }
     else {
	Cap=Previous[Ptr[i]].Intron_List;
	if(Cap!=NULL) printsites(Cap);
     }

     printf("\n");fflush(stdout);
   }


 printf ("\n\n\nPutative Genes:\n");fflush(stdout);

 k=0;

 for  (i = 0;  i < ID_Num;  i ++)
   if  (Previous[Ptr [i]] . Reject_Code != 'R') {

     k++;
    
     Lo =Previous[ Ptr [i]] . Lo;
     Hi = Previous[Ptr [i]] . Hi;
     if  (Previous[Ptr [i]] . Frame > 0) {
       printf ("%5ld %8ld ", k, Lo);fflush(stdout);
       
       last=Lo;no=0;

       Cap=Previous[Ptr[i]].Intron_List;
       while(Cap !=NULL) {
	 printf("%8ld %8ld ",(Cap->gt)-1,(Cap->ag)+2);fflush(stdout);

	 percent[no++]=Perc(last,(Cap->gt)-1,&t);
	 percent[no++]=Perc(Cap->gt,(Cap->ag)+1,&t);
	 last=(Cap->ag)+2;

	 Cap=Cap->leg;
       }
       printf("%8ld",Hi+3);fflush(stdout);
       
       percent[no++]=Perc(last,Hi+3,&t);
       
       printf(" AT%: ");fflush(stdout);
       for(j=0;j<no;j++) 
	 printf("%3.2f ",percent[j]); fflush(stdout);

     }
     else {

       last=Hi;no=0;

       printf("%5ld %8ld ",k,Hi);fflush(stdout);
       Cap=Previous[Ptr[i]].Intron_List;

       if(Cap!=NULL) printreverse(Cap);
       printf ("%8ld", Lo-3);fflush(stdout);
       percent[no++]=Perc(Lo-3,last,&t);

       printf(" AT%: ");fflush(stdout);
       for(j=0;j<no;j++) 
	 printf("%3.2f ",percent[j]); fflush(stdout);


     }

     j=Previous[Ptr[i]].Status; 

     if(j) {

       if(Previous[j].newno) { printf(" Bad overlap with gene %ld.",Previous[j].newno);fflush(stdout); }

       //printf("\nGene is: %ld, former %ld\n",k,i);
       //printf("Rej of %ld: %c\n",k,Previous[i].Reject_Code);
       //printf("\nGene is: %ld, former %ld\n",Previous[j].newno,j);
       //printf("Rej of %ld: %c\n",Previous[j].newno,Previous[j].Reject_Code);
     }
     
     putchar ('\n');
     
   }
 
 return  0;
}

void printsites(intron *Cap)
{
  double poz[100];
  int l,i;
  
  l=-1;

  while(Cap!=NULL) {
    poz[++l]=IR[Cap->ag].score;
    poz[++l]=IR[Cap->gt].score;
    Cap=Cap->leg;
  }

  printf("          start ");fflush(stdout);

  for(i=l;i>=0;i--)
    { printf("%f ",poz[i]);fflush(stdout);}
 
  printf("stop");fflush(stdout);
}

int Perc_OK(int frame,long int lo,long int hi,intron *cap)
{
  int first;
  float p, genescore;
  long int l;
  long int poz[100];
  int lr,i,exonframe;
  long int length,at,sum;
  char temp[125000];
  long int len,min_first,min_exon;
  int Score[7];

  length=0;

  sum=0;

  min_first=MIN_FIRST_EXON_LENGTH;
  min_exon=MIN_EXON_LENGTH;

  first=1;

  genescore=0;

  if(frame>0) {
    l=lo;

    while(cap!=NULL){
      p=Perc(l,(cap->gt)-1,&at);  // exon

      len=(cap->gt)-l;

      if(first) {
	if(len<min_first) {
	  return(0);
	}
      }
      else {
	if(len<min_exon) return(0);
      }

      strncpy(temp,Data+l,len);
      temp[len]='\0';
      if(strstr(temp,"TATTATTAT") != NULL) {
	return(0);
      }
      sum+=at;

      exonframe=(3-length%3)%3;
      Score_String(Data,l+exonframe,(cap->gt)-1,Data_Len,Context_Delta,Ch_Ct,Score);

      genescore+=Score[0]*len;

#ifdef MESSAGE
      printf("Exon %ld %ld: %f Score: %d Length: %ld\n",l,(cap->gt)-1,p,Score[0],len);
      fflush(stdout);

#endif MESSAGE

      if(Score[0]<Threshold_Score-5 && len >=250) return(0);
      length+=len;

      //if(p>80) return(0);
      if(first) {
	if(p>=86 && len>15) {
	  return(0); 
	}
	first=0;
      }
      else if(p>80) {
	return(0);
      }

      p=Perc(cap->gt,(cap->ag)+1,&at); //intron

#ifdef MESSAGE
      printf("Intron %ld %ld: %f\n",cap->gt,(cap->ag)+1,p);fflush(stdout);
#endif MESSAGE

      if(p<=65) {
	return(0);
      }
      l=(cap->ag)+2;
      cap=cap->leg;
    }

    p=Perc(l,hi+3,&at); //exon

    len=hi+4-l;

    strncpy(temp,Data+l,len);
    temp[len]='\0';
    if(strstr(temp,"TATTATTAT") != NULL) {
      return(0);
    }

    sum+=at;

    exonframe=(3-length%3)%3;
    Score_String(Data,l+exonframe,hi,Data_Len,Context_Delta,Ch_Ct,Score);

    genescore+=Score[0]*len;

#ifdef MESSAGE
    printf("Exon %ld %ld: %f Score: %d Length: %ld\n",l,hi+3,p,Score[0],len);
    fflush(stdout);    
#endif MESSAGE

    if(Score[0]<Threshold_Score-5 && len>=250) return(0);


    //    if(p>80) return(0);
    if(first) {
      if(p>=83) {
	return(0);
      }
    }
    else if(p>80) {
      return(0);
    }

    length+=len;

    genescore/=length;

    if(length<Min_Gene_Len || genescore<=Threshold_Score-15) return(0);

    p=sum*100;
    p=p/length;
    if(p>=83) return(0);

    return(1);
  }
  else {
    l=hi;

    lr=-1;

    while(cap!=NULL) {
      poz[++lr]=(cap->ag)-2;
      poz[++lr]=(cap->gt)+1;
      cap=cap->leg;
    }  

    for(i=lr;i>=0;i--) {
      if(i%2) {
	p=Perc(poz[i],l,&at); //exon

	len=l-poz[i]+1;


	if(first) {
	  if(len<min_first) return(0);
	}
	else {
	  if(len<min_exon) return(0);
	}

	strncpy(temp,Data+poz[i],len);
	temp[len]='\0';
	if(strstr(temp,"ATAATAATA") != NULL) return(0);

       	sum+=at;

	exonframe=(3-length%3)%3;
	Score_String(Data,l-exonframe,poz[i],Data_Len,Context_Delta,Ch_Ct,Score);

	genescore+=Score[0]*len;

#ifdef MESSAGE
	printf("Exon %ld %ld: %f Score: %d Length: %ld\n",poz[i],l,p,Score[0],len);
	fflush(stdout);
#endif MESSAGE

	if(Score[0]<Threshold_Score-5 && len>=250) return(0);

	//	if(p>80) return(0);
	if(first) {if(p>=86 && len>15) return(0); first=0;}
	else if(p>80) return(0);

	length+=len;

      }
      else {
	p=Perc(poz[i]+1,l-1,&at);

#ifdef MESSAGE
	printf("Intron %ld %ld: %f\n",poz[i]+1,l-1,p);fflush(stdout);
#endif MESSAGE

	if(p<=65) return(0); //intron

      }
      l=poz[i]; 
       
    } 

    len=l-lo+4;
    strncpy(temp,Data+lo-3,len);
    temp[len]='\0';
    if(strstr(temp,"ATAATAATA") != NULL) return(0);   

    p=Perc(lo-3,l,&at); //exon

    exonframe=(3-length%3)%3;
    Score_String(Data,l-exonframe,lo,Data_Len,Context_Delta,Ch_Ct,Score);

    genescore+=Score[0]*len;

#ifdef MESSAGE
    printf("Exon %ld %ld: %f Score: %d Length: %ld\n",lo-3,l,p,Score[0],len);
    fflush(stdout);
#endif MESSAGE

    if(Score[0]<Threshold_Score-5 && len >=250) return(0);

    length+=len;
    //if(p>80) return(0);
    if(first) {if(p>=83) return(0);}
    else if(p>80) return(0);

    sum+=at;
    p=sum*100;
    p=p/length;
    if(p>=83) return(0);

    genescore/=length;

    if(length<Min_Gene_Len || genescore<=Threshold_Score-15) return(0);
    return(1);
  }

}

void printreverse(intron *Cap)
{
  long int poz[100],t;
  int l,i;
  
  l=-1;

  while(Cap!=NULL) {
    poz[++l]=(Cap->ag)-2;
    poz[++l]=(Cap->gt)+1;
    Cap=Cap->leg;
  }

  for(i=l;i>=0;i--) {
    printf("%8ld ",poz[i]);fflush(stdout);
    if(i%2) percent[no++]=Perc(poz[i],last,&t);
    else percent[no++]=Perc(poz[i]+1,last-1,&t);
    last=poz[i];

  }
      
}


void Intron_Left(long int poz, int Frame, long int *For_Prev,
		 Gene_Ref *&Previous, long int *ID_Num,long int iL, 
		 long int Len,intron *Cap) 
// forward direction
{
  long int i,j,l,t;
  int  k,NFrame;
  intron *node;

  // look for all possible AG at left 
  for(i=poz-MIN_EXON_LENGTH-1;i>=poz-MAX_EXON_LENGTH-1;i--) { 

    if(i<=1) break;

    // if frame contains a stop codon break;
    if(((poz-i)%3==0)&& StopCod(i+3)) { 
      break;
    }

    if(IF[i].yn==2) { // found an AG

      for(j=i+1-MIN_INTRON_LENGTH;j>=i+1-MAX_INTRON_LENGTH;j--) { //look for GT
	
	if(j<=1) break;

	if(IF[j].yn==1) { // found a GT

	  k=(i-j+2)%3;   // compute new frame


	  NFrame=(3+Frame-k)%3;

	    
	  // create new node in list of introns Cap


	  node=(intron *)malloc(sizeof(intron));  
	  if (node == NULL) {
	    fprintf(stderr,"Memory allocation for node failure.\n"); 
	    abort();
	  }
	  node->gt=j;
	  node->ag=i;
	  node->leg=Cap;
	  Cap=node;

      
	  if(has_right_exon_f(Cap)) { // see if gene is a good one respecting ad

	    Score_Gene_F(NFrame,For_Prev[NFrame],Previous,ID_Num,iL,Len,Cap);

	    l= poz-i+j-2;
	    t=j-1;
	    while((l-t)%3!=0) t--;
	    
	    Intron_Left(t,NFrame,For_Prev,Previous,ID_Num,iL,Len,Cap);

	  }

	  Cap=node->leg;
	  free(node);
	  
	}
      }

    }

  }

}

int has_right_exon_r(intron *Cap)
{
  int B[40];
  long int i,next;
  intron *temp;
  double score;


  temp=Cap;
  next=0;

  while(temp != NULL) {
    
    if(next) {
      for(i=0;i<20;i++) {
	switch(Data[18-i+temp->ag]) {
	case 'A':
	case 'a': B[i]=3;break;
	case 'C':
	case 'c': B[i]=2;break;
	case 'G':
	case 'g': B[i]=1;break;
	case 'T':
	case 't': B[i]=0;break;
	default: B[i]=1;
	}
      }
      for(i=20;i<40;i++){
	switch(Data[20-i+next]) {
	case 'A':
	case 'a': B[i]=3;break;
	case 'C':
	case 'c': B[i]=2;break;
	case 'G':
	case 'g': B[i]=1;break;
	case 'T':
	case 't': B[i]=0;break;
	default: B[i]=1;
	}
      }

      //Is_AD(B,&score);

      if(!Is_AD(B,&score)) return(0);
    }
    
    next=temp->gt;
      
    temp=temp->leg;
  }
  return(1);
}


int has_right_exon_f(intron *Cap)
{
  int B[40];
  long int i,next;
  intron *temp;
  double score;



  /*printf("List of introns:\n");
  temp=Cap;
  while(temp!=NULL) {
    printf("gt=%ld - ag=%ld\n",temp->gt,temp->ag);
    temp=temp->leg;
  }*/
  

  temp=Cap;
  next=0;

  while(temp != NULL) {
    
    if(next) {
      for(i=0;i<20;i++) {
	switch(Data[i-18+next]) {
	case 'A':
	case 'a': B[i]=0;break;
	case 'C':
	case 'c': B[i]=1;break;
	case 'G':
	case 'g': B[i]=2;break;
	case 'T':
	case 't': B[i]=3;break;
	default: B[i]=1;
	}
      }
      for(i=20;i<40;i++){
	switch(Data[i-20+temp->gt]) {
	case 'A':
	case 'a': B[i]=0;break;
	case 'C':
	case 'c': B[i]=1;break;
	case 'G':
	case 'g': B[i]=2;break;
	case 'T':
	case 't': B[i]=3;break;
	default: B[i]=1;
	}
      }

      //Is_AD(B,&score);

      if(!Is_AD(B,&score)) return(0);
    }
    
    next=temp->ag;
      
    temp=temp->leg;
  }
  return(1);
}
            

void Intron_Right(long int poz, int Frame,
		 long int *Rev_Prev,Gene_Ref *&Previous,
		 long int *ID_Num,long int iL, long int Len,intron *Cap) 
// reversed direction
{
  long int i,j,l,t;
  int  k;
  intron *node,*temp;

  // look for all possible AG at right 
  for(i=poz+MIN_EXON_LENGTH;i<=poz+MAX_EXON_LENGTH+1;i++) { 
   
    if(i>=Data_Len) break;

    // if frame contains a stop codon break;
    if(((i-poz)%3==0)&& StopCod(-i)) break;

    if(IR[i].yn==2) { // found an AG


      for(j=i+1+MIN_INTRON_LENGTH;j<=i+MAX_INTRON_LENGTH;j++) { //look for GT

	if(j>=Data_Len) break;

	if(IR[j].yn==1) { // found a GT    
	  	  
	  // create new node in list of introns Cap

	  node=(intron *)malloc(sizeof(intron));  
	  if (node == NULL) {
	    fprintf(stderr,"Memory allocation for node failure.\n"); 
	    abort();
	  }
	  node->gt=j;
	  node->ag=i;
	  node->leg=NULL;
	    
	  if(Cap==NULL) Cap=node;
	  else { 
	    temp=Cap;
	    while(temp->leg != NULL) temp=temp->leg;
	    temp->leg=node;
	  }

	  if(has_right_exon_r(Cap)) {

	    Score_Gene_R(Frame,Rev_Prev[Frame],Previous, ID_Num,iL,Len,Cap);

	    l=j-i+2;
	    t=j+1-l;
	    while((t-poz)%3!=0) t++;
	    
	    Intron_Right(t+l,Frame,Rev_Prev,Previous,ID_Num,iL,Len,Cap);

	  }


	  if(Cap->leg == NULL) Cap=NULL;
	  else {
	    temp=Cap;
	    while((temp->leg)->leg !=NULL) 
	      temp=temp->leg;
	    temp->leg=NULL;
	  }

	  free(node);
	  
	}
      }
    }
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

int StartCod(long int i)
{
 if(i>0) {
    if(Data[i-2]=='a' && Data[i-1]=='t' && Data[i]=='g') return(1);
    return(0);
  }
  else {
    i=-i;
    if(Data[i]=='t' && Data[i-1]=='a' && Data[i-2]=='c') return(1);
    return(0);
  }
}

void  Add_Problem (char Code, Gene_Ref & R, long int From, long int Olap,
                   int Score, long int Delay)

/* Add a node to  R 's problem list for the problem indicated by  Code  caused
*  by gene number  From  with overlap length  Olap  and in-frame score  Score .
*  Delay  is the number of bases by which the start can be postponed to
*  eliminate a reject condition. */

  {
   Problem_Node *  P;

   P = (Problem_Node *) Safe_malloc (sizeof (Problem_Node));

   P -> Problem_Code = Code;
   P -> From = From;
   P -> Olap = Olap;
   P -> Score = Score;
   P -> Delay = Delay;

   P -> Next = P;

   /*if  (R . Problem_List == NULL)
     else
     {
     P -> Next = R . Problem_List -> Next;
     R . Problem_List -> Next = P;
     }
     R . Problem_List = P;*/

   return;
  }



double  Bulge_Cost  (int N)

/* Return the energy cost of a bulge on one side of  N  bases. */

  {
   if  (N <= 0)
       return  BIG_NEGATIVE;

   if  (N < 4)
       return  -3.3;

   return  -3.3 - (N - 3) * (15.8 - 3.3) / 27.0;
  }



unsigned  Ch_Mask  (char Ch)

/* Returns a bit mask representing character  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  0x1;
      case  'c' :
        return  0x2;
      case  'g' :
        return  0x4;
      case  't' :
        return  0x8;
      case  'r' :     // a or g
        return  0x5;
      case  'y' :     // c or t
        return  0xA;
      case  's' :     // c or g
        return  0x6;
      case  'w' :     // a or t
        return  0x9;
      case  'm' :     // a or c
        return  0x3;
      case  'k' :     // g or t
        return  0xC;
      case  'b' :     // c, g or t
        return  0xE;
      case  'd' :     // a, g or t
        return  0xD;
      case  'h' :     // a, c or t
        return  0xB;
      case  'v' :     // a, c or g
        return  0x7;
      default :       // anything
        return  0xF;
     }
  }


void Score_Gene_R(int Frame,long int Rev_Prev,
		  Gene_Ref * &Previous,long int *ID_Num,
		  long int i,long int Len,intron *Cap) 
{
  int Score[7];
  long int Gene_Len,Start,j,end,Copy_Len,k,Rev_Start;
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

  while(temp!=NULL) {

    // eliminate the portion between gt and ag;
    strncat(Copy+1,Data+j,temp->ag-j-1);
    j=temp->gt+1;
    end=temp->ag-2;
    temp=temp->leg;
  }

  strcat(Copy+1,Data+j);
  
  Copy_Len=strlen(Copy+1);

  if(j<i) end=i;
  else end=j;

  k=end+1;
  j=end+1-(Data_Len-Copy_Len);
  while((j-i)%3!=0) {j++;k++;}

  Rev_Start=0;  

  for(j=k;j<=Data_Len;j+=3) {
    if(StopCod(-j)) break;
    if(StartCod(-j)) Rev_Start=j;
  }


  if  (Rev_Start != 0) {
	 
    Gene_Len = Rev_Start - i -(Data_Len-Copy_Len) ;
    Start = Choose_Start (Rev_Start, - Gene_Len);

    Gene_Len = Start - i -(Data_Len-Copy_Len);

    if  (Gene_Len >= Min_Gene_Len && end<=Start) {
      
      Score_String (Copy, Start-(Data_Len-Copy_Len),
		    i  + 1, Copy_Len,
		    Context_Delta, Ch_Ct, Score);
      In_Frame_Score = Score [0];
      /* original variant

	 if  (In_Frame_Score >= Threshold_Score) 
	 printf ("%5ld ", ++ (*ID_Num));
	 else
	 printf ("%5s ", "");
	 printf (" R%1d %8ld %8ld %8ld %8ld %5ld   ",
	 1 + ((2 - Frame) + 1) % 3,
	 Rev_Prev+1, Start,
	 i+1, Len,
	 Gene_Len);
	 printf (" %4d   ", In_Frame_Score);
	 Permute (Score, - Frame - 1);
	 for  (j = 0;  j < 6;  j ++)
	 if  (Score [j] < 0)
	 printf ("  _");
	 else
	 printf (" %2d", Score [j]);
	 if  (Use_Independent)
	 printf ("   %2d", Score [6]);
	 putchar ('\n');*/

      // my variant

      Permute (Score, - Frame - 1);

      // end my variant

#ifdef MESSAGE

      printf("Score gene:\n");fflush(stdout);
      printgene(Start,i+1,-Frame-1,Cap);
      printf("Score is %d\n",In_Frame_Score);fflush(stdout);

#endif MESSAGE 
      
      if  (In_Frame_Score >= Threshold_Score) {
	Check_Previous (Previous, ID_Num, i + 1,
			Start, - Frame - 1,Cap,Score);
      }
    }
  }
  free(Copy);

}



void Score_Gene_F(int Frame, long int For_Prev, 
		  Gene_Ref * &Previous,	long int *ID_Num,long int i,
		  long int Len,intron *Cap) 
{
  int Score[7];
  long int Gene_Len,Start,j,end,Copy_Len,k,For_Start;
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
  For_Start=0;
 
  k=end-1;
  while((i-Data_Len+Copy_Len-k)%3!=0) k--;

  for(j=k;j>=1;j-=3) {
    if(StopCod(j)) break;
    if(StartCod(j)) For_Start=j;
  }


  if  (For_Start != 0){  // If have seen a start codon in this orf
    
    For_Start-=2;

    // see if the gene is a good one


	 
    Gene_Len = 1 + i - 3 - For_Start-(Data_Len-Copy_Len);
    
    // modify choose start in order to work for eukariots 
    Start = Choose_Start (For_Start, Gene_Len);

    // Allows use of a fancy start codon selection alg
    Gene_Len = 1 + i - 3 - Start;
    if  (Gene_Len >= Min_Gene_Len && Start<end) {
      
      // Score this orf in all 6 frames plus the random model

      Score_String (Copy, Start, i - 3-(Data_Len-Copy_Len),
		    Copy_Len, Context_Delta, Ch_Ct,
		    Score);

	
      In_Frame_Score = Score [0]; //Frame is relative to beginning of orf


      // my variant
      
      Permute (Score, Frame + 1);

      // end my variant

      
      // !!!! modifica aici sa acepti toate genele deocamdata
      //  Compare this orf to previous genes and check for overlaps
      //    and add to  Previous  array

#ifdef MESSAGE

      printf("Score gene:\n");
      printgene(i-3,Start,Frame+1,Cap);
      printf("Score is %d\n",In_Frame_Score);fflush(stdout);

#endif MESSAGE 
      
      if  (In_Frame_Score >= Threshold_Score) {
	  
	Check_Previous (Previous, ID_Num, Start, i - 3, Frame + 1,Cap,Score);
      }


      // aici^ (in check_previous trebuie adaugata si lista mea de introni 

    }
  }
  free(Copy);

}


int Compare_length(intron *Cap1,long int lo1,long int hi1,int f1,intron *Cap2,long int lo2,long int hi2,int f2)
{
  intron *temp;
  long int l1, l2;
  long int last;
  int allowed;

  // take the whole frame if possible
  if(lo1==lo2 && hi1==hi2) {
    if(Cap1==NULL) return(0);
    if(Cap2==NULL) return(1);
  }


  l1=0;
  l2=0;

  if(f1>=0) {
    last=lo1;
    temp=Cap1;
    while(temp!=NULL) {

      l1+=(temp->gt)-last;
      last=temp->ag+2;

      temp=temp->leg;
    }
    l1+=hi1-last+1;


  }
  else {
    last=lo1;
    temp=Cap1;
    while(temp!=NULL) {
      l1+=(temp->ag)-last-1;

      last=temp->gt + 1;
      temp=temp->leg;
    }
    l1+=hi1-last+1;

  } 

  if(f2>=0) {
    last=lo2;
    temp=Cap2;
    while(temp!=NULL) {
      l2+=(temp->gt)-last;

      last=temp->ag+2;
      temp=temp->leg;
    }
    l2+=hi2-last+1;
  }
  else {
    last=lo2;
    temp=Cap2;
    while(temp!=NULL) {
      l2+=(temp->ag)-last-1;

      last=temp->gt+1;
      temp=temp->leg;
    }
    l2+=hi2-last+1;

  } 

  allowed=allover;

#ifdef MESSAGE
  printf("l1=%ld l2=%ld\n",l1,l2);
#endif MESSAGE

  if(l2-l1>allowed) {
    if(l1>5*Min_Gene_Len && f1*f2<0) return(2);
    
    if(Cap1==NULL) return(2); // keep whole genes when possbile
    
    return(1); // gene2 is much longer
  }
  if(l1-l2>allowed) {
    if(l2>5*Min_Gene_Len && f1*f2<0) return(3);
    
    if(Cap2==NULL) return(3); // keep whole genes when possbile
    
    return(0); // gene1 is much longer
    
  }
  if(l1-l2<allowed && l1-l2>=0) return(3); // gene1 is longer
  if(l2-l1<allowed && l2-l1>=0) return(2); // gene2 is longer
  
}

void replace(Gene_Ref * & Previous,long int i,long int Lo,long int Hi,int Frame,intron *Cap,long int ID,int Score[7])
{
  intron *temp;
  int j;


  Previous [i] . Lo = Lo;
  Previous [i] . Hi = Hi;
  Previous [i] . Frame = Frame;
  Previous [i] . Status = OK;
  Previous [i] . Reject_Code = ' ';
  //Previous [i] . Problem_List = NULL;
  if  (Previous [ID ] . Hi < Hi)
    Previous [i] . Max_Hi = Hi;
  else
    Previous [i] . Max_Hi = Previous [ID ] . Hi;
	 
  for(j=0;j<7;j++)     Previous[i].Score[j]=Score[j];


  if(Previous[i].Intron_List!=NULL) free(Previous[i].Intron_List);

	 
  if(Cap!=NULL) {
    Previous[i].Intron_List=(intron *)malloc(sizeof(intron));  
    if (Previous[i].Intron_List == NULL) {
      fprintf(stderr,"Memory allocation for node failure.\n"); 
      abort();
    }
    temp=Previous[i].Intron_List;
    temp->gt=Cap->gt;
    temp->ag=Cap->ag;
    temp->leg=NULL;
    
    while(Cap->leg!=NULL) {
      temp->leg=(intron *)malloc(sizeof(intron));  
      if (temp == NULL) {
	fprintf(stderr,"Memory allocation for node failure.\n"); 
	abort();
      }
      temp=temp->leg;
      Cap=Cap->leg;
      temp->gt=Cap->gt;
      temp->ag=Cap->ag;
      temp->leg=NULL;
      
    }
  }
  else Previous[i].Intron_List=NULL;
}


int Simple_Overlap(long int lo1,long int hi1,long int lo2, long int hi2)
{
  long int len;

  if(lo1<=lo2 && lo2<=hi2 && hi2<=hi1) {
    len=hi2-lo2;
    if(len>Min_Olap) return(1); else return(0);
  }
  if(lo1<=lo2 && lo2<=hi1 && hi1<=hi2) {
    len=hi1-lo2;
    if(len>Min_Olap) return(1); else return(0);
  }
  if(lo2<=lo1 && lo1<=hi2 && hi2<=hi1) {
    len=hi2-lo1;
    if(len>Min_Olap) return(1); else return(0);
  }
  if(lo2<=lo1 && lo1<=hi1 && hi1<=hi2) {
    len=hi1-lo1;
    if(len>Min_Olap) return(1); else return(0);
  }
  return(0);
}

long int Overlap(long int lo1,long int hi1,intron *Cap1,int Frame1,long int lo2, long int hi2,intron *Cap2,int Frame2)
{
  Node Score[100];
  int N,i,col_e,col_p;
  intron *temp;
  long int overlap;

  Score[0].val=lo1;
  Score[1].val=lo2;
  Score[2].val=hi1;
  Score[3].val=hi2;
  Score[0].mark=1;
  Score[1].mark=0;
  Score[2].mark=1;
  Score[3].mark=0;

  N=4;

  temp=Cap1;
  if(Frame1>0) {
    while(temp !=NULL) {
      Score[N].val=(temp->gt)-1;
      Score[N++].mark=1;
      Score[N].val=(temp->ag)+2;
      Score[N++].mark=1;
      temp=temp->leg;
    }
  }
  else {
    while(temp !=NULL) {
      Score[N].val=(temp->gt)+1;
      Score[N++].mark=1;
      Score[N].val=(temp->ag)-2;
      Score[N++].mark=1;
      temp=temp->leg;
    }
  }

  temp=Cap2;
  if(Frame2>0) {
    while(temp !=NULL) {
      Score[N].val=(temp->gt)-1;
      Score[N++].mark=0;
      Score[N].val=(temp->ag)+2;
      Score[N++].mark=0;
      temp=temp->leg;
    }
  }
  else {
    while(temp !=NULL) {
      Score[N].val=(temp->gt)+1;
      Score[N++].mark=0;
      Score[N].val=(temp->ag)-2;
      Score[N++].mark=0;
      temp=temp->leg;
    }
  }

  qsort(Score,N,sizeof(Node),comp_record);

  overlap=0;

  col_e=0;
  col_p=0;

  if(Score[0].mark) col_e=1;
  else col_p=1;

  for(i=1;i<N;i++) {
    if(col_e && col_p) {
      overlap+=Score[i].val-Score[i-1].val+1;
    }

    if(Score[i].mark) col_e=1-col_e;
    else col_p=1-col_p;
  }

  return(overlap);
}
    

int comp_record(const void *a, const void *b)
{ 
  if(((Node *)a)->val > ((Node *)b)->val) return(1);
  else if (((Node *)a)->val == ((Node *)b)->val) return(0);
  else return(0);
}  
  

void printoverlap(Gene_Ref Prev,long int Hi,long int Lo,int Frame,intron *Cap)
{
  intron *temp;
 
  printf("Overlap between genes:\n"); fflush(stdout);

  if  (Prev. Frame > 0) {
    printf ("Gene 1 %8ld ", Prev.Lo);fflush(stdout);
    temp=Prev.Intron_List;
    while(temp !=NULL) {
      printf("%8ld %8ld ",(temp->gt)-1,(temp->ag)+2);
      temp=temp->leg;
    }
    printf("%8ld",Prev.Hi+3);
  }
  else {
    printf("Gene 1 %8ld ",Prev.Hi); fflush(stdout);
    temp=Prev.Intron_List;last=Hi;
    if(temp != NULL) printreverse(temp);no=0;
    printf ("%8ld", Prev.Lo-3);fflush(stdout);
    
  }
     
  putchar ('\n');

  if  (Frame > 0) {
    printf ("Gene 2 %8ld ", Lo);fflush(stdout);
    temp=Cap;
    while(temp !=NULL) {
      printf("%8ld %8ld ",(temp->gt)-1,(temp->ag)+2);fflush(stdout);
      temp=temp->leg;
    }
    printf("%8ld",Hi+3);fflush(stdout);
  }
  else {
    printf("Gene 2 %8ld ",Hi);fflush(stdout);
    temp=Cap;no=0;
    if(temp !=NULL) printreverse(temp);
    printf ("%8ld", Lo-3);fflush(stdout);
    
  }
     
  putchar ('\n');
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


int CmpGenes(long int Lo1,long int Hi1,intron *I1,int F1,long int Lo2,long int Hi2,intron *I2,int F2, double *Score1,double *Score2)
{
  char Gene1[MAX_GENE_LENGTH], Gene2[MAX_GENE_LENGTH],Copy[MAX_GENE_LENGTH];
  long int len1,len2,last,i,k;
  intron *temp;
  double S1,S2;

  len1=0;
  len2=0;

  strcpy(Gene1+1,"");
  strcpy(Gene2+2,"");

  if(F1>=0) {
    last=Lo1;
    temp=I1;

    while(temp!=NULL) {
      len1+=(temp->gt)-last;

      strncat(Gene1+1,Data+last,(temp->gt)-last);

      last=temp->ag+2;
      temp=temp->leg;
    }
    len1+=Hi1-last+1;

    strncat(Gene1+1,Data+last,Hi1-last+1);
    Gene1[len1+1]='\0';
  }
  else {
    Copy[0]='\0';
    last=Lo1;
    temp=I1;

    while(temp!=NULL) {
      
      len1+=(temp->ag)-last-1;

      strncat(Copy,Data+last,(temp->ag)-last-1);

      last=temp->gt+1;
      temp=temp->leg;
    }
    len1+=Hi1-last+1;

    strncat(Copy,Data+last,Hi1-last+1);
    Copy[len1]='\0';

    for(i=len1;i>=1;i--) {

      switch(Copy[len1-i]) {
      case 'a':
      case 'A': Gene1[i]='t';break;
      case 't':
      case 'T': Gene1[i]='a';break;
      case 'c':
      case 'C': Gene1[i]='g';break;
      case 'g':
      case 'G': Gene1[i]='c';break;
      default: Gene1[i]='c';break;
      }
    }
    Gene1[len1+1]='\0';

  }

 if(F2>=0) {
    last=Lo2;
    temp=I2;
    while(temp!=NULL) {
      len2+=(temp->gt)-last;

      strncat(Gene2+1,Data+last,(temp->gt)-last);

      last=temp->ag+2;
      temp=temp->leg;
    }
    len2+=Hi2-last+1;

    strncat(Gene2+1,Data+last,Hi2-last+1);
    Gene2[len2+1]='\0';
  }
  else {
    Copy[0]='\0';
    last=Lo2;
    temp=I2;
   
    while(temp!=NULL) {
      len2+=(temp->ag)-last-1;

      strncat(Copy,Data+last,(temp->ag)-last-1);

      last=temp->gt+1;
      temp=temp->leg;
    }
    len2+=Hi2-last+1;

    strncat(Copy,Data+last,Hi2-last+1);
    Copy[len2]='\0';

    for(i=len2;i>=1;i--) {
      switch(Copy[len2-i]) {
      case 'a':
      case 'A': Gene2[i]='t';break;
      case 't':
      case 'T': Gene2[i]='a';break;
      case 'c':
      case 'C': Gene2[i]='g';break;
      case 'g':
      case 'G': Gene2[i]='c';break;
      default: Gene2[i]='c';break;
      }
    }
    Gene2[len2+1]='\0';
  }


 //printf("Gene1=%s\n",Gene1+1);
 //printf("Gene2=%s\n",Gene2+1);fflush(stdout); 


  S1=One_Score (Gene1, len1, MODEL_LEN,Context_Delta, Ch_Ct);
  S2=One_Score (Gene2, len2, MODEL_LEN,Context_Delta, Ch_Ct);

  *Score1=S1;
  *Score2=S2;

  //printf("Score for gene %ld-%ld with len=%ld = %f\n",Lo1,Hi1,len1,S1);
  //printf("Score for gene %ld-%ld with len=%ld = %f\n",Lo2,Hi2,len2,S2);

  if(S1<S2) return(1);
  return(0);
}


double  One_Score  (char X [], int T, int Model_Len,
                     String_Array * Delta [6], double Ch_Ct [ALPHABET_SIZE])


/* Set  Score  to the probabilites of string  X [1 .. T]  being
*  generated in each of the 3 forward and 3 reverse reading frames
*  using simple nonhomogeneous Markov models in  Delta []  with
*  model length equal to  Model_Len . */

  {
   double  Max, Min, Sum, S [7];
   int  i, Has_Stop [7];

   Fast_Evaluate (X, T, Model_Len, * Delta [0], * Delta [1], * Delta [2],
                    S [0]);
   Fast_Evaluate (X, T, Model_Len, * Delta [2], * Delta [0], * Delta [1],
                    S [1]);
   Fast_Evaluate (X, T, Model_Len, * Delta [1], * Delta [2], * Delta [0],
                    S [2]);
   Fast_Evaluate (X, T, Model_Len, * Delta [3], * Delta [4], * Delta [5],
                    S [3]);
   Fast_Evaluate (X, T, Model_Len, * Delta [5], * Delta [3], * Delta [4],
                    S [4]);
   Fast_Evaluate (X, T, Model_Len, * Delta [4], * Delta [5], * Delta [3],
                    S [5]);

   if  (Use_Independent)
       Indep_Eval (X, T, Ch_Ct, S [6]);
     else
       S [6] = MIN_LOG_PROB_FACTOR * T;

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

   // S[0]=exp(S[0]);
   S[0]=-S[0];

   return(S[0]);

  }

int Delayed(Gene_Ref * & Previous,long int i,long int *Hi,long int *Lo,intron *Cap,int Frame)
{
  long int stop1,stop2,hi1,lo1,hi2,lo2;
  intron *temp;
  int f1,f2;
  int ind;

  hi1=Previous[i].Hi;
  lo1=Previous[i].Lo;
  hi2=*Hi;
  lo2=*Lo;
  f1=Previous[i].Frame;
  f2=Frame;



  // find how far we can delay the ATG
  temp=Previous[i].Intron_List;
  if(f1>=0) {
    if(temp==NULL) stop1=hi1-MIN_EXON_LENGTH;
    else stop1=(temp->gt)-1-MIN_INTRON_LENGTH;
  }
  else {
    if(temp==NULL) stop1=lo1+MIN_INTRON_LENGTH;
    else {
      while(temp->leg!=NULL) temp=temp->leg;
      stop1=temp->gt+1+MIN_EXON_LENGTH;
    }
  }

  temp=Cap;
  if(f2>=0) {
    if(temp==NULL) stop2=hi2-MIN_EXON_LENGTH;
    else stop2=(temp->gt)-1-MIN_INTRON_LENGTH;
  }
  else {
    if(temp==NULL) stop2=lo2+MIN_INTRON_LENGTH;
    else {
      while(temp->leg!=NULL) temp=temp->leg;
      stop2=temp->gt+1+MIN_EXON_LENGTH;
    }
  }


  //printf("Try to delay:\n");
  //printf("Gene1: %ld %ld %d stop=%ld\n",lo1,hi1,f1,stop1);
  //printf("Gene2: %ld %ld %d stop=%ld\n",lo2,hi2,f2,stop2);

      
  // case: A: -------->
  //       B:   ----------> 
  // delay B

  if(f1>=0 && f2 >= 0) {
    // gene1 is A and gene2 is B
    if(lo1<=lo2 && hi2-hi1>Min_Gene_Len) {
      while(lo2<=hi1-Min_Olap) lo2+=3;

      while(!StartCod(lo2+2)&&lo2<=hi2) lo2+=3;
      if(lo2<=stop2 && hi2-lo2+1>=Min_Gene_Len) {
	*Lo=lo2;
	return(1);
      }
      return(0);
    }
    // gene1 is B and gene2 is A
    if(lo2<=lo1 && hi1-hi2>Min_Gene_Len) {
      while(lo1<=hi2-Min_Olap) lo1+=3;
      while(!StartCod(lo1+2)&&lo1<=hi1) lo1+=3;
      if(lo1<=stop1 && hi1-lo1+1>=Min_Gene_Len) {
	Previous[i].Lo=lo1;
	return(1);
      }
      return(0);
    }
    return(0);
  }


  // case: A: --------->
  //       B:   <-----------
  // can't delay


  // case: A: <--------
  //       B:    -------->
  // delay both

  // gene1 is A and gene2 is B
  if(f1<0 && f2>=0){
    if(lo1<lo2 && hi2-hi1>Min_Gene_Len) {
      ind=0;
      while(!ind) {
	// try to delay gene1
	hi1-=3;
	while(!StartCod(-hi1) && hi1>=lo1) hi1-=3;
	if(StartCod(-hi1)) {
	  if(hi1>=stop1 && hi1-lo2<=Min_Olap&&hi1-lo1+1>=Min_Gene_Len) {
	    ind=1;
	    Previous[i].Hi=hi1;
	    return(1);
	  }
	  else { // try to delay gene2
	    lo2+=3;

	    while(!StartCod(lo2+2)&&lo2<=hi2) lo2+=3;
	    if(StartCod(lo2+2)) {
	      if(lo2<=stop2 && hi1-lo2<=Min_Olap && hi2-lo2+1>=Min_Gene_Len) {
		ind=1;
		*Lo=lo2;
		return(1);
	      }
	    }
	    else return(0);
	  }
	}
	else return(0);
      }
    }
    return(0);
  }

  // gene1 is B and gene2 is A
  if(f2<0 && f1>=0){
    if(lo2<lo1 && hi1-hi2>Min_Gene_Len) {

      ind=0;
      while(!ind) {

	// try to delay gene2
	hi2-=3;
	while(!StartCod(-hi2) && hi2>=lo2) hi2-=3;
	if(StartCod(-hi2)) {
	  if(hi2>=stop2 && hi2-lo1<=Min_Olap && hi2-lo2+1>=Min_Gene_Len) {
	    ind=1;
	    *Hi=hi2;
	    return(1);
	  }
	  else { // try to delay gene1
	    lo1+=3;
	    while(!StartCod(lo1+2)&&lo1<=hi1) lo1+=3;
	    if(StartCod(lo1+2)) {
	      if(lo1<=stop1 && hi2-lo1<=Min_Olap && hi1-lo1>=Min_Gene_Len) {
		ind=1;
		Previous[i].Lo=lo1;
		return(1);
	      }
	    }
	    else return(0);
	  }
	}
	else return(0);
      }
    }
    return(0);
  }

  // case: A <-------
  //       B:      <--------
  if(f1<0 && f2<0) {
    // case gene1 is A and gene2 is B
    if(lo1<lo2 && lo2<=hi1 && hi2-hi1+1>=Min_Gene_Len){
      while(hi1-lo2>Min_Olap) hi1-=3;
      while(hi1>=lo1 && !StartCod(-hi1)) hi1-=3;
      if(StartCod(-hi1) && hi1>=stop1){
	Previous[i].Hi=hi1;
	return(1);
      }
      return(0);
    }
    // case gene2 is A and gene1 is B
    if(lo2<lo1 && lo1<=hi2 && hi1-hi2+1>=Min_Gene_Len){
      while(hi2-lo1>Min_Olap) hi2-=3;
      while(hi2>=lo2 && !StartCod(-hi2)) hi2-=3;
      if(StartCod(-hi2) && hi2>=stop1){
	*Hi=hi2;
	return(1);
      }
      return(0);
    }
  }
  return(0);

}


void  Check_Previous  (Gene_Ref * & Previous, long int *ID, long int Lo, 
		       long int Hi, int Frame,intron *Cap, int Score[7])

//  Add reference for new gene  ID  at positions  Lo .. Hi  in
//  frame  Frame  to  Previous .  Then check for overlaps with other
//  genes already in  Previous  and score them and set their status.  
//  Overlaps must be
//  at least  Min_Olap  long to be reported.   Lo  and  Hi  may be
//  greater than  Data_Len .  Frame is +1, +2 or +3 for forward genes,
//  and -1, -2 or -3 for reverse genes.

  {
   long int  i, Curr_Len;
   intron *temp;
   int shadow,last,ret;
   unsigned int status;
   char reject;
   int l1,l2;

   status=OK;
   reject=' ';
  
   // see if the current gene shadow other gene or is shadowed by one

   shadow=0;last=0;

#ifdef MESSAGE
    temp=Cap;
    printf("Check gene:\n");fflush(stdout);
    printgene(Hi,Lo,Frame,temp);fflush(stdout);
#endif MESSAGE

   temp=Cap;

   if(Perc_OK(Frame,Lo,Hi,temp)) {  // if the percentages are ok!

#ifdef MESSAGE
      printf("Perc_ok!\n");fflush(stdout);
#endif MESSAGE

     //printf("Check gene %ld: %ld-%ld\n",(*ID)+1,Lo,Hi); fflush(stdout);

     // keep the longest gene
     for(i=*ID;i>0;i--) 
       if(Previous[i].Reject_Code!='R'){
	 
	 if(Simple_Overlap(Previous[i].Lo,Previous[i].Hi,Lo,Hi)) 
	   // genes overlap by more then MIN_OVERLAP_LENGTH
	
	   if(Overlap(Previous[i].Lo,Previous[i].Hi,Previous[i].Intron_List,Previous[i].Frame,Lo,Hi,Cap,Frame)){

#ifdef MESSAGE
	     printoverlap(Previous[i],Hi,Lo,Frame,Cap);
#endif MESSAGE

	     //if(!Delayed(Previous,i,&Hi,&Lo,Cap,Frame)) {
	     
	     shadow=1;		
	   
#ifdef MESSAGE  
	     //printf("ret=%d\n",ret);
	     printf("Compare length:\n");
#endif MESSAGE

	     ret=Compare_length(Previous[i].Intron_List,Previous[i].Lo,Previous[i].Hi,Previous[i].Frame,Cap,Lo,Hi,Frame);
	     
#ifdef MESSAGE
	     printf("ret=%d\n",ret);
#endif MESSAGE

	     // keep genes with introns

	     // variant 1
	     //if(ret==0 && Previous[i].Intron_List==NULL && Cap!=NULL) ret=3;
	     //if(ret==1 && Previous[i].Intron_List!=NULL && Cap==NULL) ret=2;

	     // variant 2
	     /* temp=Previous[i].Intron_List;
		l1=0;
		while(temp!=NULL) { l1++; temp=temp->leg;}
	      
		temp=Cap;
		l2=0;
		while(temp!=NULL) { l2++; temp=temp->leg;}
	      
		if(ret==0 && l1<l2) ret=3;
		if(ret==1 && l1>l2) ret=2;*/

	     
	     //printf("ret=%d\n",ret);

	     if(ret==0) break;


	     if(ret==1) {
	       // have to replace the shortest gene
	     
	       //if(CmpGenes(Previous[i].Lo,Previous[i].Hi,Previous[i].Intron_List,
	       //Previous[i].Frame,Lo,Hi,Cap,Frame,&S1,&S2))  {
	       // if second gene is better replace the oldest one
	       
	       if(!last) {
		 replace(Previous,i,Lo,Hi,Frame,Cap,*ID,Score);

#ifdef MESSAGE
		 printf("Replaced gene %ld: %ld-%ld\n",i,Previous[i].Lo,Previous[i].Hi); fflush(stdout);
#endif MESSAGE		
 
		 last=1;
	       }
	       else {
		 Previous[i].Reject_Code='R';

#ifdef MESSAGE
		 printf("Rejected gene %ld: %ld-%ld\n",i,Previous[i].Lo,Previous[i].Hi);fflush(stdout);
#endif MESSAGE	      
 

	       }
	     }
	     else if(ret==2) { 
	       // difference of length is very short, gene1 is shorter
	       Previous[i].Reject_Code='O';
	       Previous[i].Status=(*ID)+1;
	       reject='o';
	       if(!last) shadow=0;
	       
#ifdef MESSAGE
	       printf("Gene %ld is a little shorter\n",i);
#endif MESSAGE	     

	     }
	     
	     else if(ret==3) {
	       // second gene is shorter
	       
	       status=i;
	       shadow=0;
	       reject='O';
	       Previous[i].Reject_Code='o';

#ifdef MESSAGE
	       printf("Gene %ld is a little shorter\n",(*ID)+1);
#endif MESSAGE

	     }
	   	 
	     /* { else {
		printf("Delayed with gene %ld\n",i);
	      
		}*/
	   }
       }
   }
   else shadow=1;

   if(!shadow) {

#ifdef MESSAGE
      printf("Adding gene %ld: %ld-%ld\n",(*ID)+1,Lo,Hi); fflush(stdout);
#endif MESSAGE

     ++(*ID);

     Previous = (Gene_Ref *) Safe_realloc (Previous,
					   (1 + *ID) * sizeof (Gene_Ref));


     Previous [*ID] . Lo = Lo;
     Previous [*ID] . Hi = Hi;
     Previous [*ID] . Frame = Frame;
     Previous [*ID] . Status = status;
     Previous [*ID] . Reject_Code = reject;
     //Previous [*ID] . Problem_List = NULL;
     if  (*ID == 1 || Previous [*ID - 1] . Hi < Hi)
       Previous [*ID] . Max_Hi = Hi;
     else
       Previous [*ID] . Max_Hi = Previous [*ID - 1] . Hi;

     for(i=0;i<7;i++) Previous[*ID].Score[i]=Score[i];

     Curr_Len = 1 + Hi - Lo;


     if(Cap!=NULL) {
       Previous[*ID].Intron_List=(intron *)malloc(sizeof(intron));  
       if (Previous[*ID].Intron_List == NULL) {
	 fprintf(stderr,"Memory allocation for node failure.\n"); 
	 abort();
       }
       temp=Previous[*ID].Intron_List;
       temp->gt=Cap->gt;
       temp->ag=Cap->ag;
       temp->leg=NULL;

       while(Cap->leg!=NULL) {
	 temp->leg=(intron *)malloc(sizeof(intron));  
	 if (temp == NULL) {
	   fprintf(stderr,"Memory allocation for node failure.\n"); 
	   abort();
	 }
	 temp=temp->leg;
	 Cap=Cap->leg;
	 temp->gt=Cap->gt;
	 temp->ag=Cap->ag;
	 temp->leg=NULL;

       }
     }
     else Previous[*ID].Intron_List=NULL;

   }

   return;

  }



int  Choose_Score  (int S [7], int F)

//  Return the score in  S  corresponding to frame  F .

  {
   switch  (F)
     {
      case  1 :
        return  S [0];
      case  2 :
        return  S [1];
      case  3 :
        return  S [2];
      case  -1 :
        return  S [3];
      case  -2 :
        return  S [5];
      case  -3 :
        return  S [4];
      default :
        fprintf (stderr, "ERROR:  Bad frame value  %d  in  Choose_Score ()\n", F);
     }

   return  S [0];
  }



long int  Choose_Start  (long int Begin, long int Len)

  /* Return the position of the first base of the most likely start
   *  codon for the gene whose first start codon is at base  Begin
   *  and has length  Len , which is positive in the forward direction
   *  and negative in the reverse complement direction. */


  /* modify it such that it should choose the most likely start codon  */
   


  {
   double  Score;
   long int  Original_Start;
   char  Buffer [1 + UPSTREAM_LEN], Codon [4];
   int B[20],i,k;

   if  (Choose_First_Start_Codon)
       return  Begin;

   Original_Start = Begin;

   if  (Len > 0)
       {
        Transfer (Codon, Begin, 3,Data,Data_Len);
        do
          {
           Transfer (Buffer, Begin - UPSTREAM_LEN - UPSTREAM_OFFSET, UPSTREAM_LEN,Data,Data_Len);

	   // Here the score shoud be computed using an atg score
	   k=0;
	   for(i=Begin-12;i<Begin+6;i++) {
	     switch(Data[i]) {
	     case 'a': B[k++]=0;break; 
	     case 'c': B[k++]=1;break; 
	     case 'g': B[k++]=2;break; 
	     case 't': B[k++]=3;break; 
	     default: B[k++]=0;break; 
	     }
	   }
	   Is_Atg(B,&Score);

	   if  (Score >= 0)
                    return  Begin;
           
           do
             {
              Begin += 3;
              Len -= 3;
              if  (Begin > Data_Len)
                  Begin -= Data_Len;
              Transfer (Codon, Begin, 3,Data,Data_Len);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Gene_Len);
          }  while  (! Is_Stop (Codon) && Len > Min_Gene_Len);
       }
     else
       {
        Transfer (Codon, Begin, -3,Data,Data_Len);
        do
          {
           Transfer (Buffer, Begin + UPSTREAM_LEN + UPSTREAM_OFFSET, - UPSTREAM_LEN,Data,Data_Len);
	   // Here the score shoud be computed using an atg score
	   k=0;
	   for(i=Begin+12;i>Begin-6;i--) {
	     switch(Data[i]) {
	     case 'a': B[k++]=3;break; 
	     case 'c': B[k++]=2;break; 
	     case 'g': B[k++]=1;break; 
	     case 't': B[k++]=0;break; 
	     default: B[k++]=0;break; 
	     }
	   }
	   Is_Atg(B,&Score);

	   if  (Score >= ATG_THRESHOLD)
                    return  Begin;
           
           do
             {
              Begin -= 3;
              Len -= 3;
              if  (Begin < 1)
                  Begin += Data_Len;
              Transfer (Codon, Begin, -3,Data,Data_Len);
             }  while  (! Is_Start (Codon) && ! Is_Stop (Codon)
                               && Len > Min_Gene_Len);
          }  while  (! Is_Stop (Codon) && Len > Min_Gene_Len);
       }

   return  Original_Start;
  }



int  Cmp  (const void * A, const void * B)

/* Comparison function for  qsort  to sort  Gene_Refs  by
*  descending length. */

  {
    //long int  A_Len, B_Len;
   int ret;

   if ( Previous [* (long int *)A] . Hi < Previous [* (long int *)B] . Lo) return(-1);
   if ( Previous [* (long int *)B] . Hi < Previous [* (long int *)A] . Lo) return(1);

   ret=Compare_length(Previous[* (long int *)B].Intron_List,Previous[* (long int *)B].Lo,Previous[* (long int *)B].Hi,Previous[* (long int *)B].Frame,Previous[* (long int *)A].Intron_List,Previous[* (long int *)A].Lo,Previous[* (long int *)A].Hi,Previous[* (long int *)A].Frame);

   if(ret==0 || ret==3) return(1);
   if(ret==1 || ret==2) return(-1);

   /*   A_Len=1+Previous [* (long int *)A] . Hi-Previous [* (long int *)A] . Lo;
   B_Len=1+Previous [* (long int *)B] . Hi-Previous [* (long int *)B] . Lo;

   if((* ((Gene_Ref * *) A)) -> Hi < (* ((Gene_Ref * *) B)) -> Lo) return(1);
   if((* ((Gene_Ref * *) B)) -> Hi < (* ((Gene_Ref * *) A)) -> Lo) return(-1);

   A_Len = 1 + (* ((Gene_Ref * *) A)) -> Hi - (* ((Gene_Ref * *) A)) -> Lo;
   B_Len = 1 + (* ((Gene_Ref * *) B)) -> Hi - (* ((Gene_Ref * *) B)) -> Lo;

   if  (A_Len > B_Len)
       return  -1;
   else if  (A_Len < B_Len)
       return  1;
     else
     return  0;*/
  }




double  Doublet_Score  (char A, char B, char P, char Q)

/* Return the energy released by consecutive bases  AB  binding
*  to pair  PQ  where  PQ  is in  5'-3' order.   Values from
*  Lewin's Genes V, p. 115. */

  {
   B = tolower (B);
   P = tolower (P);
   Q = tolower (Q);

   switch  (tolower (A))
     {
      case  'a' :
        switch  (B)
          {
           case  'a' :
             if  (P == 't' && Q == 't')
                 return  0.9;
             goto  Error;
           case  'c' :
             if  (P == 't' && Q == 'g')
                 return  1.8;
             goto  Error;
           case  'g' :
             if  (P == 't' && Q == 'c')
                 return  2.3;
             else if  (P == 't' && Q == 't')
                 return  1.15;                    /* Guess */
             goto  Error;
           case  't' :
             if  (P == 't' && Q == 'a')
                 return  1.1;
             else if  (P == 't' && Q == 'g')
                 return  0.65;                    /* Guess */
             goto  Error;
          }

      case  'c' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'g' && Q == 't')
                 return  2.1;
             goto  Error;
           case  'c' :
             if  (P == 'g' && Q == 'g')
                 return  2.9;
             goto  Error;
           case  'g' :
             if  (P == 'g' && Q == 'c')
                 return  3.4;
             else if  (P == 'g' && Q == 't')
                 return  1.50;                    /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'g' && Q == 'a')
                 return  2.3;
             else if  (P == 'g' && Q == 'g')
                 return  1.15;                    /* Guess */
             goto  Error;
          }

      case  'g' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'c' && Q == 't')
                 return  1.7;
             else if  (P == 't' && Q == 't')
                 return  0.85;                    /* Guess */
             goto  Error;
           case  'c' :
             if  (P == 'c' && Q == 'g')
                 return  2.0;
             else if  (P == 't' && Q == 'g')
                 return  1.00;                    /* Guess */
             goto  Error;
           case  'g' :
             if  (P == 'c' && Q == 'c')
                 return  2.9;
             else if  (P == 't' && Q == 'c')
                 return  1.45;                    /* Guess */
             else if  (P == 'c' && Q == 't')
                 return  1.45;                    /* Guess */
             else if  (P == 't' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'c' && Q == 'a')
                 return  1.8;
             else if  (P == 't' && Q == 'a')
                 return  0.9;                     /* Guess */
             else if  (P == 'c' && Q == 'g')
                 return  0.9;                     /* Guess */
             else if  (P == 't' && Q == 'g')
                 return  0.5;                     /* Guess */
             goto  Error;
          }

      case  't' :
        switch  (B)
          {
           case  'a' :
             if  (P == 'a' && Q == 't')
                 return  0.9;
             else if  (P == 'g' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  'c' :
             if  (P == 'a' && Q == 'g')
                 return  1.7;
             else if  (P == 'g' && Q == 'g')
                 return  0.85;                    /* Guess */
             goto  Error;
           case  'g' :
             if  (P == 'a' && Q == 'c')
                 return  2.1;
             else if  (P == 'g' && Q == 'c')
                 return  1.05;                    /* Guess */
             else if  (P == 'a' && Q == 't')
                 return  1.05;                    /* Guess */
             else if  (P == 'g' && Q == 't')
                 return  0.5;                     /* Guess */
             goto  Error;
           case  't' :
             if  (P == 'a' && Q == 'a')
                 return  0.9;
             else if  (P == 'g' && Q == 'a')
                 return  0.5;                     /* Guess */
             else if  (P == 'a' && Q == 'g')
                 return  0.5;                     /* Guess */
             else if  (P == 'g' && Q == 'g')
                 return  0.5;                     /* Guess */
             goto  Error;
          }

     }
   
  Error:
   fprintf (stderr, "ERROR:  Bad doublet pair  %c%c  and  %c%c\n",
              A, B, P, Q);

   return  0;
  }



double  Edit_Distance  (const char * P, const char * T)

/* Find and return the highest enery match of string  P [1 ...]  within
*  string  T [1 ...] . */

  {
   double  Max, Best, Final_Score, X;
   int  a, b, i, j, M, N, Len, Len1, First = 1, Max_i;

   M = strlen (P + 1);
   N = strlen (T + 1);

   for  (i = 1;  i <= M;  i ++)
     {
      ED_Score [i] [0] . Free_i = 0.0;
      ED_Score [i] [0] . Free_j = BIG_NEGATIVE;
      ED_Score [i] [0] . Both_Free = BIG_NEGATIVE;
      ED_Score [i] [0] . Match = BIG_NEGATIVE;
      ED_Score [i] [0] . Free_i_Len = i;
      ED_Score [i] [0] . Free_j_Len = 0;
      ED_Score [i] [0] . Both_i_Len = 0;
      ED_Score [i] [0] . Both_j_Len = 0;
     }
   for  (j = 1;  j <= N;  j ++)
     {
      ED_Score [0] [j] . Free_i = BIG_NEGATIVE;
      ED_Score [0] [j] . Free_j = 0.0;
      ED_Score [0] [j] . Both_Free = BIG_NEGATIVE;
      ED_Score [0] [j] . Match = BIG_NEGATIVE;
      ED_Score [0] [j] . Free_i_Len = 0;
      ED_Score [0] [j] . Free_j_Len = j;
      ED_Score [0] [j] . Both_i_Len = 0;
      ED_Score [0] [j] . Both_j_Len = 0;
     }
   ED_Score [0] [0] . Free_i = 0.0;
   ED_Score [0] [0] . Free_j = 0.0;
   ED_Score [0] [0] . Both_Free = 0.0;
   ED_Score [0] [0] . Match = BIG_NEGATIVE;
   ED_Score [0] [0] . Free_i_Len = 0;
   ED_Score [0] [0] . Free_j_Len = 0;
   ED_Score [0] [0] . Both_i_Len = 0;
   ED_Score [0] [0] . Both_j_Len = 0;

   Final_Score = BIG_NEGATIVE;
   for  (j = 1;  j <= N;  j ++)
     {
      Max = BIG_NEGATIVE;
      for  (i = 1;  i <= M;  i ++)
        {
         Best = BIG_NEGATIVE;
         Len = 0;
         X = ED_Score [i - 1] [j] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
             }
         X = ED_Score [i - 1] [j] . Free_i;
         a = ED_Score [i - 1] [j] . Free_i_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = 1 + a;
             }
         ED_Score [i] [j] . Free_i = Best;
         ED_Score [i] [j] . Free_i_Len = Len;
         if  (Best > Max)
             {
              Max = Best;
              Max_i = i;
             }
         
         Best = BIG_NEGATIVE;
         Len = 0;
         X = ED_Score [i] [j - 1] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
             }
         X = ED_Score [i] [j - 1] . Free_j;
         a = ED_Score [i] [j - 1] . Free_j_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = 1 + a;
             }
         ED_Score [i] [j] . Free_j = Best;
         if  (Best > BIG_NEGATIVE)
             ED_Score [i] [j] . Free_j_Len = Len;
           else
             ED_Score [i] [j] . Free_j_Len = 0;

         Best = BIG_NEGATIVE;
         Len = Len1 = 0;
         X = ED_Score [i - 1] [j - 1] . Match;
         if  (X > Best)
             {
              Best = X;
              Len = Len1 = 1;
             }
         X = ED_Score [i - 1] [j] . Free_j;
         if  (X > Best)
             {
              Best = X;
              Len = 1;
              Len1 = ED_Score [i - 1] [j] . Free_j_Len;
             }
         X = ED_Score [i] [j - 1] . Free_i;
         if  (X > Best)
             {
              Best = X;
              Len = ED_Score [i] [j - 1] . Free_i_Len;
              Len1 = 1;
             }
         X = ED_Score [i - 1] [j] . Both_Free;
         a = ED_Score [i - 1] [j] . Both_i_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = a + 1;
              Len1 = ED_Score [i - 1] [j] . Both_j_Len;
             }
         X = ED_Score [i] [j - 1] . Both_Free;
         a = ED_Score [i] [j - 1] . Both_j_Len;
         if  (X > Best && a < MAX_FREE_LEN)
             {
              Best = X;
              Len = ED_Score [i] [j - 1] . Both_i_Len;
              Len1 = a + 1;
             }
         X = ED_Score [i - 1] [j - 1] . Both_Free;
         a = ED_Score [i - 1] [j - 1] . Both_i_Len;
         b = ED_Score [i - 1] [j - 1] . Both_j_Len;
         if  (X > Best && a < MAX_FREE_LEN && b < MAX_FREE_LEN)
             {
              Best = X;
              Len = a + 1;
              Len1 = b + 1;
             }
         ED_Score [i] [j] . Both_Free = Best;
         ED_Score [i] [j] . Both_i_Len = Len;
         ED_Score [i] [j] . Both_j_Len = Len1;

         if  (! Match (P [i], T [j]))
             ED_Score [i] [j] . Match = BIG_NEGATIVE;
           else
             {
              Best = 0.0;
              X = ED_Score [i - 1] [j - 1] . Free_i
                      + Bulge_Cost (ED_Score [i - 1] [j - 1] . Free_i_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Free_j
                      + Bulge_Cost (ED_Score [i - 1] [j - 1] . Free_j_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Both_Free
                      + Loop_Cost (ED_Score [i - 1] [j - 1] . Both_i_Len,
                                      ED_Score [i - 1] [j - 1] . Both_j_Len);
              if  (X > Best)
                  Best = X;
              X = ED_Score [i - 1] [j - 1] . Match;
              if  (X != BIG_NEGATIVE)
                  {
                   X += Doublet_Score (P [i - 1], P [i], T [j - 1], T [j]);
                   if  (X > Best)
                       Best = X;
                  }
              ED_Score [i] [j] . Match = Best;
              if  (Best > Max)
                  {
                   Max = Best;
                   Max_i = i;
                  }
             }
        }
      if  (Max > Final_Score)
          Final_Score = Max;
     }

   return  Final_Score;
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



int  Is_Forward_Start  (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  forward direction.

  {
    // initial version:
    // const unsigned  Mask = ATG_MASK | GTG_MASK | TTG_MASK;

    // my version:
    const unsigned  Mask = ATG_MASK;

   return  ((Codon & Mask) == Codon);
  }



int  Is_Forward_Stop  (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  forward direction.

  {
   return  ((Codon & TAR_MASK) == Codon || (Codon & TRA_MASK) == Codon);
  }



int  Is_Reverse_Start  (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a start codon in the
//  reverse direction.

  {
    //   const unsigned  Mask = CAT_MASK | CAC_MASK | CAA_MASK;

    // my version
    const unsigned  Mask = CAT_MASK;

    return  ((Codon & Mask) == Codon);
  }



int  Is_Reverse_Stop  (unsigned Codon)

//  Return  TRUE  iff bit pattern  Codon  represents a stop codon in the
//  reverse direction.

  {
   return  ((Codon & YTA_MASK) == Codon || (Codon & TYA_MASK) == Codon);
  }



int  Is_Start  (char * S)

/* Return  TRUE  iff  S  is a start codon. */

  {
   return  (strncmp (S, "atg", 3) == 0);
	    //              || strncmp (S, "ctg", 3) == 0
	    //              || strncmp (S, "gtg", 3) == 0
	    //|| strncmp (S, "ttg", 3) == 0);
  }



int  Is_Stop  (char * S)

/* Return  TRUE  iff  S  is a stop codon. */

  {
   return  (strncmp (S, "taa", 3) == 0
              || strncmp (S, "tag", 3) == 0
              || strncmp (S, "tga", 3) == 0);
  }



double  Loop_Cost  (int M, int N)

/* Return the energy cost of a loop with  M  bases on one side and
*  N  bases on the other. */

  {
   double  Cost;
   int  Min;

   if  (M <= 0 || N <= 0)
       return  BIG_NEGATIVE;

   if  (M < N)
       Min = M;
     else
       Min = N;

   if  (Min < 4)
       Cost = -0.8;
     else
       Cost = -0.8 - (Min - 3) * (8.4 - 0.8) / 27.0;

   if  (M == N)
       return  Cost;

   return  Cost - abs (M - N) * 1.1;
  }



int  Match  (char P, char Q)

/* Return true iff bases  P  and  Q  bind together. */

  {
   Q = tolower (Q);
   switch  (tolower (P))
     {
      case  'a' :
        return  (Q == 't');
      case  'c' :
        return  (Q == 'g');
      case  'g' :
        return  (Q == 'c' || Q == 't');
      case  't' :
        return  (Q == 'a' || Q == 'g');
     }

   return  0;
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
//    -f     Use ribosome-binding energy to choose start codon
//    +f     Use first codon in orf as start codon
//    -g n   Set minimum gene length to n
//    -o n   Set minimum overlap length to n.  Overlaps shorter than this
//           are ignored.
//    -p n   Set minimum overlap percentage to n%.  Overlaps shorter than
//           this percentage of *both* strings are ignored.
//    -r     Don't use independent probability score column
//    +r     Use independent probability score column
//    -s s   Use string s as the ribosome binding pattern to find start codons.
//    -t n   Set threshold score for calling as gene to n.  If the in-frame
//           score >= n, then the region is given a number and considered
//           a potential gene.

  {
   char  * P;
   long int  W;
   double  X;
   int  i, L;

   for  (i = 3;  i < argc;  i ++)
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
              case  'f' :       // use function to choose start codon in gene
                Choose_First_Start_Codon = FALSE;
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
              case  'o' :       // minimum overlap length
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap length %s\n", P);
                  else
                    Min_Olap = W;
                assert (Min_Olap > 0);
                break;
              case  'p' :       // minimum overlap percent
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                X = strtod (P, NULL);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad minimum overlap percent %s\n", P);
                  else
                    Min_Olap_Percent = X / 100.0;
                assert (Min_Olap_Percent > 0.0 && Min_Olap_Percent < 100.0);
                break;
              case  'r' :       // don't use random/independent score column
                Use_Independent = FALSE;
                break;
              case  't' :       // threshold score for calling as gene
                errno = 0;
                if  (argv [i] [2] != '\0')
                    P = argv [i] + 2;
                  else
                    P = argv [++ i];
                W = strtol (P, NULL, 10);
                if  (errno == ERANGE)
                    fprintf (stderr, "ERROR:  Bad threshold score %s\n", P);
                  else
                    Threshold_Score = W;
                assert (Threshold_Score > 0 && Threshold_Score < 100);
                break;
              default :
                fprintf (stderr, "Unrecognized option %s\n", argv [i]);
             }
           break;
         case  '+' :
           switch  (argv [i] [1])
             {
              case  'f' :       // automatically use first start codon in gene
                Choose_First_Start_Codon = TRUE;
                break;
              case  'r' :       // use random/independent score column
                Use_Independent = TRUE;
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
                    String_Array * Delta [6], double Ch_Ct [ALPHABET_SIZE],
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
   long int  i, Len;


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

   Simple_Score (Orf_Buffer, Len, MODEL_LEN, Delta, Ch_Ct, Score);

   return;
  }



void  Simple_Score  (char X [], int T, int Model_Len,
                     String_Array * Delta [6], double Ch_Ct [ALPHABET_SIZE],
                     int Score [])

/* Set  Score  to the probabilites of string  X [1 .. T]  being
*  generated in each of the 3 forward and 3 reverse reading frames
*  using simple nonhomogeneous Markov models in  Delta []  with
*  model length equal to  Model_Len . */

  {

   double  Max, Min, Sum, S [7];
   int  i, Has_Stop [7];

   Fast_Evaluate (X, T, Model_Len, * Delta [0], * Delta [1], * Delta [2],
                    S [0]);
   Fast_Evaluate (X, T, Model_Len, * Delta [2], * Delta [0], * Delta [1],
                    S [1]);
   Fast_Evaluate (X, T, Model_Len, * Delta [1], * Delta [2], * Delta [0],
                    S [2]);
   Fast_Evaluate (X, T, Model_Len, * Delta [3], * Delta [4], * Delta [5],
                    S [3]);
   Fast_Evaluate (X, T, Model_Len, * Delta [5], * Delta [3], * Delta [4],
                    S [4]);
   Fast_Evaluate (X, T, Model_Len, * Delta [4], * Delta [5], * Delta [3],
                    S [5]);

   if  (Use_Independent)
       Indep_Eval (X, T, Ch_Ct, S [6]);
     else
       S [6] = MIN_LOG_PROB_FACTOR * T;

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

