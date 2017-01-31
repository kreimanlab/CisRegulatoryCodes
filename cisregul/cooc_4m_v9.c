/* cooc_4m_v9.c 
 * compute co-occurrences of 4 motifs 
 */

/* new in version 9
 * 03-02-2004: call new versions of programs such that memory is allocated here and not in the sub-functions
 */

/* new in version 8
 * 02-20-2004: call include files cooc_dist_lib.c and cooc_file_lib.c
 * 02-20-2004: call dist_2v_allt_v3
 * 02-20-2004: changes in cooc_dist_lib.c allow order_constraint to be -1 (in which case strand information is not used)
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "../lib/nr.c"
#include "../lib/my_math_package.c"
#include "../lib/vec_methods.c"
#include "../lib/my_text_package.c"
#include "../lib/my_file_io.c"
#include "cooc_dist_lib.c"
#include "cooc_file_lib.c"

int max_n_ovdi12=1000000;       // maximum number of entries for ovdi12
int max_n_ovdi123=1000000;      // maximum number of entries for ovdi123
int max_n_ovdi1234=1000000;     // maximum number of entries for ovdi1234
int max_n_ovdi_single_transcript=100000; // maximum number of co-occurrences for a single transcript
int max_transcript_occurrences=50000;    // maximum number of occurrences per transcript
int max_occur=100000;           // maximum total number of occurrences for a given motif
int max_dist=100;               // maximum allowed distance between motifs [this should be a parameter!]
int min_dist;                   // this depends on the motif lengths
float min_dist_factor=0.5;      // regulates amount of overlap between factors (any overlap allowed if min_dist_factor=0, no overlap if min_dist_factor>=1)
FILE *log_file;                 // log file
int verbose=0;                  // debug mode, 1 to print out the transcript positions, 2 to print debug info, 3 for full debug info
int n_motifs;                   // number of motifs
float default_min=0.05;         // if there are no occurrences in the background, assume this minimum
float r_transcripts_cl_bck;     // n_transcripts_bck/n_transcripts_cl
int n_transcripts_bck;          // number of transcripts in background set
int n_transcripts_cl;           // number of transcripts in cluster set
int errorcode;                  // >0 if there was a memory overflow error during the processing of the current combination
int prog_version=9;             // program version
char last_modified[200];        // date of latest modification of code
char file_label[200];           // used in the output file names (to be able to run multiple copies in parallel)
int filteron=0;                 // 1 to filter the positions
int pos_format=1;               // 1 for positions with respect to TSS and -1 for positions with respect to 5' end
int debug_mode=0;

/* total memory requirements 
 * data_t_cl    n_motifs x max_occur   4 x 300 x 2 x 10^5   =  240 x 10^6 (for 300 motifs)
 * data_p_cl    n_motifs x max_occur   
 * data_t_bck   n_motifs x max_occur
 * data_p_bck   n_motifs x max_occur
 * data_s_cl    n_motifs x max_occur   2 x 300 x 2 x 10^5   =  120 x 10^6 (for 300 motifs)   
 * data_s_bck   n_motifs x max_occur

 * ovdi12_cl    max_n_ovdi12 x 7       2 x 7 x 2 x 10^6     = 28 x 10^6
 * ovdi12_bck   max_n_ovdi12 x 7       
 * ovdi123_cl   max_n_ovdi123 x 10     2 x 10 x 2 x 10^6    = 40 x 10^6      
 * ovdi123_bck  max_n_ovdi123 x 10     
 * ovdi1234_cl  max_n_ovdi1234 x 13    2 x 13 x 2 x 10^6    = 52 x 10^6

 * *m1_t_cl;    max_n_occur x 1        16 x 2 x 10^5        = 3.2 x 10^6
 * *m2_t_cl;
 * *m3_t_cl;
 * *m4_t_cl;
 * *m1_p_cl;
 * *m2_p_cl;
 * *m3_p_cl;
 * *m4_p_cl;
 * *m1_t_bck;
 * *m2_t_bck;
 * *m3_t_bck;
 * *m4_t_bck;
 * *m1_p_bck;
 * *m2_p_bck;
 * *m3_p_bck;
 * *m4_p_bck;
 */

void compute_pr(int n_cl,int n_tr_cl,int n_bck,int n_tr_bck,float *r,double *p);
int compute_min_dist(int *ml,int n_motifs,float min_dist_factor);
void usage(void);

main(int argc,char *argv[])
{
  FILE *input_file;                  // used to read several input files
  FILE *output_file;                 // output file 
  FILE *positions_file;              // positions file
  char filename_pos_bck[200];        // filename where the occurrences of each motif are stored (background)
  char filename_pos_cl[200];         // filename where the occurrences of each motif are stored (cluster)
  char filename_strands_bck[200];    // filename where the strands are stored, background
  char filename_strands_cl[200];     // filename where the strands are stored, cluster
  char infotext[10000];
  char log_filename[200];            // log file name
  char output_filename[200];         // output file name
  char positions_filename[200];      // positions file name
  char temptext2[1000];
  char temptext[10000];
  char transcript_lengths_filename_bck[200];
  char transcript_lengths_filename_cl[200];
  double c1,c2,c3,c4;                    // number of possible combinations of 1/2/3/4 motifs
  //double cminian_p_tr_1m;
  double cminian_p_tr_2m;            // 1.0-minian_p_tr_2m (used in get_max_n_transcripts)
  double cminian_p_tr_3m;
  double cminian_p_tr_4m;
  //double fudge_threshold_1m=10.0;        // fudge factor for multiple comparisons
  double fudge_threshold_2m=10.0;
  double fudge_threshold_3m=10.0;
  double fudge_threshold_4m=10.0;
  //double minian_p_tr_1m;
  double minian_p_tr_2m;             // maximum p value for transcripts binomial computation for 2 motifs
  double minian_p_tr_3m;
  double minian_p_tr_4m;
  double p[2];                       // to retrieve results from compute_pr
  double p_bck_oc;                   // occurrences per total number of transcripts in the background set
  double p_bck_tr;                   // proportion of transcripts with co-occurring motifs in the background set
  double p_binom_oc;                 // binomial probability, with occurrences
  double p_binom_tr;                 // binomial probability, with transcripts
  float cumulative_memory_req=0.0;
  float memory_req;
  float minian_p_tr=0.01;            // probability threshold (experiment wise)
  float r[2];                        // to retrieve results from compute_pr
  float r_oc;                        // n123_cl/n123_bck (or minian_r_tr*10 if the denominator is 0)
  float r_tr;                        // n123_tr_cl/n123_tr_bck (or minian_r_tr*10 if the denominator is 0)
  int **data_p_bck;                  // n_motifs x max_occur, positions, background
  int **data_p_cl;                   // n_motifs x max_occur, positions, cluster
  int **data_t_bck;                  // n_motifs x max_occur, transcripts, background
  int **data_s_cl;                   // n_motifs x max_occur, strands, cluster
  int **data_s_bck;                  // n_motifs x max_occur, strands, bakcground
  int **data_t_cl;                   // n_motifs x max_occur, transcripts, cluster
  int **ovdi1234_bck;                // transcripts and positions, 4 motifs, background set
  int **ovdi1234_cl;                 // "                                  , cluster set
  int **ovdi123_bck;                 // "                        , 3 motifs, background set
  int **ovdi123_cl;
  int **ovdi12_bck;
  int **ovdi12_cl;
  int **transcript_lengths_bck;
  int **transcript_lengths_cl;
  int *m1_p_bck;
  int *m1_p_cl;
  int *m1_t_bck;
  int *m1_t_cl;
  int *m2_p_bck;
  int *m2_p_cl;
  int *m2_t_bck;
  int *m2_t_cl;
  int *m3_p_bck;
  int *m3_p_cl;
  int *m3_t_bck;
  int *m3_t_cl;
  int *m4_p_bck;
  int *m4_p_cl;
  int *m4_t_bck;
  int *m4_t_cl;
  char *m1_s_bck;
  char *m1_s_cl;
  char *m2_s_bck;
  char *m2_s_cl;
  char *m3_s_bck;
  char *m3_s_cl;
  char *m4_s_bck;
  char *m4_s_cl;
  int *max_n_transcripts_2m;
  int *max_n_transcripts_3m;
  int *max_n_transcripts_4m;
  int *motif_lengths;
  int curr_arg;                          // for parameter input
  int curr_tr;                       // temporary transcript number
  int finit_i;                       // mostly for debuggin, final conditions
  int finit_j;
  int finit_k;
  int finit_l;
  int i,j,k,ll,mm;
  int index2;                        // auxiliary index variable
  int index3;
  int init_i=1;                        // mostly for debugging, initial conditions
  int init_j=1;
  int init_k=1;
  int init_l=1;
  int itemp;
  int max_el=5000;
  int max_il=5000;
  int max_n_motifs=2;                // if max_n_motifs=2, process only pairs, if max_n_motifs=3, then process only triplets, if max_n_motifs=4, process quadruplets
  int max_n_transcripts;             // maximum number of transcripts so as to remain within threshold (use in call to cooc for background
  int max_ul=5000;
  int minian_tr=4;
  int ml[5];
  int ml_1,ml_2,ml_3,ml_4;       // motif lengths
  int ms1,ms2;                      // computation time
  int n1234_bck;
  int n1234_cl;
  int n1234_stop_bck;
  int n1234_stop_cl;
  int n1234_tr_bck;
  int n1234_tr_cl;
  int n123_bck;     
  int n123_cl;             // total number of co-occurrences (note that A-A-B counts as 2 co-occurrences here
  int n123_stop_bck;
  int n123_stop_cl;
  int n123_tr_bck;        
  int n123_tr_cl;          // number of transcripts in which the 3 motifs co-occur within the appropriate distance constraints (cluster)
  int n12_bck;
  int n12_cl;
  int n12_stop_bck;        
  int n12_stop_cl;         // stop point in call to dist_2v
  int n12_tr_bck;
  int n12_tr_cl; 
  int n1_bck,n2_bck,n3_bck,n4_bck;
  int n1_cl,n2_cl,n3_cl,n4_cl;
  int n[4];                         // to retrieve results from dist2v [0]=n_cooc, [1]=n_tr_cooc, [2]=errorcode, [3]=stop point
  int order_constraint=1;            // if order_constraint=1, then order matters
  int progress_report_int=1;
  int report_pos=0;
  int curr_report_pos;               // current report_pos, to allow report_pos to change throughout the code
  int sortpos=1;
  int temp_n_motifs;                 // used to compute the threshold statistics
  long difftime;                    // time difference
  long tdif;                         // to compute time differences
  struct timeb partial_ti1,partial_tf1;
  struct timeb partial_ti2,partial_tf2;
  struct timeb partial_ti3,partial_tf3;
  struct timeb partial_ti4,partial_tf4;
  struct timeb t1,t2;               // computation time
  struct timeb ti,tf;               // to compute partial times
  struct timeb ti1m,ti2m,ti3m,ti4m,tf1m,tf2m,tf3m,tf4m;
  time_t    time_now;
  int curr_motif1,curr_motif2,curr_motif3,curr_motif4;
  int data_format=1;                   // if pos_format=-1 then data_format=2; this is used in load_and_sort_scan_data_v3
 
  int *m1p_single_transcript;          // used in call to dist_2v_allt_v4 (max_transcript_occurrences x 1)
  int *m2p_single_transcript;          // used in call to dist_2v_allt_v4 (max_transcript_occurrences x 1)
  char *m1s_single_transcript;         // used in call to dist_2v_allt_v4 (max_transcript_occurrences x 1)
  char *m2s_single_transcript;         // used in call to dist_2v_allt_v4 (max_transcript_occurrences x 1)
  int *m3p_single_transcript;          // used in call to dist_3v_allt_v4 (max_transcript_occurrences x 1)
  char *m3s_single_transcript;         // used in call to dist_3v_allt_v4 (max_transcript_occurrences x 1)
  int *m12_p_single_transcript;        // used in call to dist_3v_allt_v4 (max_transcript_occurrences x 1)
  char *m12_s_single_transcript;       // used in call to dist_3v_allt_v4 (max_transcript_occurrences x 1)
  char *m4s_single_transcript;         // used within dist_4v_allt_v4 (max_transcript_occurrences x 1)
  int *m4p_single_transcript;          // used within dist_4v_allt_v4 (max_transcript_occurrences x 1)
  int *m123_p_single_transcript;       // "                           (max_transcript_occurrences x 1)
  char *m123_s_single_transcript;      // "                           (max_transcript_occurrences x 1)
  int **ovdi_single_transcript6;       // used in call to dist_2v_allt_v4 (max_n_ovdi_single_transcripts x 6)
  int **temp_ovdi6;                    // used in call to dist_2/3/4v_allt_v4 (max_n_ovdi_single_transcripts x 6)
  int **ovdi_single_transcript9;       // used in call to dist_3v_allt_v4 (max_n_ovdi_single_transcript x 9)
  int **temp_ovdi9;                    // used in call to dist_3v_allt_v4 (max_n_ovdi_single_transcript x 9)
  int **ovdi_single_transcript12;      // used within dist_4v_allt_v4 (max_n_ovdi_single_transcript x 12)
  int **temp_ovdi12;                   // "                           (max_n_ovdi_single_transcript x 12)
  int **all_ovdi6;                     // additional matrix required by dist_3v_ordered_v4

  time(&time_now);     /* get time in seconds */
  sprintf(last_modified,"03_02_2004");
  
  if (argc<4) 
    usage();
  curr_arg=0;
  curr_arg++;sprintf(filename_pos_cl,argv[curr_arg]);
  curr_arg++;sprintf(filename_pos_bck,argv[curr_arg]);
  curr_arg++;sprintf(filename_strands_cl,argv[curr_arg]);
  curr_arg++;sprintf(filename_strands_bck,argv[curr_arg]);
  curr_arg++;n_motifs=atoi(argv[curr_arg]);
  curr_arg++;n_transcripts_cl=atoi(argv[curr_arg]);
  curr_arg++;n_transcripts_bck=atoi(argv[curr_arg]);
  if (n_transcripts_bck>0) {
    r_transcripts_cl_bck=(float)n_transcripts_bck/(float)n_transcripts_cl;
  } else {
    printf("error! n_transcripts_bck must be > 0 (currently n_transcripts_bck=%d\n",n_transcripts_bck);
    exit(1);
  }
  curr_arg++;
  if (argc>curr_arg) 
    max_dist=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    minian_tr=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    minian_p_tr=atof(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    min_dist_factor=atof(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_ul=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_el=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_il=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    init_i=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    init_j=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    init_k=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    init_l=atoi(argv[curr_arg]);
  finit_i=n_motifs;finit_j=n_motifs;finit_k=n_motifs;finit_l=n_motifs;
  curr_arg++;
  if (argc>curr_arg) 
    finit_i=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    finit_j=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    finit_k=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    finit_l=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    sprintf(file_label,argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    report_pos=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    order_constraint=atoi(argv[curr_arg]);        // default = 1
  curr_arg++;
  if (argc>curr_arg)
    sortpos=atoi(argv[curr_arg]);                 // default = 1
  curr_arg++;
  if (argc>curr_arg)
    max_n_motifs=atoi(argv[curr_arg]);            // default = 2
  curr_arg++;sprintf(transcript_lengths_filename_cl,argv[curr_arg]);
  curr_arg++;sprintf(transcript_lengths_filename_bck,argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    filteron=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    pos_format=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    verbose=atoi(argv[curr_arg]);                 // default = 0
  curr_arg++;
  if (argc>curr_arg) 
    debug_mode=atoi(argv[curr_arg]);

  if ( (finit_i>n_motifs) | (finit_j>n_motifs) | (finit_k>n_motifs) | (finit_l>n_motifs) ) {
    printf("error! finit_i=%d, finit_j=%d, finit_k=%d finit_l=%d should be less than n_motifs (%d)\n",finit_i,finit_j,finit_k,finit_l,n_motifs);
    exit(1);
  }
  if (finit_i<=0)
    finit_i=n_motifs;
  if (finit_j<=0)
    finit_j=n_motifs;
  if (finit_k<=0)
    finit_k=n_motifs;
  if (finit_l<=0)
    finit_l=n_motifs;
  if (max_n_motifs == 2) 
    finit_k=init_k-1;
  if (max_n_motifs == 3) 
    finit_l=init_l-1;

  printf("retrieved all parameters... (%d) \n",argc);

  ftime(&t1);
  sprintf(log_filename,"c.%d.%d.%d.%d.cooc.%s.log.txt",max_dist,minian_tr,max_n_motifs,order_constraint,file_label);
  log_file=fopen(log_filename,"w");
  if (!log_file) {
    printf("i could not open the file %s for writing\n",log_filename);
    exit(1);
  }

  sprintf(output_filename,"c.%d.%d.%d.%d.cooc.%s.out.txt",max_dist,minian_tr,max_n_motifs,order_constraint,file_label);
  output_file=fopen(output_filename,"w");
  if (!output_file) {
    printf("i could not open the file %s for writing\n",output_filename);
    exit(1);
  }

  sprintf(positions_filename,"c.%d.%d.%d.%d.cooc.%s.pos.txt",max_dist,minian_tr,max_n_motifs,order_constraint,file_label);
  positions_file=fopen(positions_filename,"w");
  if (!positions_file) {
    printf("i could not open the file %s for writing\n",positions_filename);
    exit(1);
  }

  /***************************************/
  /* estimate probability cut-off values */
  /***************************************/
  temp_n_motifs=finit_i-init_i+1;
  if (order_constraint<=0) {
    c1=temp_n_motifs;                        
    c2=dbico(temp_n_motifs,2);
    c3=dbico(temp_n_motifs,3);
    c4=dbico(temp_n_motifs,4);
  } else {
    c1=temp_n_motifs;
    c2=temp_n_motifs*temp_n_motifs;
    c3=c2*temp_n_motifs;
    c4=c3*temp_n_motifs;
  }
  //minian_p_tr_1m=fudge_threshold_1m*(double)minian_p_tr/c1;
  minian_p_tr_2m=fudge_threshold_2m*(double)minian_p_tr/c2;
  minian_p_tr_3m=fudge_threshold_3m*(double)minian_p_tr/c3;
  minian_p_tr_4m=fudge_threshold_4m*(double)minian_p_tr/c4;

  //cminian_p_tr_1m=1.0-minian_p_tr_1m;
  cminian_p_tr_2m=1.0-minian_p_tr_2m;
  cminian_p_tr_3m=1.0-minian_p_tr_3m;
  cminian_p_tr_4m=1.0-minian_p_tr_4m;

  sprintf(infotext,"%% cooc_4m_v%d.c\n",prog_version);
  sprintf(temptext,"%% last modified = %s\n",last_modified);strcat(infotext,temptext);
  sprintf(temptext,"%% %s", asctime(localtime(&time_now)));strcat(infotext,temptext);
  sprintf(temptext,"%% ${CODEDIR}/cpp/executables/cooc_4m_v%d.exe ",prog_version);
  for (i=1;i<argc;i++) {
    sprintf(temptext2,"%s ",argv[i]);
    strcat(temptext,temptext2);
  }
  sprintf(temptext2,"\n");strcat(temptext,temptext2);
  strcat(infotext,temptext);
  sprintf(temptext,"%% filename_pos_cl=\t%s\n",filename_pos_cl);strcat(infotext,temptext);
  sprintf(temptext,"%% filename_pos_bck=\t%s\n",filename_pos_bck);strcat(infotext,temptext);
  sprintf(temptext,"%% filename_strands_cl=\t%s\n",filename_strands_cl);strcat(infotext,temptext);
  sprintf(temptext,"%% filename_strands_bck=\t%s\n",filename_strands_bck);strcat(infotext,temptext);
  sprintf(temptext,"%% n_motifs=\t%d\n",n_motifs);strcat(infotext,temptext);
  sprintf(temptext,"%% n_transcripts_cl=\t%d\n",n_transcripts_cl);strcat(infotext,temptext);
  sprintf(temptext,"%% n_transcirpts_bck=\t%d\n",n_transcripts_bck);strcat(infotext,temptext);
  sprintf(temptext,"%% r_transcripts_cl_bck=\t%.2f\n",r_transcripts_cl_bck);strcat(infotext,temptext);
  sprintf(temptext,"%% max_dist=\t%d\n",max_dist);strcat(infotext,temptext);
  sprintf(temptext,"%% minian_tr=\t%d\n",minian_tr);strcat(infotext,temptext);
  sprintf(temptext,"%% minian_p_tr=\t%1.2g\n",minian_p_tr);strcat(infotext,temptext);
  sprintf(temptext,"%% min_dist_factor=\t%.2f\n",min_dist_factor);strcat(infotext,temptext);
  sprintf(temptext,"%% max_ul=\t%d\tmax_el=\t%d\tmax_il=\t%d\n",max_ul,max_el,max_il);strcat(infotext,temptext);
  sprintf(temptext,"%% initial conditions=\ti=%d j=%d k=%d l=%d\n",init_i,init_j,init_k,init_l);strcat(infotext,temptext);
  sprintf(temptext,"%% final conditions=\ti=%d j=%d k=%d l=%d\n",finit_i,finit_j,finit_k,finit_l);strcat(infotext,temptext);
  sprintf(temptext,"%% c1=%1.2g\tc2=%1.2g\tc3=%1.2g\tc4=%1.2g\n",c1,c2,c3,c4);strcat(infotext,temptext);
  sprintf(temptext,"%% minian_p_tr_2m=%1.2g\tminian_p_tr_3m=%1.2g\tminian_p_tr_4m=%1.2g\n",minian_p_tr_2m,minian_p_tr_3m,minian_p_tr_4m);strcat(infotext,temptext);
  sprintf(temptext,"%% order_constraint=\t%d\n",order_constraint);strcat(infotext,temptext);
  sprintf(temptext,"%% report_pos=\t%d\n",report_pos);strcat(infotext,temptext);
  sprintf(temptext,"%% sortpos=\t%d\n",sortpos);strcat(infotext,temptext);
  sprintf(temptext,"%% file_label=%s\n",file_label);strcat(infotext,temptext);
  sprintf(temptext,"%% filteron=%d\n",filteron);strcat(infotext,temptext);
  sprintf(temptext,"%% pos_format=%d\n",pos_format);strcat(infotext,temptext);
  sprintf(temptext,"%% verbose=%d\n",verbose);strcat(infotext,temptext);
  sprintf(temptext,"%% max_n_motifs=%d\n",max_n_motifs);strcat(infotext,temptext);
  sprintf(temptext,"%% transcript_lengths_filename_cl=%s\n",transcript_lengths_filename_cl);strcat(infotext,temptext);
  sprintf(temptext,"%% transcript_lengths_filename_bck=%s\n",transcript_lengths_filename_bck);strcat(infotext,temptext);
  sprintf(temptext,"%%");strcat(infotext,temptext);

  printinfo(infotext,verbose,log_file);
  fprintf(output_file,"%s\n",infotext);
  if (report_pos==1)
    fprintf(positions_file,"%s\n",infotext);
  fprintf(output_file,"%% i\tj\tk\tl\tml_1\tml_2\tml_3\tml_4\tn_cl\tn_tr_cl\tn_bck\tn_tr_bck\tr_oc\tr_tr\tp_binom_oc\tp_binom_tr\n");

  /*********************/
  /* memory allocation */
  /*********************/
  printf("memory allocation\n");
  memory_req=n_motifs*max_occur*4/1000;
  cumulative_memory_req+=(memory_req*6);
  printf("data_t_cl,data_p_cl,data_t_bck,data_p_bck,data_s_cl,data_s_bck [int, %d x %d, %.1f kb]\n",n_motifs,max_occur,memory_req);
  data_t_cl=imatrix(1,n_motifs,0,max_occur);
  data_p_cl=imatrix(1,n_motifs,0,max_occur);
  data_t_bck=imatrix(1,n_motifs,0,max_occur);
  data_p_bck=imatrix(1,n_motifs,0,max_occur);
  data_s_cl=imatrix(1,n_motifs,0,max_occur);
  data_s_bck=imatrix(1,n_motifs,0,max_occur);
  for (i=1;i<=n_motifs;i++) {
    data_t_cl[i][0]=0;
    data_p_cl[i][0]=0;
    data_t_bck[i][0]=0;
    data_p_bck[i][0]=0;
    data_s_cl[i][0]=0;
    data_s_bck[i][0]=0;
    for (j=1;j<=max_occur;j++) {
      data_t_cl[i][j]=-1;
      data_p_cl[i][j]=-1;
      data_s_cl[i][j]=-1;
      data_t_bck[i][j]=-1;
      data_p_bck[i][j]=-1;
      data_s_bck[i][j]=-1;
    }
  }

  memory_req=max_n_ovdi12*5*4/1000;
  cumulative_memory_req+=(memory_req*2);
  printf("ovdi12_cl,ovdi12_bck [int, %d x %d, %.1f kb]\n",max_n_ovdi12,5,memory_req);
  ovdi12_cl=imatrix(1,max_n_ovdi12,1,7);
  ovdi12_bck=imatrix(1,max_n_ovdi12,1,7);
  if (max_n_motifs>2) {
    memory_req=max_n_ovdi123*7*4/1000;
    cumulative_memory_req+=(memory_req*2);
    printf("ovdi123_cl,ovdi123_bck [int, %d x %d, %.1f kb]\n",max_n_ovdi123,7,memory_req);
    ovdi123_cl=imatrix(1,max_n_ovdi123,1,10);
    ovdi123_bck=imatrix(1,max_n_ovdi123,1,10);
  }
  if (max_n_motifs>3) {
    memory_req=max_n_ovdi1234*9*4/1000;
    cumulative_memory_req+=(memory_req*2);
    printf("ovdi1234_cl,ovdi1234_bck [int, %d x %d, %.1f kb]\n",max_n_ovdi1234,9,memory_req);
    ovdi1234_cl=imatrix(1,max_n_ovdi1234,1,13);
    ovdi1234_bck=imatrix(1,max_n_ovdi1234,1,13);
  }
  motif_lengths=ivector(1,n_motifs);

  memory_req=max_occur*4/1000;
  cumulative_memory_req+=memory_req*12;
  printf("cluster:m1/2_t,m1/2_p,m1/2_s [int, %d x %d, %.1f kb]\n",max_occur,1,memory_req);
  m1_t_cl=ivector(1,max_occur);
  m1_p_cl=ivector(1,max_occur);
  m1_s_cl=cvector(1,max_occur);
  m2_t_cl=ivector(1,max_occur);
  m2_p_cl=ivector(1,max_occur);
  m2_s_cl=cvector(1,max_occur);
  m1_t_bck=ivector(1,max_occur);
  m1_p_bck=ivector(1,max_occur);
  m1_s_bck=cvector(1,max_occur);
  m2_t_bck=ivector(1,max_occur);
  m2_p_bck=ivector(1,max_occur);
  m2_s_bck=cvector(1,max_occur);
  if (max_n_motifs>2) {
    memory_req=max_occur*4/1000;
    cumulative_memory_req+=memory_req*6;
    printf("cluster:m3_t,m3_p,m3_s [int, %d x %d, %.1f kb]\n",max_occur,1,memory_req);
    m3_t_cl=ivector(1,max_occur);
    m3_p_cl=ivector(1,max_occur);
    m3_s_cl=cvector(1,max_occur);
    m3_t_bck=ivector(1,max_occur);
    m3_p_bck=ivector(1,max_occur);
    m3_s_bck=cvector(1,max_occur);
  }
  if (max_n_motifs>3) {
    memory_req=max_occur*4/1000;
    cumulative_memory_req+=memory_req*6;
    printf("cluster:m4_t,m4_p,m4_s [int, %d x %d, %.1f kb]\n",max_occur,1,memory_req);
    m4_t_cl=ivector(1,max_occur);
    m4_p_cl=ivector(1,max_occur);
    m4_s_cl=cvector(1,max_occur);
    m4_t_bck=ivector(1,max_occur);
    m4_p_bck=ivector(1,max_occur);
    m4_s_bck=cvector(1,max_occur);
  }

  memory_req=n_transcripts_cl*3*4/1000;
  printf("transcript_lengths_cl [int, %d x 3, %.1f kb]\n",n_transcripts_cl,memory_req);
  transcript_lengths_cl=imatrix(1,n_transcripts_cl,1,3);
  memory_req=n_transcripts_bck*3*4/1000;
  printf("transcript_lengths_bck [int, %d x 3, %.1f kb]\n",n_transcripts_bck,memory_req);
  transcript_lengths_bck=imatrix(1,n_transcripts_bck,1,3);

  memory_req=n_transcripts_cl*4/1000;
  printf("max_n_transcripts_2/3/4m [int, %d x %d, %.1f kb]\n",n_transcripts_cl,1,memory_req);
  max_n_transcripts_2m=ivector(1,n_transcripts_cl);
  max_n_transcripts_3m=ivector(1,n_transcripts_cl);
  max_n_transcripts_4m=ivector(1,n_transcripts_cl);
  cumulative_memory_req+=memory_req;
  for (i=1;i<=n_transcripts_cl;i++) { 
    max_n_transcripts_2m[i]=0;
    max_n_transcripts_3m[i]=0;
    max_n_transcripts_4m[i]=0;
  }

  /* memory allocation for internal variables in the dist_2/3/4v functions */
  printf("allocating memory for interval variables in the dist_2/3/4v functions...\n");
  memory_req=max_transcript_occurrences*4/1000;cumulative_memory_req+=(2*memory_req);                       // 0.8 Mb with default parameters
  printf("m1/2p_single_transcript [int, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m1p_single_transcript=ivector(1,max_transcript_occurrences);m2p_single_transcript=ivector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*2/1000;cumulative_memory_req+=(2*memory_req);                       // 0.4 Mb with default parameters
  printf("m1/2s_single_transcript [char, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m1s_single_transcript=cvector(1,max_transcript_occurrences);m2s_single_transcript=cvector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*4/1000;cumulative_memory_req+=memory_req;                           // 0.4 Mb with default parameters
  printf("m3p_single_transcript [int, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m3p_single_transcript=ivector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*2/1000;cumulative_memory_req+=memory_req;                           // 0.2 Mb with default parameters
  printf("m3s_single_transcript [char, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m3s_single_transcript=cvector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*4/1000;cumulative_memory_req+=memory_req;                           // 0.4 Mb with default parameters
  printf("m4p_single_transcript [int, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m4p_single_transcript=ivector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*2/1000;cumulative_memory_req+=memory_req;                           // 0.2 Mb with default parameters
  printf("m4s_single_transcript [char, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m4s_single_transcript=cvector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*4/1000;cumulative_memory_req+=memory_req;                           // 0.4 Mb with default parameters
  printf("m12_p_single_transcript [int, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m12_p_single_transcript=ivector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*2/1000;cumulative_memory_req+=memory_req;                           // 0.2 Mb with default parameters
  printf("m12_s_single_transcript [char, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m12_s_single_transcript=cvector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*4/1000;cumulative_memory_req+=memory_req;                           // 0.4 Mb with default parameters
  printf("m123_p_single_transcript [int, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m123_p_single_transcript=ivector(1,max_transcript_occurrences);
  memory_req=max_transcript_occurrences*2/1000;cumulative_memory_req+=memory_req;                           // 0.2 Mb with default parameters
  printf("m123_s_single_transcript [char, %d x %d, %.1f kb]\n",max_transcript_occurrences,1,memory_req);
  m123_s_single_transcript=cvector(1,max_transcript_occurrences);

  memory_req=max_n_ovdi_single_transcript*6*4/1000;cumulative_memory_req+=memory_req;                       // 24 Mb
  printf("ovdi_single_transcript6 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,6,memory_req);
  ovdi_single_transcript6=imatrix(1,max_n_ovdi_single_transcript,1,6);
  memory_req=max_n_ovdi_single_transcript*6*4/1000;cumulative_memory_req+=memory_req;                       // 24 Mb
  printf("temp_ovdi6 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,6,memory_req);
  temp_ovdi6=imatrix(1,max_n_ovdi_single_transcript,1,6);
  memory_req=max_n_ovdi_single_transcript*6*4/1000;cumulative_memory_req+=memory_req;                       // 24 Mb
  printf("all_ovdi6 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,6,memory_req);
  all_ovdi6=imatrix(1,max_n_ovdi_single_transcript,1,6);
  memory_req=max_n_ovdi_single_transcript*9*4/1000;cumulative_memory_req+=memory_req;                       // 36 Mb
  printf("ovdi_single_transcript9 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,9,memory_req);
  ovdi_single_transcript9=imatrix(1,max_n_ovdi_single_transcript,1,9);
  memory_req=max_n_ovdi_single_transcript*9*4/1000;cumulative_memory_req+=memory_req;                       // 36 Mb
  printf("temp_ovdi9 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,9,memory_req);
  temp_ovdi9=imatrix(1,max_n_ovdi_single_transcript,1,9);
  memory_req=max_n_ovdi_single_transcript*12*4/1000;cumulative_memory_req+=memory_req;                      // 48 Mb
  printf("ovdi_single_transcript12 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,12,memory_req);
  ovdi_single_transcript12=imatrix(1,max_n_ovdi_single_transcript,1,12);
  memory_req=max_n_ovdi_single_transcript*12*4/1000;cumulative_memory_req+=memory_req;                      // 48 Mb
  printf("temp_ovdi12 [int, %d x %d, %.1f kb]\n",max_n_ovdi_single_transcript,12,memory_req);
  temp_ovdi12=imatrix(1,max_n_ovdi_single_transcript,1,12);

  sprintf(infotext,"cumulative_memory_req=%.0f kb\n",cumulative_memory_req);
  printinfo(infotext,verbose,log_file);

  /***************************/
  /* load transcript lengths */
  /***************************/
  if (verbose) 
    fprintf(log_file,"loading transcript lengths from %s\n",transcript_lengths_filename_cl);
  input_file=fopen(transcript_lengths_filename_cl,"r");
  if (!input_file) {
    printf("error! i could not open the transcript_lengths_filename_cl = %s; exiting...\n",transcript_lengths_filename_cl);
    exit(1);
  } 
  i=0;
  while (feof(input_file)==0) {
    i++;
    fscanf(input_file,"%d\n",&k);                // read sequence number
    for (j=1;j<=3;j++) {
      fscanf(input_file,"%d\n",&k);
      transcript_lengths_cl[i][j]=(int)k;
    }
  }
  if (i != n_transcripts_cl) 
    printf("\n\n ------   WARNING!!! n_transcripts_cl=%d and i=%d in cooc_4m_v%d.c, load transcript initial lengths ----- \n\n",n_transcripts_cl,i,prog_version);
  fclose(input_file);
  if (verbose) {
    fprintf(log_file,"%% transcript_lengths_cl\n");
    fprintf(log_file,"%% ul\tel\til\n");
    printf("%% transcript_lengths_cl\n");
    printf("%% ul\tel\til\n");
    for (i=1;i<=n_transcripts_cl;i++) {
      fprintf(log_file,"%d\t%d\t%d\t%d\n",i,transcript_lengths_cl[i][1],transcript_lengths_cl[i][2],transcript_lengths_cl[i][3]);
      printf("%d\t%d\t%d\t%d\n",i,transcript_lengths_cl[i][1],transcript_lengths_cl[i][2],transcript_lengths_cl[i][3]);
    }
  }

  input_file=fopen(transcript_lengths_filename_bck,"r");
  if (!input_file) {
    printf("error! i could not open the transcript_lengths_filename_bck = %s; exiting...\n",transcript_lengths_filename_bck);
    exit(1);
  } 
  i=0;
  while (feof(input_file)==0) {
    i++;
    fscanf(input_file,"%d\n",&k);                // read sequence number
    for (j=1;j<=3;j++) {
      fscanf(input_file,"%d\n",&k);
      transcript_lengths_bck[i][j]=(int)k;
    }
  }
  if (i != n_transcripts_bck) 
    printf("\n\n ------   WARNING!!! n_transcripts_bck=%d and i=%d in cooc_4m_v%d.c, load transcript initial lengths ----- \n\n",n_transcripts_bck,i,prog_version);
  fclose(input_file);

  /******************/
  /* load scan data */
  /******************/
  if (pos_format==-1)
    data_format=2;
  sprintf(infotext,"%% loading cluster data from filename_pos_cl=%s",filename_pos_cl);
  printinfo(infotext,verbose,log_file);
  //errorcode=load_and_sort_scan_data(filename_cl,data_t_cl,data_p_cl,motif_lengths,neg2pos,sortpos,max_occur,n_motifs,max_ul,max_el,max_il,transcript_lengths_cl);
  errorcode=load_and_sort_scan_data_v3(filename_pos_cl,filename_strands_cl,data_t_cl,data_p_cl,data_s_cl,motif_lengths,sortpos,max_occur,n_motifs,max_ul,max_el,max_il,transcript_lengths_cl,data_format,filteron,0);           // non-verbose
  //errorcode=load_and_sort_scan_data_v2(filename_pos_cl,filename_strands_cl,data_t_cl,data_p_cl,data_s_cl,motif_lengths,neg2pos,sortpos,max_occur,n_motifs,max_ul,max_el,max_il,transcript_lengths_cl);
  //errorcode=load_and_sort_scan_data_v2(filename_pos_cl,filename_strands_cl,data_t_cl,data_p_cl,data_s_cl,motif_lengths,neg2pos,sortpos,max_occur,n_motifs,-max_ul,max_el,max_il,transcript_lengths_cl);
  sprintf(infotext,"%% \tready.");
  printinfo(infotext,verbose,log_file);
  if (errorcode != 0) {
    printf("\nERROR!\n");
    printf("infotext=%s\n",infotext);
    printf("call to load_and_sort_scan_data_v2 (cluster) returned error code = %d, exiting... \n",errorcode);
    exit(1);
  }

  sprintf(infotext,"%% loading background data from filename_pos_bck=%s",filename_pos_bck);
  printinfo(infotext,verbose,log_file);
  errorcode=load_and_sort_scan_data_v3(filename_pos_bck,filename_strands_bck,data_t_bck,data_p_bck,data_s_bck,motif_lengths,sortpos,max_occur,n_motifs,max_ul,max_el,max_il,transcript_lengths_bck,data_format,filteron,0);
  /* if (filteron==1) {
     errorcode=load_and_sort_scan_data_v2(filename_pos_bck,filename_strands_bck,data_t_bck,data_p_bck,data_s_bck,motif_lengths,neg2pos,sortpos,max_occur,n_motifs,max_ul,max_el,max_il,transcript_lengths_bck);
     } else {
     errorcode=load_and_sort_scan_data_v2(filename_pos_bck,filename_strands_bck,data_t_bck,data_p_bck,data_s_bck,motif_lengths,neg2pos,sortpos,max_occur,n_motifs,-max_ul,max_el,max_il,transcript_lengths_bck);
     }
  */
  sprintf(infotext,"%% \tready.");
  printinfo(infotext,verbose,log_file);
  if (errorcode != 0) {
    printf("\nERROR!\n");
    printf("infotext=%s\n",infotext);
    printf("call to load_and_sort_scan_data_v2 (background) returned error code = %d, exiting... \n",errorcode);
    exit(1);
  }
  
  /*******************/
  /* begin the begin */
  /*******************/
  if (verbose>=1) 
    fprintf(log_file,"%% i\tj\tk\tl\tcurr_tr\tp1\tp2\tp3\tp4\tr_oc\tr_tr\tp_binom_oc\tp_binom_tr\n");

  printf("starting the fun...\n");
  for (curr_motif1=init_i;curr_motif1<=finit_i;curr_motif1++) {
    ml_1=motif_lengths[curr_motif1];
    ftime(&ti1m);
    if ( (curr_motif1%progress_report_int) == 0) {
      if (verbose>=0) {
	sprintf(infotext,"%%processing i=%d...",curr_motif1);
	printinfo(infotext,verbose,log_file);
      }
    }
    /* get line curr_motif1, i.e. occurrences of motif number curr_motif1 */
    n1_cl=data_t_cl[curr_motif1][0];                             // total number of occurrences of motif 1
    if (verbose>=2) {
      sprintf(infotext,"processing curr_motif1=%d\tn1_cl=%d",curr_motif1,n1_cl);
      printinfo(infotext,verbose,log_file);
      if (debug_mode) printf("%s\n",infotext);
    }    
    for (ll=1;ll<=n1_cl;ll++) {
      m1_t_cl[ll]=data_t_cl[curr_motif1][ll];                    // transcripts
      m1_p_cl[ll]=data_p_cl[curr_motif1][ll];                    // positions
      m1_s_cl[ll]=data_s_cl[curr_motif1][ll];                    // strands
    }
    if (n1_cl>=minian_tr) {
      /***********/
      /* motif 2 */
      /**********/
      if (order_constraint <= 0) 
	init_j=curr_motif1;                                      // if order does not matter then ij is equivalent to ji, therefore, start at j=i
      for (curr_motif2=init_j;curr_motif2<=finit_j;curr_motif2++) {
	ml_2=motif_lengths[curr_motif2];                         // length of motif 2
	ftime(&ti2m);
	/***************/
	/* cluster set */
	/***************/
	n2_cl=data_t_cl[curr_motif2][0];                         // total number of occurrences of motif 2
	if (verbose>=2) {
	  sprintf(infotext,"\tprocessing curr_motif1=%d\tcurr_motif2=%d\tn2_cl=%d",curr_motif1,curr_motif2,n2_cl);
	  printinfo(infotext,verbose,log_file);
	  if (debug_mode) printf("%s\n",infotext);
	}	
	if (n2_cl>=minian_tr) {
	  for (ll=1;ll<=n2_cl;ll++) {
	    m2_t_cl[ll]=data_t_cl[curr_motif2][ll];              // transcripts
	    m2_p_cl[ll]=data_p_cl[curr_motif2][ll];              // positions
	    m2_s_cl[ll]=data_s_cl[curr_motif2][ll];              // strands
	  }
	  if (verbose>=4) {
	    fprintf(log_file,"\t\tmm m1_t m1_p m1_s m2_t m2_p m2_s\n");
	    for (mm=1;mm<=5;mm++) 
	      fprintf(log_file,"\t\t%d %d %d %d %d %d %d\n",mm,m1_t_cl[mm],m1_p_cl[mm],m1_s_cl[mm],m2_t_cl[mm],m2_p_cl[mm],m2_s_cl[mm]);
	  }
	  //min_dist=(int)ceil(ml_1*min_dist_factor);
	  errorcode=0;
	  ml[1]=ml_1;
	  ml[2]=ml_2;
	  min_dist=compute_min_dist(ml,2,min_dist_factor);
	  max_n_transcripts=-1;    // that is, compute for all transcripts
	  curr_report_pos=1;    // report positions by default (because they are needed at later stages)
	  if (max_n_motifs==2) 
	    curr_report_pos=report_pos;
	  //dist_2v_allt_v3(m1_t_cl,m2_t_cl,m1_p_cl,m2_p_cl,m1_s_cl,m2_s_cl,n1_cl,n2_cl,ovdi12_cl,n,min_dist,max_dist,n_transcripts_cl,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi12,pos_format,debug_mode);
	  dist_2v_allt_v4(m1_t_cl,m2_t_cl,m1_p_cl,m2_p_cl,m1_s_cl,m2_s_cl,n1_cl,n2_cl,ovdi12_cl,n,min_dist,max_dist,n_transcripts_cl,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi12,pos_format,ovdi_single_transcript6,temp_ovdi6,m1p_single_transcript,m1s_single_transcript,m2p_single_transcript,m2s_single_transcript,debug_mode);
	  //void dist_2v_allt_v4(int *m1_t,int *m2_t,int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int **ovdi_single_transcript,int **temp_ovdi,int *m1p,char *m1s,int *m2p,char *m2s,int verbose);
	  if (n[2]<0) {
	    n12_cl=n[3];        // get the stop value if there was an error
	    errorcode=errorcode+1;
	    if (verbose>0)
	      fprintf(log_file,"%% WARNING!\tcurr_motif1=%d\tcurr_motif2=%d\tcl, n12_cl=%d,n12_stop_cl=%d, memory overload, results may be inaccurate\n",curr_motif1,curr_motif2,n[0],n[3]);
	  } else {
	    n12_cl=n[0];
	  }
	  n12_tr_cl=n[1];
	  if (verbose>=2) 
	    fprintf(log_file,"\tn12_cl=%d\tn12_tr_cl=%d\tn12_stop_cl=%d\n",n[0],n[1],n[3]);
	  if (n12_tr_cl>=minian_tr) {
	    if (verbose>=3) {
	      if ( (max_n_motifs>2) | (report_pos==1) ) {
		fprintf(log_file,"\tovdi12_cl\n");
		fprintf(log_file,"\ttr\tp1\tp2\ti1\ti2\ts1\ts2\n");
		mm=0;
		while ( (mm<n12_cl) & (mm<max_n_ovdi12) ) {
		  mm++;
		  fprintf(log_file,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ovdi12_cl[mm][1],ovdi12_cl[mm][2],ovdi12_cl[mm][3],ovdi12_cl[mm][4],ovdi12_cl[mm][5],ovdi12_cl[mm][6],ovdi12_cl[mm][7]);
		}
	      }
	    }
	    /************************/
	    /* background, 2 motifs */
	    /************************/
	    n12_bck=0;
	    n12_tr_bck=0;
	    // compute statistics for background
	    n1_bck=data_t_bck[curr_motif1][0]; // number of transcripts in which motif 1 is present in the background set  
	    /* assign transcripts and positions for those occurrences of motif 1 in background */
	    for (mm=1;mm<=n1_bck;mm++) {
	      m1_t_bck[mm]=data_t_bck[curr_motif1][mm];
	      m1_p_bck[mm]=data_p_bck[curr_motif1][mm];
	      m1_s_bck[mm]=data_s_bck[curr_motif1][mm];
	    }
	    n2_bck=data_t_bck[curr_motif2][0]; // number of transcripts in which motif 2 is present in the background set
	    if (verbose>=2) {
	      fprintf(log_file,"\t\t\t\tprocessing curr_motif1=%d\tcurr_motif2=%d\tn1_bck=%d\tn2_bck=%d\n",curr_motif1,curr_motif2,n1_bck,n2_bck);  
	      if (debug_mode) printf("%s\n",infotext);
	    }
	    if (n2_bck>0) {	      // assign transcripts and positions for those occurrences of motif 2 in background 
	      for (mm=1;mm<=n2_bck;mm++) {
		m2_t_bck[mm]=data_t_bck[curr_motif2][mm];
		m2_p_bck[mm]=data_p_bck[curr_motif2][mm];
		m2_s_bck[mm]=data_s_bck[curr_motif2][mm];
	      }	      
	      //min_dist=(int)ceil(ml_1*min_dist_factor);
	      if (max_n_motifs==2) {
		if (max_n_transcripts_2m[n12_tr_cl]==0) {
		  max_n_transcripts=get_max_n_transcripts(n_transcripts_cl,n12_tr_cl,cminian_p_tr_2m,n_transcripts_bck);	/* compute max_n_transcripts */
		  max_n_transcripts_2m[n12_tr_cl]=max_n_transcripts;
		} else {
		  max_n_transcripts=max_n_transcripts_2m[n12_tr_cl];
		} 
		curr_report_pos=0;                  // we do not need all the positions in the background if we are going to stop here
	      } else {
		max_n_transcripts=-1;           // because we will need the information in ovdi12_bck when adding more motifs
		curr_report_pos=1;                  // "
	      }
	      //dist_2v_allt_v3(m1_t_bck,m2_t_bck,m1_p_bck,m2_p_bck,m1_s_bck,m2_s_bck,n1_bck,n2_bck,ovdi12_bck,n,min_dist,max_dist,n_transcripts_bck,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi12,pos_format,0);
	      dist_2v_allt_v4(m1_t_bck,m2_t_bck,m1_p_bck,m2_p_bck,m1_s_bck,m2_s_bck,n1_bck,n2_bck,ovdi12_bck,n,min_dist,max_dist,n_transcripts_bck,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi12,pos_format,ovdi_single_transcript6,temp_ovdi6,m1p_single_transcript,m1s_single_transcript,m2p_single_transcript,m2s_single_transcript,0);
	      if (n[2]<0) {
		n12_bck=n[3];
		if (verbose>0)
		  fprintf(log_file,"%% WARNING!\tcurr_motif1=%d\tcurr_motif2=%d\tbck, n12_bck=%d,n12_stop_bck=%d, memory overload, results may be inaccurate\n",curr_motif1,curr_motif2,n[0],n[3]);	    
		errorcode=errorcode+10;
	      } else {
		n12_bck=n[0];
	      }
	      n12_tr_bck=n[1];
	      if (verbose>=2) 
		fprintf(log_file,"\t\t\t\t\tn12_bck=%d\tn12_tr_bck=%d\n",n12_bck,n12_tr_bck);
	    } else {                // if n2_bck==0
	      n12_bck=0;
	      n12_tr_bck=0;
	    }                       // close else to check on n2_bck>0
	    compute_pr(n12_cl,n12_tr_cl,n12_bck,n12_tr_bck,r,p);
	    r_oc=r[0];
	    r_tr=r[1];
	    p_binom_oc=p[0];
	    p_binom_tr=p[1];
	    
	    curr_motif3=-1;curr_motif4=-1;ml_3=-1;ml_4=-1;
	    if (p_binom_tr>minian_p_tr_2m) {
	      if (verbose>=2) 
	       fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n12_cl,n12_tr_cl,n12_bck,n12_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
	    } else {
	      fprintf(output_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n12_cl,n12_tr_cl,n12_bck,n12_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
	      if (verbose>=1) {
		mm=0;
		fprintf(log_file,"m1\tm2\ttr\tp1\tp2\ts1\ts2\tr_oc\tr_tr\tp_oc\tp_tr\tec\n");
		while ( (mm<n12_cl) & (mm<max_n_ovdi12) ) {
		  mm++;
		  curr_tr=ovdi12_cl[mm][1];
		  fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_tr,ovdi12_cl[mm][2],ovdi12_cl[mm][3],ovdi12_cl[mm][6],ovdi12_cl[mm][7],r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
		} // close mm loop to report results
	      } // close check on verbose>=1

	      if (report_pos==1) {
		fprintf(positions_file,"%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);	      // print positions to file 
		mm=0;
		while ( (mm<n12_cl) & (mm<max_n_ovdi12) ) {
		  mm++;
		  curr_tr=ovdi12_cl[mm][1];
		  fprintf(positions_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n",curr_tr,ovdi12_cl[mm][2],ovdi12_cl[mm][3],-1,-1,ovdi12_cl[mm][6],ovdi12_cl[mm][7]);
		}     // close mm loop to report results			  
	      }       // close check on report_pos
	    }         // close check on p_binom_tr > minian_p_tr_2m 
	    /*********************/
	    /* motif 3           */
	    /*********************/
	    if (order_constraint <= 0) 
	      init_k=curr_motif2;
	    for (curr_motif3=init_k;curr_motif3<=finit_k;curr_motif3++) {
	      ml_3=motif_lengths[curr_motif3];
	      /*********************/
	      /* cluster, 3 motifs */
	      /*********************/
	      n3_cl=data_t_cl[curr_motif3][0];
	      if (verbose>=2) {
		sprintf(infotext,"\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tn3_cl=%d",curr_motif1,curr_motif2,curr_motif3,n3_cl);
		printinfo(infotext,verbose,log_file);
		if (debug_mode) printf("%s\n",infotext);
	      }
	      //printf("\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tn3_cl=%d",curr_motif1,curr_motif2,curr_motif3,n3_cl);
	      if (n3_cl>=minian_tr) {
		for (ll=1;ll<=n3_cl;ll++) {
		  m3_t_cl[ll]=data_t_cl[curr_motif3][ll];
		  m3_p_cl[ll]=data_p_cl[curr_motif3][ll];
		  m3_s_cl[ll]=data_s_cl[curr_motif3][ll];
		}
		if (verbose>=3) {
		  fprintf(log_file,"mm\tm1_t\tm1_p\tm1_s\tm2_t\tm2_p\tm2_s\tm3_t\tm3_p\tm3_s\n");
		  for (mm=1;mm<=5;mm++) 
		    fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",mm,m1_t_cl[mm],m1_p_cl[mm],m1_s_cl[mm],m2_t_cl[mm],m2_p_cl[mm],m2_s_cl[mm],m3_t_cl[mm],m3_p_cl[mm],m3_s_cl[mm]);
		  fprintf(log_file,"\n");
		}
		//min_dist=(int)ceil(ml_2*min_dist_factor);
		errorcode=0;
		ml[3]=ml_3;
		min_dist=compute_min_dist(ml,3,min_dist_factor); 
		max_n_transcripts=-1;              // therefore compute for all transcripts
		curr_report_pos=1;
		if (max_n_motifs==3) 
		  curr_report_pos=report_pos;
		//dist_3v_allt_v3(ovdi12_cl,m3_t_cl,m3_p_cl,m3_s_cl,n12_cl,n3_cl,ovdi123_cl,n,min_dist,max_dist,n_transcripts_cl,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi123,pos_format,debug_mode);
		dist_3v_allt_v4(ovdi12_cl,m3_t_cl,m3_p_cl,m3_s_cl,n12_cl,n3_cl,ovdi123_cl,n,min_dist,max_dist,n_transcripts_cl,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi123,pos_format,ovdi_single_transcript9,temp_ovdi9,ovdi_single_transcript6,temp_ovdi6,all_ovdi6,m3p_single_transcript,m3s_single_transcript,m12_p_single_transcript,m12_s_single_transcript,debug_mode);
		//void dist_3v_allt_v4(int **ovdi12,int *m3_t,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int **ovdi_single_transcript,int **temp_ovdi,int **ovdi12_p,int *ovdi123,int **all_ovdi,int *m3p,char *m3s,int *m12_p,char *m12_s,int verbose)
		if (n[2]<0) { 
		  n123_cl=n[3];
		  if (verbose>0)
		    fprintf(log_file,"%% WARNING!\tcurr_motif1=%dcurr_motif2=%d\tcurr_motif3=%d\tcl, n123_cl=%d,n123_stop_cl=%d, memory overload, results may be inaccurate\n",curr_motif1,curr_motif2,curr_motif3,n[0],n[3]);
		  errorcode=errorcode+100;
		} else {
		  n123_cl=n[0];
		}
		n123_tr_cl=n[1];
		if (verbose>=2) 
		  fprintf(log_file,"\tn123_cl=%d\tn123_tr_cl=%d\n",n123_cl,n123_tr_cl);
		if (n123_tr_cl>=minian_tr) {
		  if (verbose>=3) {
		    if ( (max_n_motifs>3) | (report_pos==1) ) {
		      fprintf(log_file,"\tovdi123_cl\n");
		      fprintf(log_file,"\ttr\tp1\tp2\tp3\n");
		      mm=0;
		      while ( (mm<n123_cl) & (mm<max_n_ovdi12) ) {
			mm++;
			fprintf(log_file,"\t\t\t%d\t%d\t%d\t%d\n",ovdi123_cl[mm][1],ovdi123_cl[mm][2],ovdi123_cl[mm][3],ovdi123_cl[mm][4]);
			printf("\t\t\t%d\t%d\t%d\t%d\n",ovdi123_cl[mm][1],ovdi123_cl[mm][2],ovdi123_cl[mm][3],ovdi123_cl[mm][4]);
		      }
		    }
		  }
		  /************************/
		  /* background, 3 motifs */
		  /************************/
		  n123_bck=0;
		  n123_tr_bck=0;
		  n3_bck=data_t_bck[curr_motif3][0];
		  if (verbose>=2) {
		    sprintf(infotext,"\t\t\t\t\t\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tn3_bck=%d",curr_motif1,curr_motif2,curr_motif3,n3_bck);
		    printinfo(infotext,verbose,log_file);
		    if (debug_mode) printf("%s\n",infotext);
		  }
		  if (n3_bck>0) {
		    for (mm=1;mm<=n3_bck;mm++) {
		      m3_t_bck[mm]=data_t_bck[curr_motif3][mm];
		      m3_p_bck[mm]=data_p_bck[curr_motif3][mm];
		      m3_s_bck[mm]=data_s_bck[curr_motif3][mm];
		    }
		    if (max_n_motifs == 3) {
		      max_n_transcripts=max_n_transcripts_3m[n123_tr_cl];
		      if (max_n_transcripts_3m[n123_tr_cl]==0) {
			max_n_transcripts=get_max_n_transcripts(n_transcripts_cl,n123_tr_cl,cminian_p_tr_3m,n_transcripts_bck);	    /* compute max_n_transcripts */
			max_n_transcripts_3m[n123_tr_cl]=max_n_transcripts;
		      } 
		      curr_report_pos=0;
		    } else {
		      max_n_transcripts=-1;           /* we need the position information for later */
		      curr_report_pos=1;
		    }
		    if (n12_bck>0) {
		      //dist_3v_allt_v3(ovdi12_bck,m3_t_bck,m3_p_bck,m3_s_bck,n12_bck,n3_bck,ovdi123_bck,n,min_dist,max_dist,n_transcripts_bck,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi123,pos_format,0);
		      dist_3v_allt_v4(ovdi12_bck,m3_t_bck,m3_p_bck,m3_s_bck,n12_bck,n3_bck,ovdi123_bck,n,min_dist,max_dist,n_transcripts_bck,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi123,pos_format,ovdi_single_transcript9,temp_ovdi9,ovdi_single_transcript6,temp_ovdi6,all_ovdi6,m3p_single_transcript,m3s_single_transcript,m12_p_single_transcript,m12_s_single_transcript,0);
		      if (n[2]<0) {
			n123_bck=n[3];
			if (verbose>0)
			  fprintf(log_file,"%% WARNING!\ti=%dj=%d\tk=%d\tbck, n123_bck=%d, n123_stop_bck=%d, memory overload, results may be inaccurate\n",i,j,k,n[0],n[3]);
			errorcode=errorcode+1000;
		      } else {
			n123_bck=n[0];
		      }
		      n123_tr_bck=n[1];
		      if (verbose>=2) 
			fprintf(log_file,"\t\t\t\t\t\tn123_bck=%d\tn123_tr_bck=%d\tn123_stop_bck=%d\n",n[0],n[1],n[3]);
		    } else {
		      n123_bck=0;
		      n123_tr_bck=0;
		    }            // close check on n12_bck>0
		  } else {
		    n123_bck=0;
		    n123_tr_bck=0;
		  }              // close check on n3_bck > 0
		  compute_pr(n123_cl,n123_tr_cl,n123_bck,n123_tr_bck,r,p);
		  r_oc=r[0];
		  r_tr=r[1];
		  p_binom_oc=p[0];
		  p_binom_tr=p[1];
		  
		  curr_motif4=-1;ml_4=-1;
		  if (p_binom_tr>minian_p_tr_3m) {
		    if (verbose>=2) 
		      fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n123_cl,n123_tr_cl,n123_bck,n123_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
		  } else  {
		    fprintf(output_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n123_cl,n123_tr_cl,n123_bck,n123_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
		    if (verbose>=1) {
		      mm=0;
		      fprintf(log_file,"m1\tm2\tm3\ttr\tp1\tp2\tp3\ts1\ts2\ts3\tr_oc\tr_tr\tp_oc\tp_tr\tec\n");
		      while ( (mm<n123_cl) & (mm<max_n_ovdi12) ) {
			mm++;
			curr_tr=ovdi123_cl[mm][1];
			fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_tr,ovdi123_cl[mm][2],ovdi123_cl[mm][3],ovdi123_cl[mm][4],ovdi123_cl[mm][8],ovdi123_cl[mm][9],ovdi123_cl[mm][10],r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
		      }     // close mm loop to report results
		    }       // close check on verbose>=1
		    if (report_pos==1) {
		      fprintf(positions_file,"%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
		      mm=0;
		      while ( (mm<n123_cl) & (mm<max_n_ovdi12) ) {
			mm++;
			curr_tr=ovdi123_cl[mm][1];
			fprintf(positions_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",curr_tr,ovdi123_cl[mm][2],ovdi123_cl[mm][3],ovdi123_cl[mm][4],-1,ovdi123_cl[mm][8],ovdi123_cl[mm][9],ovdi123_cl[mm][10]);
		      }        // close mm loop to report results			  
		    }          // close check on report_pos==1
		  }            // close check on p_binom_tr > minian_p_tr_3m
		  
		  /**************************/
		  /* motif 4                */
		  /**************************/
		  if (order_constraint <= 0) 
		    init_l=curr_motif3;
		  for (curr_motif4=init_l;curr_motif4<=finit_l;curr_motif4++) {
		    ml_4=motif_lengths[curr_motif4];
		    n4_cl=data_t_cl[curr_motif4][0];
		    if (verbose>=2) {
		      sprintf(infotext,"\t\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tcurr_motif4=%d\tn4_cl=%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,n4_cl);
		      printinfo(infotext,verbose,log_file);
		      if (debug_mode) printf("%s\n",infotext);
		    }
		    //printf("\t\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tcurr_motif4=%d\tn4_cl=%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,n4_cl);
		    if (n4_cl>=minian_tr) {
		      for (mm=1;mm<=n4_cl;mm++) {
			m4_t_cl[mm]=data_t_cl[curr_motif4][mm];
			m4_p_cl[mm]=data_p_cl[curr_motif4][mm];  
			m4_s_cl[mm]=data_s_cl[curr_motif4][mm];
		      }
		      //min_dist=(int)ceil(ml_3*min_dist_factor);
		      ml[4]=ml_4;
		      min_dist=compute_min_dist(ml,4,min_dist_factor);
		      max_n_transcripts=-1;            // therefore all transcripts are included
		      //dist_4v_allt_v3(ovdi123_cl,m4_t_cl,m4_p_cl,m4_s_cl,n123_cl,n4_cl,ovdi1234_cl,n,min_dist,max_dist,n_transcripts_cl,ml,report_pos,order_constraint,max_n_transcripts,max_n_ovdi1234,pos_format,debug_mode);
		      dist_4v_allt_v4(ovdi123_cl,m4_t_cl,m4_p_cl,m4_s_cl,n123_cl,n4_cl,ovdi1234_cl,n,min_dist,max_dist,n_transcripts_cl,ml,report_pos,order_constraint,max_n_transcripts,max_n_ovdi1234,pos_format,ovdi_single_transcript12,temp_ovdi12,ovdi_single_transcript9,temp_ovdi6,ovdi_single_transcript6,m4p_single_transcript,m4s_single_transcript,m123_p_single_transcript,m123_s_single_transcript,debug_mode);
		      //void dist_4v_allt_v4(int **ovdi123,int *m4_t,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int **ovdi_single_transcript,int **temp_ovdi,int **ovdi123_p,int **ovdi1234,int **all_ovdi,int *m4p,char *m4s,int *m123_p,char *m123_s,int verbose)
		      if (n[2]<0) {
			n1234_cl=n[3];
			if (verbose>0)
			  fprintf(log_file,"%% WARNING!\tcurr_motif1=%dcurr_motif2=%d\tcurr_motif3=%d\tcurr_motif4=%d\tcl, n1234_cl=%d, n1234_stop_cl=%d, memory overload, results may be inaccurate\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,n[0],n[3]);
			errorcode=errorcode+10000;
		      } else {
			n1234_cl=n[0];
		      }
		      n1234_tr_cl=n[1];
		      if (verbose>=2) {
			fprintf(log_file,"\t\t\tn1234_cl=%d\tn1234_tr_cl=%d\n",n1234_cl,n1234_tr_cl);
			if (verbose>=3) {
			  if (report_pos==1) {
			    fprintf(log_file,"\t\tovdi1234_cl\n");
			    fprintf(log_file,"\t\t\ttr\tp3\tp4\n");
			    mm=0;
			    while ( (mm<n1234_cl) & (mm<max_n_ovdi12) ) {
			      mm++;
			      fprintf(log_file,"\t\t\t%d\t%d\t%d\n",ovdi1234_cl[mm][1],ovdi1234_cl[mm][2],ovdi1234_cl[mm][3]);
			    }
			  }
			}
		      }
		      if (n1234_tr_cl>=minian_tr) {
			n4_bck=data_t_bck[curr_motif4][0];
			if (verbose>=2) {
			  fprintf(log_file,"\t\t\t\t\t\tprocessing curr_motif1=%d\tcurr_motif2=%d\tcurr_motif3=%d\tcurr_motif4=%d\tn4_bck=%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,n4_bck);
			  if (debug_mode) printf("%s\n",infotext);
			}
			if (n4_bck>0) {
			  for (mm=1;mm<=n4_bck;mm++) {
			    m4_t_bck[mm]=data_t_bck[curr_motif4][mm];
			    m4_p_bck[mm]=data_p_bck[curr_motif4][mm]; 
			    m4_s_bck[mm]=data_s_bck[curr_motif4][mm];
			  }			  
			  max_n_transcripts=max_n_transcripts_4m[n1234_tr_cl];
			  if (max_n_transcripts==0) {
			    max_n_transcripts=get_max_n_transcripts(n_transcripts_cl,n1234_tr_cl,cminian_p_tr_4m,n_transcripts_bck);
			    max_n_transcripts_4m[n1234_tr_cl]=max_n_transcripts;
			  } 
			  curr_report_pos=0;
			  if (n123_bck>0) {
			    //dist_4v_allt_v3(ovdi123_bck,m4_t_bck,m4_p_bck,m4_s_bck,n123_bck,n4_bck,ovdi1234_bck,n,min_dist,max_dist,n_transcripts_bck,ml,curr_report_pos,order_constraint,max_n_transcripts,max_n_ovdi1234,pos_format,debug_mode);
			    dist_4v_allt_v4(ovdi123_bck,m4_t_bck,m4_p_bck,m4_s_bck,n123_bck,n4_bck,ovdi1234_bck,n,min_dist,max_dist,n_transcripts_bck,ml,report_pos,order_constraint,max_n_transcripts,max_n_ovdi1234,pos_format,ovdi_single_transcript12,temp_ovdi12,ovdi_single_transcript9,temp_ovdi6,ovdi_single_transcript6,m4p_single_transcript,m4s_single_transcript,m123_p_single_transcript,m123_s_single_transcript,0);
			    n1234_bck=n[0];
			    n1234_tr_bck=n[1];
			    if (n[2]<0) { 
			      if (verbose>=2) 
				fprintf(log_file,"\t\t\t\t\t\tn1234_bck=%d\tn1234_tr_bck=%d\n",n1234_bck,n1234_tr_bck);
			    } 
			  } else {
			    n1234_bck=0;n1234_tr_bck=0;
			  }             // close check on n123_bck>0
			} else {
			  fprintf(log_file,"n4_bck=%d\n",n4_bck);
			  n1234_bck=0;n1234_tr_bck=0;
			}               // close check on n4_bck>0
    			compute_pr(n1234_cl,n1234_tr_cl,n1234_bck,n1234_tr_bck,r,p);
			r_oc=r[0];
			r_tr=r[1];
			p_binom_oc=p[0];
			p_binom_tr=p[1];
			  
			if (p_binom_tr > minian_p_tr_4m) {
			  if (verbose>=2) 
			      fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n1234_cl,n1234_tr_cl,n1234_bck,n1234_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
			} else {
			  fprintf(output_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,ml_1,ml_2,ml_3,ml_4,n1234_cl,n1234_tr_cl,n1234_bck,n1234_tr_bck,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
			  if (verbose>=1) {
			    mm=0;
			    fprintf(log_file,"m1\tm2\tm3\tm4\ttr\tp1\tp2\tp3\tp4\ts1\ts2\ts3\ts4\tr_oc\tr_tr\tp_oc\tp_tr\tec\n");
			    while ( (mm<n1234_cl) & (mm<max_n_ovdi12) ) {
			      mm++;
			      curr_tr=ovdi1234_cl[mm][1];
			      fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,curr_tr,ovdi1234_cl[mm][2],ovdi1234_cl[mm][3],ovdi1234_cl[mm][4],ovdi1234_cl[mm][5],ovdi1234_cl[mm][10],ovdi1234_cl[mm][11],ovdi1234_cl[mm][12],ovdi1234_cl[mm][13],r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
			    } // close mm loop to report results
			  } // close check on verbose>=1

			  if (report_pos==1) {
			    fprintf(positions_file,"%d\t%d\t%d\t%d\t%.2f\t%.2f\t%1.2g\t%1.2g\t%d\n",curr_motif1,curr_motif2,curr_motif3,curr_motif4,r_oc,r_tr,p_binom_oc,p_binom_tr,errorcode);
			    mm=0;
			    while ( (mm<n1234_cl) & (mm<max_n_ovdi12) ) {
			      mm++;
			      curr_tr=ovdi1234_cl[mm][1];
			      fprintf(positions_file,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",curr_tr,ovdi1234_cl[mm][2],ovdi1234_cl[mm][3],ovdi1234_cl[mm][4],ovdi1234_cl[mm][5],ovdi1234_cl[mm][10],ovdi1234_cl[mm][11],ovdi1234_cl[mm][12],ovdi1234_cl[mm][13]);
			    }    // close mm loop to report results	
			  }      // close check on report_pos==1		  
			}        // close check on p_binom_tr threshold
		      } else {
			if (verbose>=2) 
			  fprintf(log_file,"n4_cl<minian_tr\n");
		      }          // close check on n1234_cl>minian_tr
		    }            // close check on n4_cl>minian
		  }              // close curr_motif4 loop
		} else {
		  if (verbose>=2)
		    fprintf(log_file,"n3_cl<minian_tr\n");
		}                 // close check on n123_cl>minian_tr
	      }                   // close check on n3_cl>minian_tr 
	    }   // close k loop
	  } else {
	    if (verbose>=2) 
	      fprintf(log_file,"n12_tr_cl<minian_tr\n");
	  }     // close check on n12_tr_cl>minian_tr
	} else {
	  if (verbose>=2) 
	    fprintf(log_file,"n2_cl<minian_tr\n");
	} // close check on n2_cl>minian_tr
      }   // close j loop
    } else {
      if (verbose>=2)
	fprintf(log_file,"n1_cl<minian_tr\n");
    }      // close check on on n1_cl>minian_tr
  }        // close i loop

  ftime(&t2);
  difftime=t2.time-t1.time;
  printf("%% t1=%d sec\t%d msec\n",t1.time,t1.millitm);
  printf("%% t2=%d sec\t%d msec\n",t2.time,t2.millitm);
  printf("%% differences (sec) = %ld\n",difftime);
  if (verbose>0) {
    fprintf(log_file,"%% t1=%d sec\t%d msec\n",t1.time,t1.millitm);
    fprintf(log_file,"%% t2=%d sec\t%d msec\n",t2.time,t2.millitm);
    fprintf(log_file,"%% difference (sec) = %ld\n",difftime);
  }
  fclose(output_file);
  fclose(positions_file);

  /* free memory */
  if (verbose>0)
    fprintf(log_file,"%% freeing memory...");
  printf("freeing memory...");
  printf("m1_t_cl...");fprintf(log_file,"m1_t_cl...");free_ivector(m1_t_cl,1,max_occur);
  printf("m1_p_cl...");fprintf(log_file,"m1_p_cl...");free_ivector(m1_p_cl,1,max_occur);
  printf("m1_s_cl...");fprintf(log_file,"m1_s_cl...");free_cvector(m1_s_cl,1,max_occur);
  printf("m2_t_cl...");fprintf(log_file,"m2_t_cl...");free_ivector(m2_t_cl,1,max_occur);
  printf("m2_p_cl...");fprintf(log_file,"m2_p_cl...");free_ivector(m2_p_cl,1,max_occur);
  printf("m2_s_cl...");fprintf(log_file,"m2_s_cl...");free_cvector(m2_s_cl,1,max_occur);
  if (max_n_motifs>2) {
    printf("m3_t_cl...");fprintf(log_file,"m3_t_cl...");free_ivector(m3_t_cl,1,max_occur);
    printf("m3_p_cl...");fprintf(log_file,"m3_p_cl...");free_ivector(m3_p_cl,1,max_occur);
    printf("m3_s_cl...");fprintf(log_file,"m3_s_cl...");free_cvector(m3_s_cl,1,max_occur);
  }
  if (max_n_motifs>3) {
    printf("m4_t_cl...");fprintf(log_file,"m4_t_cl...");free_ivector(m4_t_cl,1,max_occur);
    printf("m4_p_cl...");fprintf(log_file,"m4_p_cl...");free_ivector(m4_p_cl,1,max_occur);
    printf("m4_s_cl...");fprintf(log_file,"m4_s_cl...");free_cvector(m4_s_cl,1,max_occur);
  }
  printf("m1_t_bck...");fprintf(log_file,"m1_t_bck...");free_ivector(m1_t_bck,1,max_occur);
  printf("m1_p_bck...");fprintf(log_file,"m1_p_bck...");free_ivector(m1_p_bck,1,max_occur);
  printf("m1_s_bck...");fprintf(log_file,"m1_s_bck...");free_cvector(m1_s_bck,1,max_occur);
  printf("m2_t_bck...");fprintf(log_file,"m2_t_bck...");free_ivector(m2_t_bck,1,max_occur);
  printf("m2_p_bck...");fprintf(log_file,"m2_p_bck...");free_ivector(m2_p_bck,1,max_occur);
  printf("m2_s_bck...");fprintf(log_file,"m2_s_bck...");free_cvector(m2_s_bck,1,max_occur);
  if (max_n_motifs>2) {
    printf("m3_t_bck...");fprintf(log_file,"m3_t_bck...");free_ivector(m3_t_bck,1,max_occur);
    printf("m3_p_bck...");fprintf(log_file,"m3_p_bck...");free_ivector(m3_p_bck,1,max_occur);
    printf("m3_s_bck...");fprintf(log_file,"m3_s_bck...");free_cvector(m3_s_bck,1,max_occur);
  }
  if (max_n_motifs>3) {
    printf("m4_t_bck...");fprintf(log_file,"m4_t_bck...");free_ivector(m4_t_bck,1,max_occur);
    printf("m4_p_bck...");fprintf(log_file,"m4_p_bck...");free_ivector(m4_p_bck,1,max_occur);
    printf("m4_s_bck...");fprintf(log_file,"m4_s_bck...");free_cvector(m4_s_bck,1,max_occur);
  }

  printf("motif_lengths");fprintf(log_file,"motif_lengths");free_ivector(motif_lengths,1,n_motifs);
  printf("data_p_cl...");fprintf(log_file,"data_p_cl...");free_imatrix(data_p_cl,1,n_motifs,0,max_occur);
  printf("data_t_cl...");fprintf(log_file,"data_t_cl...");free_imatrix(data_t_cl,1,n_motifs,0,max_occur);
  printf("data_p_bck...");fprintf(log_file,"data_p_bck...");free_imatrix(data_p_bck,1,n_motifs,0,max_occur);
  printf("data_t_bck...");fprintf(log_file,"data_t_bck...");free_imatrix(data_t_bck,1,n_motifs,0,max_occur);
  printf("data_s_cl...");fprintf(log_file,"data_s_cl...");free_imatrix(data_s_cl,1,n_motifs,0,max_occur);
  printf("data_s_bck...");fprintf(log_file,"data_s_bck...");free_imatrix(data_s_bck,1,n_motifs,0,max_occur);

  printf("ovdi12_cl...");fprintf(log_file,"ovdi12_cl...");free_imatrix(ovdi12_cl,1,max_n_ovdi12,1,7);
  printf("ovdi12_bck...");fprintf(log_file,"ovdi12_bck...");free_imatrix(ovdi12_bck,1,max_n_ovdi12,1,7);
  if (max_n_motifs>2) {
    printf("ovdi123_cl...");fprintf(log_file,"ovdi123_cl...");free_imatrix(ovdi123_cl,1,max_n_ovdi123,1,10);
    printf("ovdi123_bck...");fprintf(log_file,"ovdi123_bck...");free_imatrix(ovdi123_bck,1,max_n_ovdi123,1,10);
  } 
  if (max_n_motifs>3) {
    printf("ovdi1234_cl...");fprintf(log_file,"ovdi1234_cl...");free_imatrix(ovdi1234_cl,1,max_n_ovdi1234,1,13);
    printf("ovdi1234_bck...");fprintf(log_file,"ovdi1234_bck...");free_imatrix(ovdi1234_bck,1,max_n_ovdi1234,1,13);
  }
  printf("transcript_lengths_cl/bck...,");free_imatrix(transcript_lengths_cl,1,n_transcripts_cl,1,3);free_imatrix(transcript_lengths_bck,1,n_transcripts_bck,1,3);
  printf("max_n_transcripts_2m....");free_ivector(max_n_transcripts_2m,1,n_transcripts_cl);
  printf("max_n_transcripts_3m....");free_ivector(max_n_transcripts_3m,1,n_transcripts_cl);
  printf("max_n_transcripts_4m....");free_ivector(max_n_transcripts_4m,1,n_transcripts_cl);

  printf("m1/2/3/4p_single_transcript...,");
  free_ivector(m1p_single_transcript,1,max_transcript_occurrences);free_ivector(m2p_single_transcript,1,max_transcript_occurrences);
  free_ivector(m3p_single_transcript,1,max_transcript_occurrences);free_ivector(m4p_single_transcript,1,max_transcript_occurrences);
  printf("m1/2/3/4s_single_transcript...,");
  free_cvector(m1s_single_transcript,1,max_transcript_occurrences);free_cvector(m2s_single_transcript,1,max_transcript_occurrences);
  free_cvector(m3s_single_transcript,1,max_transcript_occurrences);free_cvector(m4s_single_transcript,1,max_transcript_occurrences);
  printf("m12_p/s_single_transcript...,");
  free_ivector(m12_p_single_transcript,1,max_transcript_occurrences);free_cvector(m12_s_single_transcript,1,max_transcript_occurrences);
  printf("m123_p/s_single_transcript...,");
  free_ivector(m123_p_single_transcript,1,max_transcript_occurrences);free_cvector(m123_s_single_transcript,1,max_transcript_occurrences);

  printf("ovdi_single_transcript6/9/12...,");
  free_imatrix(ovdi_single_transcript6,1,max_n_ovdi_single_transcript,1,6);
  free_imatrix(ovdi_single_transcript9,1,max_n_ovdi_single_transcript,1,9);
  free_imatrix(ovdi_single_transcript12,1,max_n_ovdi_single_transcript,1,12);
  printf("temp_ovdi6/9/12...,");
  free_imatrix(temp_ovdi6,1,max_n_ovdi_single_transcript,1,6);
  free_imatrix(temp_ovdi9,1,max_n_ovdi_single_transcript,1,9);
  free_imatrix(temp_ovdi12,1,max_n_ovdi_single_transcript,1,12);
  printf("all_ovdi6\n");
  free_imatrix(all_ovdi6,1,max_n_ovdi_single_transcript,1,6);

  fprintf(log_file,"\n");
  printf("\n\nready. thanks for flying the friendly skies!\n");

  fclose(log_file);
  printf("system_status=0\n");

  return 0;
}

void compute_pr(int n_cl,int n_tr_cl,int n_bck,int n_tr_bck,float *r,double *p)
{
  // compute p_oc, p_tr, r_oc, r_tr
  float r_oc;       // ratio of co-occurrences
  float r_tr;       // ratio of transcripts with co-occurrences
  double p_oc;      // binomial cumulative probability of obtaining x>=n_cl given p_bck_oc
  double p_tr;      // binomial cumulative probability of obtainint x>=n_tr_cl given p_bck_tr;
  double p_bck_oc;     // probability of co-occurrence in background
  double p_bck_tr;     // probability of co-occurrence per transcript in background
  int itemp;

  if (n_bck>0) {
    r_oc=( (float)n_cl / (float)n_bck ) *r_transcripts_cl_bck;
    p_bck_oc=(double)n_bck/(double)n_transcripts_bck;
  } else {
    r_oc=( (float) n_cl / (float)default_min ) *r_transcripts_cl_bck;
    p_bck_oc=(double)default_min/(double)n_transcripts_bck;
  }

  if (n_tr_bck>0) {
    r_tr=( (float)n_tr_cl/ (float)n_tr_bck) *r_transcripts_cl_bck;
    p_bck_tr=(double)n_tr_bck/(double)n_transcripts_bck;
  } else {
    r_tr=( (float) n_tr_cl/ (float)default_min ) *r_transcripts_cl_bck;
    p_bck_tr=(double)default_min/(double)n_transcripts_bck;
  }
	
  if (verbose>=2) 
    fprintf(log_file,"calling binocdf with parameters %d %d %1.2g (occurrences)\n",n_transcripts_cl,n_cl,p_bck_oc);
		      
  if (p_bck_oc>=1.0) { 
    p_oc=0.0;    // for sure there cannot be less than n_cl in the cluster (all should show up) 
  } else {
    if (p_bck_oc==0.0) {
      p_oc=1.0;  // for sure there should be less than n_cl in the cluster (none should show up)
    } else {
      if (n_cl<n_transcripts_cl) {
	itemp=n_cl-1;  // because we want the strict x<n_cl so that we can consider n_cl within the complement
	p_oc=binocdf_mat(n_transcripts_cl,itemp,p_bck_oc);
      } else {
	p_oc=1.0;  // for sure there should be less than n_cl
      }
    }
  }

  if (verbose>=2) 
    fprintf(log_file,"calling binocdf with parameters %d %d %1.2g (transcripts)\n",n_transcripts_cl,n_tr_cl,p_bck_tr);
		      
  if (p_bck_tr>=1.0) {
    p_tr=0.0;
  } else {
    if (p_bck_tr==0.0) {
      p_tr=1.0;
    } else {
      itemp=n_tr_cl-1;    // because we want the strict x<n_cl so that we can consider the actual value n_tr_cl within the complement
      p_tr=binocdf_mat(n_transcripts_cl,itemp,p_bck_tr);
    }
  }

  p_oc=1.0-p_oc;
  p_tr=1.0-p_tr;
  if (p_oc<0.0) 
    p_oc=0.0;
  if (p_tr<0.0) 
    p_tr=0.0;

  r[0]=r_oc;
  r[1]=r_tr;
  p[0]=p_oc;
  p[1]=p_tr;
}

int compute_min_dist(int *ml,int n_motifs,float min_dist_factor) 
{
  
  int i;
  int md=1000000;
  int max_overlap;

  i=0;
  while (i<n_motifs) {
    i++;
    max_overlap=(int)floor(ml[i]*min_dist_factor);
    if (max_overlap<md)
      md=max_overlap;
  }

  return md;
}

void usage(void)
{
  printf("cooc_4m_v%d.c\n",prog_version);
  printf("study co-occurrences of up to 3 motifs\n");
  printf("usage:\n");
  printf("cooc_4m_v%d <filename_pos_cl> <filename_pos_bck> <filename_strands_cl> <filename_strands_bck> <n_motifs> <n_transcripts_cl> <n_transcripts_bck> (<max_dist> <minian_tr> <minian_p_tr> <min_dist_factor> <max_ul> <max_el> <max_il> <init_i,j,k> <finit_i,j,k> <file_label> <report_pos> <order_constraint> <sortpos> <max_n_motifs> <transcript_lengths_filename_cl> <transcript_lengths_filename_bck> <filteron> <pos_format> <verbose>)\n",prog_version);
  printf("\n\n<max_dist>\tmaximum distance between adjacent motifs [default=100]\n");
  printf("<minian_tr>\t(optional, default=4) minimum number of transcripts in which there is co-occurrence in the cluster set\n");
  printf("<minian_p_tr>\t(optional, default=0.01) maximum p values (this is then corrected by the number of combinations\n");
  printf("<min_dist_factor>\t(optional, default=0.5) separation between motifs = min_dist_factor * motif_length\n");
  printf("<max_n_motifs>\tif max_n_motifs=2, then process only pairs, if max_n_motifs=3, then process only triplets, if max_n_motifs=4, process quadruplets [default=2]\n");
  printf("<filteron>\t1 to filter the positions during reading and 0 otherwise [default=0]\n");
  printf("<pos_format>\t1 for positions with respect to TSS and -1 for positions with respect to 5' end [default=+1]\n");
  printf("\n\noutput is stored in log.txt (log file), out.txt (list of motifs and statistical significance) and pos.txt (position in each transcript) \n");

  gk();
  exit(1);
}

/* new in version 7
 * 01-31-2004: call load_and_sort_scan_data_v3
 * 01-31-2004: added data_format and verbose as input to load_and_sort_scan_data_v3
 * 01-31-2004: convert positions to TSS format
 * 02-05-2004: added a debug_mode
 */

/* new in version 6
 * 01-01-2004: correct bugs to be able to set order_constraint=1
 * 01-30-2004: added input filteron to decide whether to filter the positions during reading
 * 01-31-2004: added pos_format +1 for positions with respect to TSS (default), -1 for positions with respect to 5' end
 * see below code for full history information
 */

/* new in version 5
 * trying to study where the time committment is
 */

/* new in version 4
 * 11-04-2003: add system status at the end
 * 11-04-2003: added neg2pos,sortpos
 * 11-04-2003: include report of memory requirements
 * 11-04-2003: add report_pos to optionally report positions
 * 11-04-2003: use infotext,temptext for output
 * 11-04-2003: add processing of 1 motif clusters also
 * 11-04-2003: tough on the thresholds
 * 11-04-2003: accept file label to be able to run multiple copies in parallel
 * 11-04-2003: use more recent version with curr_arg for input of parameters (add max_ul,max_el,max_il)
 * 11-01-2003: allow the usage of specific sequence segments
 * 11-01-2003: store the max_n_transcripts
 * 11-01-2003: use load_and_sort_scan_data from crm_lib.c
 */

/* project goals:
 * extend up to 4 motifs
 * to save time, compute also 2m and 3m co-occurrences
 * allow overlap of motifs
 * change the output format to remove all this nonsense tabs
 * change the selection for output, have threshold dependent on the binomial probability, adjust threshold depending on the number of motifs
 */

/* new in version 3:
 * 09-14-2003: use of n[3] to indicate early stop in reporting occurrences due to overflow
 * 09-13-2003: added max_n_ovdi in call to dist_2v_allt, dist_3v_allt, dist_4v_allt
 * 08-29-2003: more realistic file names
 * call latest release of cooc_single_lib.c where we store the intermediate computations to make it faster (as we did in the original cooc_4m.c!)
 */

/* new in version 2:
   process on a transcript by transcript basis for each motif combination
 */

/* 02-28-2003: changed float to double variables for the computation of binocdf */
/* 03-03-2003: subtracted 1 from n in the call to binocdf_mat to compute the cumulative accurately */
/* 03-03-2003: last call to dist2v in bckgrnd modified dist2v_noovdi to just compute the n */
/* 03-14-2003: add positions file as output */
/* 03-19-2003: make it >= minian_tr */
/* 06-25-2003: accept negative positions (corresponding to minus strand). here we do not take strand into account for computations, therefore, we just convert these to positive positions */

