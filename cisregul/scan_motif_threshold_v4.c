/* scan_motif_threshold.c
 * version 4
 * given a motif weight matrix and a background model, compute the threshold by estimating the false positive percentiles
 */

/* new in version 4
 * 07-12-2003 do not compute logarithm when evaluating score (this is now done when estimating the weight matrix 
 * 07-12-2003 do not "fudge" motif to avoid 0s. this is all done now when estimating the weight matrix
 */

/* new in version 3
 * allow scanning of forward as well as reverse strand
 */

/* new in version 2 
 * 05-06-2003 print out 4 significan digits 
 * 05-06-2003 modify the weight matrix if there are zeros according to the number of sequences used to compute it
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "my_config.c"
#include "../lib/nr.c"
#include "../lib/my_math_package.c"
#include "regul_search_methods_v7.c"

float max_wm_score(float **weight_matrix,int motif_length);
float min_wm_score(float **weight_matrix,int motif_length);
void scan_motif_dist(float **weight_matrix,char *sequence,int motif_length,int seq_length,float *scores);
void usage(void);

int prog_version=4;

main(int argc, char *argv[])
{
  int i,j;
  float **weight_matrix;                     // matrix with frequencies in each position
  //float **n_matrix;                          // matrix with number of entries in each position
  FILE *weight_matrix_file;
  char weight_matrix_filename[200];
  char background_sequence[11000];
  //char *background_sequence;
  char *num_background_sequence;
  char *temp_num_background_sequence;
  FILE *background_sequence_file;
  char background_sequence_filename[200];
  int motif_length;
  int sequence_length=0;
  float temp_float;
  float *scores;
  float *rc_scores;
  float *all_scores;
  int n_scores;
  int n_scores_all;
  float p1=0.99;
  float p2=0.999;
  float p3=0.9999;
  float p4=0.99999;
  float p5=0.999999;
  int p1_index;
  int p2_index;
  int p3_index;
  int p4_index;
  int p5_index;
  float mean_score;
  float min_score;
  float max_score;
  float std_score;
  float threshold1;
  float threshold2;
  float threshold3;
  float threshold4;
  float threshold5;
  int temp_sequence_length;
  int verbose;
  FILE *output_file;
  char output_filename[200];
  float max_possible_score;
  float min_possible_score;
  int n_seqs4motif;
  //int anyzeros=0;                     // if anyzeros==1 then there is at least one zero in the weight matrix and we recompute it
  int n;                              // temporary integer for the number of sequences used in the weight matrix computation
  char *revcomp_upstream_sequence;          // reverse complementary sequence
  int rc=0;                           // rc=1 to analyse the reverse complementary sequence as well

  if (argc<5) {
    usage();
  }

  sprintf(background_sequence_filename,argv[1]);
  sprintf(weight_matrix_filename,argv[2]);
  motif_length=atoi(argv[3]);
  n_seqs4motif=atoi(argv[4]);
  if (argc>5) 
    rc=atoi(argv[5]);
  if (argc>6) {
    verbose=atoi(argv[6]);
  } else {
    verbose=0;
  }

  printf("background_sequence_filename=%s\n",background_sequence_filename);
  printf("weight_matrix_filename=%s\n",weight_matrix_filename);
  printf("motif_length=%d\n",motif_length);
  printf("n_seqs4motif=%d\n",n_seqs4motif);
  printf("rc=%d\n",rc);
  printf("verbose=%d\n",verbose);

  num_background_sequence=cvector(1,10000000);
  temp_num_background_sequence=cvector(1,1000000);
  printf("reading background sequence %s\n",background_sequence_filename);
  background_sequence_file=fopen(background_sequence_filename,"r");
  if (!background_sequence_file) {
    printf("ERROR! I could not find the file %s\n",background_sequence_filename);
    exit(1);
  }
  j=0;
  //printf("weight_matrix_filename=%s\n",weight_matrix_filename);
  while (feof(background_sequence_file)==0) {
    j++;
    fscanf(background_sequence_file,"%s\n",&background_sequence);
    temp_sequence_length=strlen(background_sequence);
    seq2num(background_sequence,temp_sequence_length,temp_num_background_sequence,0);
    /* begin debug
       for (i=1;i<=10;i++) 
       printf("%d ",temp_num_background_sequence[i]);
       printf("\n");
       exit(1);
       end debug 
    */
    /* begin debug
       for (i=1;i<=temp_sequence_length;i++) {
       if (temp_num_background_sequence[i]<=0) {
       printf("i=%d\tseq=%d\n",i,temp_num_background_sequence[i]);
       exit(1);
       }
       } 
       end debug
    */
    for (i=1;i<=temp_sequence_length;i++) {
      num_background_sequence[i+sequence_length]=temp_num_background_sequence[i];
    }
    sequence_length+=temp_sequence_length;
    //printf("j=%d\tweight_matrix=%s\n",j,weight_matrix_filename);
  }
  printf("read %d background sequences, length %d\n",j,sequence_length);
  fclose(background_sequence_file);
  //printf("weight_matrix_filename=%s\n",weight_matrix_filename);
 
  weight_matrix=matrix(1,4,1,motif_length);
  //n_matrix=matrix(1,4,1,motif_length);
  printf("reading the weight matrix file %s\n",weight_matrix_filename);  
  weight_matrix_file=fopen(weight_matrix_filename,"r");
  if (!weight_matrix_file) {
    printf("ERROR! I could not find the file %s\n",weight_matrix_filename);
    exit(1);
  }
  for (i=1;i<=4;i++) {
    for (j=1;j<=motif_length;j++) {
      fscanf(weight_matrix_file,"%f\n",&temp_float);
      weight_matrix[i][j]=temp_float;
      //n_matrix[i][j]=n_seqs4motif*temp_float;
      //if (temp_float==0) 
      //anyzeros=1;
    }
  }
  fclose(weight_matrix_file);

  //if (anyzeros==1) {
  // there is at least one zero in the weight matrix, recompute
  //  for (j=1;j<=motif_length;j++) {
  //  n=n_seqs4motif;
  //  for (i=1;i<=4;i++) {
  //if (n_matrix[i][j]==0) {
  //n++;
  //  n_matrix[i][j]=1;
  //}
  //  }
  // if (n_seqs4motif != n) {
  //for (i=1;i<=4;i++) 
  //  weight_matrix[i][j]=n_matrix[i][j]/(float)n;
  //  }
  //}
  //}

  max_possible_score=max_wm_score(weight_matrix,motif_length);
  min_possible_score=min_wm_score(weight_matrix,motif_length);
  printf("max_possible_score=%.4f\n",max_possible_score);
  printf("min_possible_score=%.4f\n",min_possible_score);

  n_scores=sequence_length-motif_length+1;
  scores=vector(1,n_scores);
  for (i=1;i<=n_scores;i++)
    scores[i]=1;
  scan_motif_dist(weight_matrix,num_background_sequence,motif_length,sequence_length,scores);

  if (rc==1) {
    rc_scores=vector(1,n_scores);
    for (i=1;i<=n_scores;i++) 
      rc_scores[i]=1;
    revcomp_upstream_sequence=cvector(1,sequence_length);
    revcomp(num_background_sequence,revcomp_upstream_sequence,sequence_length);

    /* 
       start debug
       for (i=1;i<=10;i++)
       printf("%d\t",revcomp_upstream_sequence[i]);
       printf("\n");
       
       free_matrix(weight_matrix,1,4,1,6);
       free_cvector(num_background_sequence,1,10000000);
       free_vector(scores,1,n_scores);
       free_cvector(temp_num_background_sequence,1,1000000);
       free_vector(rc_scores,1,n_scores);
       free_cvector(revcomp_upstream_sequence,1,sequence_length);
       exit(1);
       end debug
    */

    scan_motif_dist(weight_matrix,revcomp_upstream_sequence,motif_length,sequence_length,rc_scores);
    n_scores_all=2*n_scores;
    all_scores=vector(1,n_scores_all);
    for (i=1;i<=n_scores;i++) {
      all_scores[i]=scores[i];
      j=n_scores+i;
      all_scores[j]=rc_scores[i];
    }
  } else {
    n_scores_all=n_scores;
    all_scores=vector(1,n_scores_all);
    for (i=1;i<=n_scores;i++)
      all_scores[i]=scores[i];
  }

  mean_score=meanf(all_scores,n_scores_all);
  std_score=stdf(all_scores,mean_score,n_scores_all);
  min_score=minf(all_scores,n_scores_all);
  max_score=maxf(all_scores,n_scores_all);

  printf("n_scores=%d\n",n_scores_all);
  printf("mean_score=%.4f\n",mean_score);
  printf("std_score=%.4f\n",std_score);
  printf("min_score=%.4f\n",min_score);
  printf("max_score=%.4f\n",max_score);

  p1_index=(int)ceil(p1*n_scores_all);
  p2_index=(int)ceil(p2*n_scores_all);
  p3_index=(int)ceil(p3*n_scores_all);
  p4_index=(int)ceil(p4*n_scores_all);
  p5_index=(int)ceil(p5*n_scores_all);

  threshold1=select_klargest(p1_index,n_scores_all,all_scores);
  threshold2=select_klargest(p2_index,n_scores_all,all_scores);
  threshold3=select_klargest(p3_index,n_scores_all,all_scores);
  threshold4=select_klargest(p4_index,n_scores_all,all_scores);
  threshold5=select_klargest(p5_index,n_scores_all,all_scores);

  printf("thresholds\n");
  printf("p0.99=%.4f\n",threshold1);
  printf("p0.999=%.4f\n",threshold2);
  printf("p0.9999=%.4f\n",threshold3);
  printf("p0.99999=%.4f\n",threshold4);
  printf("p0.999999=%.4f\n",threshold5);

  if (verbose) {
    sprintf(output_filename,"%s.out",weight_matrix_filename);
    output_file=fopen(output_filename,"w");
    for (i=1;i<=n_scores;i++) {
      fprintf(output_file,"%.4f\n",scores[i]);
    }
    fclose(output_file);
  }

  free_matrix(weight_matrix,1,4,1,6);
  free_cvector(num_background_sequence,1,10000000);
  free_vector(scores,1,n_scores);
  free_cvector(temp_num_background_sequence,1,1000000);
  free_vector(all_scores,1,n_scores_all);
  if (rc==1) {
    free_vector(rc_scores,1,n_scores);
    free_cvector(revcomp_upstream_sequence,1,sequence_length);
  }
  return 0;
}

void scan_motif_dist(float **weight_matrix,char *sequence,int motif_length,int seq_length,float *scores)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report the whole distribution of scores
   * this is used to derive an empirical threshold based on this distribution
   */

  /* 07-12-2003: do not take the logarithm, this was already done when estimating the weight matrix */

  int i,j;
  int i_stop;
  char c;
  float p;
  float temp_s;
  int n;

  i_stop=seq_length-motif_length;

  /* begin debug
     for (i=1;i<=4;i++) {
     for (j=1;j<=motif_length;j++) 
     printf("%.4f\t",weight_matrix[i][j]);
     printf("\n");
     } 
     end debug */
  
  for (i=1;i<=i_stop;i++) {
    //for (i=1;i<100;i++) { // debug
    temp_s=0;
    for (j=1;j<=motif_length;j++) {
      c=sequence[i+j-1];
      if ( (c>-1) | (c<5) ) {
	p=weight_matrix[c][j];
	//temp_s+=log(p);
	temp_s+=p;
      } else {
	temp_s+=-7.6;     // this correspond to a p value of 5x10^-4, penalize for non a-c-g-t nucleotides
      }
    }
    
    scores[i]=temp_s;
  }
}

float max_wm_score(float **weight_matrix,int motif_length)
{
  /* given a weight matrix, compute the maximum possible score */
  
  /* 07-12-2003: do not compute log, this is done outside of this program */

  int i,j;
  float score;
  float col_max;

  score=0;
  for (i=1;i<=motif_length;i++) {
    col_max=-1;
    for (j=1;j<=4;j++) {
      if (weight_matrix[j][i]>col_max)
	col_max=weight_matrix[j][i];
    }
    //score+=log(col_max);
    score+=col_max;
  }

  return score;
}

float min_wm_score(float **weight_matrix,int motif_length)
{
  /* given a weight matrix, compute the maximum possible score */
  
  /* 07-12-2003: do not compute log, this is done outside of this program */
  int i,j;
  float score;
  float col_min;

  score=0;
  for (i=1;i<=motif_length;i++) {
    col_min=1000;
    for (j=1;j<=4;j++) {
      if (weight_matrix[j][i]<col_min)
	col_min=weight_matrix[j][i];
    }
    //score+=log(col_min);
    score+=col_min;
  }

  return score;
}

void usage(void)
{
  printf("scan_motif_threshold_v%d.c\n",prog_version);
  printf("usage\n");
  printf("scan_motif_threshold_v%d <background_sequence> <weight_matrix> <motif_length> <n_seqs4motif> <rc> <verbose>\n",prog_version);
  printf("\nn_seqs4motif is the number of sequences used to compute the motif\n");
  printf("\nverbose=1 to save all the distribution of scores\n");

  printf("\nnew in version 3: allow scanning of forward as well as reverse strand\n");
  printf("new in version 4:\n07-12-2003 do not compute logarithm when evaluating score (this is now done when estimating the weight matrix\n");
  printf("07-12-2003 do not 'fudge' motif to avoid 0s. this is all done now when estimating the weight matrix\n");

  exit(1);
}
