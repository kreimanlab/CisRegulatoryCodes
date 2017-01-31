/* regul_search_methods_v6.c
 * library of methods used in the search for regulatory regions 
 */

/* new in version 2
 * to avoid memory problems, vectors are passed as arguments to the functions defined here
 * wrote new function count_perfect_matches and incorporated this within count_motif_occurrences
 */

/* new in version 3
 * attempt to optimize count_motif_occurrences_prob
 */

/* new in version 4
 * use binary files in the computation of statistical values
 */

/* new in version 5
 * add a search function
 * get_index_from_textmotif
 */

/* new in version 6
 * add simple motif comparison
 * removed motif_length parameter from quick_2mismatch and quick_mismatch_v2
 * removed mism1,mism2 parameters from motif_stats_bin_v2
 * removed verbose parameter from count_motif_occurrences
 */

/* new in version 7
 * add count_perfect_matches_v2 (returning the positions of the matches)
 * 070102 add count_perfect_matches_v3 (returning the positions as unsigned short)
 * 02-19-2003: added push_double_vector
 */

/* lower case ascii codes */
#define AA 97
#define CC 99
#define GG 103
#define TT 116
#define NN 110

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
//#include "../lib/my_math_package.c"

/* function definitions */
void seq2num(char *sequence,int sequence_length,char *numsequence,int init_i);
char char2num(char c);
void revcomp(char *numseq,char *revcomp,int n);
float compare_aligned_seqs(char *s1,char *s2);
float count_motif_occurrences_prob(char *motif,char *upstream_sequence,float threshold);
int compare_2motifs_v1(char *motif1,char *motif2,int motif_length,int max_overlap,int mism1_m1,int mism2_m1,int mism1_m2,int mism2_m2);
int count_matches_v1(char *upstream_sequence,char *motif,int *pos,int mism_pos);
int count_perfect_matches(char *upstream_sequence,char *motif);
int count_perfect_matches_v2(char *upstream_sequence,char *motif,int *pos);
int count_perfect_matches_v2c(char *upstream_sequence,char *motif,char *motifc,int *pos);
int count_perfect_matches_v3(char *upstream_sequence,char *motif,unsigned short *pos);
int faster_find_pos(float *vec,float value,int vec_length);
int find_pos(float *vec,float value,int vec_length);
int get_complement_from_index(int index);
int get_index_from_nummotif(char *motif);
int get_index_from_textmotif(char *motif,int motif_length);
int get_motiflength_from_index(int index);
int total_motifs(int length);
void compute_residuals(float *residuals,float top_f,float *top_motif_count,int n_transcripts);
void convert_nummotif2text(char *nummotif,char *textmotif);
void convert_text2nummotif(char *motif,char *nummotif);
void get_complement_from_textmotif(char *motif,char *motifc,int motif_length);
void get_motif(int index,char *motif);
void get_nummotif_from_index(int index,char *motif);
void get_textmotif_from_index(int index,char *motif);
void get_textmotif_from_index_v2(int index,char *motif,int mism1,int mism2);
void mismatch(int index,int motif_length,int *mismatch_indices);
void motif2_stats(char *motif1,char* motif2,float *expression,char *input_filename_upstream,FILE *output_file,int n_transcripts);
void motif2_stats_withcounts(char *motif1,char* motif2,float *expression,int *motif_count1,int *motif_count2,FILE *output_file,int n_transcripts);
void motif_initialization(char *motif,int motif_length);
void motif_stats(char *motif,float *expression,char *input_filename_upstream,FILE *output_file,int n_transcripts);
void motif_stats_bin(int motif_index,float *expression,char *filename_motifcount_match,char *filename_motifcount_mismatch,FILE *output_file,int n_transcripts);
void motif_stats_bin_v2(int motif_index,float *expression,char *filename_motifcount_match,int *motif_count,FILE *output_file,int n_transcripts);
void motif_stats_withcounts(char *motif,float *expression,int *motif_count,FILE *output_file,int n_transcripts);
void print_motif(char *motif,float *expression,char *input_filename_upstream,int n_around,FILE *output_file);
void quick_2mismatch(int index,int *mismatch_indices);
void quick_mismatch(int index,int *mismatch_indices);
void quick_mismatch_v2(int index,int **mismatch_indices);
void push_int_vector(int new_value,int *ivec,int value_pos,int vec_length);
void push_float_vector(float new_value,float *ivec,int value_pos,int vec_length);
void push_double_vector(double new_value,double *ivec,int value_pos,int vec_length);
void read_fasta_seq(FILE *file_handle,char *info,char *seq);
void read_fasta_seq2(FILE *file_handle,char *info,char *seq);
int faster_find_dpos(double *vec,double value,int vec_length);
int find_dpos(double *vec,double value,int vec_length);
char *StrnCpy(char *dest,const char *src,size_t n);
char *remove_carriage_return(char *str_noret,char *str);
void num2seq(char *numsequence,int sequence_length,char *sequence,int init_i);
char num2char(char nc);
float score_motif(float **weight_matrix,char *sequence,int motif_length);
void load_weight_matrix(char *weight_matrix_filename,int motif_length,float **weight_matrix);
int scan_motif(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,int rc);
int scan_motif2(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc);
int scan_motif3(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,FILE *log_file,int seq_nr);
int scan_motif4(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,FILE *log_file,int seq_nr);
int scan_motif5(float **weight_matrix,float *mean_weight_matrix,char *sequence,char *revcomp_seq,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,FILE *log_file,int seq_nr);
int scan_motif6(float **weight_matrix,float *mean_weight_matrix,char *sequence,char *revcomp_seq,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,float *cum_threshold,FILE *log_file,int seq_nr);
void temp_threshold(float *max_weight_matrix,int motif_length,float score_threshold,float *thresh);
void weight_matrix_max(float **weight_matrix,int motif_length,float *max_weight_matrix);
void weight_matrix_mean(float **weight_matrix,int motif_length,float *mean_weight_matrix);

int scan_motif6(float **weight_matrix,float *mean_weight_matrix,char *sequence,char *revcomp_seq,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,float *cum_threshold,FILE *log_file,int seq_nr)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * new in version 2: 08-16-2003 return all positions as positive, return strand information (+ -> 1, - -> 0) in a separate vector
   * new in version 3: also print information to log file 
   * new in version 4: 
   * 12-28-2003 only report if distance between adjacent occurrences is > min_dist (set to 0 to report all)
   * 12-28-2003 if seq_nr>0, then report also the actual sequences to the log_file
   * new in version 5: 
   * 02-24-2004: do not allocate space for revcomp here, do it outside this code
   * 02-24-2004: mean_weight_matrix as input (mean for each nucleotide, use this when there is "n")
   * new in version 6:
   * 03-10-2004: added cum_threshold to try to speed things up
   */

  int i,j,k;                    // loop counters
  int i_stop;                   // stop processing when this point in the sequence is reached
  char c;                       // single nucleotide (between 0 and 4, penalize otherwise)
  float p;                      // current entry in the weigh matrix
  int n;                        // total number of occurrences (this is the output of this function)
  float current_score;          // current score
  //char *revcomp_seq;            // reverse complement of the sequence
  int curr_strand;              // strand (+1 or 0)
  char *word;                   // used to report the actual sequence if seq_nr>0
  char *numword;                // used to report the actual sequence if seq_nr>0
  int previous_position;        // previous position

  i_stop=seq_length-motif_length;
  n=0;
  if (seq_nr>0) {
    numword=cvector(1,motif_length);
    word=cvector(0,motif_length-1);
  }

  previous_position=-motif_length;     // to process from the very beginning
  curr_strand=1;
  for (i=1;i<=i_stop;i++) {
    if ( (i-previous_position)>min_dist ) {
      current_score=0.0;
      j=1;
      while ( (j<=motif_length) & (current_score>=cum_threshold[j]) ) {
	c=sequence[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=mean_weight_matrix[j];     
	}
	j++;
      }
      if (current_score>=threshold) {
	previous_position=i;
	n++;
	pos[n]=seq_length-i;                       // measure all positions with respect to 3' end
	strand[n]=curr_strand;
	
	if (seq_nr>0) {
	  for (j=1;j<=motif_length;j++) 
	    numword[j]=sequence[i+j-1];
	  sprintf(word,"");
	  num2seq(numword,motif_length,word,1);
	  fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	}      // close check on seq_nr>0
      }        // close check on score
    }          // close check on distance to previous occurrence
  }            // close i loop

  curr_strand=0;
  previous_position=-motif_length;     // to process from the very beginning
  if (rc==1) {
    //revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    for (i=1;i<=i_stop;i++) {
      if ( (i-previous_position)>min_dist ) {      
	current_score=0.0;
	j=1;
	while ( (j<=motif_length) & (current_score>=cum_threshold[j]) ) {
	  c=revcomp_seq[i+j-1];
	  if ( (c>0) & (c<5) ) {
	    p=weight_matrix[c][j];
	    current_score+=p;
	  } else {
	    current_score+=mean_weight_matrix[j];     
	  }
	  j++;
	}
	if (current_score>=threshold) {
	  previous_position=i;
	  n++;
	  pos[n]=i;                          // positions in the minus strand are indicated as negative
	  strand[n]=0;
	  
	  if (seq_nr>0) {
	    for (j=1;j<=motif_length;j++) 
	      numword[j]=revcomp_seq[i+j-1];
	    sprintf(word,"");
	    num2seq(numword,motif_length,word,1);
	    fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	  }       // close check on seq_nr>0
	}         // close check on score
      }           // close check on distance to previous occurrence
    }             // close i loop
    //free_cvector(revcomp_seq,1,seq_length);
  }               // close check on rc==1

  if (seq_nr>0) {
    free_cvector(numword,1,motif_length);
    free_cvector(word,0,motif_length-1);
  }

  return n;
}

void weight_matrix_max(float **weight_matrix,int motif_length,float *max_weight_matrix) 
{
  /* given a weight matrix, compute the maximum at each position */
  int i;
  int j;
  float current_max;

  for (i=1;i<=motif_length;i++) {
    current_max=-1000000;
    for (j=1;j<=4;j++) {
      if (weight_matrix[j][i]>current_max)
	current_max=weight_matrix[j][i];
    }
    max_weight_matrix[i]=current_max;
  }
}

void weight_matrix_mean(float **weight_matrix,int motif_length,float *mean_weight_matrix) 
{
  /* given a weight matrix, compute the mean at each position */
  int i;
  int j;
  float current_mean;

  for (i=1;i<=motif_length;i++) {
    current_mean=0;
    for (j=1;j<=4;j++) 
      current_mean+=weight_matrix[j][i];
    mean_weight_matrix[i]=current_mean;
  }
}

void temp_threshold(float *max_weight_matrix,int motif_length,float score_threshold,float *thresh)
{
  /* given a vector max_weight_matrix and a score threshold, compute the cumulative value that must be reached at each position to be able to pass the score */
  int i,j;
  float cumsum;

  thresh[1]=-1;                       // so that initializing score=0 will enter the loop because score>thresh[1]
  for (i=2;i<=motif_length;i++) {
    cumsum=0.0;
    for (j=i;j<=motif_length;j++) 
      cumsum+=max_weight_matrix[j];
    thresh[i]=score_threshold-cumsum;
  }
}

int scan_motif5(float **weight_matrix,float *mean_weight_matrix,char *sequence,char *revcomp_seq,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,FILE *log_file,int seq_nr)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * new in version 2: 08-16-2003 return all positions as positive, return strand information (+ -> 1, - -> 0) in a separate vector
   * new in version 3: also print information to log file 
   * new in version 4: 
   * 12-28-2003 only report if distance between adjacent occurrences is > min_dist (set to 0 to report all)
   * 12-28-2003 if seq_nr>0, then report also the actual sequences to the log_file
   * new in version 5: 
   * 02-24-2004: do not allocate space for revcomp here, do it outside this code
   * 02-24-2004: mean_weight_matrix as input (mean for each nucleotide, use this when there is "n")
   */

  int i,j,k;                    // loop counters
  int i_stop;                   // stop processing when this point in the sequence is reached
  char c;                       // single nucleotide (between 0 and 4, penalize otherwise)
  float p;                      // current entry in the weigh matrix
  int n;                        // total number of occurrences (this is the output of this function)
  float current_score;          // current score
  //char *revcomp_seq;            // reverse complement of the sequence
  int curr_strand;              // strand (+1 or 0)
  char *word;                   // used to report the actual sequence if seq_nr>0
  char *numword;                // used to report the actual sequence if seq_nr>0
  int previous_position;        // previous position

  i_stop=seq_length-motif_length;
  n=0;
  if (seq_nr>0) {
    numword=cvector(1,motif_length);
    word=cvector(0,motif_length-1);
  }

  previous_position=-motif_length;     // to process from the very beginning
  curr_strand=1;
  for (i=1;i<=i_stop;i++) {
    if ( (i-previous_position)>min_dist ) {
      current_score=0.0;
      for (j=1;j<=motif_length;j++) {
	c=sequence[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=mean_weight_matrix[j];     
	}
      }
      if (current_score>=threshold) {
	previous_position=i;
	n++;
	pos[n]=seq_length-i;                       // measure all positions with respect to 3' end
	strand[n]=curr_strand;
	
	if (seq_nr>0) {
	  for (j=1;j<=motif_length;j++) 
	    numword[j]=sequence[i+j-1];
	  sprintf(word,"");
	  num2seq(numword,motif_length,word,1);
	  fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	}      // close check on seq_nr>0
      }        // close check on score
    }          // close check on distance to previous occurrence
  }            // close i loop

  curr_strand=0;
  previous_position=-motif_length;     // to process from the very beginning
  if (rc==1) {
    //revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    for (i=1;i<=i_stop;i++) {
      if ( (i-previous_position)>min_dist ) {      
	current_score=0.0;
	for (j=1;j<=motif_length;j++) {
	  c=revcomp_seq[i+j-1];
	  if ( (c>0) & (c<5) ) {
	    p=weight_matrix[c][j];
	    current_score+=p;
	  } else {
	    current_score+=mean_weight_matrix[j];     
	  }
	}
	if (current_score>=threshold) {
	  previous_position=i;
	  n++;
	  pos[n]=i;                          // positions in the minus strand are indicated as negative
	  strand[n]=0;
	  
	  if (seq_nr>0) {
	    for (j=1;j<=motif_length;j++) 
	      numword[j]=revcomp_seq[i+j-1];
	    sprintf(word,"");
	    num2seq(numword,motif_length,word,1);
	    fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	  }       // close check on seq_nr>0
	}         // close check on score
      }           // close check on distance to previous occurrence
    }             // close i loop
    //free_cvector(revcomp_seq,1,seq_length);
  }               // close check on rc==1

  if (seq_nr>0) {
    free_cvector(numword,1,motif_length);
    free_cvector(word,0,motif_length-1);
  }

  return n;
}

int scan_motif4(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,int min_dist,FILE *log_file,int seq_nr)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * new in version 2: 08-16-2003 return all positions as positive, return strand information (+ -> 1, - -> 0) in a separate vector
   * new in version 3: also print information to log file 
   * new in version 4: 12-28-2003 only report if distance between adjacent occurrences is > min_dist (set to 0 to report all)
   * 12-28-2003 if seq_nr>0, then report also the actual sequences to the log_file
   */

  int i,j,k;                    // loop counters
  int i_stop;                   // stop processing when this point in the sequence is reached
  char c;                       // single nucleotide (between 0 and 4, penalize otherwise)
  float p;                      // current entry in the weigh matrix
  int n;                        // total number of occurrences (this is the output of this function)
  float current_score;          // current score
  char *revcomp_seq;            // reverse complement of the sequence
  int curr_strand;              // strand (+1 or 0)
  char *word;                   // used to report the actual sequence if seq_nr>0
  char *numword;                // used to report the actual sequence if seq_nr>0
  int previous_position;        // previous position

  i_stop=seq_length-motif_length;
  n=0;
  if (seq_nr>0) {
    numword=cvector(1,motif_length);
    word=cvector(0,motif_length-1);
  }

  previous_position=-motif_length;     // to process from the very beginning
  curr_strand=1;
  //printf("\ti_stop=%d\tmin_dist=%d\tmotif_length=%d\tthreshold=%.2f\n",i_stop,min_dist,motif_length,threshold);
  for (i=1;i<=i_stop;i++) {
    if ( (i-previous_position)>min_dist ) {
      current_score=0.0;
      //printf("%d\t",i);
      for (j=1;j<=motif_length;j++) {
	c=sequence[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=-7.6;     
	}
	//printf("%d",c);
      }
      //printf("\t%.2f\n",current_score);
      if (current_score>=threshold) {
	previous_position=i;
	n++;
	pos[n]=seq_length-i;                       // measure all positions with respect to 3' end
	strand[n]=curr_strand;
	//printf("HIT %d\t%d\t%d\n",n,pos[n],strand[n]);
	
	if (seq_nr>0) {
	  for (j=1;j<=motif_length;j++) 
	    numword[j]=sequence[i+j-1];
	  sprintf(word,"");
	  num2seq(numword,motif_length,word,1);
	  fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	}      // close check on seq_nr>0
      }        // close check on score
    }          // close check on distance to previous occurrence
  }            // close i loop

  curr_strand=0;
  previous_position=-motif_length;     // to process from the very beginning
  if (rc==1) {
    revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    for (i=1;i<=i_stop;i++) {
      if ( (i-previous_position)>min_dist ) {      
	current_score=0.0;
	for (j=1;j<=motif_length;j++) {
	  c=revcomp_seq[i+j-1];
	  if ( (c>0) & (c<5) ) {
	    p=weight_matrix[c][j];
	    current_score+=p;
	  } else {
	    current_score+=-7.6;     
	  }
	}
	if (current_score>=threshold) {
	  previous_position=i;
	  n++;
	  pos[n]=i;                          // positions in the minus strand are indicated as negative
	  strand[n]=0;
	  
	  if (seq_nr>0) {
	    for (j=1;j<=motif_length;j++) 
	      numword[j]=revcomp_seq[i+j-1];
	    sprintf(word,"");
	    num2seq(numword,motif_length,word,1);
	    fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
	  }       // close check on seq_nr>0
	}         // close check on score
      }           // close check on distance to previous occurrence
    }             // close i loop
    free_cvector(revcomp_seq,1,seq_length);
  }               // close check on rc==1

  if (seq_nr>0) {
    free_cvector(numword,1,motif_length);
    free_cvector(word,0,motif_length-1);
  }

  return n;
}

int scan_motif3(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc,FILE *log_file,int seq_nr)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * new in version 2: 08-16-2003 return all positions as positive, return strand information (+ -> 1, - -> 0) in a separate vector
   * new in version 3: also print information to log file 
   */
  int i,j,k;
  int i_stop;
  char c;
  float p;
  int n;
  float current_score;
  char *revcomp_seq;
  int curr_strand;
  char *word;
  char *numword;

  i_stop=seq_length-motif_length;
  n=0;

  numword=cvector(1,motif_length);
  word=cvector(0,motif_length-1);

  curr_strand=1;
  for (i=1;i<=i_stop;i++) {
    current_score=0.0;
    for (j=1;j<=motif_length;j++) {
      c=sequence[i+j-1];
      if ( (c>0) & (c<5) ) {
	p=weight_matrix[c][j];
	current_score+=p;
      } else {
	// penalise the cases where it is "n" or something strange happens
	current_score+=-7.6;     
      }
    }
    //printf("%d %d %.2f\n",i,c,current_score);
    if (current_score>=threshold) {
      n++;
      pos[n]=seq_length-i;                       // measure all positions with respect to 3' end
      strand[n]=curr_strand;
 
      for (j=1;j<=motif_length;j++) 
	numword[j]=sequence[i+j-1];
      sprintf(word,"");
      num2seq(numword,motif_length,word,1);
      fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
    }
  }

  curr_strand=0;
  if (rc==1) {
    revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    //for (i=1;i<=20;i++)
    //printf("%d",revcomp_seq[i]);
    //printf("\n");
    for (i=1;i<=i_stop;i++) {
      current_score=0;
      for (j=1;j<=motif_length;j++) {
	c=revcomp_seq[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=-7.6;     
	}
      }
      if (current_score>=threshold) {
	n++;
	pos[n]=i;                          // positions in the minus strand are indicated as negative
	strand[n]=0;
	/*for (j=1;j<=motif_length;j++)
	  printf("%d",revcomp_seq[i+j-1]);
	  printf("\n");
	  printf("%d\n",i);
	  exit(1);
	*/

	for (j=1;j<=motif_length;j++) 
	  numword[j]=revcomp_seq[i+j-1];
	sprintf(word,"");
	num2seq(numword,motif_length,word,1);
	fprintf(log_file,"%d\t%d\t%.4f\t%s\t%d\n",seq_nr,pos[n],current_score,word,curr_strand);
      }
    }                                       // close i loop
    free_cvector(revcomp_seq,1,seq_length);
  }                                         // close check on rc==1

  free_cvector(numword,1,motif_length);
  free_cvector(word,0,motif_length-1);

  return n;
}

int scan_motif2(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,char *strand,int rc)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * new in version 2: 08-16-2003 return all positions as positive, return strand information (+ -> 1, - -> 0) in a separate vector
   */
  int i,j,k;
  int i_stop;
  char c;
  float p;
  int n;
  float current_score;
  char *revcomp_seq;

  i_stop=seq_length-motif_length;
  n=0;

  for (i=1;i<=i_stop;i++) {
    current_score=0.0;
    for (j=1;j<=motif_length;j++) {
      c=sequence[i+j-1];
      if ( (c>0) & (c<5) ) {
	p=weight_matrix[c][j];
	current_score+=p;
      } else {
	// penalise the cases where it is "n" or something strange happens
	current_score+=-7.6;     
      }
    }
    if (current_score>=threshold) {
      n++;
      pos[n]=seq_length-i;                       // measure all positions with respect to 3' end
      strand[n]=1;
    }
  }

  if (rc==1) {
    revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    for (i=1;i<=i_stop;i++) {
      current_score=0;
      for (j=1;j<=motif_length;j++) {
	c=revcomp_seq[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=-7.6;     
	}
      }
      if (current_score>=threshold) {
	n++;
	pos[n]=i;                          // positions in the minus strand are indicated as negative
	strand[n]=0;
      }
    }                                       // close i loop
    free_cvector(revcomp_seq,1,seq_length);
  }                                         // close check on rc==1

  return n;
}

int scan_motif(float **weight_matrix,char *sequence,int motif_length,int seq_length,float threshold,int *pos,int rc)
{
  /* given a sequence (in numeric format), scan for the occurrence of a motif given by a weight matrix 
   * report all hits above a certain threshold
   * 08-09-2003 moved this function to regul_search_methods_v7.c (from scan_motif_v7.c)
   * 07-13-2003 do not compute the log (this is done when computing the weight matrix now
   * 04-29-2003 strongly penalize non a,c,g,t
   * 04-29-2003 since we already took care of the null entries in the weight matrix above, do not check again here
   * 08-13-2003: do not use global variables: strand, log_file
   * 08-13-2003: do not print out sequences at least in this version
   * 08-13-2003: scan reverse complementary sequence also if rc is 1
   */
  int i,j,k;
  int i_stop;
  char c;
  float p;
  int n;
  float current_score;
  char *revcomp_seq;

  i_stop=seq_length-motif_length;
  n=0;

  /* printf("scan_motif\n");
     printf("seq_length=%d\n",seq_length);
     printf("motif_length=%d\n",motif_length);
     printf("threshold=%.2f\n",threshold);
  */

  for (i=1;i<=i_stop;i++) {
    current_score=0.0;
    for (j=1;j<=motif_length;j++) {
      c=sequence[i+j-1];
      if ( (c>0) & (c<5) ) {
	p=weight_matrix[c][j];
	current_score+=p;
      } else {
	// penalise the cases where it is "n" or something strange happens
	current_score+=-7.6;     
      }
    }
    if (current_score>=threshold) {
      n++;
      pos[n]=seq_length-i;                        // measure all positions with respect to 3' end
      /* for (j=1;j<=motif_length;j++) {
	 numword[j]=sequence[i+j-1];
	 } */
    }
  }
  /* printf("n=%d\n",n);
     printf("current_score=%.2f\n",current_score);
     for (j=1;j<=motif_length;j++) 
     printf("%d",sequence[i-motif_length+j]);
     printf("\n");
  */

  if (rc==1) {
    revcomp_seq=cvector(1,seq_length);
    revcomp(sequence,revcomp_seq,seq_length);
    for (i=1;i<=i_stop;i++) {
      current_score=0;
      for (j=1;j<=motif_length;j++) {
	c=revcomp_seq[i+j-1];
	if ( (c>0) & (c<5) ) {
	  p=weight_matrix[c][j];
	  current_score+=p;
	} else {
	  current_score+=-7.6;     
	}
      }
      if (current_score>=threshold) {
	n++;
	pos[n]=-i;                          // positions in the minus strand are indicated as negative
      }
    }                                       // close i loop
    free_cvector(revcomp_seq,1,seq_length);
  }                                         // close check on rc==1

  return n;
}


void load_weight_matrix(char *weight_matrix_filename,int motif_length,float **weight_matrix) {
  /* load_weight_matrix
   * load a weight matrix into the variable weight_matrix
   * assumes that weight_matrix is 4 x motif_length
   * assumes that there are no comments lines in the weight matrix file
   */

  FILE *weight_matrix_file;
  int i,j;
  float temp_float;

  weight_matrix_file=fopen(weight_matrix_filename,"r");
  check_file_handle(weight_matrix_file,weight_matrix_filename);

  for (i=1;i<=4;i++) {
    for (j=1;j<=motif_length;j++) {
      fscanf(weight_matrix_file,"%f\n",&temp_float);
      weight_matrix[i][j]=temp_float;
    }
  }
  fclose(weight_matrix_file);
}

float score_motif(float **weight_matrix,char *sequence,int motif_length)
{
  /* given a sequence (in numeric format), compute its score given the weight_matrix
   */

  int i,j;
  int i_stop;
  char c;
  float p;
  float score=0;
  float penalty=-7.6;
  
  for (j=1;j<=motif_length;j++) {
    c=sequence[j];
    //printf("c=%d\t",c);
    if ( (c>0) & (c<5) ) {
      p=weight_matrix[c][j];
      //printf("p=%.4f\n",p);
      score+=log(p);
    } else {
      // penalise the cases where it is "n" or something strange happens
      score+=penalty;     // this correspond to a p value of 5x10^-4
    }
  }

  return score;
}

void read_fasta_seq(FILE *file_handle,char *info,char *seq)
{
  char txtline[10000];        // one line of text
  int n_chars=10000;          // maximum number of chars per line of text
  char c[10];                 // first character in the line
  char linfo[10000];          // temporary information line 
  char lseq[1000000];          // cumulative sequence [NOTE: maximum length=1,000,000 nucleotides]
  //int n;
  char temp[10000];           // current line without the carriage return
  char temp2[1000000];          // temporary entry holding the cumulative sequence [NOTE: maximum length=1,000,000 nucleotides]
  //int i=0;
  
  sprintf(lseq,"");
  sprintf(linfo,"");
  sprintf(info,"");
  sprintf(c,"");
  sprintf(txtline,"");
  while ( (strcmp(c,">") != 0) & (feof(file_handle)==0) ){
    //i++;
    //printf("i=%d\n",i);
  
    fgets(txtline,n_chars,file_handle);
    remove_carriage_return(temp,txtline);
    StrnCpy(c,txtline,1);
    if (strcmp(c,">")==0) {
      strcpy(linfo,temp);
    } else {
      sprintf(temp2,"%s%s",lseq,temp);
      strcpy(lseq,temp2);
      sprintf(temp,"");
      sprintf(temp2,"");
    }
  }

  sprintf(info,"%s",linfo);
  strcpy(seq,lseq);

  /* printf("regul_search_methods,read_fasta_seq\n");
  printf("\tinfo=%s\n",info);
  printf("\tseq=%s\n",seq); */
}

void read_fasta_seq2(FILE *file_handle,char *info,char *seq)
{
  char txtline[10000];        // one line of text
  int n_chars=10000;          // maximum number of chars per line of text
  char c[10];                 // first character in the line
  char temp[10000];           // current line without the carriage return

  sprintf(info,"");
  sprintf(c,"");
  sprintf(txtline,"");
  sprintf(seq,"");
  
  while ( (strcmp(c,">") != 0) & (feof(file_handle)==0) ){
    fgets(txtline,n_chars,file_handle);
    remove_carriage_return(temp,txtline);
    StrnCpy(c,txtline,1);
    if (strcmp(c,">")==0) {
      strcpy(info,temp);      
    } else {
      strcat(seq,temp);
    }
  }

  /* printf("regul_search_methods,read_fasta_seq\n");
     printf("\tinfo=%s\n",info);
     printf("\tseq=%s\n",seq); */
}

void read_fasta_seq3(FILE *file_handle,char *info,char *seq)
{
  /* read a fasta sequence
   * new in v3: attempt to solve a bug in reading the last line of the last entry
   */

  char txtline[10000];        // one line of text
  int n_chars=10000;          // maximum number of chars per line of text
  char c[10];                 // first character in the line
  char temp[10000];           // current line without the carriage return

  sprintf(info,"");
  sprintf(c,"");
  sprintf(txtline,"");
  sprintf(seq,"");
  
  while ( (strcmp(c,">") != 0) & (feof(file_handle)==0) ){
    fgets(txtline,n_chars,file_handle);
    remove_carriage_return(temp,txtline);
    StrnCpy(c,txtline,1);
    if (strcmp(c,">")==0) {
      strcpy(info,temp);      
    } else {
      if (feof(file_handle)==0)
	strcat(seq,temp);
    }
  }
}

char *remove_carriage_return(char *str_noret,char *str) 
{
  int n;
  char *d=str_noret;

  n=strlen(str);
  //n=n-2;
  n=n-1;
  StrnCpy(d,str,n);

  return (str_noret);
}

char *StrnCpy(char *dest,const char *src,size_t n)
{
  char *d=dest;
  if (!dest) return(NULL);
  if (!src) {
    *dest=0;
    return(dest);
  }
  while (n-- && (*d++ = *src++) ) ;
  *d=0;
  return (dest);
}

//float *count_motif_occurrences(char *motif,char *input_filename_upstream,int n_transcripts,float score_threshold,int verbose,float *n)
float *count_motif_occurrences(char *motif,char *input_filename_upstream,int n_transcripts,float score_threshold,float *n)
{
	float curr_score;
	int i;
	int motif_length;
	char upstream_sequence[2000];
	FILE *input_file_upstream;
	char *search_string;
	char *current_string;
	int curr_n;
	
	//motif_length=strlen(motif);	
	for (i=1;i<=n_transcripts;i++)
		n[i]=-1;
		
	input_file_upstream=fopen(input_filename_upstream,"r");
	if (input_file_upstream==NULL)	{
		printf("ERROR! I could not open your input upstream file %s. Exiting now...\n",input_filename_upstream);
		exit(1);
	}
		
	if (score_threshold>0.99)	{
		/* search for perfect matches */
		for (i=1;i<=n_transcripts;i++)	{
			fscanf(input_file_upstream,"%s\n",&upstream_sequence);
			curr_n=count_perfect_matches(upstream_sequence,motif);
			n[i]=curr_n;
		}
	}
	else	{
		for (i=1;i<=n_transcripts;i++)	{
			fscanf(input_file_upstream,"%s\n",&upstream_sequence);
			curr_score=count_motif_occurrences_prob(motif,upstream_sequence,score_threshold);
			n[i]=n[i]+curr_score+1.0;
		}
	}

	fclose(input_file_upstream);
	return n;	
}

float compare_aligned_seqs(char* s1,char* s2)
{
	/* directly compare two aligned sequences given by s1 and s2 and return a siilarity score
   */

	int i=0;
  int seq_length;
  float score=0;

  seq_length=strlen(s1);

	while ( (i<seq_length) & (score>i-2) ) {
  	if (s1[i]==s2[i])
    	score++;
    i++;	
   }
   score=(float)score/(float)seq_length;
   return score;
}

void revcomp(char *numseq,char *revcomp,int n) {
  /* reverse complement a sequence */
  
  int i;
  int c;
  int j;

  for (i=1;i<=n;i++) {
    j=n-i+1;
    revcomp[i]=5-numseq[j];
  }
}

float count_motif_occurrences_prob(char *motif,char *upstream_sequence,float threshold)
{
	/* count the number of occurrences of motif in a given upstream sequence
   */

  int i,j,k;
	float n;
	int upstream_sequence_length;
	int motif_length;
  float score;  /* similarity measure between two strings */

	motif_length=strlen(motif);
  upstream_sequence_length=strlen(upstream_sequence);
	
  n=0;
  i=0;
  while (i<=(upstream_sequence_length-motif_length))  {
  	j=0;score=0;
		while ( (j<motif_length) & (score>j-2) ) {
  		if (motif[j]==upstream_sequence[i+j])
    		score++;
    	j++;	
   	}
   	score=(float)score/(float)motif_length;
    if (score>=threshold)  {
    	/* only consider the strongest matches after a threshold is imposed */
     	n+=score;
     	i+=motif_length-1;
    }
    i++;
	}
  return n;
}

void print_motif(char *motif,float *expression,char *input_filename_upstream,int n_around,FILE *output_file)
{
	/* search and print occurrences of a given motif */
	float curr_expression;
	int i;
	int motif_length;
	FILE *input_file_upstream;
	int n;
	char upstream_sequence[2000];
	int upstream_sequence_length;
	char *search_string;
	char current_string[2000];
	char toprint[100];
	char *temp_toprint;
	int toprint_length;
	int max_toprint;
	int total_occurrences=0;
	int min_n_per_transcript=10;
	int max_n_per_transcript=0;
	int n_transcripts_present=0;
	float mean_n_per_transcript;
	float min_expression=1000;
	float max_expression=-1000;
	float mean_expression_per_occurrence=0;
	float std_expression_per_occurrence=0;
	float mean_expression_per_transcript=0;
	float std_expression_per_transcript=0;
	int current_string_length;	
	int n_transcripts;

	motif_length=strlen(motif);
	printf("**** Occurrences of %s in %s ****\n",motif,input_filename_upstream);
	printf("transc\texpression\tupstream\n");
	
	input_file_upstream=fopen(input_filename_upstream,"r");
	if (input_file_upstream==NULL)
		printf("I could not open the file %s",input_file_upstream);

	//max_toprint=motif_length+2*n_around;

	i=0;			
	while (feof(input_file_upstream)!=1)	{
		i++;	/* one more transcript */		
		curr_expression=expression[i];
		n=0;	/* number of occurrences of motif in current transcript */
		fscanf(input_file_upstream,"%s\n",&upstream_sequence);				
		upstream_sequence_length=strlen(upstream_sequence);
		strcpy(current_string,upstream_sequence);
		search_string=strstr(current_string,"a");
		while (search_string!=NULL)	{
			search_string=strstr(current_string,motif);
			if (search_string!=NULL)	{
					n++;	/* one more occurrence of motif in current transcript */
					total_occurrences++;	/* one more occurrence of motif in upstream sequences */
					if (strcmp(current_string,upstream_sequence)==1)	{
						/* first occurrence for this transcript */
						n_transcripts_present++;
						mean_expression_per_transcript+=curr_expression;
						std_expression_per_transcript+=curr_expression*curr_expression;
					}
										
					max_toprint=motif_length+2*n_around;
					printf("max_toprint=%d\n",max_toprint);
					printf("search_string=%s\n",search_string);
					printf("n_around=%d\n",n_around);
  				temp_toprint=search_string-n_around;
  				printf("temp_toprint=%s\n",temp_toprint);
  				printf("n_around=%d\n",n_around);
  				toprint_length=strlen(temp_toprint);
  				printf("toprint_length=%d\n",toprint_length);
  				printf("upstream_sequence_length=%d\n",upstream_sequence_length);
  				if (toprint_length>upstream_sequence_length)
  					temp_toprint=temp_toprint+upstream_sequence_length-toprint_length+1;
  				if (toprint_length<max_toprint)
  					max_toprint=toprint_length-1;
  				printf("%s\n",temp_toprint);	
  				printf("max_toprint=%d\n",max_toprint);
  				strncpy(toprint,temp_toprint,max_toprint);	
  				toprint_length=strlen(toprint);
  				printf("toprint_length=%d\n",toprint_length);
  				printf("%d\t%.2f\t%s\n",i,expression[i],toprint);
		 			  				
		 			strcpy(current_string,search_string);
		 			//current_string_length=strlen(current_string);
		 			
		 			printf("\n\n");

		 			if (curr_expression<min_expression)
		 				min_expression=curr_expression;
		 			if (curr_expression>max_expression)
		 				max_expression=curr_expression;		 			
		 			mean_expression_per_occurrence+=curr_expression;
		 			std_expression_per_occurrence+=curr_expression*curr_expression;
				}
		}				
		if (n>max_n_per_transcript)
			max_n_per_transcript=n;
		if ((n>0) & (n<min_n_per_transcript))
			min_n_per_transcript=n;
	}

	n_transcripts=i;	
	mean_n_per_transcript=(float)total_occurrences/(float)n_transcripts_present;			
	mean_expression_per_transcript=mean_expression_per_transcript/(float)n_transcripts_present;
	mean_expression_per_occurrence=mean_expression_per_occurrence/(float)total_occurrences;
	std_expression_per_transcript=std_expression_per_transcript/(float)n_transcripts_present;
	std_expression_per_occurrence=std_expression_per_occurrence/(float)total_occurrences;
	std_expression_per_transcript=std_expression_per_transcript-mean_expression_per_transcript;
	std_expression_per_occurrence=std_expression_per_occurrence-mean_expression_per_occurrence;
	std_expression_per_transcript=sqrt(std_expression_per_transcript);
	std_expression_per_occurrence=sqrt(std_expression_per_occurrence);
			
	printf("***** Motif statistics: %s *****\n",motif);		
	printf("\ttotal number of transcripts searched:\t%s\n",n_transcripts);
	printf("\tnumber of transcripts the motif is present in:\t%d\n",n_transcripts_present);
	printf("\ttotal number of motif occurrences:\t%d\n",total_occurrences);
	printf("\trange of occurrences per transcript (among those where it is present):\t%d to %d\n",min_n_per_transcript,max_n_per_transcript);
	printf("\tmean number of occurrences per transcript:\t%.2f\n",mean_n_per_transcript);
	printf("\tmean expression per transcript:\t%.2f\n",mean_expression_per_transcript);
	printf("\ts.d. expression per transcript:\t%.2f\n",std_expression_per_transcript);
	printf("\tmean expression per occurrence:\t%.2f\n",mean_expression_per_occurrence);
	printf("\ts.d. expression per occurrence:\t%.2f\n",std_expression_per_occurrence);
	printf("\trange of expression values:\t%.2f to %.2f\n",min_expression,max_expression);
	
	if (output_file!=NULL)	{
		fprintf(output_file,"***** Motif statistics: %s *****\n",motif);		
		fprintf(output_file,"\ttotal number of transcripts searched:\t%s\n",n_transcripts);
		fprintf(output_file,"\tnumber of transcripts the motif is present in:\t%d\n",n_transcripts_present);
		fprintf(output_file,"\ttotal number of motif occurrences:\t%d\n",total_occurrences);
		fprintf(output_file,"\trange of occurrences per transcript (among those where it is present):\t%d to %d\n",min_n_per_transcript,max_n_per_transcript);
		fprintf(output_file,"\tmean number of occurrences per transcript:\t%.2f\n",mean_n_per_transcript);
		fprintf(output_file,"\tmean expression per transcript:\t%.2f\n",mean_expression_per_transcript);
		fprintf(output_file,"\ts.d. expression per transcript:\t%.2f\n",std_expression_per_transcript);
		fprintf(output_file,"\tmean expression per occurrence:\t%.2f\n",mean_expression_per_occurrence);
		fprintf(output_file,"\ts.d. expression per occurrence:\t%.2f\n",std_expression_per_occurrence);
		fprintf(output_file,"\trange of expression values:\t%.2f to %.2f\n",min_expression,max_expression);
	}
}

int count_perfect_matches(char *upstream_sequence,char *motif)
{
	/* count the number of perfect matches for motif within the upstream sequence */
	int i;
	char *search_string;
	char current_string[2000];
	int n=0;
	int motif_length;
			
	motif_length=strlen(motif);	
	
	strcpy(current_string,upstream_sequence);
	while (search_string!=NULL)	{				
		search_string=strstr(current_string,motif);
		if (search_string!=NULL)	{
			n++;			
			strcpy(current_string,search_string+motif_length);
		}
	}
	return n;
}


int count_perfect_matches_v2(char *upstream_sequence,char *motif,int *pos)
{
  /* count the number of perfect matches for motif within the upstream sequence 
   * return the total count as well as the corresponding positions
   */

  char *search_string;
  char current_string[10000];
  int n=0;
  int motif_length;
			
  motif_length=strlen(motif);	
  strcpy(current_string,upstream_sequence);
  search_string=strstr(current_string,motif);
  while (search_string!=NULL)	{				
    search_string=strstr(current_string,motif);
    if (search_string!=NULL)	{
      n++;			
      pos[n]=strlen(search_string);
      strcpy(current_string,search_string+motif_length);
    }
  }
  return n;
}

int count_perfect_matches_v3(char *upstream_sequence,char *motif,unsigned short *pos)
{
  /* count the number of perfect matches for motif within the upstream sequence 
   * return the total count as well as the corresponding positions
   */

  char *search_string;
  char current_string[10000];
  int n=0;
  int motif_length;
			
  motif_length=strlen(motif);	
  strcpy(current_string,upstream_sequence);
  search_string=strstr(current_string,motif);
  while (search_string!=NULL)	{				
    search_string=strstr(current_string,motif);
    if (search_string!=NULL)	{
      n++;			
      pos[n]=(unsigned short)strlen(search_string);
      strcpy(current_string,search_string+motif_length);
    }
  }
  return n;
}

int count_matches_v1(char *upstream_sequence,char *motif,int *pos,int mism_pos)
{
  /* count the number of matches for motif within the upstream sequence 
   * return the total count as well as the corresponding positions
   * new in v3: allow one mismatch for the motif
   */

  char *search_string;
  char current_string[10000];
  int n=0;
  int motif_length;
  char orig_motif[10];
  int i;
  char alphabet[5];

  alphabet[1]=AA;alphabet[2]=CC;alphabet[3]=GG;alphabet[4]=TT;
  sprintf(orig_motif,motif);
  motif_length=strlen(motif);
  if (mism_pos==0) {
    strcpy(current_string,upstream_sequence);
    search_string=strstr(current_string,motif);
    while (search_string!=NULL)	{				
      search_string=strstr(current_string,motif);
      if (search_string!=NULL)	{
	n++;			
	pos[n]=strlen(search_string);
	strcpy(current_string,search_string+motif_length);
      }
    }
    return n;
  } else {
    for (i=1;i<=4;i++) {
      sprintf(motif,orig_motif);
      motif[motif_length-mism_pos]=alphabet[i];
      strcpy(current_string,upstream_sequence);
      search_string=strstr(current_string,motif);
      while (search_string!=NULL) {	
	search_string=strstr(current_string,motif);
	if (search_string!=NULL)	{
	  n++;			
	  pos[n]=strlen(search_string);
	  strcpy(current_string,search_string+motif_length);
	}
      }
    }
    return n;
  }
}

int count_perfect_matches_v2c(char *upstream_sequence,char *motif,char *motifc,int *pos)
{
  /* count the number of perfect matches for motif within the upstream sequence 
   * return the total count as well as the corresponding positions
   * here also include the complement of the motif
   */

  char *search_string;
  char current_string[10000];
  int n=0;
  int motif_length;
			
  motif_length=strlen(motif);	
  strcpy(current_string,upstream_sequence);
  search_string=strstr(current_string,motif);
  while (search_string!=NULL)	{				
    search_string=strstr(current_string,motif);
    if (search_string!=NULL)	{
      n++;			
      pos[n]=strlen(search_string);
      strcpy(current_string,search_string+motif_length);
    }
  }

  strcpy(current_string,upstream_sequence);
  search_string=strstr(current_string,motifc);
  while (search_string!=NULL)	{				
    search_string=strstr(current_string,motifc);
    if (search_string!=NULL)	{
      n++;			
      pos[n]=strlen(search_string);
      strcpy(current_string,search_string+motif_length);
    }
  }

  return n;
}

void get_motif(int index,char *motif)
{
	/* get the corresponding motif for a given index */
	int i,j,k;
	int max_motif_length=10;
	int motif_length;
	int motif_index;      /* new index relative to corresponding length */
	int total_n_motifs;
	int curr_length;
	char motifs_filename[100];
	int *total_n;	/* total number of motifs for each length */
	FILE *motif_file;
	char curr_motif[100];

	/* memory allocation */
	total_n=ivector(1,max_motif_length);
	/* validation */
	if (index<=0)	{
		sprintf(motif,"");
	}
	
	/* first compute motif length */
	motif_length=1;	
	for (j=1;j<=max_motif_length;j++) {
		total_n_motifs=1;
		for (k=1;k<=j;k++)
		    total_n_motifs=total_n_motifs*4;
		total_n_motifs=4*(total_n_motifs-1);
		total_n_motifs=total_n_motifs/(float)3;
		total_n[j]=total_n_motifs;
		if (index>total_n_motifs)
		    motif_length++;
	}

	/* compute new indices (relative to each motif file) */
	motif_index=index;
	if (motif_length>1)
	    motif_index=index-total_n[motif_length-1];

	/* open motifs file */
	sprintf(motifs_filename,"%s/temp/motifs_length%d.txt",gen_dir,motif_length);
	motif_file=fopen(motifs_filename,"r");	
	if (!motif_file)
		printf("I could not open the file %s\n",motifs_filename);
	j=1;  /* goes over each motif in the file until the desired one is reached */
	while ( (j>=1) & (j<=motif_index) )  {
		fscanf(motif_file,"%s\n",&curr_motif);
		j++;
	}
	fclose(motif_file);
	
	for (i=0;i<motif_length;i++)	
		motif[i]=curr_motif[i];
	
	/* free memory */	
	free_ivector(total_n,1,max_motif_length);
}

void motif_stats(char *motif,float *expression,char *input_filename_upstream,FILE *output_file,int n_transcripts)
{
	/* search occurrences of motif in upstream sequence and print out some basic stats about the motif */
  int i;
  int index;
	int n_init,n_finit;
  int motif_length;
  int upstream_length;
  float curr_expression;
  int n;
  int total_occurrences=0;
  int max_occurrences_per_transcript=0;
  float mean_occurrences_per_transcript;
  int n_transcripts_present=0;
  float mean_expression=0;
  float min_expression=0;
  float max_expression=0;
  float expression_var=0;
  float expression_std;
  float overall_mean_expression;
  float ptranscripts;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;
  FILE *input_file_upstream;
  char upstream_sequence[2000];	

	/* print overall expression stats */
	oa_mean_expression=meanf(expression,n_transcripts);
	oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
	oa_min_expression=minf(expression,n_transcripts);
	oa_max_expression=maxf(expression,n_transcripts);
	
	printf("Overall expression stats\n");fprintf(output_file,"Overall expression stats\n");
	printf("\tmean expression=%.2f\n",oa_mean_expression);fprintf(output_file,"\tmean expression=%.2f\n",oa_mean_expression);
	printf("\ts.d. expression=%.2f\n",oa_std_expression);fprintf(output_file,"\ts.d. expression=%.2f\n",oa_std_expression);
	printf("\tmin expression=%.2f\n",oa_min_expression);fprintf(output_file,"\tmin expression=%.2f\n",oa_min_expression);
	printf("\tmax expression=%.2f\n",oa_max_expression);fprintf(output_file,"\tmax expression=%.2f\n",oa_max_expression);		
		
	//motif_length=strlen(motif);
  printf("\tOccurences of motif %s in %s\n",motif,input_filename_upstream);fprintf(output_file,"\tOccurences of motif %s in %s\n",motif,input_filename_upstream);
	input_file_upstream=fopen(input_filename_upstream,"r");
  i=1;
  while (feof(input_file_upstream)!=1)	{
		fscanf(input_file_upstream,"%s\n",&upstream_sequence);
		n=count_perfect_matches(upstream_sequence,motif);
		if (n>0)	{
			total_occurrences=total_occurrences+n;
			if (n>max_occurrences_per_transcript)
				max_occurrences_per_transcript=n;
			n_transcripts_present++;
			curr_expression=expression[i];
			mean_expression=mean_expression+curr_expression;
			expression_var=expression_var+curr_expression*curr_expression;
      if (curr_expression<min_expression)
				min_expression=curr_expression;
			if (curr_expression>max_expression)
      	max_expression=curr_expression;			
		}
		i++;
	}
	fclose(input_file_upstream);

  printf("***** Statistics for motif = %s *****\n",motif);
  printf("\tNumber of transcripts explored = %d\n",n_transcripts);
  ptranscripts=100.0*(float)n_transcripts_present/(float)n_transcripts;
  printf("\tNumber of transcripts motif is present in = %d (%.1f)\n",n_transcripts_present,ptranscripts);
  printf("\tTotal number of motif occurrences = %d\n",total_occurrences);
  mean_occurrences_per_transcript=(float)total_occurrences/(float)n_transcripts_present;
	printf("\t\tMean number of occurrences per transcript = %.2f\n",mean_occurrences_per_transcript);
  printf("\t\tMax number of occurrences per transcript = %d\n",max_occurrences_per_transcript);
 	printf("\n\tExpression\n");
  mean_expression=mean_expression/(float)n_transcripts_present;
  printf("\tMean expression level of transcripts with motif = %.2f\n",mean_expression);
	printf("\tMin expression level of transcripts with motif = %.2f\n",min_expression);
  printf("\tMax expression level of transcripts with motif = %.2f\n",max_expression);
  expression_var=expression_var-mean_expression*mean_expression;
  expression_var=expression_var/(float)n_transcripts_present;
  //expression_std=(float)sqrt(expression_var);
  printf("\tS.D. of expression level of transcripts with motif = %.2f\n",expression_var);

  fprintf(output_file,"***** Statistics for motif = %s *****\n",motif);
  fprintf(output_file,"\tNumber of transcripts explored = %d\n",n_transcripts);
  fprintf(output_file,"\tNumber of transcripts motif is present in = %d (%.1f)\n",n_transcripts_present,ptranscripts);
  fprintf(output_file,"\tTotal number of motif occurrences = %d\n",total_occurrences);
	fprintf(output_file,"\t\tMean number of occurrences per transcript = %.2f\n",mean_occurrences_per_transcript);
  fprintf(output_file,"\t\tMax number of occurrences per transcript = %d\n",max_occurrences_per_transcript);
 	fprintf(output_file,"\n\tExpression\n");
  fprintf(output_file,"\tMean expression level of transcripts with motif = %.2f\n",mean_expression);
	fprintf(output_file,"\tMin expression level of transcripts with motif = %.2f\n",min_expression);
  fprintf(output_file,"\tMax expression level of transcripts with motif = %.2f\n",max_expression);
  fprintf(output_file,"\tS.D. of expression level of transcripts with motif = %.2f\n",expression_var);
}

void compute_residuals(float *residuals,float top_f,float *top_motif_count,int n_transcripts)
{
	/* compute residuals of model
   * completely re-wrote this method on 11-27-2001
   * just take current residuals, f value for current top motif and
   * normalized counts for current top motif
   */
	int i,j;
	float *model;

	printf("\tComputing residuals...(n_transcripts=%d)\n",n_transcripts);
	model=vector(1,n_transcripts);
	for (i=1;i<=n_transcripts;i++)
		model[i]=0;

  /* compute model predictions */
	for (i=1;i<=n_transcripts;i++)
		model[i]=top_f*top_motif_count[i];

	/* compute residuals */
	for (i=1;i<=n_transcripts;i++)
		residuals[i]=residuals[i]-model[i];

	printf("Ready.\n");
		
	/* free memory */
	free_vector(model,1,n_transcripts);		
}  /* end of compute_residuals method */

void motif2_stats_withcounts(char *motif1,char* motif2,float *expression,int *motif_count1,int *motif_count2,FILE *output_file,int n_transcripts)
{
	/* use input number of occurrences of motif1, motif2 and co-occurrences in upstream sequence and print out some basic stats about the motifs */
  int i;
  int index;
  int motif_length1;
  int motif_length2;
  float curr_expression;
  int n1,n2,nb;
  int total_occurrences1=0;
  int total_occurrences2=0;
  int total_occurrencesb_product=0;
  int total_occurrencesb_sum=0;
  int max_occurrences_per_transcript1=0;
  int max_occurrences_per_transcript2=0;
  int max_occurrences_per_transcriptb=0;
  float mean_occurrences_per_transcript1=0;
  float mean_occurrences_per_transcript2=0;
  float mean_occurrences_per_transcriptb=0;
  int n_transcripts_present1=0;
  int n_transcripts_present2=0;
  int n_transcripts_presentb=0;
  float mean_expression1,min_expression1,max_expression1,expression_var1,expression_std1;
  float mean_expression2,min_expression2,max_expression2,expression_var2,expression_std2;
  float mean_expressionb,min_expressionb,max_expressionb,expression_varb,expression_stdb;
  float overall_mean_expression;
  float ptranscripts1,ptranscripts2,ptranscriptsb;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;

  mean_expression1=0;min_expression1=0;max_expression1=0;expression_var1=0;//expression_std1=0;
  mean_expression2=0;min_expression2=0;max_expression2=0;expression_var2=0;//expression_std2=0;
  mean_expressionb=0;min_expressionb=0;max_expressionb=0;expression_varb=0;//expression_stdb=0;

	/* print overall expression stats */
	oa_mean_expression=meanf(expression,n_transcripts);
	oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
	oa_min_expression=minf(expression,n_transcripts);
	oa_max_expression=maxf(expression,n_transcripts);
	
	printf("Overall expression stats\n");fprintf(output_file,"Overall expression stats\n");
	printf("\tmean expression=%.2f\n",oa_mean_expression);fprintf(output_file,"\tmean expression=%.2f\n",oa_mean_expression);
	printf("\ts.d. expression=%.2f\n",oa_std_expression);fprintf(output_file,"\ts.d. expression=%.2f\n",oa_std_expression);
	printf("\tmin expression=%.2f\n",oa_min_expression);fprintf(output_file,"\tmin expression=%.2f\n",oa_min_expression);
	printf("\tmax expression=%.2f\n",oa_max_expression);fprintf(output_file,"\tmax expression=%.2f\n",oa_max_expression);		
		
  printf("\tOccurences of motifs %s and %s\n",motif1,motif2);fprintf(output_file,"\tOccurences of motif %s and %s\n",motif1,motif2);
  for (i=1;i<=n_transcripts;i++)	{
		n1=motif_count1[i];
		n2=motif_count2[i];
		nb=n1*n2;
		curr_expression=expression[i];
		if (n1>0)	{
			total_occurrences1=total_occurrences1+n1;
			if (n1>max_occurrences_per_transcript1)
				max_occurrences_per_transcript1=n1;
			n_transcripts_present1++;
			mean_expression1=mean_expression1+curr_expression;
			expression_var1=expression_var1+curr_expression*curr_expression;
      if (curr_expression<min_expression1)
				min_expression1=curr_expression;
			if (curr_expression>max_expression1)
      	max_expression1=curr_expression;			
		}
		if (n2>0)	{
			total_occurrences2=total_occurrences2+n2;
			if (n2>max_occurrences_per_transcript2)
				max_occurrences_per_transcript2=n2;
			n_transcripts_present2++;
			mean_expression2=mean_expression2+curr_expression;
			expression_var2=expression_var2+curr_expression*curr_expression;
      if (curr_expression<min_expression2)
				min_expression2=curr_expression;
			if (curr_expression>max_expression2)
      	max_expression2=curr_expression;			
		}
		if (nb>0)	{
			total_occurrencesb_product=total_occurrencesb_product+n1*n2;
			total_occurrencesb_sum=total_occurrencesb_sum+n1+n2;
			if (nb>max_occurrences_per_transcriptb)
				max_occurrences_per_transcriptb=nb;
			n_transcripts_presentb++;
			mean_expressionb=mean_expressionb+curr_expression;
			expression_varb=expression_varb+curr_expression*curr_expression;
      if (curr_expression<min_expressionb)
				min_expressionb=curr_expression;
			if (curr_expression>max_expressionb)
      	max_expressionb=curr_expression;			
		}	
		i++;
	}

  printf("***** Statistics for motifs = %s and %s *****\n",motif1,motif2);
  printf("\tNumber of transcripts explored = %d\n",n_transcripts);
  ptranscripts1=100.0*(float)n_transcripts_present1/(float)n_transcripts;
  ptranscripts2=100.0*(float)n_transcripts_present2/(float)n_transcripts;
  ptranscriptsb=100.0*(float)n_transcripts_presentb/(float)n_transcripts;
  printf("\tNumber of transcripts motif is present in =\t%d (%.1f)\t%d (%.1f)\t%d (%.1f)\n",n_transcripts_present1,ptranscripts1,n_transcripts_present2,ptranscripts2,n_transcripts_presentb,ptranscriptsb);
  printf("\tTotal number of motif occurrences         =\t%d\t%d\t%d\t(sum=%d)\n",total_occurrences1,total_occurrences2,total_occurrencesb_product,total_occurrencesb_sum);
  mean_occurrences_per_transcript1=(float)total_occurrences1/(float)n_transcripts_present1;
  mean_occurrences_per_transcript2=(float)total_occurrences2/(float)n_transcripts_present2;
  mean_occurrences_per_transcriptb=(float)total_occurrencesb_product/(float)n_transcripts_presentb;
	printf("\tMean number of occurrences per transcript =\t%.2f\t%.2f\t%.2f\n",mean_occurrences_per_transcript1,mean_occurrences_per_transcript2,mean_occurrences_per_transcriptb);
  printf("\tMax number of occurrences per transcript  =\t%d\t%d\t%d\n",max_occurrences_per_transcript1,max_occurrences_per_transcript2,max_occurrences_per_transcriptb);
 	printf("\n\tExpression\n");
  mean_expression1=mean_expression1/(float)n_transcripts_present1;
  mean_expression2=mean_expression2/(float)n_transcripts_present2;
  mean_expressionb=mean_expressionb/(float)n_transcripts_presentb;
  printf("\tMean expr level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",mean_expression1,mean_expression2,mean_expressionb);
	printf("\tMin expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",min_expression1,min_expression2,min_expressionb);
  printf("\tMax expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",max_expression1,max_expression2,max_expressionb);
  expression_var1=expression_var1-mean_expression1*mean_expression1;expression_var1=expression_var1/(float)n_transcripts_present1;//expression_std1=(float)sqrt(expression_var1);
  expression_var2=expression_var2-mean_expression2*mean_expression2;expression_var2=expression_var2/(float)n_transcripts_present2;//expression_std2=(float)sqrt(expression_var2);
  expression_varb=expression_varb-mean_expressionb*mean_expressionb;expression_varb=expression_varb/(float)n_transcripts_presentb;//expression_stdb=(float)sqrt(expression_varb);
  printf("\tS.D. of expr level of transcripts w motif =\t%.2f\t%.2f\t%.2f\n",expression_var1,expression_var2,expression_varb);

  fprintf(output_file,"***** Statistics for motifs = %s and %s *****\n",motif1,motif2);
  fprintf(output_file,"\tNumber of transcripts explored = %d\n",n_transcripts);
  fprintf(output_file,"\tNumber of transcripts motif is present in =\t%d (%.1f)\t%d (%.1f)\t%d (%.1f)\n",n_transcripts_present1,ptranscripts1,n_transcripts_present2,ptranscripts2,n_transcripts_presentb,ptranscriptsb);
  fprintf(output_file,"\tTotal number of motif occurrences         =\t%d\t%d\t%d\t(sum=%d)\n",total_occurrences1,total_occurrences2,total_occurrencesb_product,total_occurrencesb_sum);
	fprintf(output_file,"\tMean number of occurrences per transcript =\t%.2f\t%.2f\t%.2f\n",mean_occurrences_per_transcript1,mean_occurrences_per_transcript2,mean_occurrences_per_transcriptb);
  fprintf(output_file,"\tMax number of occurrences per transcript  =\t%d\t%d\t%d\n",max_occurrences_per_transcript1,max_occurrences_per_transcript2,max_occurrences_per_transcriptb);
 	fprintf(output_file,"\n\tExpression\n");
  fprintf(output_file,"\tMean expr level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",mean_expression1,mean_expression2,mean_expressionb);
	fprintf(output_file,"\tMin expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",min_expression1,min_expression2,min_expressionb);
  fprintf(output_file,"\tMax expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",max_expression1,max_expression2,max_expressionb);
  fprintf(output_file,"\tS.D. of expr level of transcripts w motif =\t%.2f\t%.2f\t%.2f\n",expression_var1,expression_var2,expression_varb);
}

void motif_stats_withcounts(char *motif,float *expression,int *motif_count,FILE *output_file,int n_transcripts)
{
	/* use input number of occurrences of motif in upstream sequence and print out some basic stats about the motif */
  int i;
  int index;
  float curr_expression;
  int n;
  int total_occurrences=0;
  int max_occurrences_per_transcript=0;
  float mean_occurrences_per_transcript;
  int n_transcripts_present=0;
  float mean_expression=0;
  float min_expression=0;
  float max_expression=0;
  float expression_var=0;
  float expression_std;
  float overall_mean_expression;
  float ptranscripts;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;

	/* print overall expression stats */
	oa_mean_expression=meanf(expression,n_transcripts);
	oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
	oa_min_expression=minf(expression,n_transcripts);
	oa_max_expression=maxf(expression,n_transcripts);
	
	printf("Overall expression stats\n");fprintf(output_file,"Overall expression stats\n");
	printf("\tmean expression=%.2f\n",oa_mean_expression);fprintf(output_file,"\tmean expression=%.2f\n",oa_mean_expression);
	printf("\ts.d. expression=%.2f\n",oa_std_expression);fprintf(output_file,"\ts.d. expression=%.2f\n",oa_std_expression);
	printf("\tmin expression=%.2f\n",oa_min_expression);fprintf(output_file,"\tmin expression=%.2f\n",oa_min_expression);
	printf("\tmax expression=%.2f\n",oa_max_expression);fprintf(output_file,"\tmax expression=%.2f\n",oa_max_expression);		
		
  printf("\tOccurences of motif %s\n",motif);fprintf(output_file,"\tOccurences of motif %s\n",motif);
  for (i=1;i<=n_transcripts;i++)	{
		n=motif_count[i];
		if (n>0)	{
			total_occurrences=total_occurrences+n;
			if (n>max_occurrences_per_transcript)
				max_occurrences_per_transcript=n;
			n_transcripts_present++;
			curr_expression=expression[i];
			mean_expression=mean_expression+curr_expression;
			expression_var=expression_var+curr_expression*curr_expression;
      if (curr_expression<min_expression)
				min_expression=curr_expression;
			if (curr_expression>max_expression)
      	max_expression=curr_expression;			
		}
	}

  printf("***** Statistics for motif = %s *****\n",motif);
  printf("\tNumber of transcripts explored = %d\n",n_transcripts);
  ptranscripts=100.0*(float)n_transcripts_present/(float)n_transcripts;
  printf("\tNumber of transcripts motif is present in = %d (%.1f)\n",n_transcripts_present,ptranscripts);
  printf("\tTotal number of motif occurrences = %d\n",total_occurrences);
  mean_occurrences_per_transcript=(float)total_occurrences/(float)n_transcripts_present;
	printf("\t\tMean number of occurrences per transcript = %.2f\n",mean_occurrences_per_transcript);
  printf("\t\tMax number of occurrences per transcript = %d\n",max_occurrences_per_transcript);
 	printf("\n\tExpression\n");
  mean_expression=mean_expression/(float)n_transcripts_present;
  printf("\tMean expression level of transcripts with motif = %.2f\n",mean_expression);
	printf("\tMin expression level of transcripts with motif = %.2f\n",min_expression);
  printf("\tMax expression level of transcripts with motif = %.2f\n",max_expression);
  expression_var=expression_var-mean_expression*mean_expression;
  expression_var=expression_var/(float)n_transcripts_present;
  //expression_std=(float)sqrt(expression_var);
  printf("\tS.D. of expression level of transcripts with motif = %.2f\n",expression_var);

  fprintf(output_file,"***** Statistics for motif = %s *****\n",motif);
  fprintf(output_file,"\tNumber of transcripts explored = %d\n",n_transcripts);
  fprintf(output_file,"\tNumber of transcripts motif is present in = %d (%.1f)\n",n_transcripts_present,ptranscripts);
  fprintf(output_file,"\tTotal number of motif occurrences = %d\n",total_occurrences);
	fprintf(output_file,"\t\tMean number of occurrences per transcript = %.2f\n",mean_occurrences_per_transcript);
  fprintf(output_file,"\t\tMax number of occurrences per transcript = %d\n",max_occurrences_per_transcript);
 	fprintf(output_file,"\n\tExpression\n");
  fprintf(output_file,"\tMean expression level of transcripts with motif = %.2f\n",mean_expression);
	fprintf(output_file,"\tMin expression level of transcripts with motif = %.2f\n",min_expression);
  fprintf(output_file,"\tMax expression level of transcripts with motif = %.2f\n",max_expression);
  fprintf(output_file,"\tS.D. of expression level of transcripts with motif = %.2f\n",expression_var);
}

void motif2_stats(char *motif1,char* motif2,float *expression,char *input_filename_upstream,FILE *output_file,int n_transcripts)
{
	/* use input number of occurrences of motif1, motif2 and co-occurrences in upstream sequence and print out some basic stats about the motifs */
  int i;
  int index;
	int n_init,n_finit;
  float curr_expression;
  int n1,n2,nb;
  int motif_length1,motif_length2;
  int total_occurrences1=0;
  int total_occurrences2=0;
  int total_occurrencesb_product=0;
  int total_occurrencesb_sum=0;
  int max_occurrences_per_transcript1=0;
  int max_occurrences_per_transcript2=0;
  int max_occurrences_per_transcriptb=0;
  float mean_occurrences_per_transcript1=0;
  float mean_occurrences_per_transcript2=0;
  float mean_occurrences_per_transcriptb=0;
  int n_transcripts_present1=0;
  int n_transcripts_present2=0;
  int n_transcripts_presentb=0;
  float mean_expression1,min_expression1,max_expression1,expression_var1,expression_std1;
  float mean_expression2,min_expression2,max_expression2,expression_var2,expression_std2;
  float mean_expressionb,min_expressionb,max_expressionb,expression_varb,expression_stdb;
  float overall_mean_expression;
  float ptranscripts1,ptranscripts2,ptranscriptsb;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;
  FILE *input_file_upstream;
  char upstream_sequence[2000];	

  mean_expression1=0;min_expression1=0;max_expression1=0;expression_var1=0;//expression_std1=0;
  mean_expression2=0;min_expression2=0;max_expression2=0;expression_var2=0;//expression_std2=0;
  mean_expressionb=0;min_expressionb=0;max_expressionb=0;expression_varb=0;//expression_stdb=0;

	/* print overall expression stats */
	oa_mean_expression=meanf(expression,n_transcripts);
	oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
	oa_min_expression=minf(expression,n_transcripts);
	oa_max_expression=maxf(expression,n_transcripts);
	
	printf("Overall expression stats\n");fprintf(output_file,"Overall expression stats\n");
	printf("\tmean expression=%.2f\n",oa_mean_expression);fprintf(output_file,"\tmean expression=%.2f\n",oa_mean_expression);
	printf("\ts.d. expression=%.2f\n",oa_std_expression);fprintf(output_file,"\ts.d. expression=%.2f\n",oa_std_expression);
	printf("\tmin expression=%.2f\n",oa_min_expression);fprintf(output_file,"\tmin expression=%.2f\n",oa_min_expression);
	printf("\tmax expression=%.2f\n",oa_max_expression);fprintf(output_file,"\tmax expression=%.2f\n",oa_max_expression);		
		
	//motif_length1=strlen(motif1);
	//motif_length2=strlen(motif2);
  printf("\tOccurences of motifs %s and %s in %s\n",motif1,motif2,input_filename_upstream);fprintf(output_file,"\tOccurences of motif %s and %s in %s\n",motif1,motif2,input_filename_upstream);
	input_file_upstream=fopen(input_filename_upstream,"r");
  i=1;
  while (feof(input_file_upstream)!=1)	{
		fscanf(input_file_upstream,"%s\n",&upstream_sequence);
		n1=count_perfect_matches(upstream_sequence,motif1);
		n2=count_perfect_matches(upstream_sequence,motif2);
		nb=n1*n2;
		curr_expression=expression[i];
		if (n1>0)	{
			total_occurrences1=total_occurrences1+n1;
			if (n1>max_occurrences_per_transcript1)
				max_occurrences_per_transcript1=n1;
			n_transcripts_present1++;
			mean_expression1=mean_expression1+curr_expression;
			expression_var1=expression_var1+curr_expression*curr_expression;
      if (curr_expression<min_expression1)
				min_expression1=curr_expression;
			if (curr_expression>max_expression1)
      	max_expression1=curr_expression;			
		}
		if (n2>0)	{
			total_occurrences2=total_occurrences2+n2;
			if (n2>max_occurrences_per_transcript2)
				max_occurrences_per_transcript2=n2;
			n_transcripts_present2++;
			mean_expression2=mean_expression2+curr_expression;
			expression_var2=expression_var2+curr_expression*curr_expression;
      if (curr_expression<min_expression2)
				min_expression2=curr_expression;
			if (curr_expression>max_expression2)
      	max_expression2=curr_expression;			
		}
		if (nb>0)	{
			total_occurrencesb_product=total_occurrencesb_product+n1*n2;
			total_occurrencesb_sum=total_occurrencesb_sum+n1+n2;
			if (nb>max_occurrences_per_transcriptb)
				max_occurrences_per_transcriptb=nb;
			n_transcripts_presentb++;
			mean_expressionb=mean_expressionb+curr_expression;
			expression_varb=expression_varb+curr_expression*curr_expression;
      if (curr_expression<min_expressionb)
				min_expressionb=curr_expression;
			if (curr_expression>max_expressionb)
      	max_expressionb=curr_expression;			
		}	
		i++;
	}
	fclose(input_file_upstream);

  printf("***** Statistics for motifs = %s and %s *****\n",motif1,motif2);
  printf("\tNumber of transcripts explored = %d\n",n_transcripts);
  ptranscripts1=100.0*(float)n_transcripts_present1/(float)n_transcripts;
  ptranscripts2=100.0*(float)n_transcripts_present2/(float)n_transcripts;
  ptranscriptsb=100.0*(float)n_transcripts_presentb/(float)n_transcripts;
  printf("\tNumber of transcripts motif is present in =\t%d (%.1f)\t%d (%.1f)\t%d (%.1f)\n",n_transcripts_present1,ptranscripts1,n_transcripts_present2,ptranscripts2,n_transcripts_presentb,ptranscriptsb);
  printf("\tTotal number of motif occurrences         =\t%d\t%d\t%d\t(sum=%d)\n",total_occurrences1,total_occurrences2,total_occurrencesb_product,total_occurrencesb_sum);
  mean_occurrences_per_transcript1=(float)total_occurrences1/(float)n_transcripts_present1;
  mean_occurrences_per_transcript2=(float)total_occurrences2/(float)n_transcripts_present2;
  mean_occurrences_per_transcriptb=(float)total_occurrencesb_product/(float)n_transcripts_presentb;
	printf("\tMean number of occurrences per transcript =\t%.2f\t%.2f\t%.2f\n",mean_occurrences_per_transcript1,mean_occurrences_per_transcript2,mean_occurrences_per_transcriptb);
  printf("\tMax number of occurrences per transcript  =\t%d\t%d\t%d\n",max_occurrences_per_transcript1,max_occurrences_per_transcript2,max_occurrences_per_transcriptb);
 	printf("\n\tExpression\n");
  mean_expression1=mean_expression1/(float)n_transcripts_present1;
  mean_expression2=mean_expression2/(float)n_transcripts_present2;
  mean_expressionb=mean_expressionb/(float)n_transcripts_presentb;
  printf("\tMean expr level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",mean_expression1,mean_expression2,mean_expressionb);
	printf("\tMin expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",min_expression1,min_expression2,min_expressionb);
  printf("\tMax expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",max_expression1,max_expression2,max_expressionb);
  expression_var1=expression_var1-mean_expression1*mean_expression1;expression_var1=expression_var1/(float)n_transcripts_present1;//expression_std1=(float)sqrt(expression_var1);
  expression_var2=expression_var2-mean_expression2*mean_expression2;expression_var2=expression_var2/(float)n_transcripts_present2;//expression_std2=(float)sqrt(expression_var2);
  expression_varb=expression_varb-mean_expressionb*mean_expressionb;expression_varb=expression_varb/(float)n_transcripts_presentb;//expression_stdb=(float)sqrt(expression_varb);
  printf("\tS.D. of expr level of transcripts w motif =\t%.2f\t%.2f\t%.2f\n",expression_var1,expression_var2,expression_varb);

  fprintf(output_file,"***** Statistics for motifs = %s and %s *****\n",motif1,motif2);
  fprintf(output_file,"\tNumber of transcripts explored = %d\n",n_transcripts);
  fprintf(output_file,"\tNumber of transcripts motif is present in =\t%d (%.1f)\t%d (%.1f)\t%d (%.1f)\n",n_transcripts_present1,ptranscripts1,n_transcripts_present2,ptranscripts2,n_transcripts_presentb,ptranscriptsb);
  fprintf(output_file,"\tTotal number of motif occurrences         =\t%d\t%d\t%d\t(sum=%d)\n",total_occurrences1,total_occurrences2,total_occurrencesb_product,total_occurrencesb_sum);
	fprintf(output_file,"\tMean number of occurrences per transcript =\t%.2f\t%.2f\t%.2f\n",mean_occurrences_per_transcript1,mean_occurrences_per_transcript2,mean_occurrences_per_transcriptb);
  fprintf(output_file,"\tMax number of occurrences per transcript  =\t%d\t%d\t%d\n",max_occurrences_per_transcript1,max_occurrences_per_transcript2,max_occurrences_per_transcriptb);
 	fprintf(output_file,"\n\tExpression\n");
  fprintf(output_file,"\tMean expr level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",mean_expression1,mean_expression2,mean_expressionb);
	fprintf(output_file,"\tMin expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",min_expression1,min_expression2,min_expressionb);
  fprintf(output_file,"\tMax expr. level of transcripts with motif =\t%.2f\t%.2f\t%.2f\n",max_expression1,max_expression2,max_expressionb);
  fprintf(output_file,"\tS.D. of expr level of transcripts w motif =\t%.2f\t%.2f\t%.2f\n",expression_var1,expression_var2,expression_varb);
}

int get_index_from_nummotif(char *motif)
{
  int motif_length;
  int index;
  int i;

  motif_length=strlen(motif);
  index=0;
  for (i=0;i<motif_length;i++) {
    index+=powers_of_four[i]*(int)motif[motif_length-i-1];
  }
  return index;
}

void get_nummotif_from_index(int index,char *motif)
{
  int temp_index;
  int rem;
  int i;
  int motif_length;

  motif_length=strlen(motif);
  i=1;
  temp_index=index;
  while (temp_index>0) {
    rem = temp_index % 4;
    if (rem==0) {
      temp_index-=4;
      motif[motif_length-i]=4;
    }
    else {
      temp_index-=rem;
      motif[motif_length-i]=rem;
    }
    temp_index/=4;
    i++;
  }
}

void convert_nummotif2text(char *nummotif,char *textmotif)
{
  int motif_length;
  int i;

  motif_length=strlen(nummotif);
  for (i=0;i<motif_length;i++) {
    switch (nummotif[i]) {
    case 1: textmotif[i]=AA;break;
    case 2: textmotif[i]=CC;break;
    case 3: textmotif[i]=GG;break;
    case 4: textmotif[i]=TT;break;
    }
  }
}

void convert_text2nummotif(char *motif,char *nummotif)
{
  int motif_length;
  int i;

  motif_length=strlen(motif);
  for (i=0;i<motif_length;i++) {
    switch (motif[i]) {
    case AA: nummotif[i]=1;break; // 97
    case CC: nummotif[i]=2;break; // 99
    case GG: nummotif[i]=3;break; // 103
    case TT: nummotif[i]=4;break; // 116
    }
  }
}

void mismatch(int index,int motif_length,int *mismatch_indices)
{
  char motif[10];
  int i,j,k,ll;
  char curr_char;
  char other_char[3];
  char motif_mismatch[10];
  int mismatch_index;

  switch (motif_length) {
  case 1: sprintf(motif,"n");break;
  case 2: sprintf(motif,"nn");break;
  case 3: sprintf(motif,"nnn");break;
  case 4: sprintf(motif,"nnnn");break;
  case 5: sprintf(motif,"nnnnn");break;
  case 6: sprintf(motif,"nnnnnn");break;
  case 7: sprintf(motif,"nnnnnnn");break;
  case 8: sprintf(motif,"nnnnnnnn");break;
  case 9: sprintf(motif,"nnnnnnnnn");break;
  case 10: sprintf(motif,"nnnnnnnnnn");break;
  }

  get_nummotif_from_index(index,motif);

  /* total_n_mismatches=3*motif_length; */
  k=1;
  for (i=0;i<motif_length;i++) {
    curr_char=motif[i];
    switch (curr_char) {
    case 1: other_char[0]=2;other_char[1]=3;other_char[2]=4;break;
    case 2: other_char[0]=1;other_char[1]=3;other_char[2]=4;break;
    case 3: other_char[0]=1;other_char[1]=2;other_char[2]=4;break;
    case 4: other_char[0]=1;other_char[1]=2;other_char[2]=3;break;
    }
    for (j=0;j<3;j++) {
      strcpy(motif_mismatch,motif);
      motif_mismatch[i]=other_char[j];
      mismatch_indices[k]=get_index_from_nummotif(motif_mismatch);
      k++;
    }
  }
}

int total_motifs(int length)
{
	/* returns the total number of motifs of a given length */
	int i;
	int n_motifs;

	n_motifs=1;
	for (i=1;i<=length;i++)
		n_motifs=n_motifs*4;
	
	return n_motifs;	
}

float bootstrap_stats(float *bootstrap,int n,float actual_value)
{
	/* give the proportion of bootstrapped entries that are larger than the actual value */
	int n_larger=0;
	float p_larger;
	int i;
	
	for (i=1;i<=n;i++) {
		if (bootstrap[i]>actual_value)
			n_larger++;
	}
	p_larger=(float)n_larger/(float)n;
	return p_larger;
}

void motif_stats_bin(int motif_index,float *expression,char *filename_motifcount_match,char *filename_motifcount_mismatch,FILE *output_file,int n_transcripts)
{
	/* display statistics using the binary files as input to quickly retrieve the
	/* use input number of occurrences of motif in upstream sequence and print out some basic stats about the motif */
  int i;
  int index;
  float curr_expression;
  char n_match;
  char n_mismatch;
  int total_occurrences_match=0;
  int total_occurrences_mismatch=0;
  int max_occurrences_per_transcript_match=0;
  int max_occurrences_per_transcript_mismatch=0;
  float mean_occurrences_per_transcript_match;
  float mean_occurrences_per_transcript_mismatch;
  int n_transcripts_present_match=0;
  int n_transcripts_present_mismatch=0;
  float mean_expression_match=0;
  float min_expression_match=0;
  float max_expression_match=0;
  float expression_var_match=0;
  float expression_std_match;
  float mean_expression_mismatch=0;
  float min_expression_mismatch=0;
  float max_expression_mismatch=0;
  float expression_var_mismatch=0;
  float expression_std_mismatch;
  float overall_mean_expression;
  float ptranscripts_match;
  float ptranscripts_mismatch;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;
  FILE *file_motifcount_match;
  FILE *file_motifcount_mismatch;

	/* print overall expression stats */
	oa_mean_expression=meanf(expression,n_transcripts);
	oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
	oa_min_expression=minf(expression,n_transcripts);
	oa_max_expression=maxf(expression,n_transcripts);
	
	printf("Overall expression stats\n");
	printf("\tmean expression=%.2f\n",oa_mean_expression);fprintf(output_file,"\t%.2f",oa_mean_expression);
	printf("\ts.d. expression=%.2f\n",oa_std_expression);fprintf(output_file,"\t%.2f",oa_std_expression);
	printf("\tmin expression=%.2f\n",oa_min_expression);fprintf(output_file,"\t%.2f",oa_min_expression);
	printf("\tmax expression=%.2f\n",oa_max_expression);fprintf(output_file,"\t%.2f",oa_max_expression);		
	
	file_motifcount_match=fopen(filename_motifcount_match,"r");
	if (!file_motifcount_match) {
		printf("ERROR!!! I could not open the file %s for reading\n",filename_motifcount_match);
		exit(1);
	}
	file_motifcount_mismatch=fopen(filename_motifcount_mismatch,"r");
	if (!file_motifcount_mismatch) {
		printf("ERROR!!! I could not open the file %s for reading\n",filename_motifcount_mismatch);
		exit(1);
	}		
  for (i=1;i<=n_transcripts;i++)	{
		index=n_transcripts*(motif_index-1)+i-1; /* index of current entry in motif count files */
  	fseek(file_motifcount_match,sizeof(char)*index,SEEK_SET);
  	fread(&n_match,sizeof(char),1,file_motifcount_match);										
  	fseek(file_motifcount_mismatch,sizeof(char)*index,SEEK_SET);
  	fread(&n_mismatch,sizeof(char),1,file_motifcount_mismatch);												
		if (n_match>0)	{
			total_occurrences_match=total_occurrences_match+n_match;
			if (n_match>max_occurrences_per_transcript_match)
				max_occurrences_per_transcript_match=n_match;
			n_transcripts_present_match++;
			curr_expression=expression[i];
			mean_expression_match=mean_expression_match+curr_expression;
			expression_var_match=expression_var_match+curr_expression*curr_expression;
      if (curr_expression<min_expression_match)
				min_expression_match=curr_expression;
			if (curr_expression>max_expression_match)
      	max_expression_match=curr_expression;			
		}
		if (n_mismatch>0)	{
			total_occurrences_mismatch=total_occurrences_mismatch+n_mismatch;
			if (n_mismatch>max_occurrences_per_transcript_mismatch)
				max_occurrences_per_transcript_mismatch=n_mismatch;
			n_transcripts_present_mismatch++;
			curr_expression=expression[i];
			mean_expression_mismatch=mean_expression_mismatch+curr_expression;
			expression_var_mismatch=expression_var_mismatch+curr_expression*curr_expression;
      if (curr_expression<min_expression_mismatch)
				min_expression_mismatch=curr_expression;
			if (curr_expression>max_expression_mismatch)
      	max_expression_mismatch=curr_expression;			
		}		
		i++;
	}

  printf("\tNumber of transcripts explored = %d\n",n_transcripts);
  ptranscripts_match=100.0*(float)n_transcripts_present_match/(float)n_transcripts;
  ptranscripts_mismatch=100.0*(float)n_transcripts_present_mismatch/(float)n_transcripts;
  printf("\tNumber of transcripts motif is present in (match)= %d (%.1f)\n",n_transcripts_present_match,ptranscripts_match);
  printf("\tNumber of transcripts motif is present in (mismatch)= %d (%.1f)\n",n_transcripts_present_mismatch,ptranscripts_mismatch);
  printf("\tTotal number of motif occurrences (match)= %d\n",total_occurrences_match);
  printf("\tTotal number of motif occurrences (mismatch)= %d\n",total_occurrences_mismatch);
  mean_occurrences_per_transcript_match=(float)total_occurrences_match/(float)n_transcripts_present_match;
  mean_occurrences_per_transcript_mismatch=(float)total_occurrences_mismatch/(float)n_transcripts_present_mismatch;
	printf("\t\tMean number of occurrences per transcript (match) = %.2f\n",mean_occurrences_per_transcript_match);
	printf("\t\tMean number of occurrences per transcript (mismatch) = %.2f\n",mean_occurrences_per_transcript_mismatch);
  printf("\t\tMax number of occurrences per transcript (match)= %d\n",max_occurrences_per_transcript_match);
  printf("\t\tMax number of occurrences per transcript (mismatch)= %d\n",max_occurrences_per_transcript_mismatch);
 	printf("\n\tExpression\n");
  mean_expression_match=mean_expression_match/(float)n_transcripts_present_match;
  mean_expression_mismatch=mean_expression_mismatch/(float)n_transcripts_present_mismatch;
  printf("\tMean expression level of transcripts with motif (match)= %.2f\n",mean_expression_match);
  printf("\tMean expression level of transcripts with motif (mismatch)= %.2f\n",mean_expression_mismatch);
	printf("\tMin expression level of transcripts with motif (match) = %.2f\n",min_expression_match);
	printf("\tMin expression level of transcripts with motif (mismatch) = %.2f\n",min_expression_mismatch);
  printf("\tMax expression level of transcripts with motif (match)= %.2f\n",max_expression_match);
  printf("\tMax expression level of transcripts with motif (mismatch)= %.2f\n",max_expression_mismatch);
  expression_var_match=expression_var_match/(float)n_transcripts_present_match;
  expression_var_mismatch=expression_var_mismatch/(float)n_transcripts_present_mismatch;
  expression_var_match=expression_var_match-mean_expression_match*mean_expression_match;
  expression_var_mismatch=expression_var_mismatch-mean_expression_mismatch*mean_expression_mismatch;
  //expression_std_match=(float)sqrt(expression_var_match);
  //expression_std_mismatch=(float)sqrt(expression_var_mismatch);
  printf("\tS.D. of expression level of transcripts with motif (match)= %.2f\n",expression_var_match);
  printf("\tS.D. of expression level of transcripts with motif (mismatch)= %.2f\n",expression_var_mismatch);

  fprintf(output_file,"\t%d\t%.1f",n_transcripts_present_match,ptranscripts_match);
  fprintf(output_file,"\t%d",total_occurrences_match);
	fprintf(output_file,"\t%.2f",mean_occurrences_per_transcript_match);
  fprintf(output_file,"\t%d",max_occurrences_per_transcript_match);
  fprintf(output_file,"\t%.2f",mean_expression_match);
	fprintf(output_file,"\t%.2f",min_expression_match);
  fprintf(output_file,"\t%.2f",max_expression_match);
  fprintf(output_file,"\t%.2f",expression_var_match);

  fprintf(output_file,"\t%d\t%.1f",n_transcripts_present_mismatch,ptranscripts_mismatch);
  fprintf(output_file,"\t%d",total_occurrences_mismatch);
	fprintf(output_file,"\t%.2f",mean_occurrences_per_transcript_mismatch);
  fprintf(output_file,"\t%d",max_occurrences_per_transcript_mismatch);
  fprintf(output_file,"\t%.2f",mean_expression_mismatch);
	fprintf(output_file,"\t%.2f",min_expression_mismatch);
  fprintf(output_file,"\t%.2f",max_expression_mismatch);
  fprintf(output_file,"\t%.2f",expression_var_mismatch);
}

//void motif_stats_bin_v2(int motif_index,int mism1,int mism2,float *expression,char *filename_motifcount_match,int *motif_count,FILE *output_file,int n_transcripts);
void motif_stats_bin_v2(int motif_index,float *expression,char *filename_motifcount_match,int *motif_count,FILE *output_file,int n_transcripts)
{
  /* display statistics using the binary files as input to quickly retrieve the motif count numbers 
   * new in version 2: only match file is used as input and we also take mism1 and mism2 as input and 
   *                   do not print info to screen
   */

  /* note: motif count refers to the motif count including the corresponding mismatch counts, etc */

  int i;
  int index;
  float curr_expression;
  char n_match;
  char n_mismatch;
  int total_occurrences_match=0;
  int total_occurrences_mismatch=0;
  int max_occurrences_per_transcript_match=0;
  int max_occurrences_per_transcript_mismatch=0;
  float mean_occurrences_per_transcript_match;
  float mean_occurrences_per_transcript_mismatch;
  int n_transcripts_present_match=0;
  int n_transcripts_present_mismatch=0;
  float mean_expression_match=0;
  float min_expression_match=0;
  float max_expression_match=0;
  float expression_var_match=0;
  float expression_std_match;
  float mean_expression_mismatch=0;
  float min_expression_mismatch=0;
  float max_expression_mismatch=0;
  float expression_var_mismatch=0;
  float expression_std_mismatch;
  float overall_mean_expression;
  float ptranscripts_match;
  float ptranscripts_mismatch;
  float oa_std_expression;
  float oa_mean_expression;
  float oa_min_expression;
  float oa_max_expression;
  FILE *file_motifcount_match;

  /* print overall expression stats */
  oa_mean_expression=meanf(expression,n_transcripts);
  oa_std_expression=stdf(expression,oa_mean_expression,n_transcripts);
  oa_min_expression=minf(expression,n_transcripts);
  oa_max_expression=maxf(expression,n_transcripts);
	
  fprintf(output_file,"\t%.2f",oa_mean_expression);
  fprintf(output_file,"\t%.2f",oa_std_expression);
  fprintf(output_file,"\t%.2f",oa_min_expression);
  fprintf(output_file,"\t%.2f",oa_max_expression);		
	
  file_motifcount_match=fopen(filename_motifcount_match,"r");
  if (!file_motifcount_match) {
    printf("ERROR!!! I could not open the file %s for reading\n",filename_motifcount_match);exit(1);
  }
  index=n_transcripts*(motif_index-1);
  fseek(file_motifcount_match,sizeof(char)*index,SEEK_SET);
  for (i=1;i<=n_transcripts;i++)	{
    fread(&n_match,sizeof(char),1,file_motifcount_match);
    n_mismatch=motif_count[i];

    if (n_match>0)	{
      total_occurrences_match=total_occurrences_match+n_match;
      if (n_match>max_occurrences_per_transcript_match)
	max_occurrences_per_transcript_match=n_match;
      n_transcripts_present_match++;
      curr_expression=expression[i];
      mean_expression_match=mean_expression_match+curr_expression;
      expression_var_match=expression_var_match+curr_expression*curr_expression;
      if (curr_expression<min_expression_match)
	min_expression_match=curr_expression;
      if (curr_expression>max_expression_match)
      	max_expression_match=curr_expression;			
    }
    if (n_mismatch>0)	{
      total_occurrences_mismatch=total_occurrences_mismatch+n_mismatch;
      if (n_mismatch>max_occurrences_per_transcript_mismatch)
	max_occurrences_per_transcript_mismatch=n_mismatch;
      n_transcripts_present_mismatch++;
      curr_expression=expression[i];
      mean_expression_mismatch=mean_expression_mismatch+curr_expression;
      expression_var_mismatch=expression_var_mismatch+curr_expression*curr_expression;
      if (curr_expression<min_expression_mismatch)
	min_expression_mismatch=curr_expression;
      if (curr_expression>max_expression_mismatch)
      	max_expression_mismatch=curr_expression;			
    }		
  }

  ptranscripts_match=100.0*(float)n_transcripts_present_match/(float)n_transcripts;
  ptranscripts_mismatch=100.0*(float)n_transcripts_present_mismatch/(float)n_transcripts;
  mean_occurrences_per_transcript_match=(float)total_occurrences_match/(float)n_transcripts_present_match;
  mean_occurrences_per_transcript_mismatch=(float)total_occurrences_mismatch/(float)n_transcripts_present_mismatch;
  mean_expression_match=mean_expression_match/(float)n_transcripts_present_match;
  mean_expression_mismatch=mean_expression_mismatch/(float)n_transcripts_present_mismatch;
  expression_var_match=expression_var_match/(float)n_transcripts_present_match;
  expression_var_mismatch=expression_var_mismatch/(float)n_transcripts_present_mismatch;
  expression_var_match=expression_var_match-mean_expression_match*mean_expression_match;
  expression_var_mismatch=expression_var_mismatch-mean_expression_mismatch*mean_expression_mismatch;
  //expression_std_match=(float)sqrt(expression_var_match);
  //expression_std_mismatch=(float)sqrt(expression_var_mismatch);

  fprintf(output_file,"\t%d\t%.1f",n_transcripts_present_match,ptranscripts_match);
  fprintf(output_file,"\t%d",total_occurrences_match);
  fprintf(output_file,"\t%.2f",mean_occurrences_per_transcript_match);
  fprintf(output_file,"\t%d",max_occurrences_per_transcript_match);
  fprintf(output_file,"\t%.2f",mean_expression_match);
  fprintf(output_file,"\t%.2f",min_expression_match);
  fprintf(output_file,"\t%.2f",max_expression_match);
  fprintf(output_file,"\t%.2f",expression_var_match);
  
  fprintf(output_file,"\t%d\t%.1f",n_transcripts_present_mismatch,ptranscripts_mismatch);
  fprintf(output_file,"\t%d",total_occurrences_mismatch);
  fprintf(output_file,"\t%.2f",mean_occurrences_per_transcript_mismatch);
  fprintf(output_file,"\t%d",max_occurrences_per_transcript_mismatch);
  fprintf(output_file,"\t%.2f",mean_expression_mismatch);
  fprintf(output_file,"\t%.2f",min_expression_mismatch);
  fprintf(output_file,"\t%.2f",max_expression_mismatch);
  fprintf(output_file,"\t%.2f",expression_var_mismatch);
}

//void quick_mismatch(int index,int motif_length,int *mismatch_indices)
void quick_mismatch(int index,int *mismatch_indices)
{
  /* given an index number, compute the indices of all possible mismatches */
  int temp_index;
  int rem;
  int i,j,k;
  int delta_index;
  char other_char[3];
  int temp;

  i=1;	/* position starting from the rightmost (3') position */
	k=0;	/* running number of mismatches */
	temp_index=index;
  while (temp_index>0) {
    rem = temp_index % 4;
    switch (rem) {
    case 1: other_char[0]=2;other_char[1]=3;other_char[2]=4;break;
    case 2: other_char[0]=1;other_char[1]=3;other_char[2]=4;break;
    case 3: other_char[0]=1;other_char[1]=2;other_char[2]=4;break;
    case 0: other_char[0]=1;other_char[1]=2;other_char[2]=3;rem=4;break;
    }
    temp_index-=rem;

    for (j=0;j<3;j++) {
    	k++;
      delta_index=(other_char[j]-rem)*powers_of_four[i-1];
			mismatch_indices[k]=delta_index+index;
    }
    temp_index/=4;
    i++;
  }
}

//void quick_mismatch_v2(int index,int motif_length,int **mismatch_indices)
void quick_mismatch_v2(int index,int **mismatch_indices)
{
  /* given an index number, compute the indices of all possible mismatches */
  /* new in v2
   * here the indices are returned in a matrix indicating the mismatches for each position in
   * separate rows
   */

  int temp_index;
  int rem;
  int i,j,k;
  int delta_index;
  char other_char[3];
  int temp;

  i=1;	/* position starting from the rightmost (3') position */
  temp_index=index;
  while (temp_index>0) {
    rem = temp_index % 4;
    switch (rem) {
    case 1: other_char[0]=2;other_char[1]=3;other_char[2]=4;break;
    case 2: other_char[0]=1;other_char[1]=3;other_char[2]=4;break;
    case 3: other_char[0]=1;other_char[1]=2;other_char[2]=4;break;
    case 0: other_char[0]=1;other_char[1]=2;other_char[2]=3;rem=4;break;
    }
    temp_index-=rem;

    for (j=0;j<3;j++) {
      delta_index=(other_char[j]-rem)*powers_of_four[i-1];
      mismatch_indices[i][j+1]=delta_index+index;
    }
    temp_index/=4;
    i++;
  }
}

//void quick_2mismatch(int index,int motif_length,int *mismatch_indices)
void quick_2mismatch(int index,int *mismatch_indices)
{
  /* given an index number, compute the indices of all possible motifs with 2 mismatches */
  int temp_index1,temp_index2;
  int rem1,rem2;
  int i,j,k,ll,mm;
  int delta_index;
  char other_char1[3];
  char other_char2[3];
  int temp;

  i=1;	/* position starting from the rightmost (3') position */
  k=0;	/* running number of mismatches */
  temp_index1=index;
  while (temp_index1>0) {
    rem1 = temp_index1 % 4;
    switch (rem1) {
    case 1: other_char1[0]=2;other_char1[1]=3;other_char1[2]=4;break;
    case 2: other_char1[0]=1;other_char1[1]=3;other_char1[2]=4;break;
    case 3: other_char1[0]=1;other_char1[1]=2;other_char1[2]=4;break;
    case 0: other_char1[0]=1;other_char1[1]=2;other_char1[2]=3;rem1=4;break;
    }
    temp_index1-=rem1;
    temp_index1/=4;

    for (j=0;j<3;j++) {
      ll=1;  		
      temp_index2=temp_index1;
      while (temp_index2>0) {
	rem2 = temp_index2 % 4;
	switch (rem2) {
	case 1: other_char2[0]=2;other_char2[1]=3;other_char2[2]=4;break;
	case 2: other_char2[0]=1;other_char2[1]=3;other_char2[2]=4;break;
	case 3: other_char2[0]=1;other_char2[1]=2;other_char2[2]=4;break;
	case 0: other_char2[0]=1;other_char2[1]=2;other_char2[2]=3;rem2=4;break;
	}
	temp_index2-=rem2;

	for (mm=0;mm<3;mm++) {
	  k++;
	  delta_index=(other_char1[j]-rem1)*powers_of_four[i-1]+(other_char2[mm]-rem2)*powers_of_four[ll+i-1];
	  mismatch_indices[k]=delta_index+index;
	}
	temp_index2/=4;
	ll++;
      }
    }
    i++;
  }
}

void get_textmotif_from_index(int index,char *motif)
{
  int temp_index;
  int rem;
  int i;
  int motif_length;

  if ( index == 0 ) {
    sprintf(motif,"n");
  } else {
    motif_length=get_motiflength_from_index(index);
    switch (motif_length) {
    case 1: strcpy(motif,"n");break;
    case 2: strcpy(motif,"nn");break;
    case 3: strcpy(motif,"nnn");break;
    case 4: strcpy(motif,"nnnn");break;
    case 5: strcpy(motif,"nnnnn");break;
    case 6: strcpy(motif,"nnnnnn");break;
    case 7: strcpy(motif,"nnnnnnn");break;
    case 8: strcpy(motif,"nnnnnnnn");break;
    case 9: strcpy(motif,"nnnnnnnnn");break;
    case 10: strcpy(motif,"nnnnnnnnnn");break;
    }
    
    i=1;
    temp_index=index;
    while (temp_index>0) {
      rem = temp_index % 4;
      switch (rem) {
      case 0:	temp_index-=4;motif[motif_length-i]=116;break;
      case 1:	temp_index-=rem;motif[motif_length-i]=97;break;
      case 2:	temp_index-=rem;motif[motif_length-i]=99;break;
      case 3:	temp_index-=rem;motif[motif_length-i]=103;break;
      }
      temp_index/=4;
      i++;
    }
  }
}

void get_textmotif_from_index_v2(int index,char *motif,int mism1,int mism2)
{
  /* get text motif from index
   * in this version, allow also for mismatches: "n" in positions mism1 and mism2
   * note: mism1=1 for the rightmost 3' position whereas this same position corresponds to 
   *       motif[motif_length-1]. Be careful about this!
   */

  int temp_index;
  int rem;
  int i;
  int motif_length;

  if ( index == 0 ) {
    sprintf(motif,"n");
  } else {
    motif_length=get_motiflength_from_index(index);
    motif_initialization(motif,motif_length);
    
    i=1;
    temp_index=index;
    while (temp_index>0) {
      rem = temp_index % 4;
      switch (rem) {
      case 0:	temp_index-=4;motif[motif_length-i]=116;break;
      case 1:	temp_index-=rem;motif[motif_length-i]=97;break;
      case 2:	temp_index-=rem;motif[motif_length-i]=99;break;
      case 3:	temp_index-=rem;motif[motif_length-i]=103;break;
      }
      temp_index/=4;
      i++;
    }
    
    if (mism1>0) {
      if (mism1<=motif_length) {
	motif[motif_length-mism1]=110;
      }
      if (mism2>0) {
	if (mism2<=motif_length) {
	  motif[motif_length-mism2]=110;
	}
      }
    }
  }
}

int get_motiflength_from_index(int index)
{
  int temp_index;
  int rem;
  int i;

  i=0;
  temp_index=index;
  while (temp_index>0) {
    rem = temp_index % 4;
    if (rem==0)
    	rem=4;
    temp_index-=rem;
    temp_index/=4;
    i++;
  }

  return i;
}

void get_complement_from_textmotif(char *motif,char *motifc,int motif_length)
{
  int i;

  for (i=0;i<motif_length;i++) {
    switch (motif[i]) {
    case AA: motifc[i]=TT;break;
    case CC: motifc[i]=GG;break;
    case GG: motifc[i]=CC;break;
    case TT: motifc[i]=AA;break;
    }
  }
}

int get_complement_from_index(int index)
{
  /* given a motif index, return the index of the complementary motif */
  int temp_index;
  int rem;
  int i;
  int comp_index=0;

  i=0;
  temp_index=index;
  /* this code goes from right-most to left-most */
  while (temp_index>0) {
    rem = temp_index % 4;
    if (rem==0) {
      temp_index-=4;
      comp_index+=(5-4)*powers_of_four[i];
    }
    else {
      temp_index-=rem;
      comp_index+=(5-rem)*powers_of_four[i];
    }
    temp_index/=4;
    i++;
  }
  return comp_index;
}

int isearch_element(int *list,int element,int n)
{
	/* search for the first occurrence of element in list
	 * return index if found and -1 otherwise
	 */
	
	int i=1;
	int index=-1;
	
	while ( (i>0) & (i<n) ) {
		if (list[i]==element) {
			index=i;
			i=-1;
		}		
		i++;	 	
	}			
	return index;
}

int get_index_from_textmotif(char *motif,int motif_length)
{
  /* get the motif index from the text motif
   * faster version than get_index
   */

  /* motif goes from index=0 to index=motif_length-1 */
  int index;
  int i;
  int nummotif;

  index=0;
  i=0;
  while (i<motif_length) {
    // i++;                              // 05-23-2003 i moved this to the end of this loop
    //for (i=0;i<motif_length;i++) {
    nummotif=-1;                         // default in case it is not a,c,g,t
    //printf("%d\n",motif[motif_length-i-1]);
    switch (motif[motif_length-i-1]) {
    case AA: nummotif=1;break; // 97
    case CC: nummotif=2;break; // 99
    case GG: nummotif=3;break; // 103
    case TT: nummotif=4;break; // 116
    otherwise: nummotif=-1;break;
    }
    if (nummotif<0) {
      i=motif_length;
      index=-1;
    } else {
      index+=powers_of_four[i]*nummotif;
    }
    i++;
  }

  return index;
}

void push_ivector(int new_value,int* ivec,int n)
{
	/* push the contents of the vector up to add a new value at index 1 */
	int i;
	for (i=n;i>=2;i--)
		ivec[i]=ivec[i-1];
	ivec[1]=new_value;		
}

void push_int_vector(int new_value,int *ivec,int value_pos,int vec_length)
{
  /* push the contents of the vector to add a new value at a specific location
   * note: index 2 --goes to--> index 1
   *       index 3 --goes to--> index 2
   *       etc
   */
  int i;

  if ( (value_pos>0) & (value_pos<=vec_length) ) {
    for (i=1;i<=(value_pos-1);i++) 
      ivec[i]=ivec[i+1];
    ivec[value_pos]=new_value;
  }

}

void push_float_vector(float new_value,float *ivec,int value_pos,int vec_length)
{
  /* push the contents of the vector to add a new value at a specific location
   * note: index 2 --goes to--> index 1
   *       index 3 --goes to--> index 2
   *       etc
   */
  int i;

  if ( (value_pos>0) & (value_pos<=vec_length) ) {
    for (i=1;i<=(value_pos-1);i++) 
      ivec[i]=ivec[i+1];
    ivec[value_pos]=new_value;
  }
}

void push_double_vector(double new_value,double *dvec,int value_pos,int vec_length)
{
  /* push the contents of the vector to add a new value at a specific location
   * note: index 2 --goes to--> index 1
   *       index 3 --goes to--> index 2
   *       etc
   */
  int i;

  if ( (value_pos>0) & (value_pos<=vec_length) ) {
    for (i=1;i<=(value_pos-1);i++) 
      dvec[i]=dvec[i+1];
    dvec[value_pos]=new_value;
  }
}

void indices_addnt(int index,int motif_length,int *indices)
{
  /* given an index number, compute the indices of all possible motifs that have one additional nucleotide at the ex */
  int temp_index;
  int rem;
  int i;
  int right_index=0;

	/* add to the right */
  i=1;
  temp_index=index;
  while (temp_index>0) {
    rem = temp_index % 4;
    if (rem==0) {
      temp_index-=4;
			right_index+=4*powers_of_four[i];
    }
    else {
      temp_index-=rem;
			right_index+=rem*powers_of_four[i];
    }
    temp_index/=4;
    i++;
  }
  for (i=1;i<=4;i++)
		indices[i]=right_index+i;

	/* add to the left */
	for (i=1;i<=4;i++)
		indices[4+i]=powers_of_four[motif_length]*i+index;
			
}

void motif_initialization(char *motif,int motif_length)
{
  switch (motif_length) {
  case 1: strcpy(motif,"n");break;
  case 2: strcpy(motif,"nn");break;
  case 3: strcpy(motif,"nnn");break;
  case 4: strcpy(motif,"nnnn");break;
  case 5: strcpy(motif,"nnnnn");break;
  case 6: strcpy(motif,"nnnnnn");break;
  case 7: strcpy(motif,"nnnnnnn");break;
  case 8: strcpy(motif,"nnnnnnnn");break;
  case 9: strcpy(motif,"nnnnnnnnn");break;
  case 10: strcpy(motif,"nnnnnnnnnn");break;
  case 11: strcpy(motif,"nnnnnnnnnnn");break;
  case 12: strcpy(motif,"nnnnnnnnnnnn");break;
  case 13: strcpy(motif,"nnnnnnnnnnnnn");break;
  case 14: strcpy(motif,"nnnnnnnnnnnnnn");break;
  case 15: strcpy(motif,"nnnnnnnnnnnnnnn");break;
  }
}

int find_pos(float *vec,float value,int vec_length)
{
  /* find the position of value within vec assuming that vec is sorted in ascending order */
  /* 12-11-2002: if it is equal to the maximum, set the value position to be that of the maximum */
  int value_pos;
  int i;

  if (value>=vec[vec_length]) {
    value_pos=vec_length;
  } else {
    i=1;
    value_pos=vec_length;
    while (i<=vec_length) {
      if (vec[i]>=value) {
	value_pos=i-1;
	i=vec_length+1;
      }
      i++;
    }    
  }
  return value_pos;

}

int find_dpos(double *vec,double value,int vec_length)
{
  /* find the position of value within vec assuming that vec is sorted in ascending order */
  /* 12-11-2002: if it is equal to the maximum, set the value position to be that of the maximum */
  
  int value_pos;
  int i;

  if (value>=vec[vec_length]) {
    value_pos=vec_length;
  } else {
    i=1;
    value_pos=vec_length;
    while (i<=vec_length) {
      if (vec[i]>=value) {
	value_pos=i-1;
	i=vec_length+1;
      }
      i++;
    }    
  }
  return value_pos;

}

int faster_find_pos(float *vec,float value,int vec_length)
{
  /* find the position of value within vec assuming that vec is sorted in ascending order */
  int value_pos;
  int i;
  int i_inf,i_sup;

  i_inf=1;
  i_sup=vec_length;
  while ( (i_sup-i_inf)>1 ) {
    i=(i_sup+i_inf)/2;
    if (vec[i]>=value) {
      i_sup=i;
    }
    else {
      i_inf=i;
    }
  }
  value_pos=i_inf;

  return value_pos;
}

int faster_find_dpos(double *vec,double value,int vec_length)
{
  /* find the position of value within vec assuming that vec is sorted in ascending order */
  /* here double values instead of float */
  int value_pos;
  int i;
  int i_inf,i_sup;

  i_inf=1;
  i_sup=vec_length;
  while ( (i_sup-i_inf)>1 ) {
    i=(i_sup+i_inf)/2;
    if (vec[i]>=value) {
      i_sup=i;
    }
    else {
      i_inf=i;
    }
  }
  value_pos=i_inf;

  return value_pos;
}

int compare_2motifs_v1(char *motif1,char *motif2,int motif_length,int max_overlap,int mism1_m1,int mism2_m1,int mism1_m2,int mism2_m2)
{
  int i;
  int shift;
  int max_shift;
  int max_coincidences=0;
  int curr_coincidences;

  mism1_m1=motif_length-mism1_m1;
  mism2_m1=motif_length-mism2_m1;
  mism1_m2=motif_length-mism1_m2;
  mism2_m2=motif_length-mism2_m2;

  max_shift=motif_length-max_overlap;
  if (max_shift<0) {
    max_coincidences=0;
    return max_coincidences;
    /* error: max_overlap must be <= motif_length */
  }

  shift=-max_shift;
  while ( (max_coincidences<max_overlap) & (shift<=max_shift) ) {
    curr_coincidences=0;
    if (shift>=0) {
      /* shift motif 1 to the right */
      for (i=0;i<(motif_length-shift);i++) {
	if ( (i==mism1_m1) | (i==mism2_m1) | (i==mism1_m2) | (i==mism2_m2) ) {
	  curr_coincidences++;
	}
	else {
	  if (motif1[i]==motif2[i+shift]) 
	    curr_coincidences++;
	}
      }
    }
    else {
      /* shift motif 1 to the left */
      for (i=-shift;i<motif_length;i++) {
	if ( (i==mism1_m1) | (i==mism2_m1) | (i==mism1_m2) | (i==mism2_m2) ) {
	  curr_coincidences++;
	}
	else {
	  if (motif1[i]==motif2[i+shift])
	    curr_coincidences++;
	}
      }
    }
    if (curr_coincidences>max_coincidences)
      max_coincidences=curr_coincidences;

    shift++;
  }
  
  return max_coincidences;
}

int count_ns(char *sequence)
{
  /* count the number of "n" in a given string */
  int n;
  int sequence_length;
  int i;

  n=0;
  sequence_length=strlen(sequence);
  for (i=0;i<sequence_length;i++) {
    if (sequence[i]==NN) {
      n++;
    }
  }
  return n;
}

void seq2num(char *sequence,int sequence_length,char *numsequence,int init_i)
{
  int i,j;
  char c;
  char numc;

  j=0;
  for (i=init_i;i<sequence_length;i++) {
    c=sequence[i];
    numc=char2num(c);
    j++;
    numsequence[j]=numc;
    //printf("%d\t%d\t%d\n",c,numc,numsequence[j]);
  }
}

void num2seq(char *numsequence,int sequence_length,char *sequence,int init_i)
{
  int i,j;
  char c;
  char numc;

  for (i=init_i;i<=sequence_length;i++) {
    numc=numsequence[i];
    c=num2char(numc);
    //printf("%d %d %d\n",i,numc,c);
    j=i-1;
    sequence[j]=c;
    //printf("%d %d %d\n",i,numc,c);
  }
  //printf("sequence=%s\n",sequence);
}

char char2num(char c)
{
  /* nucleotide to number conversion 
   * a -> 1
   * c -> 2
   * g -> 3
   * t -> 4
   */
  char numchar;

  numchar=0;
  switch (c) {
  case 97: numchar=1;break;
  case 99: numchar=2;break;
  case 103: numchar=3;break;
  case 116: numchar=4;break;
  }

  return numchar;
}

char num2char(char nc)
{
  /* number to nucleotide conversion
   * 1 -> a
   * 2 -> c
   * 3 -> g
   * 4 -> t
   */
  char c;

  c='n';
  switch (nc) {
  case 1: c='a';break;
  case 2: c='c';break;
  case 3: c='g';break;
  case 4: c='t';break;
  }
  
  return c;
}
