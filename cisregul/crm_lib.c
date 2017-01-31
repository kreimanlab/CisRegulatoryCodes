/* crm_lib.c 
 * functions to search for a module of transcription factor binding sites (tfbs)
 * last modified 01-02-2004
 */ 

/* 01-31-2004: added load_and_sort_scan_data_v3
 * 01-31-2004: added filter_positions_v5
 * 01-02-2004: added load_and_sort_scan_data_v2
 * 12-30-2004: added filter_positions_v3
 * 12-30-2004: added filter_positions_v4
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/************************/
/* function definitions */
/************************/
int filter_positions(int *pos,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos);
int filter_positions_v2(int *pos,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos);
int filter_positions_v3(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos,char *filtered_strand);
int filter_positions_v4(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos,char *filtered_strand);
int filter_positions_v5(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int el,int *filtered_pos,char *filtered_strand);
int load_and_sort_scan_data(char *filename,int **data_t,int **data_p,int *motif_lengths,int neg2pos,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths);
int load_and_sort_scan_data_v2(char *filename_positions,char *filename_strands,int **data_t,int **data_p,int **data_s,int *motif_lengths,int neg2pos,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths);
int load_and_sort_scan_data_v3(char *filename_positions,char *filename_strands,int **data_t,int **data_p,int **data_s,int *motif_lengths,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths,int data_format,int filteron,int verbose);

void crm_faster_1motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output);
void crm_faster_2motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output);
void crm_faster_3motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int m3,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output);
void crm_faster_4motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int m3,int m4,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output);
void count_motif(int *m_p,int n_occurrences,int init_index,int lower_border,int upper_border,int *out);

void crm_1motifs(int *m1_p,int n1,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions);
void crm_1motifs_allt(int *m1_t,int *m1_p,int n1,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions);
void crm_2motifs(int *m1_p,int *m2_p,int n1,int n2,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions);
void crm_2motifs_allt(int *m1_t,int *m2_t,int *m1_p,int *m2_p,int n1,int n2,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions);
void crm_3motifs(int *m1_p,int *m2_p,int *m3_p,int n1,int n2,int n3,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions);
void crm_3motifs_allt(int *m1_t,int *m2_t,int *m3_t,int *m1_p,int *m2_p,int *m3_p,int n1,int n2,int n3,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int report_all,int n_transcripts,int *n,int **crm_positions);
void crm_4motifs(int *m1_p,int *m2_p,int *m3_p,int *m4_p,int n1,int n2,int n3,int n4,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions); 
void crm_4motifs_allt(int *m1_t,int *m2_t,int *m3_t,int *m4_t,int *m1_p,int *m2_p,int *m3_p,int *m4_p,int n1,int n2,int n3,int n4,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions);

int vector_overlap_v2(int *original_vector,int n,int max_overlap,int *output_vector,char *original_indices_vector,char *output_indices_vector);
int vector_overlap_v3(int *original_pos,char *original_strand,int *output_pos,char *output_strand,int n,int max_overlap);

void convert_positions2tss(int *pos,int n,int el,int il);
int convert_position2tss(int pos,int el,int il);
void reflect_pos(int *pos,int n);

/*************/
/* functions */
/*************/

int filter_positions(int *pos,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos)
{
  /* filter the positions in *pos
   * so as to restrict the search to within max_ul, max_el and max_il
   * positions in pos are measured with respect to the 3' end
   * we want to restrict to within max_ul of the tss
   *                     to within max_el of the tss
   *                     to within max_il of the 1st exon - 1st intron boundary
   * ul,el,il are the actual lengths in the current sequence
   * the filtered positions are returned in filtered_pos
   * the output number is the output of this function
   */

  int i;
  int in;
  int nout=0;             // number of ouptut entries
  int current_pos;

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    if (current_pos>il) {
      current_pos=current_pos-il;
      if (current_pos>el) {
	/* occurrence within the upstream sequence */
	/* debug here 
	   printf("occurrence in the upstram sequence %d --> %d ",pos[i],current_pos); */
	current_pos=current_pos-el;
	//printf("--> %d\n",current_pos);
	if (current_pos<=max_ul) {
	  nout++;
	  //filtered_pos[nout]=-current_pos;
	  filtered_pos[nout]=pos[i];
	}   // close check on ul restriction
      } else {
	/* occurrence within the first exon */
	/* debug here 
	   printf("occurrence in the exon %d --> %d ",pos[i],current_pos); */
	current_pos=el-current_pos;                    // measure with respect to the tss
	//printf("--> %d\n",current_pos);
	if (current_pos<=max_el) {
	  nout++;
	  //filtered_pos[nout]=current_pos;
	  //filtered_pos[nout]=pos[1];             // NOTE: BIG BUG! solved on 01-01-2004
	  filtered_pos[nout]=pos[i];      
	}   // close check on el restriction
      }     // close check on whether it is on the first exon
    } else {
      /* occurrence within the first intron */
      /* debug here  
	 printf("occurrence in the intron %d --> %d ",pos[i],current_pos); */
      current_pos=il-current_pos;                      // measure with respect to the exon-intron boundary
      //printf("--> %d\n",current_pos);
      if (current_pos<=max_il) {
	nout++;

	filtered_pos[nout]=pos[i];
      }      // close check on il restriction
    }        // close check whether it is in the first intron
  }          // close i loop

  return nout;
}

int filter_positions_v2(int *pos,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos)
{
  /* filter the positions in *pos
   * so as to restrict the search to within max_ul, max_el and max_il
   * positions in pos are measured with respect to the 3' end
   * we want to restrict to within max_ul of the tss
   *                     to within max_el of the tss
   *                     to within max_il of the 1st exon - 1st intron boundary
   * ul,el,il are the actual lengths in the current sequence
   * the filtered positions are returned in filtered_pos
   * the output number is the output of this function
   *
   * new in version 2:
   * use scan output used for all sequences (measured with respect to TSS, upstream sequences are > 0, exon and intron are < 0)
   */

  int i;
  int in;
  int nout=0;             // number of ouptut entries
  int current_pos;

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    if (current_pos>0) {
      /* occurrence within the upstream region */
      if (current_pos<=max_ul) {
	nout++;
	filtered_pos[nout]=el+il+current_pos;           // convert to positions with respect to the 3' end of the sequence
      }   // close check on ul restriction
    } else {
      current_pos=-current_pos;
      if (current_pos<=el) {
	/* occurrence within the first exon */
	if (current_pos<=max_el) {
	  current_pos=el-current_pos;
	  nout++;
	  filtered_pos[nout]=current_pos+il;            // convert to positions with respect to the 3' end of the sequence
	}   // close check on el restriction
      } else {     // close check on whether it is on the first exon
      /* occurrence within the first intron */
	current_pos=current_pos-el;
	if (current_pos<=max_il) {
	  nout++;
	  filtered_pos[nout]=il-current_pos;            // convert to positions with respect to the 3' end of the sequence
	}
      }      // close check on il restriction
    }        // close check whether it is in the upstream region or not
  }          // close i loop

  return nout;
}

int filter_positions_v3(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos,char *filtered_strand)
{
  /* filter the positions in *pos so as to restrict the search to within max_ul, max_el and max_il
   * positions in pos are measured with respect to the 3' end
   * we want to restrict to within max_ul of the tss
   *                     to within max_el of the tss
   *                     to within max_il of the 1st exon - 1st intron boundary
   * ul,el,il are the actual lengths in the current sequence
   * the filtered positions are returned in filtered_pos
   * the output number is the output of this function
   *
   * new in version 3: this actually comes from v1, here we add the strand information
   */

  int i;
  int in;
  int nout=0;             // number of ouptut entries
  int current_pos;

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    if (current_pos>il) {
      current_pos=current_pos-il;
      if (current_pos>el) {	                /* occurrence within the upstream sequence */
	current_pos=current_pos-el;
	if (current_pos<=max_ul) {
	  nout++;
	  filtered_pos[nout]=pos[i];
	  filtered_strand[nout]=strand[i];
	}                                       // close check on ul restriction
      } else {                          	/* occurrence within the first exon */
	current_pos=el-current_pos;             // measure with respect to the tss
	if (current_pos<=max_el) {
	  nout++;
	  filtered_pos[nout]=pos[i];
	  filtered_strand[nout]=strand[i];
	}                                       // close check on el restriction
      }                                         // close check on whether it is on the first exon
    } else {                                   /* occurrence within the first intron */
      current_pos=il-current_pos;                      // measure with respect to the exon-intron boundary
      if (current_pos<=max_il) {
	nout++;
	filtered_pos[nout]=pos[i];
	filtered_strand[nout]=strand[i];
      }      // close check on il restriction
    }        // close check whether it is in the first intron
  }          // close i loop

  return nout;
}

int filter_positions_v4(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int ul,int el,int il,int *filtered_pos,char *filtered_strand)
{
  /* filter the positions in *pos
   * so as to restrict the search to within max_ul, max_el and max_il
   * positions in pos are measured with respect to the 3' end
   * we want to restrict to within max_ul of the tss
   *                     to within max_el of the tss
   *                     to within max_il of the 1st exon - 1st intron boundary
   * ul,el,il are the actual lengths in the current sequence
   * the filtered positions are returned in filtered_pos
   * the output number is the output of this function
   *
   * new in version 2:
   * use scan output used for all sequences (measured with respect to TSS, upstream sequences are > 0, exon and intron are < 0)
   * new in version 4: this comes from v2, use strand information
   */

  int i;
  int in;
  int nout=0;             // number of ouptut entries
  int current_pos;

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    if (current_pos>0) {                                /* occurrence within the upstream region */
      if (current_pos<=max_ul) {
	nout++;
	// convert to positions with respect to the 3' end of the sequence
	if (max_il==0) {
	  if (max_el==0) {
	    filtered_pos[nout]=current_pos;
	  } else {
	    filtered_pos[nout]=current_pos+el;
	  }
	} else {
	  filtered_pos[nout]=el+il+current_pos;           
	}
	filtered_strand[nout]=strand[i];
      }                                                 // close check on ul restriction
    } else {
      current_pos=-current_pos;
      if (current_pos<=el) {                    	/* occurrence within the first exon */
	if (current_pos<=max_el) {
	  current_pos=el-current_pos;
	  nout++;
	  filtered_pos[nout]=current_pos+il;            // convert to positions with respect to the 3' end of the sequence
	  filtered_strand[nout]=strand[i];    
	}                                               // close check on el restriction
      } else {                                          // close check on whether it is on the first exon
	/* occurrence within the first intron */
	current_pos=current_pos-el;
	if (current_pos<=max_il) {
	  nout++;
	  filtered_pos[nout]=il-current_pos;            // convert to positions with respect to the 3' end of the sequence
	  filtered_strand[nout]=strand[i];
	}
      }      // close check on il restriction
    }        // close check whether it is in the upstream region or not
  }          // close i loop

  return nout;
}

int filter_positions_v5(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int el,int *filtered_pos,char *filtered_strand)
{
  /* filter the positions in *pos
   * so as to restrict the search to within max_ul, max_el and max_il
   * positions in pos are measured with respect to the TSS
   * we want to restrict to within max_ul of the tss
   *                     to within max_el of the tss
   *                     to within max_il of the 1st exon - 1st intron boundary
   * ul,el,il are the actual lengths in the current sequence
   * the filtered positions are returned in filtered_pos
   * the total number of entries after filtering is the output of this function
   *
   * new in version 2:
   * use scan output used for all sequences (measured with respect to TSS, upstream sequences are > 0, exon and intron are < 0)
   * new in version 4: this comes from v2, use strand information
   * new in version 5: 
   * 01-31-2004 this comes from v4, here positions are measured with respect to the TSS, upstream is <0, a larger number indicates more 3'
   * 01-31-2004 return positions in the same format (with respect to TSS)
   */

  int i;                  // loop counter
  int nout=0;             // number of ouptut entries
  int current_pos;        // current position
  int current_strand;     // current strand

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    current_strand=strand[i];
    if (current_pos<0) {                                /* occurrence within the upstream region */
      current_pos=-current_pos;
      if (current_pos<=max_ul) {
	nout++;
	//if (max_il==0) {
	//if (max_el==0) {
	filtered_pos[nout]=-current_pos;
	//} else {
	//filtered_pos[nout]=current_pos+el;
	//}
	//} else {
	//filtered_pos[nout]=el+il+current_pos;           
	//}
	filtered_strand[nout]=current_strand;
      }                                                 // close check on ul restriction
    } else {                           
      /* occurrence within the exon or the intron */
      //current_pos=-current_pos;
      if (current_pos<=el) {                    	/* occurrence within the first exon */
	if (current_pos<=max_el) {
	  //current_pos=el-current_pos;
	  nout++;
	  //filtered_pos[nout]=current_pos+il;            // convert to positions with respect to the 3' end of the sequence
	  filtered_pos[nout]=current_pos;
	  filtered_strand[nout]=current_strand;    
	}                                               // close check on el restriction
      } else {                                          // close check on whether it is on the first exon
	/* occurrence within the first intron */
	current_pos=current_pos-el;
	if (current_pos<=max_il) {
	  nout++;
	  //filtered_pos[nout]=il-current_pos;            // convert to positions with respect to the 3' end of the sequence
	  filtered_pos[nout]=current_pos+el;
	  filtered_strand[nout]=current_strand;
	}    // close check on il restriction
      }      // close check on current_pos<=max_el
    }        // close check whether it is in the upstream region or not
  }          // close i loop

  return nout;
}

int load_and_sort_scan_data(char *filename,int **data_t,int **data_p,int *motif_lengths,int neg2pos,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths)
{
  /* this function loads output data from a motif scan and sorts the positions in each transcript
   * the format for the motif scan is the following:
   * 1             2       3          4        5          6
   * motif_length  n_seqs  transcript position transcript position ...
   * one line per motif
   * all positions measured with respect to the 3' end (the larger the value, the more 5' the position)
   *
   * input
   * filename: scan file name
   * data_t: where the output transcript numbers will be stored, one row per motif, one column per occurrence (n_motifs x max_occur)
   * data_p: where the output positions numbers will be stored (n_motifs x max_occur)
   * motif_lengths: lenght of all motifs (n_motifs x 1)
   * neg2pos: 1 to convert negative positions to positive
   * sortpos: 1 to sort positions
   * max_occur: maximum number of occurrences
   * n_motifs: number of motifs
   * max_ul: maximum upstream length now
   * max_el: maximum exon length now
   * max_il: maximum intron length now
   * transcript_lengths: for each transcript, actual upstream, exon and intron lengths (n_transcripts x 3)
   *
   * 10-27-2003: avoid using global variables
   * 10-27-2003: allow to study only specific sequence intervals
   * 10-16-2003: avoid overlapping occurrences by considering the motif length
   * 10-17-2003: add neg2pos
   * 10-17-2003: add sortpos
   * 10-17-2003: make it shorter by incorporating the header stuff and first line directly into the loop
   * return an error code if the number of motifs read does not match the actual number of motifs
   */

  FILE *input_file;       // pointer to input file
  char txtline[1000000];  // read one line of the input file
  char *gotit;            // used in the search process
  int header;           // 1 while reading the header
  int i,j,k,ll,m;         // loop counters
  int *v;                 // vector holding the transcripts and positions for running entry
  int start=3;            // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[2];              // number of entries
  int n;                  // number of entries
  int linenr=0;           // current line number (each line corresponds to a separate motif)
  char lineret;           // carriage return;
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int twice_max_occur;                  // 2 * max_occur
  int *temp_positions;                  // temporary positions, used to sort the positions in current transcript (in ascending order)
  float *dummy;                         // used in call to sorting function
  int previous_transcript;              // previous transcript
  int curr_transcript;                  // current transcript
  int n_occurrences;                    // number of occurrences in current transcript
  int max_motif_overlap;                // only allow occurrences that overlap by less than this value
  float max_motif_overlap_factor=0.5;   // max_motif_overlap = (ceil) max_motif_overlap_factor * motif_length;
  int n_no_overlap;                     // number of entries in current transcript after removing overlaps
  int *positions_no_overlap;            // entries in current transcript after removing overlaps
  long cumulative_nofilter=0;
  long cumulative_occurrences=0;        // cumulative sum of n_occurrences
  long cumulative_no_overlap=0;         // cumulative sum of n_no_overlap
  int lines_read=0;                     // total number of lines read in the file (including header lines)
  int lines_overflow=0;                 // total number of lines with a memory overflow
  char infotext[10000];                 // log information
  char temptext[1000];                  // log information
  int errorcode=0;                      // return errorcode = 1 if linenr does not match n_motifs at the end
  int *filtered_positions;              // filter for length
  int n_after_filter;                   // number of entries in current transcript after filtering
  int curr_ul,curr_el,curr_il;          // lengths for current transcript
  int n_filtered_occurrences;           // in each transcript, number of occurrences after filtering and removing overlaps

  lineret='\n';
  
  /* memory allocation */
  twice_max_occur=2*max_occur;
  v=ivector(1,twice_max_occur);
  temp_positions=ivector(1,max_occur);
  positions_no_overlap=ivector(1,max_occur);
  dummy=vector(1,max_occur);
  filtered_positions=ivector(1,max_occur);
  
  /* open file */
  input_file=fopen(filename,"r");
  if (!input_file) {
    printf("error! i could not open the file %s",filename);
    exit(1);
  } 

  while ( (feof(input_file) == 0) & (linenr<n_motifs) ) {
    lines_read++;
    fgets(txtline,n_chars,input_file);
    header=is_header(txtline);
    if (header) {
      /* this is a header line, ignore it */
    } else {
      linenr++;      /* debug here      printf("-------------------  linenr=%d  -------------------------- \n",linenr); */
      gotit=strchr(txtline,lineret);
      if (!gotit) {
	//fprintf(log_file,"ERROR!!!\tfilename=%s\nlines_read=%d\tlinenr=%d\nthere has been an overflow,strlen(txtline)=%d\n",filename,lines_read,linenr,strlen(txtline));
	printf("ERROR!!!\tfilename=%s\nlines_read=%d\tlinenr=%d\nthere has been an overflow,strlen(txtline)=%d\n",filename,lines_read,linenr,strlen(txtline));
	lines_overflow++;
      }
      splitline2int(txtline,v,start,np,twice_max_occur);        // separate the text line txtline into integers stored in v
      n=np[0];                                                  // total number of entries in this line
      if (n>(twice_max_occur)) {
	//fprintf(log_file,"ERROR!!!\nfilename=%s\nlinenr=%d\nn=%d,max_occur=%d\nexceeded the limit of number of occurrences, please increase max_occur\n",filename,linenr,n,max_occur);
	printf("ERROR!!!\nfilename=%s\nlinenr=%d\nn=%d,max_occur=%d\nexceeded the limit of number of occurrences, please increase max_occur\n",filename,linenr,n,max_occur);
	n=twice_max_occur;
      }
      if (linenr>n_motifs) {
	printf("ERROR! linenr=%d should be <= n_motifs=%d\n");
	exit(1);
      }
      motif_lengths[linenr]=v[1];      /* debug here       printf("linenr=%d\tmotif length=%d\tn=%d\n",linenr,v[1],n); */
      ll=0;                            // running number of occurrences for current motif
      if (n>0) {                       // if there are any occurrences for this motif in any transcript
	
	//DEBUG HERE, SHOULDN'T WE STORE THE NUMBER OF ENTRIES AFTER WE FILTERED THEM BY THE LENGTHS?
	data_t[linenr][0]=n/2-1;       // store the number of entries in position 0
	data_p[linenr][0]=n/2-1;
	previous_transcript=-1;        // initial previous transcript
	n_occurrences=0;               // number of occurrences in current transcript
	n_filtered_occurrences=0;      // number of occurrences after filtering and removing overlaps for current motif
	for (j=3;j<n;j+=2) {	       // for each occurrence (starts at 3 because v[1] is the motif length and v[2] the number of sequences  
	  curr_transcript=v[j];        // current transcript
	  k=j+1;                       // index for position
	  if (neg2pos==1) {            // convert negative entries to positive ones if and only if requested to do so
	    if (v[k]<0)
	      v[k]=-v[k];
	  }
	
	  if (curr_transcript == previous_transcript) {
	    n_occurrences++;
	    temp_positions[n_occurrences]=v[k];
	  } else {
	    if (previous_transcript<0) {	         /* this was the first entry, just add it */
	      n_occurrences++;
	      temp_positions[n_occurrences]=v[k];
	    } else {                    	        /* ready with this transcript */
	      if (n_occurrences>0) {                    /* there were any occurrences in the current transcript */
		cumulative_nofilter+=n_occurrences;
		n_after_filter=filter_positions(temp_positions,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
		if (n_after_filter>0) {
		  if (sortpos)                          /* if we need to sort the positions */
		    quicksorti(n_after_filter,filtered_positions,dummy);
		  //quicksorti(n_occurrences,temp_positions,dummy);        // sort the positions in current transcript
		  //cumulative_occurrences+=n_occurrences;
		  cumulative_occurrences+=n_after_filter;
		  max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
		  //n_no_overlap=vector_overlap(temp_positions,n_occurrences,max_motif_overlap,positions_no_overlap);
		  n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap);
		  cumulative_no_overlap+=n_no_overlap;
		  //for (m=1;m<=n_occurrences;m++) {
		  for (m=1;m<=n_no_overlap;m++) {
		    ll++;
		    data_t[linenr][ll]=previous_transcript;
		    //data_p[linenr][ll]=temp_positions[m];
		    data_p[linenr][ll]=positions_no_overlap[m];
		  }
		  n_filtered_occurrences+=n_no_overlap;
		}                                        // close check on n_after_filter>0
	      }                                          // close check on n_occurrences>0
	      n_occurrences=1;                           // reset the number of occurrences to 1
	      temp_positions[n_occurrences]=v[k];        // and write the first position for the new transcript
	    }                                            // close check on previous_transcript<0
	    previous_transcript=curr_transcript;         // set the previous_transcript to the current transcript
	    curr_ul=transcript_lengths[curr_transcript][1];
	    curr_el=transcript_lengths[curr_transcript][2];
	    curr_il=transcript_lengths[curr_transcript][3];
	  }                                              // close check on curr_transcript == previous_transcript
	}                                                // close j loop
     
	/* process the last transcript */
	if (n_occurrences>0) {
	  n_after_filter=filter_positions(temp_positions,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	  cumulative_nofilter+=n_occurrences;
	  if (n_after_filter>0) {
	    if (sortpos)
	      quicksorti(n_after_filter,filtered_positions,dummy);
	    cumulative_occurrences+=n_after_filter;
	    max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
	    n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap);
	    cumulative_no_overlap+=n_no_overlap;
	    for (m=1;m<=n_no_overlap;m++) {
	      ll++;
	      data_t[linenr][ll]=previous_transcript;
	      data_p[linenr][ll]=positions_no_overlap[m];
	    }
	    n_filtered_occurrences+=n_no_overlap;
	  }      // close check on n_after_filter>0
	}        // close check on n_occurrences>0 for the last transcript
	/* update the actual number of occurrences */
	data_t[linenr][0]=n_filtered_occurrences;       // store the number of entries in position 0
	data_p[linenr][0]=n_filtered_occurrences;
      }    // if n>0 (i.e. if there was any occurrence of this motif in any transcript
    }      // close check on wether this is a header line
  }        // while there is still information in the file
  fclose(input_file);

  /* print out some information about the reading process */
  sprintf(infotext,"\nnumber of lines read = %d\n",lines_read);
  sprintf(temptext,"number of lines with overflow = %d\n",lines_overflow);strcat(infotext,temptext);
  sprintf(temptext,"number of lines with scan data = %d\n",linenr);strcat(infotext,temptext);
  if (linenr != n_motifs) {
    sprintf(temptext,"\n\nERROR! linenr=%d and n_motifs=%d\n\n",linenr,n_motifs);strcat(infotext,temptext);
    errorcode=1;
  }
  sprintf(temptext,"cumulative number of occurrences (before filter)= %ld\n",cumulative_nofilter);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences (after filter)= %ld\n",cumulative_occurrences);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences with no overlap = %ld",cumulative_no_overlap);strcat(infotext,temptext);
  //printinfo(infotext,1,log_file);
  printf("%s\n",infotext);

  free_ivector(v,1,twice_max_occur);
  free_ivector(temp_positions,1,max_occur);
  free_vector(dummy,1,max_occur);
  free_ivector(positions_no_overlap,1,max_occur);
  free_ivector(filtered_positions,1,max_occur);

  return errorcode;
}

int load_and_sort_scan_data_v2(char *filename_positions,char *filename_strands,int **data_t,int **data_p,int **data_s,int *motif_lengths,int neg2pos,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths)
{
  /* this function loads output data from a motif scan and sorts the positions in each transcript
   * the format for the motif scan is the following:
   * 1             2       3          4        5          6
   * motif_length  n_seqs  transcript position transcript position ...
   * one line per motif
   * all positions measured with respect to the 3' end (the larger the value, the more 5' the position)
   *
   * input
   * filename: scan file name
   * data_t: where the output transcript numbers will be stored, one row per motif, one column per occurrence (n_motifs x max_occur)
   * data_p: where the output positions numbers will be stored (n_motifs x max_occur)
   * motif_lengths: lenght of all motifs (n_motifs x 1)
   * neg2pos: 1 to convert negative positions to positive
   * sortpos: 1 to sort positions
   * max_occur: maximum number of occurrences
   * n_motifs: number of motifs
   * max_ul: maximum upstream length now
   * max_el: maximum exon length now
   * max_il: maximum intron length now
   * transcript_lengths: for each transcript, actual upstream, exon and intron lengths (n_transcripts x 3)
   *
   * new in version 2:
   * 01-01-2004: load also strand information
   * return an error code if the number of motifs read does not match the actual number of motifs
   * 01-30-2004: call vector_overlap_v3
   * 01-30-2004: do not filter if max_ul<0 (allow to skip the filtering step if we do not want to filter)
   */

  FILE *input_file_positions;       // pointer to input file for positions
  FILE *input_file_strands;         // pointer to input file for strands
  char txtline_pos[1000000];        // read one line of the input file
  char txtline_strand[1000000];     // read one line of the input file
  char *gotit;                      // used in the search process
  int header;                       // 1 while reading the header
  int i,j,k,ll,m;                   // loop counters
  int *v;                           // vector holding the transcripts and positions for running entry
  int *w;                           // vector holding the transcripts and strands for running entry
  int start=3;                      // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[2];                        // number of entries
  int n_pos;                        // number of entries in current line of pos file
  int n_strand;                     // number of entries in current line of strand file
  int linenr=0;                     // current line number (each line corresponds to a separate motif)
  char lineret;                     // carriage return;
  int n_chars=1000000;              // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int twice_max_occur;              // 2 * max_occur
  int *temp_positions;              // temporary positions, (raw positions, before filtering and sorting)
  char *temp_strands;               // temporary strands, (raw strands, before filtering and sorting)
  float *indices;                   // used in call to sorting function
  int previous_transcript;          // previous transcript
  int curr_transcript;              // current transcript
  int curr_transcript_strand;       // current transcript from the strand file
  int n_occurrences;                // number of occurrences in current transcript
  int max_motif_overlap;            // only allow occurrences that overlap by less than this value
  float max_motif_overlap_factor=0.5;   // max_motif_overlap = (ceil) max_motif_overlap_factor * motif_length;
  int n_no_overlap;                 // number of entries in current transcript after removing overlaps
  int *positions_no_overlap;        // entries in current transcript after removing overlaps
  char *strands_no_overlap;         // strands in current transcript after removing overlap
  long cumulative_nofilter=0;       
  long cumulative_occurrences=0;        // cumulative sum of n_occurrences
  long cumulative_no_overlap=0;         // cumulative sum of n_no_overlap
  int lines_read=0;                     // total number of lines read in the file (including header lines)
  int lines_overflow=0;                 // total number of lines with a memory overflow
  char infotext[10000];                 // log information
  char temptext[1000];                  // log information
  int errorcode=0;                      // return errorcode = 1 if linenr does not match n_motifs at the end
  int *filtered_positions;              // filter for length, positions
  char *filtered_strands;               // filter for length, strands
  int n_after_filter;                   // number of entries in current transcript after filtering
  int curr_ul,curr_el,curr_il;          // lengths for current transcript
  int n_filtered_occurrences;           // in each transcript, number of occurrences after filtering and removing overlaps
  int filteron=1;                       // 1 to filter and 0 to skip filtering

  lineret='\n';
  if (max_ul<0)
    filteron=0;
  
  /* memory allocation */
  twice_max_occur=2*max_occur;                    // total maximum allowed number of entries, including the transcript and position
  v=ivector(1,twice_max_occur);                   // allocate space to store all these numbers
  w=ivector(1,twice_max_occur);                   // allocate space for strands (this contains both transcripts and strands, therefore it must be an integer)
  temp_positions=ivector(1,max_occur);            // temporary positions 
  temp_strands=cvector(1,max_occur);              // temporary strands
  positions_no_overlap=ivector(1,max_occur);      // positions after removing overlaps
  strands_no_overlap=cvector(1,max_occur);        // strands after removing overlaps
  indices=vector(1,max_occur);                    // indices when sorting the entries
  filtered_positions=ivector(1,max_occur);        // positions after filtering
  filtered_strands=cvector(1,max_occur);          // strands after filtering
  
  /* open files */
  input_file_positions=fopen(filename_positions,"r");
  if (!input_file_positions) {
    printf("error! i could not open the file %s",filename_positions);
    exit(1);
  } 
  input_file_strands=fopen(filename_strands,"r");
  if (!input_file_strands) {
    printf("error! i could not open the file %s",filename_strands);
    exit(1);
  }

  while ( (feof(input_file_positions) == 0) & (feof(input_file_strands)==0) & (linenr<n_motifs) ) {
    lines_read++;
    fgets(txtline_pos,n_chars,input_file_positions);
    fgets(txtline_strand,n_chars,input_file_strands);
    header=is_header(txtline_pos);
    if (!header) {      /* this is not a header line, ignore it */
      linenr++;         /* debug here      printf("-------------------  linenr=%d  -------------------------- \n",linenr); */
      gotit=strchr(txtline_pos,lineret);
      if (!gotit) {
	printf("ERROR!!!\tfilename_positions=%s\nlines_read=%d\tlinenr=%d\nthere has been an overflow,strlen(txtline)=%d\n",filename_positions,lines_read,linenr,strlen(txtline_pos));
	lines_overflow++;
      }
      splitline2int(txtline_pos,v,start,np,twice_max_occur);        // separate the text line txtline into integers stored in v
      n_pos=np[0];                                                  // total number of entries in this line
      splitline2int(txtline_strand,w,start,np,twice_max_occur);     // separate the text line txtline into integers stored in v
      n_strand=np[0];                                               // total number of entries in this line
      //printf("lines_read=%d\tn_pos=%d\tn_strand=%d\n",lines_read,n_pos,n_strand);       // debug here
      if (n_pos != n_strand) {
	printf("crm_lib.c\nload_and_sort_scan_data_v2\nerror! n_strand=%d\nn_pos=%d\npos=%s\nstrand=%s\n",n_strand,n_pos,txtline_pos,txtline_strand);
	exit(1);
      }
      if (n_pos>(twice_max_occur)) {
	printf("ERROR!!!\nfilename=%s\nlinenr=%d\nn=%d,max_occur=%d\nexceeded the limit of number of occurrences, please increase max_occur\n",filename_positions,linenr,n_pos,max_occur);
	n_pos=twice_max_occur;
      }
      if (linenr>n_motifs) {
	printf("ERROR! linenr=%d should be <= n_motifs=%d\n");
	exit(1);
      }
      motif_lengths[linenr]=v[1];      /* debug here       printf("linenr=%d\tmotif length=%d\tn=%d\n",linenr,v[1],n); */
      ll=0;                            // running number of occurrences for current motif
      if (n_pos>0) {                   // if there are any occurrences for this motif in any transcript	
	data_t[linenr][0]=n_pos/2-1;   // store the number of entries in position 0
	data_p[linenr][0]=n_pos/2-1;
	data_s[linenr][0]=n_pos/2-1;
	previous_transcript=-1;        // initial previous transcript
	n_occurrences=0;               // number of occurrences in current transcript
	n_filtered_occurrences=0;      // number of occurrences after filtering and removing overlaps for current motif
	for (j=3;j<n_pos;j+=2) {       // for each occurrence (starts at 3 because v[1] is the motif length and v[2] the number of sequences  
	  curr_transcript=v[j];        // current transcript
	  curr_transcript_strand=w[j];
	  if (curr_transcript != curr_transcript_strand) {
	    printf("linenr=%d\tcurr_transcript=%d\tcurr_transcript_strand=%d\n",linenr,curr_transcript,curr_transcript_strand);
	    exit(1);
	  }
	  k=j+1;                       // index for position
	  if (neg2pos==1) {            // convert negative entries to positive ones if and only if requested to do so
	    if (v[k]<0)
	      v[k]=-v[k];
	  }	
	  if (curr_transcript == previous_transcript) {
	    n_occurrences++;
	    temp_positions[n_occurrences]=v[k];
	    temp_strands[n_occurrences]=w[k];
	  } else {
	    if (previous_transcript<0) {	         /* this was the first entry, just add it */
	      n_occurrences++;
	      temp_positions[n_occurrences]=v[k];
	      temp_strands[n_occurrences]=w[k];
	    } else {                    	        /* ready with this transcript */
	      if (n_occurrences>0) {                    /* there were any occurrences in the current transcript */
		cumulative_nofilter+=n_occurrences;
		if (filteron==1) {
		  n_after_filter=filter_positions_v3(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
		} else {
		  n_after_filter=n_occurrences;
		  for (i=1;i<=n_after_filter;i++) {
		    filtered_positions[i]=temp_positions[i];
		    filtered_strands[i]=temp_strands[i];
		  }
		}
		if (n_after_filter>0) {
		  if (sortpos)      {                    /* if we need to sort the positions */
		    for (i=1;i<=n_after_filter;i++)
		      indices[i]=i;
		    quicksorti(n_after_filter,filtered_positions,indices);
		    sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
		  }
		  cumulative_occurrences+=n_after_filter;
		  max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
		  //n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap);
		  //n_no_overlap=vector_overlap_v2(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap,filtered_strands,strands_no_overlap);
		  n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_motif_overlap);
		  cumulative_no_overlap+=n_no_overlap;
		  for (m=1;m<=n_no_overlap;m++) {
		    ll++;
		    data_t[linenr][ll]=previous_transcript;
		    data_p[linenr][ll]=positions_no_overlap[m];
		    data_s[linenr][ll]=strands_no_overlap[m];
		  }
		  n_filtered_occurrences+=n_no_overlap;
		}                                        // close check on n_after_filter>0
	      }                                          // close check on n_occurrences>0
	      n_occurrences=1;                           // reset the number of occurrences to 1
	      temp_positions[n_occurrences]=v[k];        // and write the first position for the new transcript
	      temp_strands[n_occurrences]=w[k];
	    }                                            // close check on previous_transcript<0
	    previous_transcript=curr_transcript;         // set the previous_transcript to the current transcript
	    curr_ul=transcript_lengths[curr_transcript][1];
	    curr_el=transcript_lengths[curr_transcript][2];
	    curr_il=transcript_lengths[curr_transcript][3];
	  }                                              // close check on curr_transcript == previous_transcript
	}                                                // close j loop
     
	/* process the last transcript */
	if (n_occurrences>0) {
	  cumulative_nofilter+=n_occurrences;
	  if (filteron==1) {
	    n_after_filter=filter_positions_v3(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
	  } else {
	    n_after_filter=n_occurrences;
	    for (i=1;i<=n_after_filter;i++) {
	      filtered_positions[i]=temp_positions[i];
	      filtered_strands[i]=temp_strands[i];
	    }
	  }
	  if (n_after_filter>0) {
	    if (sortpos) {
	      for (i=1;i<=n_after_filter;i++)                   
		indices[i]=i;
	      quicksorti(n_after_filter,filtered_positions,indices);
	      quicksorti(n_after_filter,filtered_positions,indices);
	      sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
	    }
	    cumulative_occurrences+=n_after_filter;
	    max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
	    //n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap);
	    //n_no_overlap=vector_overlap_v2(filtered_positions,n_after_filter,max_motif_overlap,positions_no_overlap,filtered_strands,strands_no_overlap);
	    n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_motif_overlap);
	    cumulative_no_overlap+=n_no_overlap;
	    for (m=1;m<=n_no_overlap;m++) {
	      ll++;
	      data_t[linenr][ll]=previous_transcript;
	      data_p[linenr][ll]=positions_no_overlap[m];
	      data_s[linenr][ll]=strands_no_overlap[m];
	    }
	    n_filtered_occurrences+=n_no_overlap;
	  }      // close check on n_after_filter>0
	}        // close check on n_occurrences>0 for the last transcript
	/* update the actual number of occurrences */
	data_t[linenr][0]=n_filtered_occurrences;       // store the number of entries in position 0
	data_p[linenr][0]=n_filtered_occurrences;
	data_s[linenr][0]=n_filtered_occurrences;
      }    // if n>0 (i.e. if there was any occurrence of this motif in any transcript
    }      // close check on wether this is a header line
  }        // while there is still information in the file
  fclose(input_file_positions);
  fclose(input_file_strands);

  /* print out some information about the reading process */
  sprintf(infotext,"\nnumber of lines read = %d\n",lines_read);
  sprintf(temptext,"number of lines with overflow = %d\n",lines_overflow);strcat(infotext,temptext);
  sprintf(temptext,"number of lines with scan data = %d\n",linenr);strcat(infotext,temptext);
  if (linenr != n_motifs) {
    sprintf(temptext,"\n\nERROR! linenr=%d and n_motifs=%d\n\n",linenr,n_motifs);strcat(infotext,temptext);
    errorcode=1;
  }
  sprintf(temptext,"cumulative number of occurrences (before filter)= %ld\n",cumulative_nofilter);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences (after filter)= %ld\n",cumulative_occurrences);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences with no overlap = %ld",cumulative_no_overlap);strcat(infotext,temptext);
  //printinfo(infotext,1,log_file);
  printf("%s\n",infotext);

  free_ivector(v,1,twice_max_occur);
  free_ivector(w,1,twice_max_occur);
  free_ivector(temp_positions,1,max_occur);
  free_cvector(temp_strands,1,max_occur);
  free_vector(indices,1,max_occur);
  free_ivector(positions_no_overlap,1,max_occur);
  free_cvector(strands_no_overlap,1,max_occur);
  free_ivector(filtered_positions,1,max_occur);
  free_cvector(filtered_strands,1,max_occur);

  return errorcode;
}

int load_and_sort_scan_data_v3(char *filename_positions,char *filename_strands,int **data_t,int **data_p,int **data_s,int *motif_lengths,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths,int data_format,int filteron,int verbose)
{
  /* this function loads output data from a motif scan and sorts the positions in each transcript
   * the format for the motif scan is the following:
   * 1             2       3          4        5          6
   * motif_length  n_seqs  transcript position transcript position ...
   * one line per motif
   * all positions measured with respect to the 3' end (the larger the value, the more 5' the position)
   *
   * input
   * filename: scan file name
   * data_t: where the output transcript numbers will be stored, one row per motif, one column per occurrence (n_motifs x max_occur)
   * data_p: where the output positions numbers will be stored (n_motifs x max_occur)
   * motif_lengths: lenght of all motifs (n_motifs x 1)
   * sortpos: 1 to sort positions
   * max_occur: maximum number of occurrences
   * n_motifs: number of motifs
   * max_ul: maximum upstream length now
   * max_el: maximum exon length now
   * max_il: maximum intron length now
   * transcript_lengths: for each transcript, actual upstream, exon and intron lengths (n_transcripts x 3)
   *
   * new in version 3:
   * 01-31-2004: convert positions to TSS
   * 01-31-2004: added data_format
   * 01-31-2004: added verbose
   * 01-31-2004: removed neg2pos stuff
   * 01-31-2004: added filteron
   *
   * new in version 2:
   * 01-01-2004: load also strand information
   * return an error code if the number of motifs read does not match the actual number of motifs
   * 01-30-2004: call vector_overlap_v3
   * 01-30-2004: do not filter if max_ul<0 (allow to skip the filtering step if we do not want to filter)
   */

  FILE *input_file_positions;       // pointer to input file for positions
  FILE *input_file_strands;         // pointer to input file for strands
  char txtline_pos[1000000];        // read one line of the input file
  char txtline_strand[1000000];     // read one line of the input file
  char *gotit;                      // used in the search process
  int header;                       // 1 while reading the header
  int i,j,k,ll,m;                   // loop counters
  int *v;                           // vector holding the transcripts and positions for running entry
  int *w;                           // vector holding the transcripts and strands for running entry
  int start=3;                      // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[2];                        // number of entries
  int n_pos;                        // number of entries in current line of pos file
  int n_strand;                     // number of entries in current line of strand file
  int linenr=0;                     // current line number (each line corresponds to a separate motif)
  char lineret;                     // carriage return;
  int n_chars=1000000;              // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int twice_max_occur;              // 2 * max_occur
  int *temp_positions;              // temporary positions, (raw positions, before filtering and sorting)
  char *temp_strands;               // temporary strands, (raw strands, before filtering and sorting)
  float *indices;                   // used in call to sorting function
  int previous_transcript;          // previous transcript
  int curr_transcript;              // current transcript
  int curr_transcript_strand;       // current transcript from the strand file
  int n_occurrences;                // number of occurrences in current transcript
  int max_motif_overlap;            // only allow occurrences that overlap by less than this value
  float max_motif_overlap_factor=0.5;   // max_motif_overlap = (ceil) max_motif_overlap_factor * motif_length;
  int n_no_overlap;                 // number of entries in current transcript after removing overlaps
  int *positions_no_overlap;        // entries in current transcript after removing overlaps
  char *strands_no_overlap;         // strands in current transcript after removing overlap
  long cumulative_nofilter=0;       
  long cumulative_occurrences=0;        // cumulative sum of n_occurrences
  long cumulative_no_overlap=0;         // cumulative sum of n_no_overlap
  int lines_read=0;                     // total number of lines read in the file (including header lines)
  int lines_overflow=0;                 // total number of lines with a memory overflow
  char infotext[10000];                 // log information
  char temptext[1000];                  // log information
  int errorcode=0;                      // return errorcode = 1 if linenr does not match n_motifs at the end
  int *filtered_positions;              // filter for length, positions
  char *filtered_strands;               // filter for length, strands
  int n_after_filter;                   // number of entries in current transcript after filtering
  int curr_ul,curr_el,curr_il;          // lengths for current transcript
  int n_filtered_occurrences;           // in each transcript, number of occurrences after filtering and removing overlaps

  lineret='\n';
  
  /* memory allocation */
  twice_max_occur=2*max_occur;                    // total maximum allowed number of entries, including the transcript and position
  v=ivector(1,twice_max_occur);                   // allocate space to store all these numbers
  w=ivector(1,twice_max_occur);                   // allocate space for strands (this contains both transcripts and strands, therefore it must be an integer)
  temp_positions=ivector(1,max_occur);            // temporary positions 
  temp_strands=cvector(1,max_occur);              // temporary strands
  positions_no_overlap=ivector(1,max_occur);      // positions after removing overlaps
  strands_no_overlap=cvector(1,max_occur);        // strands after removing overlaps
  indices=vector(1,max_occur);                    // indices when sorting the entries
  filtered_positions=ivector(1,max_occur);        // positions after filtering
  filtered_strands=cvector(1,max_occur);          // strands after filtering
  
  /* open files */
  input_file_positions=fopen(filename_positions,"r");
  if (!input_file_positions) {
    printf("error! i could not open the file %s",filename_positions);
    exit(1);
  } 
  input_file_strands=fopen(filename_strands,"r");
  if (!input_file_strands) {
    printf("error! i could not open the file %s",filename_strands);
    exit(1);
  }

  while ( (feof(input_file_positions) == 0) & (feof(input_file_strands)==0) & (linenr<n_motifs) ) {
    lines_read++;
    fgets(txtline_pos,n_chars,input_file_positions);
    fgets(txtline_strand,n_chars,input_file_strands);
    header=is_header(txtline_pos);
    if (!header) {      /* this is not a header line, ignore it */
      linenr++;         /* debug here      printf("-------------------  linenr=%d  -------------------------- \n",linenr); */
      gotit=strchr(txtline_pos,lineret);
      if (!gotit) {
	printf("ERROR!!!\tfilename_positions=%s\nlines_read=%d\tlinenr=%d\nthere has been an overflow,strlen(txtline)=%d\n",filename_positions,lines_read,linenr,strlen(txtline_pos));
	lines_overflow++;
      }
      splitline2int(txtline_pos,v,start,np,twice_max_occur);        // separate the text line txtline into integers stored in v
      n_pos=np[0];                                                  // total number of entries in this line
      splitline2int(txtline_strand,w,start,np,twice_max_occur);     // separate the text line txtline into integers stored in v
      n_strand=np[0];                                               // total number of entries in this line
      if (verbose)
	printf("lines_read=%d\tn_pos=%d\tn_strand=%d\n",lines_read,n_pos,n_strand);       
      if (n_pos != n_strand) {
	printf("crm_lib.c\nload_and_sort_scan_data_v2\nerror! n_strand=%d\nn_pos=%d\npos=%s\nstrand=%s\n",n_strand,n_pos,txtline_pos,txtline_strand);
	exit(1);
      }
      if (n_pos>(twice_max_occur)) {
	printf("ERROR!!!\nfilename=%s\nlinenr=%d\nn=%d,max_occur=%d\nexceeded the limit of number of occurrences, please increase max_occur\n",filename_positions,linenr,n_pos,max_occur);
	n_pos=twice_max_occur;
      }
      if (linenr>n_motifs) {
	printf("ERROR! linenr=%d should be <= n_motifs=%d\n");
	exit(1);
      }
      motif_lengths[linenr]=v[1];      /* debug here       printf("linenr=%d\tmotif length=%d\tn=%d\n",linenr,v[1],n); */
      ll=0;                            // running number of occurrences for current motif
      if (n_pos>0) {                   // if there are any occurrences for this motif in any transcript	
	data_t[linenr][0]=n_pos/2-1;   // store the number of entries in position 0
	data_p[linenr][0]=n_pos/2-1;
	data_s[linenr][0]=n_pos/2-1;
	previous_transcript=-1;        // initial previous transcript
	n_occurrences=0;               // number of occurrences in current transcript
	n_filtered_occurrences=0;      // number of occurrences after filtering and removing overlaps for current motif
	for (j=3;j<n_pos;j+=2) {       // for each occurrence (starts at 3 because v[1] is the motif length and v[2] the number of sequences  
	  curr_transcript=v[j];        // current transcript
	  curr_transcript_strand=w[j];
	  if (curr_transcript != curr_transcript_strand) {
	    printf("linenr=%d\tcurr_transcript=%d\tcurr_transcript_strand=%d\n",linenr,curr_transcript,curr_transcript_strand);
	    exit(1);
	  }
	  k=j+1;                       // index for position
	  if (curr_transcript == previous_transcript) {
	    n_occurrences++;
	    temp_positions[n_occurrences]=v[k];
	    temp_strands[n_occurrences]=w[k];
	  } else {
	    if (previous_transcript<0) {	         /* this was the first entry, just add it */
	      n_occurrences++;
	      temp_positions[n_occurrences]=v[k];
	      temp_strands[n_occurrences]=w[k];
	    } else {                    	        /* ready with this transcript */
	      if (n_occurrences>0) {                    /* there were any occurrences in the current transcript */
		cumulative_nofilter+=n_occurrences;
		if (verbose) {
		  printf("transcript=%d\n",previous_transcript);
		  printf("\ttemp_positions (%d): ",n_occurrences);
		  for (i=1;i<=n_occurrences;i++)
		    printf("%d %d,",temp_positions[i],temp_strands[i]);
		  printf("\n");
		}

		/* format conversion */
		if (data_format==2) {
		  /* positions with respect to 3' end */
		  convert_positions2tss(temp_positions,n_occurrences,curr_el,curr_il);
		  if (verbose) {
		    printf("\ttemp_positions after conversion (%d): ",n_occurrences);
		    for (i=1;i<=n_occurrences;i++)
		      printf("%d %d,",temp_positions[i],temp_strands[i]);
		    printf("\n");
		  }
		}
		if (filteron==1) {
		  //n_after_filter=filter_positions_v3(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
		  n_after_filter=filter_positions_v5(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_el,filtered_positions,filtered_strands);
		} else {
		  n_after_filter=n_occurrences;
		  for (i=1;i<=n_after_filter;i++) {
		    filtered_positions[i]=temp_positions[i];
		    filtered_strands[i]=temp_strands[i];
		  }
		}
		if (verbose) {
		  printf("\tafter filtering (%d): ",n_after_filter);
		  for (i=1;i<=n_after_filter;i++)
		    printf("%d %d,",filtered_positions[i],filtered_strands[i]);
		  printf("\n");
		}
		if (n_after_filter>0) {
		  if (sortpos)      {                    /* if we need to sort the positions */
		    for (i=1;i<=n_after_filter;i++)
		      indices[i]=i;
		    quicksorti(n_after_filter,filtered_positions,indices);
		    sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
		    if (verbose) {
		      printf("\tafter sorting (%d): ",n_after_filter);
		      for (i=1;i<=n_after_filter;i++)
			printf("%d %d,",filtered_positions[i],filtered_strands[i]);
		      printf("\n");
		    }
		  }
		  cumulative_occurrences+=n_after_filter;
		  max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
		  n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_motif_overlap);
		  cumulative_no_overlap+=n_no_overlap;
		  if (verbose) {
		    printf("\tno overlap (%d): ",n_no_overlap);
		    for (i=1;i<=n_no_overlap;i++)
		      printf("%d %d,",positions_no_overlap[i],strands_no_overlap[i]);
		    printf("\n");
		  }
		  for (m=1;m<=n_no_overlap;m++) {
		    ll++;
		    data_t[linenr][ll]=previous_transcript;
		    data_p[linenr][ll]=positions_no_overlap[m];
		    data_s[linenr][ll]=strands_no_overlap[m];
		  }
		  n_filtered_occurrences+=n_no_overlap;
		}                                        // close check on n_after_filter>0
	      }                                          // close check on n_occurrences>0
	      n_occurrences=1;                           // reset the number of occurrences to 1
	      temp_positions[n_occurrences]=v[k];        // and write the first position for the new transcript
	      temp_strands[n_occurrences]=w[k];
	    }                                            // close check on previous_transcript<0
	    previous_transcript=curr_transcript;         // set the previous_transcript to the current transcript
	    curr_ul=transcript_lengths[curr_transcript][1];
	    curr_el=transcript_lengths[curr_transcript][2];
	    curr_il=transcript_lengths[curr_transcript][3];
	  }                                              // close check on curr_transcript == previous_transcript
	}                                                // close j loop
     
	/* process the last transcript */
	if (n_occurrences>0) {
	  cumulative_nofilter+=n_occurrences;
	  if (verbose) {
	    printf("transcript=%d\n",previous_transcript);
	    printf("\ttemp_positions (%d): ",n_occurrences);
	    for (i=1;i<=n_occurrences;i++)
	      printf("%d %d,",temp_positions[i],temp_strands[i]);
	    printf("\n");
	  }
	  if (data_format==2) {	  /* format conversion for positions with respect to 3' end to TSS coordinates */
	    convert_positions2tss(temp_positions,n_occurrences,curr_el,curr_il);
	    if (verbose) {
	      printf("\ttemp_positions after conversion (%d): ",n_occurrences);
	      for (i=1;i<=n_occurrences;i++)
		printf("%d %d,",temp_positions[i],temp_strands[i]);
	      printf("\n");
	    }
	  }
	  
	  if (filteron==1) {
	    //n_after_filter=filter_positions_v3(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
	    n_after_filter=filter_positions_v5(temp_positions,temp_strands,n_occurrences,max_ul,max_el,max_il,curr_el,filtered_positions,filtered_strands);
	  } else {
	    n_after_filter=n_occurrences;
	    for (i=1;i<=n_after_filter;i++) {
	      filtered_positions[i]=temp_positions[i];
	      filtered_strands[i]=temp_strands[i];
	    }
	  }
	  if (n_after_filter>0) {
	    if (sortpos) {
	      for (i=1;i<=n_after_filter;i++)                   
		indices[i]=i;
	      quicksorti(n_after_filter,filtered_positions,indices);
	      quicksorti(n_after_filter,filtered_positions,indices);
	      sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
	    }
	    cumulative_occurrences+=n_after_filter;
	    max_motif_overlap=(int)floor(max_motif_overlap_factor*v[1]); // maximum allowed overlap between motifs
	    n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_motif_overlap);
	    cumulative_no_overlap+=n_no_overlap;
	    for (m=1;m<=n_no_overlap;m++) {
	      ll++;
	      data_t[linenr][ll]=previous_transcript;
	      data_p[linenr][ll]=positions_no_overlap[m];
	      data_s[linenr][ll]=strands_no_overlap[m];
	    }
	    n_filtered_occurrences+=n_no_overlap;
	  }      // close check on n_after_filter>0
	}        // close check on n_occurrences>0 for the last transcript
	/* update the actual number of occurrences */
	data_t[linenr][0]=n_filtered_occurrences;       // store the number of entries in position 0
	data_p[linenr][0]=n_filtered_occurrences;
	data_s[linenr][0]=n_filtered_occurrences;
      }    // if n>0 (i.e. if there was any occurrence of this motif in any transcript
    }      // close check on wether this is a header line
  }        // while there is still information in the file
  fclose(input_file_positions);
  fclose(input_file_strands);

  /* print out some information about the reading process */
  sprintf(infotext,"\nnumber of lines read = %d\n",lines_read);
  sprintf(temptext,"number of lines with overflow = %d\n",lines_overflow);strcat(infotext,temptext);
  sprintf(temptext,"number of lines with scan data = %d\n",linenr);strcat(infotext,temptext);
  if (linenr != n_motifs) {
    sprintf(temptext,"\n\nERROR! linenr=%d and n_motifs=%d\n\n",linenr,n_motifs);strcat(infotext,temptext);
    errorcode=1;
  }
  sprintf(temptext,"cumulative number of occurrences (before filter)= %ld\n",cumulative_nofilter);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences (after filter)= %ld\n",cumulative_occurrences);strcat(infotext,temptext);
  sprintf(temptext,"cumulative number of occurrences with no overlap = %ld",cumulative_no_overlap);strcat(infotext,temptext);
  printf("%s\n",infotext);

  free_ivector(v,1,twice_max_occur);
  free_ivector(w,1,twice_max_occur);
  free_ivector(temp_positions,1,max_occur);
  free_cvector(temp_strands,1,max_occur);
  free_vector(indices,1,max_occur);
  free_ivector(positions_no_overlap,1,max_occur);
  free_cvector(strands_no_overlap,1,max_occur);
  free_ivector(filtered_positions,1,max_occur);
  free_cvector(filtered_strands,1,max_occur);

  return errorcode;
}

void crm_faster_1motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output)
{
  /* search for windows with at least N occurrences of motif m1 in all transcripts 
   * faster version using the stored counts of each motif
   *
   * input
   * crm_positions = [1] = transcript [2] = positions of the putative crms (bin), [2] = number of motifs in each position
   * max_crm_positions = maximum number of rows allocated in crm_positions
   * m1 = index for motif 1
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * window_bins = number of bins in a window (used to make a larger step when there was a hit to avoid overlapping counts)
   * n_transcripts = number of transcripts
   * max_n_transcripts = stop computation if there are more than max_n_transcripts with module (set to -1 to compute for all transcritps)
   * occurrences_per_transcript = total number of occurrences per transcript
   * report_all = set to 1 to report positions of all modules found, set to 0 to report only the number
   * 
   * output
   * crm_output[0] = number of putative crms, crm_output[1] = number of transcripts with putative crms, crm_output[2] = errorcode (overflow of crm_positions)
   *
   * 10-20-2003: added report_all
   */

  int errorcode=0;
  int transcript;                    // running transcript
  int n_tr=0;                        // number of transcripts with hits
  int n_crm_init;                    // n_crm before counting in current transcript (used to determine whether there was a crm in the current transcript)
  int bin_init;
  int bin_finit;
  int bin;                           // current bin
  int n_crm=0;                       // number of crm
  int total_n1;                      // total number of occurrences of the motif in the current transcript
  int n_cooc;

  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts;

  transcript=0;
  while ( (transcript<n_transcripts) & (n_tr<=max_n_transcripts) ) {
    transcript++;
    total_n1=occurrences_per_transcript[m1][transcript];

    if (total_n1>=N) {
      bin_init=init_bin_boundary[m1][transcript];    /* get initial bin */
      bin_finit=finit_bin_boundary[m1][transcript];  /* get final bin   */
    
      bin=bin_init;
      n_crm_init=n_crm;
      if (report_all == 1) {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    if (n_crm > max_crm_positions) {
	      errorcode++;
	    } else {
	      crm_positions[n_crm][1]=transcript;
	      crm_positions[n_crm][2]=bin;
	      crm_positions[n_crm][3]=n_cooc;
	    }
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      } else {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      }

      if (n_crm>n_crm_init)
	n_tr++;
    }       // close check on total_n1>N
  }         // close transcript loop

  crm_output[0]=n_crm;
  crm_output[1]=n_tr;
  crm_output[2]=errorcode;
}

void crm_faster_2motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output)
{
  /* overlap of four motifs in all transcripts 
   * faster version using the stored counts of each motif
   *
   * input
   * crm_positions = [1] = transcript [2] = positions of the putative crms (bin), [2] = number of motifs in each position
   * max_crm_positions = maximum number of rows allocated in crm_positions
   * m1,m2 = indices for motifs 1 through 4
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * window_bins = number of bins in a window (used to make a larger step when there was a hit to avoid overlapping counts)
   * n_transcripts = number of transcripts
   * 
   * output
   * crm_output[0] = number of putative crms, crm_output[1] = number of transcripts with putative crms, crm_output[2] = errorcode (overflow of crm_positions)
   *
   * 10-20-2003: added report_all
   */

  int errorcode=0;
  int transcript;                    // running transcript
  int n_tr=0;                        // number of transcripts with hits
  int n_crm_init;                    // n_crm before counting in current transcript (used to determine whether there was a crm in the current transcript)
  int bin_init1,bin_init2,bin_init3,bin_init4,bin_init;
  int bin_finit1,bin_finit2,bin_finit3,bin_finit4,bin_finit;
  int bin;                           // current bin
  int n_crm=0;                       // number of crm
  int n1,n2;                         // count of occurrences in a specific bin for a specific transcript for each motif
  int n_cooc;
  int total_n1,total_n2;             // total number of occurrences of motifs m1 and m2 in the current transcript

  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts;

  transcript=0;
  while ( (transcript<n_transcripts) & (n_tr<=max_n_transcripts) ) {
    transcript++;
    total_n1=occurrences_per_transcript[m1][transcript];
    total_n2=occurrences_per_transcript[m2][transcript];

    if ( (total_n1+total_n2)>=N) {

      /* get initial bin */
      bin_init1=init_bin_boundary[m1][transcript];
      bin_init2=init_bin_boundary[m2][transcript];
      bin_init=min2(bin_init1,bin_init2);
      /* get final bin */
      bin_finit1=finit_bin_boundary[m1][transcript];
      bin_finit2=finit_bin_boundary[m2][transcript];
      bin_finit=max2(bin_finit1,bin_finit2);
    
      bin=bin_init;
      n_crm_init=n_crm;

      if (report_all == 1) {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    if (n_crm > max_crm_positions) {
	      errorcode++;
	    } else {
	      crm_positions[n_crm][1]=transcript;
	      crm_positions[n_crm][2]=bin;
	      crm_positions[n_crm][3]=n_cooc;
	    }
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      } else {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      }
      if (n_crm>n_crm_init)
	n_tr++;
    }       // close check on (total_n1+total_n2)>=N
  }         // close transcript loop

  crm_output[0]=n_crm;
  crm_output[1]=n_tr;
  crm_output[2]=errorcode;
}

void crm_faster_3motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int m3,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output)
{
  /* search for windows with at least N occurrences of any of the three motifs m1,m2,m3
   * faster version using the stored counts of each motif
   *
   * input
   * crm_positions = [1] = transcript [2] = positions of the putative crms (bin), [2] = number of motifs in each position
   * max_crm_positions = maximum number of rows allocated in crm_positions
   * m1,m2,m3 = indices for motifs 1 through 3
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * window_bins = number of bins in a window (used to make a larger step when there was a hit to avoid overlapping counts)
   * n_transcripts = number of transcripts
   * 
   * output
   * crm_output[0] = number of putative crms, crm_output[1] = number of transcripts with putative crms, crm_output[2] = errorcode (overflow of crm_positions)
   *
   * 10-20-2003: added report_all
   */

  int errorcode=0;
  int transcript;                    // running transcript
  int n_tr=0;                        // number of transcripts with hits
  int n_crm_init;                    // n_crm before counting in current transcript (used to determine whether there was a crm in the current transcript)
  int bin_init1,bin_init2,bin_init3,bin_init;
  int bin_finit1,bin_finit2,bin_finit3,bin_finit;
  int bin;                           // current bin
  int n_crm=0;                       // number of crm
  int n1,n2,n3;                      // count of occurrences in a specific bin for a specific transcript for each motif
  int n_cooc;                        // number of co-occurrences
  int total_n1,total_n2,total_n3;    // total number of occurrences of motifs m1,m2,m3 in the current transcript
  int function_verbose=0;            // to output performance time plus parameters
  struct timeb ti,tf;                // used to compute performance time
  long tdif;                         // used to compute performance time

  if (function_verbose==1) 
    ftime(&ti);

  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts;

  transcript=0;
  while ( (transcript<n_transcripts) & (n_tr<=max_n_transcripts) ) {
    transcript++;
    total_n1=occurrences_per_transcript[m1][transcript];
    total_n2=occurrences_per_transcript[m2][transcript];
    total_n3=occurrences_per_transcript[m3][transcript];

    if ( (total_n1+total_n2+total_n3) >= N) {

      /* get initial bin */
      bin_init1=init_bin_boundary[m1][transcript];
      bin_init2=init_bin_boundary[m2][transcript];
      bin_init3=init_bin_boundary[m3][transcript];
      bin_init=min3(bin_init1,bin_init2,bin_init3);
      /* get final bin */
      bin_finit1=finit_bin_boundary[m1][transcript];
      bin_finit2=finit_bin_boundary[m2][transcript];
      bin_finit3=finit_bin_boundary[m3][transcript];
      bin_finit=max3(bin_finit1,bin_finit2,bin_finit3);
      
      bin=bin_init;
      n_crm_init=n_crm;

      if (report_all == 1) {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript]+motif_counts[m3][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    if (n_crm > max_crm_positions) {
	      errorcode++;
	    } else {
	      crm_positions[n_crm][1]=transcript;
	      crm_positions[n_crm][2]=bin;
	      crm_positions[n_crm][3]=n_cooc;
	    }
	    bin+=window_bins;
	  } else {
	    bin++;
	  }    // close check on n_cooc>=N
	}        // close bin loop
      } else {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript]+motif_counts[m3][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    bin+=window_bins;
	  } else {
	    bin++;
	  }    // close check on n_cooc>=N
	}        // close bin loop
      }
      if (n_crm>n_crm_init)
	n_tr++;
    }      // close check on (total_n1+total_n2+total_n3)>=N
  }       // close transcript loop

  crm_output[0]=n_crm;
  crm_output[1]=n_tr;
  crm_output[2]=errorcode;

  if (function_verbose==1) {
    ftime(&tf);
    tdif=time_difference(ti,tf);
    printf("crm_faster_3motifs_allt\t%ld msec cpu time\n",tdif);
  }

}

void crm_faster_4motifs_allt(char ***motif_counts,int **init_bin_boundary,int **finit_bin_boundary,int **crm_positions,int max_crm_positions,int m1,int m2,int m3,int m4,int N,int window_bins,int n_transcripts,int max_n_transcripts,int **occurrences_per_transcript,int report_all,int *crm_output)
{
  /* overlap of four motifs in all transcripts 
   * faster version using the stored counts of each motif
   *
   * input
   * crm_positions = [1] = transcript [2] = positions of the putative crms (bin), [2] = number of motifs in each position
   * max_crm_positions = maximum number of rows allocated in crm_positions
   * m1,m2,m3,m4 = indices for motifs 1 through 4
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * window_bins = number of bins in a window (used to make a larger step when there was a hit to avoid overlapping counts)
   * n_transcripts = number of transcripts
   * 
   * output
   * crm_output[0] = number of putative crms, crm_output[1] = number of transcripts with putative crms, crm_output[2] = errorcode (overflow of crm_positions)
   *
   * 10-20-2003: added report_all
   */

  int errorcode=0;
  int transcript;                    // running transcript
  int n_tr=0;                        // number of transcripts with hits
  int n_crm_init;                    // n_crm before counting in current transcript (used to determine whether there was a crm in the current transcript)
  int bin_init1,bin_init2,bin_init3,bin_init4,bin_init;
  int bin_finit1,bin_finit2,bin_finit3,bin_finit4,bin_finit;
  int bin;                           // current bin
  int n_crm=0;                       // number of crm
  int n1,n2,n3,n4;                   // count of occurrences in a specific bin for a specific transcript for each motif
  int n_cooc;
  int total_n1,total_n2,total_n3,total_n4;

  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts;

  transcript=0;
  while ( (transcript<n_transcripts) & (n_tr<=max_n_transcripts) ) {
    transcript++;
    total_n1=occurrences_per_transcript[m1][transcript];
    total_n2=occurrences_per_transcript[m2][transcript];
    total_n3=occurrences_per_transcript[m3][transcript];
    total_n4=occurrences_per_transcript[m4][transcript];

    if ( (total_n1+total_n2+total_n3+total_n4)>=N ) {
      /* get initial bin */
      bin_init1=init_bin_boundary[m1][transcript];
      bin_init2=init_bin_boundary[m2][transcript];
      bin_init3=init_bin_boundary[m3][transcript];
      bin_init4=init_bin_boundary[m4][transcript];
      bin_init=min4(bin_init1,bin_init2,bin_init3,bin_init4);
      /* get final bin */
      bin_finit1=finit_bin_boundary[m1][transcript];
      bin_finit2=finit_bin_boundary[m2][transcript];
      bin_finit3=finit_bin_boundary[m3][transcript];
      bin_finit4=finit_bin_boundary[m4][transcript];
      bin_finit=max4(bin_finit1,bin_finit2,bin_finit3,bin_finit4);
      
      bin=bin_init;
      n_crm_init=n_crm;

      if (report_all == 1) {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript]+motif_counts[m3][bin][transcript]+motif_counts[m4][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    if (n_crm > max_crm_positions) {
	      errorcode++;
	    } else {
	      crm_positions[n_crm][1]=transcript;
	      crm_positions[n_crm][2]=bin;
	      crm_positions[n_crm][3]=n_cooc;
	    }
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      } else {
	while (bin<=bin_finit) {
	  n_cooc=motif_counts[m1][bin][transcript]+motif_counts[m2][bin][transcript]+motif_counts[m3][bin][transcript]+motif_counts[m4][bin][transcript];
	  if (n_cooc>=N) {
	    n_crm++;
	    bin+=window_bins;
	  } else {
	    bin++;
	  }   // close check on n_cooc>=N
	}     // close bin loop
      }

      if (n_crm>n_crm_init)
	n_tr++;
    }     // close check on (total_n1+total_n2+total_n3+total_n4)>=N
  }       // close transcript loop

  crm_output[0]=n_crm;
  crm_output[1]=n_tr;
  crm_output[2]=errorcode;
}

void count_motif(int *m_p,int n_occurrences,int init_index,int lower_border,int upper_border,int *out)
{
  /* count the number of occurrences of the motif within lower_order and upper_border
   * m_p = positions of the motif
   * n_occurrences = number of occurrences of the motif
   * init_index = initial occurrence of the motif to start searching
   * lower_border = lower end of the window
   * upper_border = upper end of the window
   * out = output, out[0]=number of occurrences, out[1]=init_index
   */

  int i;                             // running index for the motif
  int p;                             // current position of the motif
  int n=0;                           // running number of occurrences of the motif

  //printf("\t\tn_occurrences=%d\tinit_index=%d\tlower_border=%d\tupper_border=%d\n",n_occurrences,init_index,lower_border,upper_border);

  i=init_index;
  p=m_p[i]         ;                 // start from init_index
  while (p<upper_border) {           // stop if we are past the border
    if (p<lower_border) {           
      init_index++;                  // if the position still has not reached the border, then next time start with a higher position
    } else {
      n++;                           // add one to the count
    }
    //printf("\t\t\ti=%d\tp=%d\tn=%d\n",i,p,n);
    i++;                             // next index
    if (i <= n_occurrences) {
      p=m_p[i];                        // position for the next index
    } else {
      p=upper_border+1;
    }
  }

  out[0]=n;
  out[1]=init_index;
}

void crm_1motifs(int *m1_p,int n1,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions)
{
  /* search for windows of size dw with at least N occurrences of the motif
   * 
   * input
   *
   * m1_p = positions of motif 1
   * n1 = number of occurrences of motif 1
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * dw = window length (nucleotides)
   * max_crm_positions = maximum number of crm positions (to avoid overflows)
   * bin_size = step size for scanning through the transcript
   * report_all = 1 to report all positions
   *
   * output = function output, output[0] = number of CRMS, output[1] = errorcode
   * errorcode = number of times where the position of a CRM could not be recorded due to overflow of max_crm_positions
   * crm_positions = positions of the putative crms plus number of occurrences in window
   *
   * this function performs the search in a single transcript
   *
   * 10-18-2003: fixed bug, for loop in previous version for pos should be a while loop
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the input
   */

  int min_pos;                                             // minimum position      
  int max_pos;                                             // maximum position
  int n_crm=0;                                             // number of putative CRMs
  int pos;                                                 // current position 
  int upper_border;                                        // window upper border
  int init_index1;                                         // initial index for motif 1 in call to count_motif
  int out[2];                                              // output of count_motif
  int n_cooc;                                              // number of co-occurrences within the current window
  int errorcode=0;                                         // error code
  int i;

  if (n1>=N) {
    min_pos=m1_p[1];                                         /* determine minimum position to start scanning */
    max_pos=m1_p[n1];                                        /* determine minimum position to stop scanning */
    max_pos=max_pos-dw+bin_size;                                      // because the last window extends from this point for dw nucleotides
    pos =min_pos;

    if (report_all==1) {
      while (pos<=max_pos) {
	upper_border=pos+dw;
	n_cooc=0;
	if (init_index1<=n1) {
	  count_motif(m1_p,n1,init_index1,pos,upper_border,out);
	  n_cooc=out[0];
	  init_index1=out[1];
	}
	if (n_cooc>=N) {
	  n_crm++;
	  if (n_crm<max_crm_positions) {
	    crm_positions[n_crm][1]=pos;
	    crm_positions[n_crm][2]=n_cooc;
	  } else {
	    errorcode++;
	  }
	  pos+=dw;
	} else {
	  pos+=bin_size;
	} // close check on n_cooc>=N
      }   // close while loop on pos
    } else {
      while (pos<=max_pos) {
	upper_border=pos+dw;
	n_cooc=0;
	if (init_index1<=n1) {
	  count_motif(m1_p,n1,init_index1,pos,upper_border,out);
	  n_cooc=out[0];
	  init_index1=out[1];
	}
	if (n_cooc>=N) {
	  n_crm++;
	  pos+=dw;
	} else {
	  pos+=bin_size;
	} // close check on n_cooc>=N
      }   // close while loop on pos
    }     // close check on report_all == 1
  }       // close check on n1>=N
  
  output[0]=n_crm;
  output[1]=errorcode;
}

void crm_1motifs_allt(int *m1_t,int *m1_p,int n1,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions)
{
  /* call crm_1motif
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the input
   */

  int **crm_positions_transcript;                 // crm positions in current transcript
  int *m1p;
  int curr_n1;                           // number of occurrences of motif in current transcript
  int curr_tr1;
  int ec=0;                                      // error code in current function
  int errorcode;
  int i,j;
  int i1;
  int max_crm_positions_transcript=1000;         // maximum number of crm positions in a transcript (memory = 1000 x 1 (col) x (2 bytes/int) = 2 kb)
  int max_transcript_occurrences=1000;           // maximum number of occurrences of one motif per transcript (memory = 1000x1(cols)x2(bytes/int)=2kb)
  int n1_tr=0;
  int n_crms;                                    // number of crms (output of crm_1motif
  int n_crms_reported;                           // number of reported crms (output of crm_1motif; this is the same as n_crms unless there was a memory overflow
  int output[2];                                 // output from the crm_2motifs program
  int running_entries=0;                         // running number of entries in ovdi
  int running_entries_reported;                  // number of reported crms; while n_crms and n_crms_reported refer to a single transcript, this refers to the hole function
  int transcript=0;

  i1=1;curr_tr1=m1_t[i1];                         // first transcript, motif 1
  
  m1p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 1 in current transcript
  if (report_all==1) 
    crm_positions_transcript=imatrix(1,max_crm_positions_transcript,1,2);
  
  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;           // n12_tr will always be less than this

  while ( (transcript < n_transcripts) & (n1_tr<=max_n_transcripts) ) {
    transcript++;
    if (curr_tr1<0) 
      transcript=n_transcripts+1;       // done with the mtoif, let's call it a day
    
    while ( (curr_tr1<transcript) & (curr_tr1>=0) & (i1<n1) ) {
      i1++;curr_tr1=m1_t[i1];
    }
    curr_n1=0;
    if (curr_tr1==transcript) {
      while ( (curr_tr1==transcript) & (i1<=n1) ) {
	curr_n1++;
	m1p[curr_n1]=m1_p[i1];
	i1++;
	if (i1<=n1) 
	  curr_tr1=m1_t[i1];
      }
      i1--;                       // set to the last entry corresponding to current transcript
      curr_tr1=transcript;        // set curr_tr1 to current transcript
    }

    if (curr_n1>=N) {
      crm_1motifs(m1p,curr_n1,N,dw,max_crm_positions_transcript,bin_size,report_all,output,crm_positions_transcript);
      n_crms=output[0];
      errorcode=output[1];
      n_crms_reported=n_crms-errorcode;
    } else {
      n_crms=0;
      n_crms_reported=0;
    }
    if (n_crms>0)
      n1_tr++;
      
    if (report_all==1) {
      if ( (running_entries+n_crms) > max_crm_positions ) {
	ec--;
	running_entries_reported=running_entries;
	printf("running_entries=%d\n",running_entries);
	printf("running_entries_reported=%d\n",running_entries_reported);
      } else {
	for (i=1;i<=n_crms_reported;i++) {
	  running_entries++;
	  crm_positions[running_entries][1]=transcript;
	  crm_positions[running_entries][2]=crm_positions_transcript[i][1];
	  crm_positions[running_entries][3]=crm_positions_transcript[i][2];
	}
      }   // close check on whether running_entries will be within max_crm_positions
    }     // close if (report_all==1)
  }       // while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n1_tr;
  n[2]=ec;
  if (ec<0) {
    n[3]=running_entries_reported;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_ivector(m1p,1,max_transcript_occurrences);     
  if (report_all==1)
    free_imatrix(crm_positions_transcript,1,max_crm_positions_transcript,1,2);
}

void crm_2motifs(int *m1_p,int *m2_p,int n1,int n2,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions)
{
  /* search for co-occurrences of at least N motifs in windows of size dw 
   * 
   * input
   *
   * m1_p = positions of motif 1
   * m2_p = positions of motif 2
   * n1 = number of occurrences of motif 1
   * n2 = number of occurrences of motif 2
   * crm_positions = [1] = positions of the putative crms, [2] = number of motifs in each position
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * dw = window length (nucleotides)
   * max_crm_positions = maximum number of crm positions (to avoid overflows)
   * bin_size = step size for scanning through the transcript
   * report_all = if 1, then return all positions in crm_positions, otherwise, just the number
   *
   * output = function output, output[0] = number of CRMS, output[1] = errorcode
   * errorcode = number of times where the position of a CRM could not be recorded due to overflow of max_crm_positions
   *
   * this function performs the search in a single transcript
   *
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the input so that input is first and output goes last
   */

  int min_pos=20000;                           // minimum position of motif 1 and motif 2      
  int max_pos=-20000;                           // maximum position of motif 1 and motif 2
  int n_crm=0;                                             // number of putative CRMs
  int pos;                                                 // current position 
  int upper_border;                                        // window upper border
  int init_index1=1;                                         // initial index for motif 1 in call to count_motif
  int init_index2=1;                                         // initial index for motif 2 in call to count_motif
  int out[2];                                              // output of count_motif
  int n_cooc;                                              // number of co-occurrences within the current window
  int errorcode=0;                                         // error code
  int i;

 /* determine minimum and position to start scanning */
  if (n1>0) {
    if (m1_p[1]<min_pos)
      min_pos=m1_p[1];
    if (m1_p[n1]>max_pos)
      max_pos=m1_p[n1];
  }
  if (n2>0) {
    if (m2_p[1]<min_pos)
      min_pos=m2_p[1];
    if (m2_p[n2]>max_pos)
      max_pos=m2_p[n2];
  }
  max_pos=max_pos-dw+bin_size;                                 // because the last window extends from this point for dw nucleotides
  pos=min_pos;

  /* debug here 
     printf("min_pos=%d\tmax_pos=%d\n",min_pos,max_pos); */
  
  if (report_all==1) {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;

      /* printf("%d %d %d %d %d %d ",n1,n2,pos,upper_border,init_index1,init_index2); */
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);
	n_cooc+=out[0];
	init_index1=out[1];
      }
      
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);
	n_cooc+=out[0];
	init_index2=out[1];
      }
      
      /* debug here 
	 printf("%d %d\n",n_cooc,n_crm); */

      if (n_cooc>=N) {
	n_crm++;
	if (n_crm<max_crm_positions) {
	  crm_positions[n_crm][1]=pos;
	  crm_positions[n_crm][2]=n_cooc;
	} else {
	  errorcode++;
	}
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }
  } else {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;
      
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);
	n_cooc+=out[0];
	init_index1=out[1];
      }
      
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);
	n_cooc+=out[0];
	init_index2=out[1];
      }
      
      if (n_cooc>=N) {
	n_crm++;
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }
  }
  
  output[0]=n_crm;
  output[1]=errorcode;
}

void crm_2motifs_allt(int *m1_t,int *m2_t,int *m1_p,int *m2_p,int n1,int n2,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions)
{
  /* call crm_2motifs
   * 
   * inputs
   * m1_t: transcripts, motif 1
   * m2_t: transcripts, motif 2
   * m1_p: positions, motif 1
   * m2_p: positions, motif 2
   * n1: total number of occurrences, motif 1
   * n2: total number of occurrences, motif 2
   * max_n_transcripts: stop computation if the number of transcripts with crm is bigger than this number
   * max_crm_positions: to avoid overflow in crm_positions
   * N: minimum number of motif occurrences in a window of dw nucleotides
   * dw: window size (nucleotides)
   * bin_size: bin size (nucleotides)
   * n_transcripts: number of transcripts 
   * report_all: if 1, then report all positions, otherwise just the numbers
   * 
   * output 
   * n
   *        n[0]=running_entries;
   *        n[1]=n12_tr;
   *        n[2]=errorcode;
   *        n[3]=running_entries_reported
   * crm_positions: for each crm, return the transcript (col1), position (col2), number of motifs (col3)
   *
   * 10-17-2003: do not require that all motifs must be present in each module
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the inputs
   */

  int **crm_positions_transcript;                 // crm positions in current transcript
  int *m1p;                                       // motif 1 occurrences in running transcript
  int *m2p;                                       // motif 2 occurrences in running transcript
  int curr_n1,curr_n2;                            // number of occurrences of motif in current transcript
  int curr_tr1,curr_tr2;                          // running transcript for motif 1 and motif 2
  int ec=0;                                       // error code in current function
  int errorcode=0;                                // error code in call to crm_2motifs
  int i,j;
  int i1,i2;
  int max_crm_positions_transcript=1000;         // maximum number of crm positions in a transcript (memory = 1000 x 1 (col) x (2 bytes/int) = 2 kb)
  int max_transcript_occurrences=1000;           // maximum number of occurrences of one motif per transcript (memory = 1000x1(cols)x2(bytes/int)=2kb)
  int n12;
  int n12_tr=0;
  int n_crms;                                    // number of crms (output of crm_1motif
  int n_crms_reported;                           // number of reported crms (output of crm_1motif; this is the same as n_crms unless there was a memory overflow
  int output[2];                                 // output from the crm_2motifs program
  int running_entries=0;                         // running number of entries in ovdi
  int running_entries_reported;                  // number of reported crms; while n_crms and n_crms_reported refer to a single transcript, this refers to the hole function
  int transcript=0;

  i1=1;curr_tr1=m1_t[i1];                         // first transcript, motif 1
  i2=1;curr_tr2=m2_t[i2];                         // first transcrpit, motif 2
  
  m1p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 1 in current transcript
  m2p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 2 in current transcript
  if (report_all==1)
    crm_positions_transcript=imatrix(1,max_crm_positions_transcript,1,2);
  
  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;           // n12_tr will always be less than this

  while ( (transcript < n_transcripts) & (n12_tr<=max_n_transcripts) ) {
    transcript++;

    while ( (curr_tr1<transcript) & (curr_tr1>=0) & (i1<n1) ) {
      i1++;curr_tr1=m1_t[i1];
    }
    curr_n1=0;
    if (curr_tr1==transcript) {
      while ( (curr_tr1==transcript) & (i1<=n1) ) {
	curr_n1++;
	m1p[curr_n1]=m1_p[i1];
	i1++;
	if (i1<=n1) 
	  curr_tr1=m1_t[i1];
      }
      i1--;                       // set to the last entry corresponding to current transcript
      curr_tr1=transcript;        // set curr_tr1 to current transcript
    }                             // note that now we do not require that both motifs should be present in every instance of the crm

    while ( (curr_tr2<transcript) & (curr_tr2>=0) & (i2<n2) ) {
      i2++;curr_tr2=m2_t[i2];
    }
    
    curr_n2=0;
    if (curr_tr2==transcript) {
      while ( (curr_tr2==transcript) & (i2<=n2) ) {
	curr_n2++;
	m2p[curr_n2]=m2_p[i2];
	i2++;
	if (i2<=n2) 
	  curr_tr2=m2_t[i2];
      }
      i2--;
      curr_tr2=transcript;
    }

    if ( (curr_n1+curr_n2) >= N) {
      /* printf("transcript=%d\n",transcript); debug here */ 
      crm_2motifs(m1p,m2p,curr_n1,curr_n2,N,dw,max_crm_positions_transcript,bin_size,report_all,output,crm_positions_transcript);
      n_crms=output[0];
      errorcode=output[1];
      n_crms_reported=n_crms-errorcode;
    } else {
      n_crms=0;
      n_crms_reported=0;
    }
    if (n_crms>0)
      n12_tr++;
    
    if (report_all==1) {
      if ( (running_entries+n_crms) > max_crm_positions ) {
	ec--;
	running_entries_reported=running_entries;
	printf("running_entries=%d\n",running_entries);
	printf("running_entries_reported=%d\n",running_entries_reported);
      } else {
	for (i=1;i<=n_crms_reported;i++) {
	  running_entries++;
	  crm_positions[running_entries][1]=transcript;
	  crm_positions[running_entries][2]=crm_positions_transcript[i][1];
	  crm_positions[running_entries][3]=crm_positions_transcript[i][2];
	}
      }      // close check on running_entries+n_crms > max_crm_positions
    }        // close check on report_all == 1
  }          // close while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n12_tr;
  n[2]=ec;
  if (errorcode<0) {
    n[3]=running_entries_reported;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_ivector(m1p,1,max_transcript_occurrences);     
  free_ivector(m2p,1,max_transcript_occurrences);     
  if (report_all==1)
    free_imatrix(crm_positions_transcript,1,max_crm_positions_transcript,1,2);
}

void crm_3motifs(int *m1_p,int *m2_p,int *m3_p,int n1,int n2,int n3,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions)
{
  /* search for co-occurrences of at least N motifs in windows of size dw 
   * 
   * input
   *
   * m1_p = positions of motif 1
   * m2_p = positions of motif 2
   * m3_p = positions of motif 3
   * n1 = number of occurrences of motif 1
   * n2 = number of occurrences of motif 2
   * n3 = number of occurrences of motif 3
   * crm_positions = [1] = positions of the putative crms, [2] = number of motifs in each position
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * dw = window length (nucleotides)
   * max_crm_positions = maximum number of crm positions (to avoid overflows)
   * bin_size = step size for scanning through the transcript
   *
   * output = function output, output[0] = number of CRMS, output[1] = errorcode
   * errorcode = number of times where the position of a CRM could not be recorded due to overflow of max_crm_positions
   *
   * this function performs the search in a single transcript
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the inputs
   */

  int min_pos=20000;                           // minimum position of motif 1, 2 and 3
  int max_pos=-20000;                          // maximum position of motif 1, 2 and 3
  int n_crm=0;                                             // number of putative CRMs
  int pos;                                                 // current position 
  int upper_border;                                        // window upper border
  int init_index1=1;                                         // initial index for motif 1 in call to count_motif
  int init_index2=1;                                         // initial index for motif 2 in call to count_motif
  int init_index3=1;                                       // initial index for motif 3 in call to count_motif
  int out[2];                                              // output of count_motif
  int n_cooc;                                              // number of co-occurrences within the current window
  int errorcode=0;                                         // error code
  int i;

  /* determine minimum position to start scanning */
  if (n1>0) {
    if (m1_p[1]<min_pos)
      min_pos=m1_p[1];
    if (m1_p[n1]>max_pos)
      max_pos=m1_p[n1];
  }
  if (n2>0) {
    if (m2_p[1]<min_pos)
      min_pos=m2_p[1];
    if (m2_p[n2]>max_pos)
      max_pos=m2_p[n2];
  }
  if (n3>0) {
    if (m3_p[1]<min_pos)
      min_pos=m3_p[1];
    if (m3_p[n3]>max_pos)
      max_pos=m3_p[n3];
  }
  max_pos=max_pos-dw+bin_size;                                 // because the last window extends from this point for dw nucleotides
  pos=min_pos;

  if (report_all==1) {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;
      
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);
	n_cooc+=out[0];
	init_index1=out[1];
      }
      
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);
	n_cooc+=out[0];
	init_index2=out[1];
      }
      
      if (init_index3<=n3) {
	count_motif(m3_p,n3,init_index3,pos,upper_border,out);
	n_cooc+=out[0];
	init_index3=out[1];
      }
      
      if (n_cooc>=N) {
	n_crm++;
	if (n_crm<max_crm_positions) {
	  crm_positions[n_crm][1]=pos;
	  crm_positions[n_crm][2]=n_cooc;
	} else {
	  errorcode++;
	}
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }       // close while loop on pos
  } else {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;      
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);n_cooc+=out[0];init_index1=out[1];
      }
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);n_cooc+=out[0];init_index2=out[1];
      }
      if (init_index3<=n3) {
	count_motif(m3_p,n3,init_index3,pos,upper_border,out);n_cooc+=out[0];init_index3=out[1];
      }
      
      if (n_cooc>=N) {
	n_crm++;
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }       // close while loop on pos
  }         // close check on report_all
  
  output[0]=n_crm;
  output[1]=errorcode;
}

void crm_3motifs_allt(int *m1_t,int *m2_t,int *m3_t,int *m1_p,int *m2_p,int *m3_p,int n1,int n2,int n3,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions)
{
  /* call crm_3motifs for each transcript
   * 
   * input parameters
   * m1/2/3_t = transcript
   * m1/2/3_p = position
   * n1/2/3 = number of occurrrences for each motif
   * max_n_transcripts = stop computing if the number of transcripts with modules exceeds this number
   * crm_positions = [1]=transcript, [2]=position for the windows with co-occurring motifs, [3]=number of motif occurrences in the window
   * max_crm_positions = maximum allowed number of rows in crm_positions
   * N = at least N motif occurrences in window
   * dw = window size (nucleotides)
   * bin_size = step size to scan the sequence
   * n_transcripts = number of transcripts
   * 
   * n = output
   *    n[0]=running_entries;
   *    n[1]=n123_tr;
   *    n[2]=ec;
   *    n[3]=running_entries_reported;
   *
   * last modified 10-10-2003
   * 10-17-2003: do not require that all motifs must be present in each module
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the inputs
   */

  int **crm_positions_transcript;                 // crm positions in current transcript
  int *m1p;                                       // occurrences of motif 1 in current transcript
  int *m2p;
  int *m3p;
  int curr_n1,curr_n2,curr_n3;                    // number of occurrences of motif in current transcript
  int curr_tr1,curr_tr2,curr_tr3;                 // current transcript for each motif
  int ec=0;                                       // error code in current function
  int errorcode=0;                                // error code in call to crm_3motifs
  int i;
  int i1,i2,i3;
  int max_crm_positions_transcript=1000;         // maximum number of crm positions in a transcript (memory = 1000 x 1 (col) x (2 bytes/int) = 2 kb)
  int max_transcript_occurrences=1000;           // maximum number of occurrences of one motif per transcript (memory = 1000x1(cols)x2(bytes/int)=2kb)
  int n123;       
  int n123_tr=0;
  int n_crms;                                    // number of crms (output of crm_1motif
  int n_crms_reported;                           // number of reported crms (output of crm_1motif; this is the same as n_crms unless there was a memory overflow
  int output[2];                                 // output from the crm_2motifs program
  int running_entries=0;                         // running number of entries in ovdi
  int running_entries_reported;                  // number of reported crms; while n_crms and n_crms_reported refer to a single transcript, this refers to the hole function
  int transcript=0;

  i1=1;curr_tr1=m1_t[i1];                         // first transcript, motif 1
  i2=1;curr_tr2=m2_t[i2];                         // first transcrpit, motif 2
  i3=1;curr_tr3=m3_t[i3];                         // first transcript, motif 3
  
  m1p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 1 in current transcript
  m2p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 2 in current transcript
  m3p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 3 in current transcript
  if (report_all==1)
    crm_positions_transcript=imatrix(1,max_crm_positions_transcript,1,2);
  
  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;           // n12_tr will always be less than this

  while ( (transcript < n_transcripts) & (n123_tr<=max_n_transcripts) ) {
    transcript++;
    
    while ( (curr_tr1<transcript) & (curr_tr1>=0) & (i1<n1) ) {
      i1++;curr_tr1=m1_t[i1];
    }
    curr_n1=0;
    if (curr_tr1==transcript) {
      while ( (curr_tr1==transcript) & (i1<=n1) ) {
	curr_n1++;
	m1p[curr_n1]=m1_p[i1];
	i1++;
	if (i1<=n1) 
	  curr_tr1=m1_t[i1];
      }
      i1--;                       // set to the last entry corresponding to current transcript
      curr_tr1=transcript;        // set curr_tr1 to current transcript
    }      

    while ( (curr_tr2<transcript) & (curr_tr2>=0) & (i2<n2) ) {
      i2++;curr_tr2=m2_t[i2];
    }      
    curr_n2=0;
    if (curr_tr2==transcript) {
      while ( (curr_tr2==transcript) & (i2<=n2) ) {
	curr_n2++;
	m2p[curr_n2]=m2_p[i2];
	i2++;
	if (i2<=n2) 
	  curr_tr2=m2_t[i2];
      }
      i2--;
      curr_tr2=transcript;
    }
	
    while ( (curr_tr3<transcript) & (curr_tr3>=0) & (i3<n3) ) {
      i3++;curr_tr3=m3_t[i3];
    }
    curr_n3=0;
    if (curr_tr3==transcript) {
      while ( (curr_tr3==transcript) & (i3<=n3) ) {
	curr_n3++;
	m3p[curr_n3]=m3_p[i3];
	i3++;
	if (i3<=n3) 
	  curr_tr3=m3_t[i3];
      }
      i3--;
      curr_tr3=transcript;
    }	
  
    if ( (curr_n1+curr_n2+curr_n3) >= N) {
      crm_3motifs(m1p,m2p,m3p,curr_n1,curr_n2,curr_n3,N,dw,max_crm_positions,bin_size,report_all,output,crm_positions_transcript);
      n_crms=output[0];
      errorcode=output[1];
      n_crms_reported=n_crms-errorcode;
    } else {
      n_crms=0;
      n_crms_reported=0;
    }
    if (n_crms>0)
      n123_tr++;
    
    if (report_all==1) {
      if ( (running_entries+n_crms) > max_crm_positions ) {
	ec--;
	running_entries_reported=running_entries;
	printf("running_entries=%d\n",running_entries);
	printf("running_entries_reported=%d\n",running_entries_reported);
      } else {
	for (i=1;i<=n_crms_reported;i++) {
	  running_entries++;
	  crm_positions[running_entries][1]=transcript;
	  crm_positions[running_entries][2]=crm_positions_transcript[i][1];
	  crm_positions[running_entries][3]=crm_positions_transcript[i][2];
	}
      }   // close check on running_entries+n_crms > max_crm_positions
    }     // close check on report_all==1
  }       // close while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n123_tr;
  n[2]=ec;
  if (errorcode<0) {
    n[3]=running_entries_reported;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_ivector(m1p,1,max_transcript_occurrences);     
  free_ivector(m2p,1,max_transcript_occurrences);     
  free_ivector(m3p,1,max_transcript_occurrences);     
  if (report_all==1)
    free_imatrix(crm_positions_transcript,1,max_crm_positions_transcript,1,2);
}

void crm_4motifs(int *m1_p,int *m2_p,int *m3_p,int *m4_p,int n1,int n2,int n3,int n4,int N,int dw,int max_crm_positions,int bin_size,int report_all,int *output,int **crm_positions)
{
  /* search for co-occurrences of at least N motifs in windows of size dw 
   * 
   * input
   *
   * m1_p = positions of motif 1
   * m2_p = positions of motif 2
   * m3_p = positions of motif 3
   * m4_p = positions of motif 4
   * n1 = number of occurrences of motif 1
   * n2 = number of occurrences of motif 2
   * n3 = number of occurrences of motif 3
   * n4 = number of occurrences of motif 4
   * crm_positions = [1] = positions of the putative crms, [2] = number of motifs in each position
   * N = minimum number of putative transcription factor binding sites to state that the window is a putative CRM
   * dw = window length (nucleotides)
   * max_crm_positions = maximum number of crm positions (to avoid overflows)
   * bin_size = step size for scanning through the transcript
   *
   * output = function output, output[0] = number of CRMS, output[1] = errorcode
   * errorcode = number of times where the position of a CRM could not be recorded due to overflow of max_crm_positions
   *
   * this function performs the search in a single transcript
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the inputs
   */

  int min_pos=20000;
  int max_pos=-20000;
  int n_crm=0;                                             // number of putative CRMs
  int pos;                                                 // current position 
  int upper_border;                                        // window upper border
  int init_index1=1;                                       // initial index for motif 1 in call to count_motif
  int init_index2=1;                                       // initial index for motif 2 in call to count_motif
  int init_index3=1;                                       // initial index for motif 3 in call to count_motif
  int init_index4=1;
  int out[2];                                              // output of count_motif
  int n_cooc;                                              // number of co-occurrences within the current window
  int errorcode=0;                                         // error code
  int i;

  /* determine minimum and position to start scanning */
  if (n1>0) {
    if (m1_p[1]<min_pos)
      min_pos=m1_p[1];
    if (m1_p[n1]>max_pos)
      max_pos=m1_p[n1];
  }
  if (n2>0) {
    if (m2_p[1]<min_pos)
      min_pos=m2_p[1];
    if (m2_p[n2]>max_pos)
      max_pos=m2_p[n2];
  }
  if (n3>0) {
    if (m3_p[1]<min_pos)
      min_pos=m3_p[1];
    if (m3_p[n3]>max_pos)
      max_pos=m3_p[n3];
  }
  if (n4>0) {
    if (m4_p[1]<min_pos)
      min_pos=m4_p[1];
    if (m4_p[n4]>max_pos)
      max_pos=m4_p[n4];
  }
  max_pos=max_pos-dw+bin_size;                                 // because the last window extends from this point for dw nucleotides
  pos=min_pos;

  if (report_all==1) {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;
    
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);n_cooc+=out[0];init_index1=out[1];
      }
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);n_cooc+=out[0];init_index2=out[1];
      }
      if (init_index3<=n3) {
	count_motif(m3_p,n3,init_index3,pos,upper_border,out);n_cooc+=out[0];init_index3=out[1];
      }
      if (init_index4<=n4) {
	count_motif(m4_p,n4,init_index4,pos,upper_border,out);n_cooc+=out[0];init_index4=out[1];
      }

      if (n_cooc>=N) {
	n_crm++;
	if (n_crm<max_crm_positions) {
	  crm_positions[n_crm][1]=pos;
	  crm_positions[n_crm][2]=n_cooc;
	} else {
	  errorcode++;
	}
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }            // close while loop on pos
  } else {
    while (pos<=max_pos) {
      upper_border=pos+dw;
      n_cooc=0;
      if (init_index1<=n1) {
	count_motif(m1_p,n1,init_index1,pos,upper_border,out);n_cooc+=out[0];init_index1=out[1];
      }
      if (init_index2<=n2) {
	count_motif(m2_p,n2,init_index2,pos,upper_border,out);n_cooc+=out[0];init_index2=out[1];
      }
      if (init_index3<=n3) {
	count_motif(m3_p,n3,init_index3,pos,upper_border,out);n_cooc+=out[0];init_index3=out[1];
      }
      if (init_index4<=n4) {
	count_motif(m4_p,n4,init_index4,pos,upper_border,out);n_cooc+=out[0];init_index4=out[1];
      }
      if (n_cooc>=N) {
	n_crm++;
	pos+=dw;
      } else {
	pos+=bin_size;
      }
    }            // close while loop on pos
  }              // close check on report_all==1
  
  output[0]=n_crm;
  output[1]=errorcode;
}

void crm_4motifs_allt(int *m1_t,int *m2_t,int *m3_t,int *m4_t,int *m1_p,int *m2_p,int *m3_p,int *m4_p,int n1,int n2,int n3,int n4,int max_n_transcripts,int max_crm_positions,int N,int dw,int bin_size,int n_transcripts,int report_all,int *n,int **crm_positions)
{
  /* call crm_4motifs for each transcript
   * 
   * input parameters
   * m1/2/3/4_t = transcript
   * m1/2/3/4_p = position
   * n1/2/3/4 = number of occurrrences for each motif
   * max_n_transcripts = stop computing if the number of transcripts with modules exceeds this number
   * crm_positions = [1]=transcript, [2]=position for the windows with co-occurring motifs, [3]=number of motif occurrences in the window
   * max_crm_positions = maximum allowed number of rows in crm_positions
   * N = at least N motif occurrences in window
   * dw = window size (nucleotides)
   * bin_size = step size to scan the sequence
   * n_transcripts = number of transcripts
   * 
   * n = output
   *    n[0]=running_entries;
   *    n[1]=n123_tr;
   *    n[2]=ec;
   *    n[3]=running_entries_reported;
   *
   * 10-17-2003: written from crm_3motifs_allt
   * 10-18-2003: added report_all
   * 10-18-2003: changed the order of the inputs
   */

  int **crm_positions_transcript;                 // crm positions in current transcript
  int *m1p;                                       // occurrences of motif 1 in current transcript
  int *m2p;
  int *m3p;
  int *m4p;
  int curr_n1,curr_n2,curr_n3,curr_n4;            // number of occurrences of motif in current transcript
  int curr_tr1,curr_tr2,curr_tr3,curr_tr4;        // current transcript for each motif
  int ec=0;                                       // error code in current function
  int errorcode=0;                                // error code in call to crm_3motifs
  int i;
  int i1,i2,i3,i4;
  int max_crm_positions_transcript=1000;         // maximum number of crm positions in a transcript (memory = 1000 x 1 (col) x (2 bytes/int) = 2 kb)
  int max_transcript_occurrences=1000;           // maximum number of occurrences of one motif per transcript (memory = 1000x1(cols)x2(bytes/int)=2kb)
  int n1234;       
  int n1234_tr=0;
  int n_crms;                                    // number of crms (output of crm_1motif
  int n_crms_reported;                           // number of reported crms (output of crm_1motif; this is the same as n_crms unless there was a memory overflow
  int output[2];                                 // output from the crm_2motifs program
  int running_entries=0;                         // running number of entries in ovdi
  int running_entries_reported;                  // number of reported crms; while n_crms and n_crms_reported refer to a single transcript, this refers to the hole function
  int transcript=0;

  i1=1;curr_tr1=m1_t[i1];                         // first transcript, motif 1
  i2=1;curr_tr2=m2_t[i2];                         // first transcrpit, motif 2
  i3=1;curr_tr3=m3_t[i3];                         // first transcript, motif 3
  i3=4;curr_tr4=m4_t[i4];                         // first transcript, motif 4
  
  m1p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 1 in current transcript
  m2p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 2 in current transcript
  m3p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 3 in current transcript
  m4p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 4 in current transcript
  if (report_all==1)
    crm_positions_transcript=imatrix(1,max_crm_positions_transcript,1,2);
  
  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;           // n1234_tr will always be less than this

  while ( (transcript < n_transcripts) & (n1234_tr<=max_n_transcripts) ) {
    transcript++;

    while ( (curr_tr1<transcript) & (curr_tr1>=0) & (i1<n1) ) {
      i1++;curr_tr1=m1_t[i1];
    }
    curr_n1=0;
    if (curr_tr1==transcript) {
      while ( (curr_tr1==transcript) & (i1<=n1) ) {
	curr_n1++;
	m1p[curr_n1]=m1_p[i1];
	i1++;
	if (i1<=n1) 
	  curr_tr1=m1_t[i1];
      }
      i1--;                       // set to the last entry corresponding to current transcript
      curr_tr1=transcript;        // set curr_tr1 to current transcript
    }      

    while ( (curr_tr2<transcript) & (curr_tr2>=0) & (i2<n2) ) {
      i2++;curr_tr2=m2_t[i2];
    }      
    curr_n2=0;
    if (curr_tr2==transcript) {
      while ( (curr_tr2==transcript) & (i2<=n2) ) {
	curr_n2++;
	m2p[curr_n2]=m2_p[i2];
	i2++;
	if (i2<=n2) 
	  curr_tr2=m2_t[i2];
      }
      i2--;
      curr_tr2=transcript;
    }
	
    while ( (curr_tr3<transcript) & (curr_tr3>=0) & (i3<n3) ) {
      i3++;curr_tr3=m3_t[i3];
    }
    curr_n3=0;
    if (curr_tr3==transcript) {
      while ( (curr_tr3==transcript) & (i3<=n3) ) {
	curr_n3++;
	m3p[curr_n3]=m3_p[i3];
	i3++;
	if (i3<=n3) 
	  curr_tr3=m3_t[i3];
      }
      i3--;
      curr_tr3=transcript;
    }	
  
     while ( (curr_tr4<transcript) & (curr_tr4>=0) & (i4<n4) ) {
      i4++;curr_tr4=m4_t[i4];
    }
    curr_n4=0;
    if (curr_tr4==transcript) {
      while ( (curr_tr4==transcript) & (i4<=n4) ) {
	curr_n4++;
	m4p[curr_n4]=m4_p[i4];
	i4++;
	if (i4<=n4) 
	  curr_tr4=m4_t[i4];
      }
      i4--;
      curr_tr4=transcript;
    }	

   if ( (curr_n1+curr_n2+curr_n3+curr_n4) >= N) {
      crm_4motifs(m1p,m2p,m3p,m4p,curr_n1,curr_n2,curr_n3,curr_n4,N,dw,max_crm_positions,bin_size,report_all,output,crm_positions_transcript);
      n_crms=output[0];
      errorcode=output[1];
      n_crms_reported=n_crms-errorcode;
    } else {
      n_crms=0;
      n_crms_reported=0;
    }
    if (n_crms>0)
      n1234_tr++;
    
    if (report_all==1) {
      if ( (running_entries+n_crms) > max_crm_positions ) {
	ec--;
	running_entries_reported=running_entries;
	printf("running_entries=%d\n",running_entries);
	printf("running_entries_reported=%d\n",running_entries_reported);
      } else {
	for (i=1;i<=n_crms_reported;i++) {
	  running_entries++;
	  crm_positions[running_entries][1]=transcript;
	  crm_positions[running_entries][2]=crm_positions_transcript[i][1];
	  crm_positions[running_entries][3]=crm_positions_transcript[i][2];
	}
      }   // close check on running_entries+n_crms > max_crm_positions
    }     // close check on report_all == 1
  }       // while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n1234_tr;
  n[2]=ec;
  if (errorcode<0) {
    n[3]=running_entries_reported;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_ivector(m1p,1,max_transcript_occurrences);     
  free_ivector(m2p,1,max_transcript_occurrences);     
  free_ivector(m3p,1,max_transcript_occurrences);     
  free_ivector(m4p,1,max_transcript_occurrences);     
  if (report_all==1)
    free_imatrix(crm_positions_transcript,1,max_crm_positions_transcript,1,2);
}

int vector_overlap_v2(int *original_vector,int n,int max_overlap,int *output_vector,char *original_indices_vector,char *output_indices_vector)
{
  /* given a vector of n integer values, original_vector keep only entries where the distance between adjacent entries is more than max_overlap
   * this can be used to remove duplicates by setting max_overlap=0
   * NOTE: original_vector is assumed to be sorted in ascending order
   * returns the number of entries kept and those entries in output_vector
   * new in version 2: 01-02-2004: also extract the corresponding entries in another integer vector
   */

  int i;
  int d;
  int n_kept=1;

  output_vector[n_kept]=original_vector[1];               // the first entry always stays
  output_indices_vector[n_kept]=original_indices_vector[1];
  for (i=2;i<=n;i++) {
    if ( (original_vector[i]-original_vector[i-1]) > max_overlap ) {
      n_kept++;
      output_vector[n_kept]=original_vector[i];
      output_indices_vector[n_kept]=original_indices_vector[i];
    }
  }

  return n_kept;
}

int vector_overlap_v3(int *original_pos,char *original_strand,int *output_pos,char *output_strand,int n,int max_overlap)
{
  /* given a vector of n integer values, original_vector keep only entries where the distance between adjacent entries is more than max_overlap
   * this can be used to remove duplicates by setting max_overlap=0
   * NOTE: original_vector is assumed to be sorted in ascending order
   * 
   * new in version 3: 
   * 01-30-2004: changed the order of the inputs
   * 01-30-2004: changed the name of the inputs
   * 01-30-2004: discard only if the two adjacent occurrence are in the same strand
   * returns the number of entries kept and those entries in output_vector
   */

  int i;
  int d;
  int n_kept=1;
  int in;

  output_pos[1]=original_pos[1];               // the first entry always stays
  output_strand[1]=original_strand[1];      
  for (i=2;i<=n;i++) {
    in=0;
    if ( (original_pos[i]-original_pos[i-1]) > max_overlap )
      in=1;
    if (original_strand[i]!=original_strand[i-1])
      in=1;
    if (in==1) {
      n_kept++;
      output_pos[n_kept]=original_pos[i];
      output_strand[n_kept]=original_strand[i];
    }
  }

  return n_kept;
}

void convert_positions2tss(int *pos,int n,int el,int il)
{
  /* given a set of positions measured with respect to the 3' end of the sequence,
   * convert them to positions measured with respect to the TSS
   * el,il are the actual lengths in the current sequence
   */

  int i;                      // loop counter
  int current_pos;            // current position

  for (i=1;i<=n;i++) {
    current_pos=pos[i];
    if (current_pos>il) {
      current_pos=current_pos-il;
      if (current_pos>el) {	                /* occurrence within the upstream sequence */
	current_pos=current_pos-el;
	pos[i]=-current_pos;
      } else {                          	/* occurrence within the first exon */
	current_pos=el-current_pos;             // measure with respect to the tss
	pos[i]=current_pos;
      }                                         // close check on whether it is on the first exon
    } else {                                   /* occurrence within the first intron */
      current_pos=il-current_pos;                      // measure with respect to the exon-intron boundary
      current_pos=current_pos+el;
      pos[i]=current_pos;
    }        // close check whether it is in the first intron
  }          // close i loop
}

int convert_position2tss(int pos,int el,int il)
{
  /* given a position measured with respect to the 3' end of the sequence,
   * convert it to position measured with respect to the TSS
   * el,il are the actual lengths in the current sequence
   */

  int current_pos;            // current position

  current_pos=pos;
  if (current_pos>il) {
    current_pos=current_pos-il;
    if (current_pos>el) {	                /* occurrence within the upstream sequence */
      current_pos=current_pos-el;
      pos=-current_pos;
    } else {                          	/* occurrence within the first exon */
      pos=el-current_pos;             // measure with respect to the tss
    }                                         // close check on whether it is on the first exon
  } else {                                   /* occurrence within the first intron */
    current_pos=il-current_pos;                      // measure with respect to the exon-intron boundary
    pos=current_pos+el;
  }        // close check whether it is in the first intron
 
  return pos;
}

void reflect_pos(int *pos,int n)
{
  /* reflect the points in pos
   * so that pos[i] becomes -pos[i]
   * this is used if the format of the positions is such that the upstream sequences are >0 and the exon and intron sequences are <0
   * and we wish to convert to the TSS format
   */

  int i;                      // loop counter

  for (i=1;i<=n;i++) 
    pos[i]=-pos[i];

}
