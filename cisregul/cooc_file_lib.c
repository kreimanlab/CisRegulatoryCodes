/* cooc_file_lib.c
 * functions used in cooc_4m and cooc_single_1234m
 * created 02-19-2004
 * gk, all rights reserved
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* function definitions */
void convert_positions2tss(int *pos,int n,int el,int il);
int filter_positions_v5(int *pos,char *strand,int n,int max_ul,int max_el,int max_il,int el,int *filtered_pos,char *filtered_strand);
int load_and_sort_scan_data_v3(char *filename_positions,char *filename_strands,int **data_t,int **data_p,int **data_s,int *motif_lengths,int sortpos,int max_occur,int n_motifs,int max_ul,int max_el,int max_il,int **transcript_lengths,int data_format,int filteron,int verbose);
void read_line_wstrand_v4(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il,int data_format,int filteron,int verbose);
int read_scan_file_v5(char *filename_positions,char *filename_strands,int *transcripts,int *positions,char *strands,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il,int *out,int data_format,int filteron,int verbose);
void reflect_pos(int *pos,int n);
int vector_overlap_v3(int *original_pos,char *original_strand,int *output_pos,char *output_strand,int n,int max_overlap);

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

int read_scan_file_v5(char *filename_positions,char *filename_strands,int *transcripts,int *positions,char *strands,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il,int *out,int data_format,int filteron,int verbose)
{
  /* here we read a whole scan file in the following format
   * header lines (start with '>' or '%')
   * transcript position position position ...
   * transcript position position position ...
   *
   * same format with strands instead of positions for the strands file
   *
   * data_format=1: positions measured with respect to the TSS but with the upstream as > 0 
   * data_format=2: positions measured with respect to the 3' end of the sequence 
   *
   * returns the total number of entries in transcripts and positions
   * returns also the transcripts in *transcripts and the positions in *positions
   *
   * new in v5
   * 01-31-2004: call read_line_wstrand_v4
   * 01-31-2004: get rid of neg2pos
   * 01-31-2004: added data_format
   * 01-31-2004: added verbose
   * 01-31-2004: added filteron (1 to filter the occurrences and 0 otherwise)
   * see full history at the end of this function
   */

  FILE *input_file_pos;     // handle to file with positions
  FILE *input_file_strand;  // handle to file with strands
  int t[5];                 // t[0]= current transcript  t[1]=number of positions in current transcript
  int *m_p;                 // positions in current transcript
  char *m_s;                // strands in current transcript
  int max_occur=10000;      // maximum number of positions per transcript
  int curr_tr;              // current transcript
  int curr_n;               // number of positions in current transcript
  int running_n=0;          // running number of entries in transcripts and positions
  int running_np=0;         // running number of entries in the plus strand
  int running_nm=0;         // running number of entries in the minus strand
  int i;                    // loop index
  int togoon=1;             // continue after loop depending on togoon
  //int call2or3=2;           // call read_line_wstrand_v2 unless there is an entry in the map; in that case, call read_line_wstrand_v3
  int n_orig;               // original number of entries
  int n_after_filter;       // number of entries after filtering
  int cumul_filter_out=0;   // cumulative number of entries filtered out
  int cumul_overlap_out=0;  // cumulative number of entries removed (including filter and overlap)

  /* i=1;
     while ( (i<10000) & (call2or3==2) ) {
     if (transcript_map[i][1]>0)
     call2or3=3;
     if (transcript_map[i][2]>0)
     call2or3=3;
     i++;
     }
  */

  input_file_pos=fopen(filename_positions,"r");
  if (!input_file_pos) {
    printf("error! i could not open %s\n",filename_positions);
    return running_n;
  }
  input_file_strand=fopen(filename_strands,"r");
  if (!input_file_strand) {
    printf("error! i could not open %s\n",filename_strands);
    return running_n;
  }
  m_p=ivector(1,max_occur);
  m_s=cvector(1,max_occur);
  while (togoon>=0) {
    read_line_wstrand_v4(input_file_pos,input_file_strand,t,m_p,m_s,max_occur,sortpos,max_ul,max_el,max_il,transcript_lengths,transcript_map,max_overlap,seq_ul,seq_el,seq_il,data_format,filteron,verbose);
    /* if (call2or3==3) {
       read_line_wstrand_v3(input_file_pos,input_file_strand,t,m_p,m_s,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,transcript_map,max_overlap,seq_ul,seq_el,seq_il);
       } else {
       read_line_wstrand_v2(input_file_pos,input_file_strand,t,m_p,m_s,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,max_overlap,seq_ul,seq_el,seq_il);
       } 
    */
    curr_tr=t[0];                                   // current transcript
    curr_n=t[1];                                    // number of occurrences in current transcript
    togoon=t[2];                                    // 1 if the call to read a line was ok
    n_orig=t[3];                                    // number of entries in the original scan output
    n_after_filter=t[4];                            // number of entries after filtering
    cumul_filter_out+=(n_orig-n_after_filter);      // cumulative number of entries filtered out
    cumul_overlap_out+=(n_orig-curr_n);             // cumulative number of entries eliminated due to overlap with adjacent occurrences
    if (curr_n>0) {
      for (i=1;i<=curr_n;i++) {
	running_n++;
	transcripts[running_n]=curr_tr;
	positions[running_n]=m_p[i];
	strands[running_n]=m_s[i];
	if (m_s[i]>0) 
	  running_np++;
	if (m_s[i]==0) 
	  running_nm++;
      }       // close i loop
    }         // close check on curr_n>0
  }           // close togoon
  fclose(input_file_pos);
  fclose(input_file_strand);
  out[0]=cumul_filter_out;
  out[1]=cumul_overlap_out;
  out[2]=running_np;
  out[3]=running_nm;
  free_ivector(m_p,1,max_occur);
  free_cvector(m_s,1,max_occur);

  return running_n;

  /* 
   * new in v4
   * 01-03-2004: add strand information
   *
   * new in v3
   * 10-30-2003: add seq_ul,seq_el,seq_il to allow a default sequence length to measure positions with respect to the 3' end
   * 10-28-2003: add transcript_map, if transcript_map[1][1]>0 then call read_line_v3, otherwise just call read_line_v2 and ignore transcript_map
   * 10-27-2003: call read_line_v2 (therefore allowing filtering of the scan output)
   * 10-27-2003: call to vector_overlap is now performed within read_line_v2
   * new in v2
   * 10-17-2003: make it simpler by including the header stuff into the main loop
   * 10-17-2003: ignore overlapping occurrences
   * 01-07-2004: added out variable out[0]=cumul_filter_out, out[1]=cumul_overlap_out, out[2]=running_np, out[3]=running_nm
   */
}

void read_line_wstrand_v4(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il,int data_format,int filteron,int verbose)
{
  /* given a file handle, read a single line
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   *   with the positions measured with respect to the 3' end of the sequence
   * data_format=1: positions measured with respect to the TSS but with the upstream as > 0 
   * data_format=2: positions measured with respect to the 3' end of the sequence 
   * variable t
   *    t[0] = transcript number
   *    t[1] = number of positions returned in data_p
   *    t[2] = -1 -> end of file   t[2]=0 -> n>max_occur  t[2]=1 -> ok  t[2]=2 -> n=0 (no positions)  t[2]=3 -> curr_transcript=0,odd
   *    t[3] = number of positions before filtering and before removing overlaps
   *    t[4] = number of positions after filtering before removing overlaps
   *
   * this comes from read_line_v3 to read a line from the scan for all files
   * 01-04-2004 here we also read the strand information
   * 01-30-2004 call vector_overlap_v3
   * new in version 4
   * 01-31-2004 attempt to make a single read_line_wstrand function including all formats
   * 01-31-2004 all positions are returned with respect to the TSS (the larger the number the more 3')
   * 01-31-2004 remove the neg2pos stuff
   * 01-31-2004 added filteron, if 1 then call the filter function
   */

  char txtline_pos[100000];              // read one line of the input positions file
  char txtline_strand[100000];           // read one line of the input strands file 
  int i,j,k,ll;                          // loop counters
  int *w;                                // vector holding the transcripts and strands for running entry
  int *v;                                // vector holding the transcripts and positions for running entry
  int start=0;                           // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];                             // number of entries
  int n;                                 // number of entries in a line
  int nm1;                               // n-1
  int n_chars=1000000;                   // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int header;                           // 1 if current line is a header line
  int *temp_positions;                  // temporary positions (contains the raw positions before filtering
  char *temp_strands;
  int *filtered_positions;              // filtered positions (contains the positions after filtering)
  char *filtered_strands;
  int *positions_no_overlap;            // contains the positions after filtering and removing overlaps
  char *strands_no_overlap;
  float *indices;                         // dummy variable used in call to sort function
  int curr_transcript;                  // current transcript
  int curr_ul,curr_el,curr_il;          // upstream length, exon length, intron length
  int n_after_filter;                   // number of occurrences after filtering
  int n_no_overlap=0;                   // number of occurrences that do not overlap
  int curr_transcript_index;

  t[0]=-1;t[1]=-1;t[2]=-1;t[3]=-1;t[4]=-1;
  if (feof(input_file_pos)!=0) {
    /* done with file, exit */
  } else {
    fgets(txtline_pos,n_chars,input_file_pos);
    fgets(txtline_strand,n_chars,input_file_strand);
    header=is_header(txtline_pos);
    if (header==1) {
      t[2]=2;
    } else {
      v=ivector(1,max_occur);
      w=ivector(1,max_occur);
      splitline2int(txtline_pos,v,start,np,max_occur);
      n=np[0]-1;
      nm1=n-1;
      splitline2int(txtline_strand,w,start,np,max_occur);
      if (n>0) {
	if (n>max_occur) {
	  n=max_occur;
	  nm1=max_occur;
	  t[2]=0;
	  printf("WARNING! n=%d and max_occur=%d\n",n,max_occur);
	}
	curr_transcript=v[1];
	if (data_format==1) {
	  /* positions measured with respect to the TSS but with the upstream as > 0 */
	  curr_transcript_index=v[1];
	  if (transcript_map[curr_transcript_index][2]>0)
	    curr_transcript=transcript_map[curr_transcript_index][2];
	} 
	// if (data_format==2) {
	/* positions measured with respect to the 3' end of the sequence */
	//curr_transcript=v[1];
	//}

	t[0]=curr_transcript;  // transcript number
	if (verbose)
	  printf("transcript=%d\n",curr_transcript);
	if (curr_transcript<=0) {
	  t[2]=3;
	  printf("WARNING!!! curr_transcript=%d\nthis is odd\ntxtline=%s\nn=%d\n",txtline_pos,n);
	  exit(1);
	} else {
	  /* get curr_ul,curr_el,curr_il */
	  curr_ul=transcript_lengths[curr_transcript][1];
	  if (curr_ul>seq_ul)
	    curr_ul=seq_ul;
	  curr_el=transcript_lengths[curr_transcript][2];
	  if (curr_el>seq_el)
	    curr_el=seq_el;
	  curr_il=transcript_lengths[curr_transcript][3];
	  if (curr_il>seq_il)
	    curr_il=seq_il;
	  /* memory allocation */
	  temp_positions=ivector(1,n);               //temp_positions=ivector(1,max_occur);
	  temp_strands=cvector(1,n);                 //temp_strands=cvector(1,max_occur);
	  filtered_positions=ivector(1,n);           //filtered_positions=ivector(1,max_occur);
	  filtered_strands=cvector(1,n);             //filtered_strands=cvector(1,max_occur);
	  indices=vector(1,max_occur);               //indices=vector(1,max_occur);
	  positions_no_overlap=ivector(1,n);         //positions_no_overlap=ivector(1,max_occur);
	  strands_no_overlap=cvector(1,n);           //strands_no_overlap=cvector(1,max_occur);

	  t[3]=n-1;                                  // number of positions before filtering and removing overlaps
	  t[2]=1;                                    // return ok signal
	  /*if (neg2pos==1) {
	    for (i=2;i<=n;i++) {
	    j=i-1;
	    temp_strands[j]=w[i];
	    if (v[i]>0) {
	    temp_positions[j]=v[i];
	    } else {
	    temp_positions[j]=-v[i];
	    }
	    }
	    } else { */
	    // for (i=1;i<=n;i++) {    01-02-2004 changed the loop to start at i=2; i=1 is reserved for the transcript number
	  for (i=2;i<=n;i++) {
	    j=i-1;
	    temp_positions[j]=v[i];
	    temp_strands[j]=w[i];
	  }
	  if (verbose) {
	    printf("\ttemp_positions (%d): ",n-1);
	    for (i=1;i<=nm1;i++) 
	      printf("%d %d,",temp_positions[i],temp_strands[i]);
	    printf("\n");
	  }
	  //}      // close check on neg2pos==1
	  //n_after_filter=filter_positions_v2(temp_positions,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	  
	  /* format conversion */
	  if (data_format==1) {
	    reflect_pos(temp_positions,nm1);
	  }
	  if (data_format==2) {
	    convert_positions2tss(temp_positions,nm1,curr_el,curr_il);
	  }
	  if (verbose) {
	    printf("\tconverted temp_positions (%d): ",nm1);
	    for (i=1;i<=nm1;i++) 
	      printf("%d %d,",temp_positions[i],temp_strands[i]);
	    printf("\n");
	  }
	  //n_after_filter=filter_positions_v4(temp_positions,temp_strands,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
	  if (filteron==1) {
	    n_after_filter=filter_positions_v5(temp_positions,temp_strands,nm1,max_ul,max_el,max_il,curr_el,filtered_positions,filtered_strands);
	  } else {
	    n_after_filter=nm1;
	    for (i=1;i<=n_after_filter;i++) {
	      filtered_positions[i]=temp_positions[i];
	      filtered_strands[i]=temp_strands[i];
	    }
	  }
	  t[4]=n_after_filter;
	  if (verbose) {
	    printf("\tfiltered_positions (%d): ",n_after_filter);
	    for (i=1;i<=n_after_filter;i++) 
	      printf("%d %d,",filtered_positions[i],filtered_strands[i]);
	    printf("\n");
	  }
	  if (n_after_filter>0) {
	    if (sortpos) {
	      for (i=1;i<=n_after_filter;i++)
		indices[i]=i;
	      quicksorti(n_after_filter,filtered_positions,indices);
	      sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
	    }
	    if (verbose) {
	      printf("\tsorted filtered_positions (%d): ",n_after_filter);
	      for (i=1;i<=n_after_filter;i++) 
		printf("%d %d,",filtered_positions[i],filtered_strands[i]);
	      printf("\n");
	    }
	    n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_overlap);
	    if (verbose) {
	      printf("\tpositions_no_overlap (%d): ",n_no_overlap);
	      for (i=1;i<=n_no_overlap;i++) 
		printf("%d %d,",positions_no_overlap[i],strands_no_overlap[i]);
	      printf("\n");
	    }
	    for (k=1;k<=n_no_overlap;k++) {
	      data_p[k]=positions_no_overlap[k];
	      data_s[k]=strands_no_overlap[k];
	    }
	  }     // close check on n_after_filter>0
	  t[1]=n_no_overlap;
	  /* free memory */
	  free_ivector(temp_positions,1,n);
	  free_cvector(temp_strands,1,n);
	  free_ivector(filtered_positions,1,n);
	  free_cvector(filtered_strands,1,n);
	  free_ivector(positions_no_overlap,1,n);
	  free_cvector(strands_no_overlap,1,n);
	  free_vector(indices,1,n);
	}     // close check on curr_transcript>0
      } else {
	t[2]=2;
      }   // close check on n>0
      free_ivector(v,1,max_occur);
      free_ivector(w,1,max_occur);
     }  // close if (is_header(txtline)==0)     
  }    // close feof check
}      // end of read_line_v2

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
