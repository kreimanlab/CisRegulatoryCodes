/* cooc_single_lib.c
 * functions used in cooc_single_2m/3m/4m.c
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* 01-31-2004: modified dist_2v_ordered_v2,dist_2v_allt_v2,dist_3v_ordered_v3,dist_3v_allt_v3,dist_3v_allt_from_m123_v1,dist_4v_ordered_v3,dist_4v_allt_v3,dist_4v_allt_from_m1234_v1 adding one input (pos_format) 
 * 02-02-2004: added verbose to dist_3v_allt_from_m123_v1,dist_2v_allt_v2,dist_3v_unordered_v3,dist_3v_allt_v3, dist_4v_allt_from_m1234_v1
 * 02-05-2004: moved all the functions that are not used to cooc_single_lib_oldversions.c
 */

/* 2 motif functions */
void dist_2v_ordered_v2(int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format,int verbose);
void dist_2v_unordered_v3(int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions);
void dist_2v_allt_v2(int *m1_t,int *m2_t,int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose);
/* 3 motif functions */
void dist_3v_ordered_v3(int **ovdi12_p,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format,int verbose);
void dist_3v_unordered_v3(int **ovdi12_p,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int verbose);
void dist_3v_allt_v3(int **ovdi12,int *m3_t,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose);
void dist_3v_allt_from_m123_v1(int *m1_t,int *m2_t,int *m3_t,int *m1_p,int *m2_p,int *m3_p,char *m1_s,char *m2_s,char *m3_s,int n1,int n2,int n3,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose);
/* 4 motif functions */
void dist_4v_ordered_v3(int **ovdi123_p,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format);
void dist_4v_unordered_v3(int **ovdi123_p,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions);
void dist_4v_allt_v3(int **ovdi123,int *m4_t,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format);
void dist_4v_allt_from_m1234_v1(int *m1_t,int *m2_t,int *m3_t,int *m4_t,int *m1_p,int *m2_p,int *m3_p,int *m4_p,char *m1_s,char *m2_s,char *m3_s,char *m4_s,int n1,int n2,int n3,int n4,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose);
/* file functions */
void read_line(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos);
void read_header(FILE *input_file,int *t,int *data_p,int max_occur);
int read_scan_file(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos);
int read_scan_file_v2(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos,int max_overlap);
void read_line_v2(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int max_overlap,int seq_ul,int seq_el,int seq_il);
void read_line_v3(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il);
int read_scan_file_v3(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il);
void read_line_wstrand_v2(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int max_overlap,int seq_ul,int seq_el,int seq_il);
void read_line_wstrand_v3(FILE *input_file_pos,FILE *input_file_strands,int *t,int *data_p,char *data_s,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il);
void read_line_wstrand_v4(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il,int data_format,int filteron,int verbose);
int read_scan_file_v4(char *filename_positions,char *filename_strands,int *transcripts,int *positions,char *strands,int neg2pos,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il,int *out);
int read_scan_file_v5(char *filename_positions,char *filename_strands,int *transcripts,int *positions,char *strands,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il,int *out,int data_format,int filteron,int verbose);
/* other functions */
int remove_overlaps(int **m_old,int **m_new,int *ml,int n_entries,int n_motifs);
int remove_overlaps_v2(int **m_old,int **m_new,int *ml,int n_entries,int n_motifs);
int range_constraint(int **m_old,int **m_new,int range_threshold,int n_entries,int n_motifs);
int get_max_n_transcripts(int n_transcripts_cl,int n_cl,double cdf_threshold,int n_transcripts_bck);

void dist_2v_ordered_v2(int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format,int verbose)
{
  /* here we only have the positions and strands for the two motifs in a given transcript  
   *
   * if get_positions==1, then fill ovdi with positions
   *  ovdi (n_ovdi x 4) p1 p2 i j   (where i and j correspond to the indices in m1_p and m2_p for the occurrences of motif 1 and motif 2 respectively)
   * 
   * n[0]=number of co-occurrences within the current transcript
   * n[1]=error code (0 if there are no errors, -1 if there is a memory overload
   *
   * ovdi (n_ovdi x 6):  position1  position2  i  j strand1 strand2 [note: in previous version, ovdi had 4 columns]
   *
   * 11-08-2003: if get_positions==0, then do not allocate memory for all_ovdi
   * 09-08-2003: try to speed up by stopping computation if we run out of motif 2
   * 09-08-2003: try to speed up by using an initial j and not starting from 0 each time
   * 07-25-2003: we assume here that the positions are sorted in increasing order
   * 07-25-2003: try to improve speed by avoiding computing distances when the second motif is past the max_dist limit
   * 07-25-2003: try to improve speed by avoiding computing distances when the first motif is past the maximum of motif 2 + max_dist
   * new in version 2
   * 12-29-2003: consider the strands (use m1_s and m2_s input with +1 or -1 for each strand)
   * 01-31-2004: added pos_format, pos_format=1 means: a larger number is more 3' (this corresponds to the position specification with respect to the TSS with upstream 
   *             positions as < 0 and exon and intron as > 0) 
   *             pos_format=-1 means: a larger number is more 5' (this corresponds to the position specificatino with respect to the 5' end of the sequence with all 
   *             positions as > 0)
   * 02-02-2004: added verbose
   */

  int i,j;                     // loop indices
  int p1;                      // position of motif 1
  int p2;                      // position of motif 2
  char s1;                      // strand of motif 1
  char s2;                      // strand of motif 2
  int d;                       // distance between motifs
  int n_ovdi=0;                // number of co-occurring motifs (before call to remove_overlaps)
  int **all_ovdi;              // positions of overlapping motifs
  int n_cooc;                  // number of co-occurring motifs
  int max_p2;                  // maximum position of motif 2
  int init_j=1;                // initial occurrence of motif 2
  int firstnear;               // used in the computation of init_j

  if ( (pos_format!=1) & (pos_format!=-1) ) {
    printf("ERROR!\ncooc_single_lib.c\ndist_2v_ordered_v2\npos_format=%d\nd should be either +1 or -1\nexiting...\n",pos_format);
    exit(1);
  }

  if (verbose)
    printf("max_dist=%d\tmin_dist=%d\n",max_dist,min_dist);

  max_p2=m2_p[n2]+max_dist;    // because we assume it is already sorted, the maximum is simply the last entry
  if (get_positions==1)
    all_ovdi=imatrix(1,max_n_ovdi,1,6);

  i=0;
  p1=m1_p[1];                  // start with the first value, do not even enter the loop is this is already beyond max_p2
  while ( (i<n1) & (n_ovdi<max_n_ovdi) & (p1<max_p2) ) {
    i++;
    j=init_j;                       // start processing motif 2 at this entry
    p1=m1_p[i];                     // position of motif 1
    s1=m1_s[i];                     // strand of motif 1
    d=0;                            // this took me a couple of hours ...
    firstnear=1;
    if (verbose)
      printf("i=%d\tp1=%d\ts1=%d\n",i,p1,s1);
    while ( (j<=n2) & (n_ovdi<max_n_ovdi) & (d<max_dist) ) {
      p2=m2_p[j];
      s2=m2_s[j];
      d=0;                          // this took me another couple of hours...
      if ( (s1==1) & (s2==1) )
	d=p2-p1;
      if ( (s1==0) & (s2==0) )
	d=p1-p2;
      d=d*pos_format;
      if (verbose)
	printf("\tj=%d\tp2=%d\ts2=%d\td=%d\n",j,p2,s2,d);
      if ( (d<max_dist) & (d>0) ) {
	if (firstnear==1) {
	  init_j=j;
	  firstnear=0;
	}
	if (d>min_dist) {
	  n_ovdi++;
	  if (get_positions==1) {
	    all_ovdi[n_ovdi][1]=p1;
	    all_ovdi[n_ovdi][2]=p2;
	    all_ovdi[n_ovdi][3]=i;
	    all_ovdi[n_ovdi][4]=j;
	    all_ovdi[n_ovdi][5]=s1;
	    all_ovdi[n_ovdi][6]=s2;
	  }
	}    // close if (d>min_dist)
      }      // close ff (d<max_dist)
      j++;
    }        // close j loop
    if (firstnear == 1) {
      if (p1>m2_p[init_j])
	init_j++;                // we could not find any occurrence of motif 2 near one of motif 1 and motif 2 is to the left of motif 1, start from the next one
      if (init_j>n2)
	p1=max_p2+1;             // we ran out of motif 2, stop computation
    }
  }          // close i loop

  if ( (n_ovdi>0) & (get_positions==1) ) {
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n_ovdi,2);
  } else {
    n_cooc=n_ovdi;
  }

  if ( (i<n1) & (p1<max_p2) ) {
    n[1]=-1;
  } else {
    n[1]=0;
  }
  n[0]=n_cooc;      // 0 by default

  if (get_positions==1)
    free_imatrix(all_ovdi,1,max_n_ovdi,1,6);
}

void dist_2v_unordered_v3(int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions)
{
  /* compute the co-occurrences of two motifs without regard to the order
   * here we only have the positions and strands of the two motifs in a given transcript  
   * 
   * if get_positions==1, then fill ovdi with positions
   *   ovdi (n_ovdi x 6) p1 p2 i j s1 s2  (where i and j correspond to the indices in m1_p and m2_p for the occurrences of motif 1 and motif 2 respectively)
   * 
   * n[0]=number of co-occurrences within the current transcript
   * n[1]=error code (0 if there are no errors, -1 if there is a memory overload
   *
   * here order does not matter
   * in contrast to v1, here we loop only once through the data and directly get the co-occurring events
   *
   * 09-08-2003: add init_j to try to make it faster
   * 01-12-2004: add strand information (new in version 3)
   * 01-12-2004: ovdi now has 6 columns
   */

  int i,j;                 // loop indices (indices of motif1 and motif2 respectively in m1_p and m2_p)
  int p1;                  // position1
  int p2;                  // position2
  char s1,s2;              // strand of motif 1, motif2
  int d;                   // distance between motifs
  int ad;                  // absolute value of d
  int n_ovdi=0;            // number of co-occurring entries before call to remove_overlaps
  int **all_ovdi;          // all co-occurring entries (before call to remove_overlaps)
  int n_cooc;              // number of co-occurring entires
  int max_p2;              // maximum position of motif 2
  int min_p2;              // minimum positions of motif 2
  int init_j=1;
  int firstnear;
  int init_i;

  max_p2=m2_p[n2]+max_dist;                     // maximum position of motif 2; we will compare the position of motif1 until it has passed this value
  min_p2=m2_p[1]-max_dist;                      // minimum position of motif 2;
  if (get_positions==1)
    all_ovdi=imatrix(1,max_n_ovdi,1,6);

  /* move fast until we reach min_p2 */
  init_i=1;p1=m1_p[init_i];
  while (p1<min_p2) {
    init_i++;
    p1=m1_p[init_i];
  }

  i=init_i-1;
  p1=m1_p[1];                                   // initial position for motif 1
  while ( (i<n1) & (n_ovdi<max_n_ovdi) & (p1<max_p2) ) {
    i++;
    j=init_j;
    p1=m1_p[i];
    s1=m1_s[i];
    d=0;                                        // this took another day of work (even after making the same mistake in dist_2v_ordered !!!)
    firstnear=1;
    while ( (j<=n2) & (n_ovdi<max_n_ovdi) & (d<max_dist) ) {
      p2=m2_p[j];
      s2=m2_s[j];
      d=p2-p1;                    // this is used to stop the loop after we are past the max_dist (and then ad is the absolute value of d)
      ad=0;                       // do not consider unless the strands are ok
      if ( ( (s1==1) & (s2==1) ) | ( (s1==0) & (s2==0) ) ) 
	ad=fabs(d);
      //d=p2-p1;ad=fabs(d);
      /* if (p2>p1) {
	 ad=p2-p1;
	 } else {
	 ad=p1-p2;
	 }
      */
      if ( (ad<max_dist) & (ad>0) ) {
	//if ( (ad<max_dist) & (ad>min_dist) ) {
	if (firstnear==1) {
	  init_j=j;
	  firstnear=0;
	}
	if (ad>min_dist) {
	  n_ovdi++;
	  if (get_positions==1) {
	    all_ovdi[n_ovdi][1]=p1;
	    all_ovdi[n_ovdi][2]=p2;
	    all_ovdi[n_ovdi][3]=i;
	    all_ovdi[n_ovdi][4]=j;
	    all_ovdi[n_ovdi][5]=s1;
	    all_ovdi[n_ovdi][6]=s2;
	  }
	}    // close if (ad>min_dist)
      }      // close if (ad<min_dist) 
      j++;
    }        // close j loop
  }          // close i loop

  if ( (n_ovdi>0) & (get_positions==1) ) {
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n_ovdi,2);
  } else {
    n_cooc=n_ovdi;
  }

  if ( (i<n1) & (p1<max_p2) ) {
    n[1]=-1;
  } else {
    n[1]=0;
  }
  n[0]=n_cooc;      // 0 by default

  if (get_positions==1)
    free_imatrix(all_ovdi,1,max_n_ovdi,1,6);
}

void dist_2v_allt_v2(int *m1_t,int *m2_t,int *m1_p,int *m2_p,char *m1_s,char *m2_s,int n1,int n2,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose)
{
  /* call dist_2v_unordered_v3 or dist_2v_ordered_v2 for all transcripts 
   * get_positions=1 --> to fill ovdi with the positions
   * order_constraint=1 --> order matters
   *
   * outputs: ovdi (running_entries x 7): 
   *        ovdi[running_entries][1]=transcript;
   *	    ovdi[running_entries][2]=ovdi_single_transcript[i][1];         position, motif 1
   *        ovdi[running_entries][3]=ovdi_single_transcript[i][2];         position, motif 2
   *	    ovdi[running_entries][4]=ovdi_single_transcript[i][3];         index, motif 1
   *        ovdi[running_entries][5]=ovdi_single_transcript[i][4];         index, motif 2
   *	    ovdi[running_entries][6]=ovdi_single_transcript[i][5];         strand, motif 1
   *        ovdi[running_entries][7]=ovdi_single_transcript[i][6];         strand, motif 2
   *
   * n (4 entries)
   *        n[0]=running_entries;
   *        n[1]=n12_tr;
   *        n[2]=errorcode;
   *        n[3]=running_entries when get_positions is stopped because of overflow
   *
   * new in version 2:
   * 01-01-2004: use strand information
   * 01-01-2004: call dist_2v_ordered_v2
   * 01-01-2004: output now has 7 columns (it includes the strand information)
   * 01-12-2004: call dist_2v_unordered_v3
   * 01-31-2004: added pos_format (see dist_2v_ordered_v2 for explanation of this parameter)
   */
 
  int **ovdi_single_transcript;                  // ovdi output for a single transcript
  int *m1p;                                      // positions of motif 1 in current transcript
  char *m1s;                                      // strands of motif 1 in current transcript
  int *m2p;                                      // positions of motif 2 in current transcript
  char *m2s;                                      // strands of motif 2 in current transcript
  int curr_n1,curr_n2;                           // number of occurrences of motif in current transcript
  int curr_tr1,curr_tr2;                         // current transcript for motif 1,2
  int errorcode=0;                               // error code
  int i,j;                                       // loop indices
  int i1,i2;                                     // loop indices
  int max_n_ovdi_single_transcript=50000;        // maximum number of co-occurrences in current transcript (memory = 10000x4(cols)x2(bytes/int)=80kb)
  int max_transcript_occurrences=10000;           // maximum number of occurrences of one motif per transcript (memory = 1000x1(cols)x2(bytes/int)=2kb)
  int n12;                                       // number of co-occurrences
  int n12_stop=0;                                // number of entries when the search is stopped 
  int n12_tr=0;                                  // number of transcripts with co-occurring motifs
  int running_entries=0;                         // running number of entries in ovdi
  int transcript=0;                              // transcript number=1;curr_tr1=m1_t[i1];

  i1=1;curr_tr1=m1_t[i1];
  i2=1;curr_tr2=m2_t[i2];
  
  m1p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 1 in current transcript
  m1s=cvector(1,max_transcript_occurrences);      // will hold the strands of motif 1 in current transcript
  m2p=ivector(1,max_transcript_occurrences);      // will hold the positions of motif 2 in current transcript
  m2s=cvector(1,max_transcript_occurrences);      // will hold the strands of motif 2 in current transcript
  ovdi_single_transcript=imatrix(1,max_n_ovdi_single_transcript,1,6);        // two positions, two indices, two strands
  
  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;           // n12_tr will always be less than this

  while ( (transcript < n_transcripts) & (n12_tr<=max_n_transcripts) ) {
    transcript++;
    if ( (curr_tr1<0) | (curr_tr2<0) ) 
      transcript=n_transcripts+1;       // at least one of the motifs is done, let's call it a day
    
    while ( (curr_tr1<transcript) & (curr_tr1>=0) & (i1<n1) ) {
      i1++;curr_tr1=m1_t[i1];
    }
    if (curr_tr1==transcript) {
      curr_n1=0;
      while ( (curr_tr1==transcript) & (i1<=n1) ) {
	curr_n1++;
	m1p[curr_n1]=m1_p[i1];
	m1s[curr_n1]=m1_s[i1];
	i1++;
	if (i1<=n1) 
	  curr_tr1=m1_t[i1];
      }
      i1--;                       // set to the last entry corresponding to current transcript
      curr_tr1=transcript;        // set curr_tr1 to current transcript
      
      while ( (curr_tr2<transcript) & (curr_tr2>=0) & (i2<n2) ) {
	i2++;curr_tr2=m2_t[i2];
      }
      
      if (curr_tr2==transcript) {
	curr_n2=0;
	while ( (curr_tr2==transcript) & (i2<=n2) ) {
	  curr_n2++;
	  m2p[curr_n2]=m2_p[i2];
	  m2s[curr_n2]=m2_s[i2];
	  i2++;
	  if (i2<=n2) 
	    curr_tr2=m2_t[i2];
	}
	i2--;
	curr_tr2=transcript;

	if (verbose) {
	  printf("transcript=%d\n",transcript);
	  printf("motif 1 (%d): ",curr_n1);
	  for (i=1;i<=curr_n1;i++)
	    printf("%d %d,",m1p[i],m1s[i]);
	  printf("\n");
	  printf("motif 2 (%d): ",curr_n2);
	  for (i=1;i<=curr_n2;i++)
	    printf("%d %d,",m2p[i],m2s[i]);
	  printf("\n");
	}
	if (order_constraint==1) {
	  dist_2v_ordered_v2(m1p,m2p,m1s,m2s,curr_n1,curr_n2,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions,pos_format,0);
	} else {
	  dist_2v_unordered_v3(m1p,m2p,m1s,m2s,curr_n1,curr_n2,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions);
	}
	n12=n[0];
	errorcode+=n[1];           // n[1]<0 if there are errors
	if (verbose) {
	  printf("ovdi_single_transcript (%d)\n",n12);
	  for (i=1;i<=n12;i++)
	    printf("%d %d %d %d\n",ovdi_single_transcript[i][1],ovdi_single_transcript[i][2],ovdi_single_transcript[i][3],ovdi_single_transcript[i][4]);
	}
	if (n12>0)
	  n12_tr++;
	
	if (get_positions==1) {
	  if ( (running_entries+n12)>max_n_ovdi ) {
	    get_positions=0;
	    errorcode--;
	    n12_stop=running_entries;
	    printf("running_entries=%d\n",running_entries);
	    printf("n12=%d\n",n12);
	    printf("setting get_positions to 0\n");
	  }
	  for (i=1;i<=n12;i++) {
	    running_entries++;
	    ovdi[running_entries][1]=transcript;
	    ovdi[running_entries][2]=ovdi_single_transcript[i][1];        // position of motif 1
	    ovdi[running_entries][3]=ovdi_single_transcript[i][2];        // position of motif 2
	    ovdi[running_entries][4]=ovdi_single_transcript[i][3];        // index of motif 1 (within the positions in m1p
	    ovdi[running_entries][5]=ovdi_single_transcript[i][4];        // index of motif 2 (within the positions in m2p
	    ovdi[running_entries][6]=ovdi_single_transcript[i][5];        // strand of motif 1
	    ovdi[running_entries][7]=ovdi_single_transcript[i][6];        // strand of motif 2
	  }
	} else {
	  running_entries+=n12;
	}
      }   // close if (curr_tr2==transcript)
    }     // close if (curr_tr1==transcript)
  }       // while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n12_tr;
  n[2]=errorcode;
  if (errorcode<0) {
    if (n12_stop<=0)
      n12_stop=running_entries;
    n[3]=n12_stop;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_ivector(m1p,1,max_transcript_occurrences);     
  free_cvector(m1s,1,max_transcript_occurrences);     
  free_ivector(m2p,1,max_transcript_occurrences);     
  free_cvector(m2s,1,max_transcript_occurrences);     
  free_imatrix(ovdi_single_transcript,1,max_n_ovdi_single_transcript,1,6);
}

void dist_3v_unordered_v3(int **ovdi12_p,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int verbose)
{
  /* co-occurrence of 3 motifs without order constraint
   * note that this is different from v1 in that here we constraint the 3 motifs to be within a range specified by max_dist and do not check on the individual distances
   * this is for a single transcript
   * 
   * ml (1 x 3) contains the motif lengths for each of the 3 motifs involved
   * if get_positions==1, then fill ovdi with positions
   * outputs:
   *   ovdi (n_ovdi x 9) p1 p2 p3 i j k s1 s2 s3  (where i,j,k correspond to the indices in m1_p, m2_p and m3_p for the occurrences of motif 1, 2 and 3 respectively)
   *   n[0]=number of co-occurrences within the current transcript
   *   n[1]=error code (0 if there are no errors, -1 if there is a memory overload
   *
   * here order does not matter
   * in contrast to v1, here we loop only once through the data and directly get the co-occurring events
   *
   * 09-01-03: assume ordered positions and consider d13 and d23 to avoid unnecessary visits to high j indices
   * 01-12-04: add strand information (new in version 3)
   * 01-12-04: ovdi now has 9 columns
   */

  int i,j;                          // loop counters
  int p1,p2,p3;                     // positions
  char s1,s2,s3;                    // strands
  int ad13,ad23;                    // absolute values of distances between motifs
  int d13,d23;                      // distances between motifs
  int index1,index2;                // indices
  int d;                            // distance
  int n_ovdi=0;
  int **all_ovdi;
  int n_cooc;
  float r;

  if (get_positions==1)
    all_ovdi=imatrix(1,max_n_ovdi,1,9);

  if (verbose)
    printf("\tdist_3v_unordered_v3\tmax_dist=%d\tmin_dist=%d\tn12=%d\tn3=%d\tmax_n_ovdi=%d\n",max_dist,min_dist,n12,n3,max_n_ovdi);

  i=0;
  while ( (i<n12) & (n_ovdi<max_n_ovdi) ) {
    i++;
    j=0;
    p1=ovdi12_p[i][1];
    p2=ovdi12_p[i][2];
    s1=ovdi12_p[i][5];
    s2=ovdi12_p[i][6];
    d13=0;                                                   // this has to be here also, otherwise we do not enter the loop!
    d23=0;                                                   // this has to be here also, otherwise we do not enter the loop!
    if (verbose)
      printf("\ni=%d\tp1=%d\tp2=%d\ts1=%d\ts2=%d\tn_ovdi=%d\n",i,p1,p2,s1,s2,n_ovdi);
    while ( (j<n3) & (n_ovdi<max_n_ovdi) & (d13<max_dist) & (d23<max_dist) ) {
      j++;
      p3=m3_p[j];
      s3=m3_s[j];
      if (verbose) 
	printf("\tj=%d\tp3=%d\ts3=%d\t",j,p3,s3);
      d13=0;                                                 // several days to realize that this should be in here!
      d23=0;
      if ( ( (s2==1) & (s3==1) ) | ( (s2==0) & (s3==0) ) ) 
	d13=p3-p1;
      ad13=fabs(d13);
      if (verbose)
	printf("d13=%d\tad13=%d\t",d13,ad13);
      if ( (ad13>min_dist) & (ad13<max_dist) ) {
	d23=p3-p2;
	ad23=fabs(d23);
	if (verbose)
	  printf("d23=%d\tad23=%d\n",d23,ad23);
	if ( (ad23>min_dist) & (ad23<max_dist) ) {
	  n_ovdi++;
	  if (get_positions==1) {
	    all_ovdi[n_ovdi][1]=p1;
	    all_ovdi[n_ovdi][2]=p2;
	    all_ovdi[n_ovdi][3]=p3;
	    all_ovdi[n_ovdi][4]=ovdi12_p[i][3];
	    all_ovdi[n_ovdi][5]=ovdi12_p[i][4];
	    all_ovdi[n_ovdi][6]=j;
	    all_ovdi[n_ovdi][7]=s1;
	    all_ovdi[n_ovdi][8]=s2;
	    all_ovdi[n_ovdi][9]=s3;
	  }            // close check on get_positions==1
	}              // close check on ad23>min_dist and ad23<max_dist
      }                // close check on ad13>min_dist and ad23<max_dist
      if ( (d13>max_dist) | (d23>max_dist) ) 
	j=n3+1;        // we are past max_dist, subsequent entries won't make it
    }                  // close j loop
  }                    // close i loop

  if ( (n_ovdi>0) & (get_positions==1) ) {
    if (verbose) 
      printf("\tn_ovdi=%d\tcalling remove_overlaps_v2\n",n_ovdi);
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n_ovdi,3);
    if (verbose)
      printf("\tn_cooc=%d\n",n_cooc);
  } else {
    n_cooc=n_ovdi;
  }

  if (i<n12) {
    n[1]=-1;
    if (verbose)
      printf("dist_3v_unordered_v3\ti=%d\tn12=%d\tn_ovdi=%d\tmax_n_ovdi=%d\n",i,n12,n_ovdi,max_n_ovdi);
  } else {
    n[1]=0;
  }
  n[0]=n_cooc;      // 0 by default

  if (get_positions==1)
    free_imatrix(all_ovdi,1,max_n_ovdi,1,9);
}

void dist_3v_ordered_v3(int **ovdi12_p,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format,int verbose)
{
  /* given the positions of three motifs in a single transcript, return the number of co-occurrences and the positions
   * 
   * in contrast to v1, here we try to avoid a call to dist_2v_ordered by retrieving the results directly from a previous run of dist_2v_ordered
   * note that the inputs therefore change:
   * ovdi12_p contains the positions of n12 entries for motif 1 and motif 2 (we keep motif 1 here to write them all into ovdi)
   * new in version 3: 
   * 01-02-2004: here we also use the strand information
   * 01-02-2004: note that ovdi12_p now contains 6 columns (where the last two are strand information) and ovdi now contains 9 columns (where the last 3 are strand)
   * 01-31-2004: added pos_format (see dist_2v_ordered_v2 for explanation of this parameter)
   * 02-02-2004: added verbose
   */

  int *m12_p;        // positions of co-occurrences of motifs 1 and 2
  char *m12_s;        // strands for the co-occurrences of motifs 1 and 2
  int n123;          // number of co-occurrences of motifs 1,2,3 before call to remove_overlaps
  int n_123[2];      // n_123[0]=n123, n_123[1]=errorcode
  int **ovdi123;     // positions and indices of co-occurring motifs before all to remove_overlaps
  int i;             // loop index
  int curr_index;
  int n_cooc;
  int **all_ovdi;

  m12_p=ivector(1,n12);
  m12_s=cvector(1,n12);
  ovdi123=imatrix(1,max_n_ovdi,1,6);           // position of motif 2, position of motif 3, index of motif 2 in m12_p, index of motif 3 in m3_p, strand2, strand3
  all_ovdi=imatrix(1,max_n_ovdi,1,9);          // motif positions, indices, strands

  for (i=1;i<=n12;i++) {
    m12_p[i]=ovdi12_p[i][2];
    m12_s[i]=ovdi12_p[i][6];
  }
  if (verbose) {
    printf("dist_3v_ordered_v3: motif 2 (%d): ",n12);
    for (i=1;i<=n12;i++) 
      printf("%d %d,",m12_p[i],m12_s[i]);
    printf("\n");
    printf("dist_3v_ordered_v3: motif 3 (%d): ",n3);
    for (i=1;i<=n3;i++) 
      printf("%d %d,",m3_p[i],m3_s[i]);
    printf("\n");
  }
  dist_2v_ordered_v2(m12_p,m3_p,m12_s,m3_s,n12,n3,ovdi123,n_123,min_dist,max_dist,max_n_ovdi,ml,get_positions,pos_format,verbose);
  n123=n_123[0];
  if (verbose)
    printf("\tn123=%d\n",n123);
  
  if (get_positions==1) {
    for (i=1;i<=n123;i++) {
      curr_index=ovdi123[i][3];                   // index in m12_p, corresponding to this index in ovdi12
      all_ovdi[i][1]=ovdi12_p[curr_index][1];     // position of motif 1
      all_ovdi[i][2]=ovdi123[i][1];               // position of motif 2
      all_ovdi[i][3]=ovdi123[i][2];               // position of motif 3
      all_ovdi[i][4]=ovdi12_p[curr_index][3];     // index in m1_p
      all_ovdi[i][5]=ovdi12_p[curr_index][4];     // index in m2_p
      all_ovdi[i][6]=ovdi123[i][4];               // index in m3_p
      all_ovdi[i][7]=ovdi12_p[curr_index][5];     // strand of motif 1
      all_ovdi[i][8]=ovdi123[i][5];               // strand of motif 2
      all_ovdi[i][9]=ovdi123[i][6];               // strand of motif 3
    }
  }

  if ( (n123>0) & (get_positions==1) ) {
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n123,3);
  } else {
    n_cooc=n123;
  }
  
  n[0]=n_cooc;
  n[1]=n_123[1];

  free_ivector(m12_p,1,n12);
  free_cvector(m12_s,1,n12);
  free_imatrix(ovdi123,1,max_n_ovdi,1,6);
  free_imatrix(all_ovdi,1,max_n_ovdi,1,9);
}

void dist_3v_allt_v3(int **ovdi12,int *m3_t,int *m3_p,char *m3_s,int n12,int n3,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose)
{
  /* call dist_3v_unordered_v3 or dist_3v_ordered_v3 for all transcripts 
   * get_positions = 1  --> to fill ovdi with the positions
   * order_constraint=1 --> order matters (call dist_3v_ordered)
   *
   * outputs: 
   *      ovdi (running_entres x 10)     transcript, 3 positions, 3 indices, 3 strands
   *      n    (n[0]=running_entries;n[1]=n123_tr;n[2]=errorcode;
   *
   * 09-13-2003: renamed max_n_ovdi to max_n_ovdi_single_transcript
   * 09-13-2003: added max_n_ovdi with the maximum allowed number of entries into ovdi
   * 07-24-2003: added max_n_transcripts, set max_n_transcripts to -1 if you do not wish to use this constraint
   * 01-02-2004: new in version 3, added strand information
   * 01-12-2004: call dist_3v_unordered_v3
   * 01-31-2004: added pos_format (see dist_2v_ordered_v2 for explanation)
   * 02-02-2004: added verbose
   */

  int i12,i3;
  int **ovdi_single_transcript;
  int running_entries=0;
  int errorcode=0;
  int n123_tr=0;                                 // number of transcripts with co-occurrences of all three motifs
  int transcript=0;
  int curr_tr12,curr_tr3;
  int **ovdi12_p;                                // positions (and indices) of cooccurring motifs 1 and 2 in running transcript (memory = 10kx6colsx2 bytes/int = 120 kb)
  int *m3p;                                      // positions of motif 3 in running transcript (memory = 10000 x 1 col x 2 bytes/int = 20 kb)
  char *m3s;                                     // strands for motif 3 in running transcript
  int max_transcript_occurrences=10000;          // maximum number of occurrences per transcript
  int max_n_ovdi_single_transcript=50000;        // maximum number of co-occurrences in current transcript
  int n123;                                      // total number of co-occurrences of all three motifs
  int n123_stop=0;                               // number of co-occurrences when it stops reporting because of an overflow
  int i,j;
  int curr_n12,curr_n3;

  if (verbose) 
    printf("dist_3v_allt_v3\tn12=%d\tn3=%d\tmin_dist=%d\tmax_dist=%d\tn_transcripts=%d\tmax_n_transcripts=%d\tmax_n_ovdi=%d\tpos_format=%d\n",n12,n3,min_dist,max_dist,n_transcripts,max_n_transcripts,max_n_ovdi,pos_format);

  i12=1;curr_tr12=ovdi12[i12][1];                  // initial transcript 1-2
  i3=1;curr_tr3=m3_t[i3];                          // initial transcript 3

  ovdi12_p=imatrix(1,max_n_ovdi_single_transcript,1,6);                    // positions, indices and strands for motif 1 and motif 2 in current transcript
  m3p=ivector(1,max_transcript_occurrences);                               // will hold the positions of motif 3 in current transcript
  m3s=cvector(1,max_transcript_occurrences);                               // strands for motif 3 in current transcript
  ovdi_single_transcript=imatrix(1,max_n_ovdi_single_transcript,1,9);      // output of dist_3v (positions, indices, strands)

  if (max_n_transcripts<0) 
    max_n_transcripts=n_transcripts+1;

  while ( (transcript < n_transcripts) & (n123_tr<=max_n_transcripts) ) {
    transcript++;    
    if ( (curr_tr12<0) | (curr_tr3<0) ) 
      transcript=n_transcripts+1;                 // at least one of the motifs is done, let's call it a day
    while ( (curr_tr12<transcript) & (curr_tr12>=0) & (i12<n12) ) {
      i12++;curr_tr12=ovdi12[i12][1];                     // get occurrences of motif 1 up to current transcript
    }
    if (curr_tr12==transcript) {                   // if motif 1 occurs in current transcript
      curr_n12=0;
      while ( (curr_tr12==transcript) & (i12<=n12) ) {
	curr_n12++;
	ovdi12_p[curr_n12][1]=ovdi12[i12][2];
	ovdi12_p[curr_n12][2]=ovdi12[i12][3];
	ovdi12_p[curr_n12][3]=ovdi12[i12][4];
	ovdi12_p[curr_n12][4]=ovdi12[i12][5];
	ovdi12_p[curr_n12][5]=ovdi12[i12][6];
	ovdi12_p[curr_n12][6]=ovdi12[i12][7];
	i12++;
	if (i12<=n12)
	  curr_tr12=ovdi12[i12][1];
      }
      i12--;                                       // so that we can then search for i1 in the above part of the code
      curr_tr12=transcript;
      
      while ( (curr_tr3<transcript) & (curr_tr3>=0) & (i3<n3) ) {
	i3++;curr_tr3=m3_t[i3];
      }      
      if (curr_tr3==transcript) {
	curr_n3=0;
	while ( (curr_tr3==transcript) & (i3<=n3) ) {
	  curr_n3++;
	  m3p[curr_n3]=m3_p[i3];
	  m3s[curr_n3]=m3_s[i3];
	  i3++;
	  if (i3<=n3) 
	    curr_tr3=m3_t[i3];
	}
	i3--;
	curr_tr3=transcript;

	if (verbose) {
	  printf("\ntranscript=%d\n",transcript);
	  printf("moitf 3 (%d): ",curr_n3);
	  for (i=1;i<=curr_n3;i++)
	    printf("%d %d,",m3p[i],m3s[i]);
	  printf("\n");
	}

	if (order_constraint==0) {
	  dist_3v_unordered_v3(ovdi12_p,m3p,m3s,curr_n12,curr_n3,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions,verbose);
	} else {
	  dist_3v_ordered_v3(ovdi12_p,m3p,m3s,curr_n12,curr_n3,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions,pos_format,verbose);
	}
	n123=n[0];
	errorcode+=n[1];            // n[1]<0 if there are errors	
	if (verbose) 
	  printf("\tn123=%d\tn[1]=%d\terrorcode=%d\trunning_entries=%d\n",n123,n[1],errorcode,running_entries);

	if (n123>0)
	  n123_tr++;
	if (get_positions==1) {
	  if ( (running_entries+n123)>max_n_ovdi ) {
	    get_positions=0;
	    errorcode--;
	    printf("running_entries=%d\n",running_entries);
	    printf("n123=%d\n",n123);
	    printf("setting get_positions to 0\n");
	    n123_stop=running_entries;
	  }
	  //for (i=1;i<=n123;i++) {
	  i=0;
	  while ( (i<n123) & (running_entries<max_n_ovdi) ) {
	    i++;
	    running_entries++;
	    ovdi[running_entries][1]=transcript;
	    ovdi[running_entries][2]=ovdi_single_transcript[i][1];      // position of motif 1
	    ovdi[running_entries][3]=ovdi_single_transcript[i][2];      // position of motif 2
	    ovdi[running_entries][4]=ovdi_single_transcript[i][3];      // position of motif 3
	    ovdi[running_entries][5]=ovdi_single_transcript[i][4];      // index of motif 1
	    ovdi[running_entries][6]=ovdi_single_transcript[i][5];      // index of motif 2
	    ovdi[running_entries][7]=ovdi_single_transcript[i][6];      // index of motif 3
	    ovdi[running_entries][8]=ovdi_single_transcript[i][7];      // strand of motif 1
	    ovdi[running_entries][9]=ovdi_single_transcript[i][8];      // strand of motif 2
	    ovdi[running_entries][10]=ovdi_single_transcript[i][9];     // strand of motif 3
	  }
	} else {
	  running_entries+=n123;
	}      // close if (get_posisions==1)
      }        // close if (curr_tr3==transcript)
    }            // close if (curr_tr12==transcript)
  }              // while (transcript < n_transcripts)
  
  n[0]=running_entries;
  n[1]=n123_tr;
  n[2]=errorcode;
  if (errorcode<0) {
    if (n123_stop<=0)                       // modified on 02-18-2004. before it said n123_stop<running_entries but this did not work when there was an overflow
      n123_stop=running_entries;
    n[3]=n123_stop;
  } else {
    n[3]=running_entries;
  }
  /* free memory */
  free_imatrix(ovdi12_p,1,max_n_ovdi_single_transcript,1,6);
  free_ivector(m3p,1,max_transcript_occurrences);     
  free_cvector(m3s,1,max_transcript_occurrences);
  free_imatrix(ovdi_single_transcript,1,max_n_ovdi_single_transcript,1,9);
}

void dist_3v_allt_from_m123_v1(int *m1_t,int *m2_t,int *m3_t,int *m1_p,int *m2_p,int *m3_p,char *m1_s,char *m2_s,char *m3_s,int n1,int n2,int n3,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose)
{
  /* compute cooccurring instances of 3 motifs
   * get_positions = 1    --> fill ovdi with the positions
   * order_constraint = 1 --> order matters
   *
   * note that in contrast to dist_3v_allt_v2 here we do not have the information from the cooccurrences of motif 1 and motif 2 to start with
   * 
   * this function comes from dist_3v_allt
   * 01-03-2004: we change the function name to further distinguish this from dist_3v_allt_v3
   * 01-03-2004: here we include the strand information
   * 01-03-2004: call dist_2v_allt_v2
   * 01-03-2004: call dist_3v_allt_v3
   * 01-31-2004: added pos_format
   * 02-02-2004: added debug mode verbose
   */

  int **ovdi12;                             // co-occurrences of motif 1 and motif 2
  int max_n_ovdi12=500000;                  // maximum number of co-occurrences of motif 1 and motif 2 (memory = 500000x5(cols)x2(bytes/int)=5000kb)
  int nout12[4];                            // [0]=n12, [1]=n12_tr, [2]=errorcode, [3]=n12_stop 
  int nout123[4];                           // [0]=n123, [1]=n123_tr, [3]=errorcode, [4]=n123_stop
  int n12;                                  // number of co-occurrences between motif 1 and motif 2 
  int get_positions12=1;                    // 1 to get the positions of co-occurring instances of motif 1 and motif 2 (this should always be 1 to go for 3 motifs)
  int max_n_transcripts12=-1;               // maximum number of transcripts (if < 0, then do not limit the computations)
  int errorcode=0;
  int i,j;

  ovdi12=imatrix(1,max_n_ovdi12,1,7);       // transcript, 2 positions, 2 indices, 2 strands

  /* start by computing the co-occurrences of motif 1 and motif 2 */
  if (verbose) {
    printf("motif 1 (%d): ",n1);
    for (i=1;i<=n1;i++) 
      printf("%d %d %d,",m1_t[i],m1_p[i],m1_s[i]);
    printf("\n");
    printf("motif 2 (%d): ",n2);
    for (i=1;i<=n2;i++) 
      printf("%d %d %d,",m2_t[i],m2_p[i],m2_s[i]);
    printf("\n");
    printf("motif 3 (%d): ",n3);
    for (i=1;i<=n3;i++) 
      printf("%d %d %d,",m3_t[i],m3_p[i],m3_s[i]);
    printf("\n");
  }

  dist_2v_allt_v2(m1_t,m2_t,m1_p,m2_p,m1_s,m2_s,n1,n2,ovdi12,nout12,min_dist,max_dist,n_transcripts,ml,get_positions12,order_constraint,max_n_transcripts12,max_n_ovdi12,pos_format,0);
  errorcode+=nout12[2];       // nout12[2]<0 if there are errors
  if (errorcode<0) {
    n12=nout12[3];            // when it stopped reporting positions because of overflow
  } else {
    n12=nout12[0];            // total number of co-occurrences (total number of rows in ovdi12)
  }
  if (verbose) {
    printf("ovdi12 (%d)\n",n12);
    for (i=1;i<=n12;i++) 
      printf("%d %d %d %d %d\n",ovdi12[i][1],ovdi12[i][2],ovdi12[i][3],ovdi12[i][6],ovdi12[i][7]);
  }

  /* now call dist_3v_allt_v2 */
  dist_3v_allt_v3(ovdi12,m3_t,m3_p,m3_s,n12,n3,ovdi,nout123,min_dist,max_dist,n_transcripts,ml,get_positions,order_constraint,max_n_transcripts,max_n_ovdi,pos_format,verbose);
  n[0]=nout123[0];
  n[1]=nout123[1];
  errorcode+=nout123[2];          // nout123[2]<0 if there are errors
  n[2]=errorcode;
  n[3]=nout123[3];
  if (verbose) {
    printf("ovdi123 (%d)\n",n[0]);
    for (i=1;i<=n[0];i++) 
      printf("%d %d %d %d %d %d %d\n",ovdi[i][1],ovdi[i][2],ovdi[i][3],ovdi[i][4],ovdi[i][8],ovdi[i][9],ovdi[i][10]);
  }

  free_imatrix(ovdi12,1,max_n_ovdi12,1,7);  
}

void dist_4v_ordered_v3(int **ovdi123_p,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions,int pos_format)
{
  /* given the positions of four motifs in a single transcript, return the number of co-occurrences and the positions
   * 
   * in contrast to v1, here we try to avoid a call to dist_2v_ordered by retrieving the results directly from a previous run of dist_2v_ordered
   * note that the inputs therefore change:
   * ovdi123_p contains the positions of n123 entries for motif 1, motif 2 and motif 3 (we keep motif 1 and 2 here to write them all into ovdi)
   *
   * 01-02-2004: new in version 3, use strand information
   * 01-31-2004: added pos_format (see dist_2v_ordered_v2 for explanation)
   */

  int *m123_p;        // positions for co-occurrences of motifs 1, 2 and 3
  char *m123_s;       // strands for co-occurrences of motifs 1, 2 and 3
  int n1234;          // number of co-occurrences of all 4 motifs       
  int n_1234[2];      // used in call to dist_2v_ordered to retrived the number of co-occurrences
  int **ovdi1234;     // used in call to dist_2v_orderedd to retrieve the positions of co-occurring motifs 
  int i,j;            // loop variable
  int curr_index;     // index
  int n_cooc;         // number of co-occurrences
  int **all_ovdi;     // all co-occurrences

  m123_p=ivector(1,n123);
  m123_s=cvector(1,n123);
  ovdi1234=imatrix(1,max_n_ovdi,1,6);   // position of motif 3, position of motif 4, indices of motifs 3 and 4, strands of motifs 3 and 4
  all_ovdi=imatrix(1,max_n_ovdi,1,12);  // 4 positions, 4 indices, 4 strands

  for (i=1;i<=n123;i++) {
    m123_p[i]=ovdi123_p[i][3];
    m123_s[i]=ovdi123_p[i][9];
  }
  dist_2v_ordered_v2(m123_p,m4_p,m123_s,m4_s,n123,n4,ovdi1234,n_1234,min_dist,max_dist,max_n_ovdi,ml,get_positions,pos_format,0);
  n1234=n_1234[0];
  if (get_positions==1) {
    for (i=1;i<=n1234;i++) {
      curr_index=ovdi1234[i][3];                   // index in m12_p, corresponding to this index in ovdi12
      all_ovdi[i][1]=ovdi123_p[curr_index][1];     // position of motif 1
      all_ovdi[i][2]=ovdi123_p[curr_index][2];     // position of motif 2
      all_ovdi[i][3]=ovdi123_p[curr_index][3];     // position of motif 3
      all_ovdi[i][4]=ovdi1234[i][2];               // position of motif 4
      all_ovdi[i][5]=ovdi123_p[curr_index][4];     // index of motif 1
      all_ovdi[i][6]=ovdi123_p[curr_index][5];     // index of motif 2
      all_ovdi[i][7]=ovdi123_p[curr_index][6];     // index of motif 3
      all_ovdi[i][8]=ovdi1234[i][4];               // index of motif 4
      all_ovdi[i][9]=ovdi123_p[curr_index][7];     // strand of motif 1
      all_ovdi[i][10]=ovdi123_p[curr_index][8];    // strand of motif 2
      all_ovdi[i][11]=ovdi123_p[curr_index][9];    // strand of motif 3
      all_ovdi[i][12]=ovdi1234[i][6];              // strand of motif 4
    }
  }
  
  if ( (n1234>0) & (get_positions==1) ) {
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n1234,4);
  } else {
    n_cooc=n1234;
  }
  
  n[0]=n_cooc;
  n[1]=n_1234[1];

  free_ivector(m123_p,1,n123);
  free_cvector(m123_s,1,n123);
  free_imatrix(ovdi1234,1,max_n_ovdi,1,4);
  free_imatrix(all_ovdi,1,max_n_ovdi,1,8);
}

void dist_4v_unordered_v3(int **ovdi123_p,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int max_n_ovdi,int *ml,int get_positions)
{
  /* co-occurrence of 4 motifs without order constraint
   * note that this is different from v1 in that here we constraint the 3 motifs to be within a range specified by max_dist and do not check on the individual distances
  /* this is for a single transcript
   * 
   * if get_positions==1, then fill ovdi with positions
   *  ovdi (n_ovdi x 12) p1 p2 p3 p4 i j k l s1 s2 s3 s4  (where i,j,k,l correspond to the indices in m1/2/3/4_p for the occurrences of motif 1/2/3/4 respectively)
   * 
   * n[0]=number of co-occurrences within the current transcript
   * n[1]=error code (0 if there are no errors, -1 if there is a memory overload
   *
   * ovdi (n_ovdi x 12):  position1  position2 position3 position4 i  j k l strand1 strand2 strand3 strand4
   *
   * here order does not matter
   * in contrast to v1, here we loop only once through the data and directly get the co-occurring events
   *
   * 09-01-03: assume ordered positions and consider d13,d24,d34 to avoid visiting unnecessary high j indices
   * 01-12-04: add strand information (new in version 3)
   * 01-12-04: ovdi now has 12 columns  
   */

  int i,j;                             // loop variables
  int p1,p2,p3,p4;                     // positions of the motifs
  char s1,s2,s3,s4;                    // strands of the motifs
  int d14,d24,d34;                     // distances between motifs
  int ad14,ad24,ad34;                  // absolute values of the distances between motifs
  int n_ovdi=0;                        // number of co-occurring entries (before call to remove_overlaps)
  int **all_ovdi;                      // positions of co-occurring entries
  int n_cooc;                          // number of co-occurring entries

  if (get_positions==1)
    all_ovdi=imatrix(1,max_n_ovdi,1,12);

  i=0;
  while ( (i<n123) & (n_ovdi<max_n_ovdi) ) {
    i++;
    j=0;
    p1=ovdi123_p[i][1];s1=ovdi123_p[i][7];
    p2=ovdi123_p[i][2];s2=ovdi123_p[i][8];
    p3=ovdi123_p[i][3];s3=ovdi123_p[i][9];
    d14=0;d24=0;d34=0;                   // this has to be here also; otherwise we do not enter the loop!
    while ( (j<n4) & (n_ovdi<max_n_ovdi) & (d14<max_dist) & (d24<max_dist) & (d34<max_dist) ) {
      j++;
      p4=m4_p[j];
      s4=m4_s[j];
      d14=0;                             // several days to realize that this should be in here!
      d24=0;
      d34=0;
      if ( ( (s3==1) & (s4==1) ) | ( (s3==0) & (s4==0) ) )
	d14=p4-p1;
      ad14=fabs(d14);
      if ( (ad14>min_dist) & (ad14<max_dist) ) {
	d24=p4-p2;ad24=fabs(d24);
	if ( (ad24>min_dist) & (ad24<max_dist) ) {
	  d34=p4-p3;ad34=fabs(d34);
	  if ( (ad34>min_dist) & (ad34<max_dist) ) {
	    n_ovdi++;
	    if (get_positions==1) {
	      all_ovdi[n_ovdi][1]=p1;
	      all_ovdi[n_ovdi][2]=p2;
	      all_ovdi[n_ovdi][3]=p3;
	      all_ovdi[n_ovdi][4]=p4;
	      all_ovdi[n_ovdi][5]=ovdi123_p[i][4];
	      all_ovdi[n_ovdi][6]=ovdi123_p[i][5];
	      all_ovdi[n_ovdi][7]=ovdi123_p[i][6];
	      all_ovdi[n_ovdi][8]=j;
	      all_ovdi[n_ovdi][9]=s1;
	      all_ovdi[n_ovdi][10]=s2;
	      all_ovdi[n_ovdi][11]=s3;
	      all_ovdi[n_ovdi][12]=s4;
	    }
	  }   // close check on ad34
	}     // close check on ad24
      }       // close check on ad14
      if ( (d14>max_dist) | (d24>max_dist) | (d34>max_dist) )
	j=n4+1;
    }         // close j loop
  }           // close i loop

  if ( (n_ovdi>0) & (get_positions==1) ) {
    n_cooc=remove_overlaps_v2(all_ovdi,ovdi,ml,n_ovdi,4);
  } else {
    n_cooc=n_ovdi;
  }

  if (i<n123) {
    n[1]=-1;
  } else {
    n[1]=0;
  }
  n[0]=n_cooc;      // 0 by default

  //printf("n_cooc=%d\n",n_cooc);
  if (get_positions==1)
    free_imatrix(all_ovdi,1,max_n_ovdi,1,12);
}

void read_header(FILE *input_file,int *t,int *data_p,int max_occur)
{
  /* given a file handle, read the header (lines starting with ">") as well as the first entry
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   * return the transcript in t[0], the number of positions in t[1] and the positions in data_p
   */
  char txtline[100000];    // read one line of the input file
  char searchchar;        // search whether it is a comment line
  char searchchar2;       // allow also matlab comments!
  char *gotit;            // used in the search process
  int header=1;           // 1 while reading the header
  int i,j,k,ll;           // loop counters
  int *v;                 // vector holding the transcripts and positions for running entry
  int start=0;            // start position for splitline2int
  int np[1];              // number of entries
  int n;                  // number of entries
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int neg2pos=1;
  
  if (feof(input_file)!=0 ) {
    /* nothing else to read */
    t[0]=-1;
    t[1]=-1;
    t[2]=0;
  } else {
    searchchar='>';
    searchchar2='\%';
  
    v=ivector(1,max_occur);
    
    while (header) {
      fgets(txtline,n_chars,input_file);
      gotit=strchr(txtline,searchchar);
      if (!gotit) {
	gotit=strchr(txtline,searchchar2);
	if (!gotit) 
	  header=0;
      } 
    }
    
    splitline2int(txtline,v,start,np,max_occur);
    n=np[0]-1;

    if (n>max_occur) {
      n=max_occur;
      t[2]=-1;
    } else {
      t[2]=0;
    }
    if (n>0) {
      t[0]=v[1];
      t[1]=n-1;
      t[2]=0;
      for (i=2;i<=n;i++) {
	j=i-1;
	if (neg2pos==1) {
	  if (v[i]>0) {
	    data_p[j]=v[i];
	  } else {
	    data_p[j]=-v[i];
	  }
	} else {
	  data_p[j]=v[i];
	}
      }
    } else {
      t[0]=-1;
      t[1]=-1;
      t[2]=0;
    }
    /* printf("n=%d\n",n);
       for (i=i;i<=n;i++)
       printf("%d ",data_p[i]);
       printf("\n"); */
    free_ivector(v,1,max_occur);
  }
}

void read_line(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos)
{
  /* given a file handle, read a single line
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   * return the transcript in t[0], the number of positions in t[1] and the positions in data_p
   * t[2]=-1 -> end of file   t[2]=0 -> n>max_occur  t[2]=1 -> ok  t[2]=2 -> n=0 (no positions) t[2]=3 -> header line
   * 05-31-2003: read files created with rc=1 (negative positions indicating the opposite strand
   *
   * if neg2pos==1, convert negative positions to positive ones 
   * if sortpos==1, sort the positions
   * if the end of file is reached, returns all t entries as -1
   */

  char txtline[100000];    // read one line of the input file
  int i,j,k,ll;           // loop counters
  int *v;                 // vector holding the transcripts and positions for running entry
  int start=0;            // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];              // number of entries
  int n;                  // number of entries
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  float *d;               // dummy variable used in sorting
  int header;

  //printf("max_occur=%d\n",max_occur);
  printf(".");
  if (feof(input_file)!=0) {
    /* nothing else to read in this file */
    t[0]=-1;
    t[1]=-1;
    t[2]=-1;
  } else {
    fgets(txtline,n_chars,input_file);
    //printf("processing txtline=%s\n",txtline);
    header=is_header(txtline);

    if (header==1) {
      t[0]=-1;
      t[1]=-1;
      t[2]=3;
    } else {
      v=ivector(1,max_occur);
      //printf("calling splitline2int with txtline=%s\n",txtline);
      splitline2int(txtline,v,start,np,max_occur);
      n=np[0]-1;
      if (n>max_occur) {
	n=max_occur;
	t[2]=0;
      }
      if (n>0) {
	t[0]=v[1];
	t[1]=n-1;
	t[2]=1;
	for (i=2;i<=n;i++) {
	  j=i-1;
	  if (neg2pos==1) {
	    if (v[i]>0) {
	      data_p[j]=v[i];
	    } else {
	      data_p[j]=-v[i];
	    }
	  } else {
	    data_p[j]=v[i];
	  }
	}
      } else {
	t[0]=-1;
	t[1]=-1;
	t[2]=2;
      }
      
      if ( (n>1) & (sortpos==1) ) {
	n--;
	
	/* sort the positions */
	d=vector(1,n);
	quicksorti(n,data_p,d);
	free_vector(d,1,n);
      }
      free_ivector(v,1,max_occur);
    }  // close if (is_header(txtline)==0)     
  }    // close feof check
}

void dist_4v_allt_v3(int **ovdi123,int *m4_t,int *m4_p,char *m4_s,int n123,int n4,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format)
{
  /* call dist_4v_unordered_v3 or dist_4v_ordered_v3 for all transcripts 
   * get_positions = 1  --> to fill ovdi with the positions
   * order_constraint=1 --> order matters (call dist_3v_ordered)
   *
   * 09-13-2003: renamed max_n_ovdi to max_n_ovdi_single_transcript
   * 09-13-2003: added max_n_ovdi with the maximum allowed number of entries into ovdi
   * 07-24-2003: added max_n_transcripts; set max_n_transcripts to -1 if you do not wish to use this constraint
   * 01-02-2004: (new in version 3) use strand information
   * 01-31-2004: added pos_format (see dist_2v_ordered_v2 for explanation)
   */

  int i123,i4;
  int **ovdi_single_transcript;
  int running_entries=0;
  int errorcode=0;
  int n1234_tr=0;                                 // number of transcripts with co-occurrences of all three motifs
  int transcript=0;
  int curr_tr123,curr_tr4;
  int *m4p;
  char *m4s;
  int **ovdi123_p;
  int max_transcript_occurrences=10000;
  int max_n_ovdi_single_transcript=50000;         // maximum number of co-occurrences in current transcript
  int n1234;                                      // total number of co-occurrences of all three motifs
  int n1234_stop;
  int i,j;
  int curr_n123,curr_n4;

  i123=1;curr_tr123=ovdi123[i123][1];                   // initial transcript 1-2-3
  i4=1;curr_tr4=m4_t[i4];                               // initial transcript 3

  ovdi123_p=imatrix(1,max_n_ovdi_single_transcript,1,9);  // positions, indices and strands for motif 1 and motif 2 and 3 in current transcript
  m4p=ivector(1,max_transcript_occurrences);            // will hold the positions of motif 4 in current transcript
  m4s=cvector(1,max_transcript_occurrences);            // will hold the strands of motif 4 in current transcript
  ovdi_single_transcript=imatrix(1,max_n_ovdi_single_transcript,1,12);  // 4 positions, 4 indices, 4 strands for current transcript

  if (max_n_transcripts<0)
    max_n_transcripts=n_transcripts+1;
  while ( (transcript < n_transcripts) & (n1234_tr<=max_n_transcripts) ) {
    transcript++;    
    if ( (curr_tr123<0) | (curr_tr4<0) ) 
      transcript=n_transcripts+1;                 // at least one of the motifs is done, let's call it a day
    while ( (curr_tr123<transcript) & (curr_tr123>=0) & (i123<n123) ) {
      i123++;curr_tr123=ovdi123[i123][1];                     // get occurrences of motif 1 up to current transcript
    }
    if (curr_tr123==transcript) {                   // if motif 1 occurs in current transcript
      curr_n123=0;
      while ( (curr_tr123==transcript) & (i123<=n123) ) {
	curr_n123++;
	ovdi123_p[curr_n123][1]=ovdi123[i123][2];        // position of motif 1
	ovdi123_p[curr_n123][2]=ovdi123[i123][3];        // motif 2
	ovdi123_p[curr_n123][3]=ovdi123[i123][4];        // motif 3
	ovdi123_p[curr_n123][4]=ovdi123[i123][5];        // index of motif 1
	ovdi123_p[curr_n123][5]=ovdi123[i123][6];        // motif 2
	ovdi123_p[curr_n123][6]=ovdi123[i123][7];        // motif 3
	ovdi123_p[curr_n123][7]=ovdi123[i123][8];        // strand of motif 1
	ovdi123_p[curr_n123][8]=ovdi123[i123][9];        // strand of motif 2
	ovdi123_p[curr_n123][9]=ovdi123[i123][10];       // strand of motif 3
	i123++;
	if (i123<=n123)
	  curr_tr123=ovdi123[i123][1];
      }
      i123--;                                       // so that we can then search for i1 in the above part of the code
      curr_tr123=transcript;
      while ( (curr_tr4<transcript) & (curr_tr4>=0) & (i4<n4) ) {
	i4++;curr_tr4=m4_t[i4];
      }
      
      if (curr_tr4==transcript) {
	curr_n4=0;
	while ( (curr_tr4==transcript) & (i4<=n4) ) {
	  curr_n4++;
	  m4p[curr_n4]=m4_p[i4];
	  m4s[curr_n4]=m4_s[i4];
	  i4++;
	  if (i4<=n4) 
	    curr_tr4=m4_t[i4];
	}
	i4--;
	curr_tr4=transcript;
	
	if (order_constraint==0) {
	  dist_4v_unordered_v3(ovdi123_p,m4p,m4s,curr_n123,curr_n4,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions);
	} else {
	  dist_4v_ordered_v3(ovdi123_p,m4p,m4s,curr_n123,curr_n4,ovdi_single_transcript,n,min_dist,max_dist,max_n_ovdi_single_transcript,ml,get_positions,pos_format);
	}
	n1234=n[0];
	errorcode+=n[1];       // n[1]<0 if there are errors	
	if (n1234>0)
	  n1234_tr++;
 	if (get_positions==1) {
	  if ( (running_entries+n1234)>max_n_ovdi ) {
	    get_positions=0;
	    errorcode--;
	    printf("running_entries=%d\n",running_entries);
	    printf("n1234=%d\n",n1234);
	    printf("setting get_positions to 0\n");
	    n1234_stop=running_entries;
	  }
	  i=0;
	  while ( (i<n1234) & (running_entries<max_n_ovdi) ) {
	    i++;
	    //for (i=1;i<=n1234;i++) {
	    running_entries++;
	    ovdi[running_entries][1]=transcript;
	    ovdi[running_entries][2]=ovdi_single_transcript[i][1];           // positionf of motif 1
	    ovdi[running_entries][3]=ovdi_single_transcript[i][2]; 
	    ovdi[running_entries][4]=ovdi_single_transcript[i][3];
	    ovdi[running_entries][5]=ovdi_single_transcript[i][4];
	    ovdi[running_entries][6]=ovdi_single_transcript[i][5];           // index of motif 1
	    ovdi[running_entries][7]=ovdi_single_transcript[i][6]; 
	    ovdi[running_entries][8]=ovdi_single_transcript[i][7];
	    ovdi[running_entries][9]=ovdi_single_transcript[i][8];
	    ovdi[running_entries][10]=ovdi_single_transcript[i][9];           // strand of motif 1
	    ovdi[running_entries][11]=ovdi_single_transcript[i][10];          
	    ovdi[running_entries][12]=ovdi_single_transcript[i][11];
	    ovdi[running_entries][13]=ovdi_single_transcript[i][12];

	  }
	} else {
	  running_entries+=n1234;
	}      // close if (get_posisions==1)
      }        // close if (curr_tr4==transcript)
    }          // close if (curr_tr123==transcript)
  }              // while (transcript < n_transcripts)

  n[0]=running_entries;
  n[1]=n1234_tr;
  n[2]=errorcode;
  if (errorcode<0) {
    //if (n1234_stop<running_entries)
    if (n1234_stop<=0)
      n1234_stop=running_entries;
    n[3]=n1234_stop;
  } else {
    n[3]=running_entries;
  }

  /* free memory */
  free_imatrix(ovdi123_p,1,max_n_ovdi_single_transcript,1,9);
  free_ivector(m4p,1,max_transcript_occurrences);     
  free_cvector(m4s,1,max_transcript_occurrences);
  free_imatrix(ovdi_single_transcript,1,max_n_ovdi_single_transcript,1,12);
}

int remove_overlaps(int **m_old,int **m_new,int *ml,int n_entries,int n_motifs)
{
  /* given a matrix with the positions of each motif for multiple occurrences, remove overlaps 
   * (situations with the same motifs in "basically" the same positions)
   * returns the remaining number of entries
   * copies the indices entries in ovdi matrices
   */

  int issimilar;                     /* 1 if two entries are similar and 0 otherwise */
  int i,j,k;                           
  int running_n=1;                   /* running number of entries */
  int tempsimilar;                   /* temporary variable used to compare two entries */
 
  /* the first entry always stays */
  for (j=1;j<=(2*n_motifs);j++) 
    m_new[1][j]=m_old[1][j];  

  for (i=2;i<=n_entries;i++) {
    issimilar=0;                     /* by default, the two entries are not similar */
    k=1;                             /* compare entry number i with all the previous ones starting with k=1 and up to k=running_n or until we find a similarity */
    while ( (issimilar == 0) & (k<=running_n) ) {
      tempsimilar=0;	  
      for (j=1;j<=n_motifs;j++) {	    
	if ( (fabs(m_new[k][j]-m_old[i][j])) <= (ml[j]/2) )
	  tempsimilar++;
      }
      if (tempsimilar == n_motifs) 
	issimilar=1;
      k++;	  
    }	
    if (issimilar==0) {
      //new_n++;
      running_n++;
      for (j=1;j<=(2*n_motifs);j++) 
	m_new[running_n][j]=m_old[i][j];
    }
  }

  return running_n;
  //return new_n;
    
}

void dist_4v_allt_from_m1234_v1(int *m1_t,int *m2_t,int *m3_t,int *m4_t,int *m1_p,int *m2_p,int *m3_p,int *m4_p,char *m1_s,char *m2_s,char *m3_s,char *m4_s,int n1,int n2,int n3,int n4,int **ovdi,int *n,int min_dist,int max_dist,int n_transcripts,int *ml,int get_positions,int order_constraint,int max_n_transcripts,int max_n_ovdi,int pos_format,int verbose)
{
  /* compute cooccurring instances of 4 motifs
   * get_positions = 1    --> fill ovdi with the positions
   * order_constraint = 1 --> order matters
   *
   * note that in contrast to dist_4v_allt_v2 here we do not have the information from the cooccurrences of motif 1, motif 2 and motif 3 to start with
   * 
   * this function comes from dist_4v_allt
   * 01-03-2004: we change the function name to further distinguish this from dist_4v_allt_v3
   * 01-03-2004: here we include the strand information
   * 01-03-2004: call dist_2v_allt_v2
   * 01-03-2004: call dist_3v_allt_v3 and dist_4v_allt_v3
   * 01-31-2004: added pos_format
   * 02-02-2004: added verbose
   */

  int **ovdi12;                          // list of cooccurring instances of motifs 1 and 2
  int **ovdi123;                         // list of cooccurring instances of motifs 1, 2 and 3
  int max_n_ovdi12=500000;               // maximum number of entries in call to dist_2v_allt (memory = 500k x 5 cols x 2 bytes/int = 5 Mb)
  int max_n_ovdi123=100000;              // maximum number of entries in call to dist_3v_allt (memory = 100k x 7 cols x 2 bytes/int = 1.4 Mb)
  int nout12[4];
  int nout123[4];
  int nout1234[4];                      // [0]=n1234, [1]=n1234_tr, [2]=errorcode
  int n12;
  int n123;
  int get_positions12=1;
  int max_n_transcripts12=-1;
  int max_n_transcripts123=-1;
  int errorcode=0;
  int i,j;

  ovdi12=imatrix(1,max_n_ovdi12,1,7);
  ovdi123=imatrix(1,max_n_ovdi123,1,10);

  if (verbose) 
    printf("dist_4v_allt_from_m1234_v1\tmax_dist=%d\tmin_dist=%d\tmax_n_transcripts=%d\n",max_dist,min_dist,max_n_transcripts);

  if (verbose) {
    printf("motif 1 (%d): ",n1);
    for (i=1;i<=n1;i++) 
      printf("%d %d %d,",m1_t[i],m1_p[i],m1_s[i]);
    printf("\n");
    printf("motif 2 (%d): ",n2);
    for (i=1;i<=n2;i++)
      printf("%d % d %d,",m2_t[i],m2_p[i],m2_s[i]);
    printf("\n");
  }

  /* start by computing the co-occurrences of motif 1 and motif 2 */
  dist_2v_allt_v2(m1_t,m2_t,m1_p,m2_p,m1_s,m2_s,n1,n2,ovdi12,nout12,min_dist,max_dist,n_transcripts,ml,get_positions12,order_constraint,max_n_transcripts12,max_n_ovdi12,pos_format,0);
  errorcode+=nout12[2];      // nout12[2]<0 if there are errors
  if (nout12[2]<0) {
    n12=nout12[3];           // if there was an overflow, then we wish to use the number of entries that have a position reported
  } else {
    n12=nout12[0];           // total number of co-occurrences (total number of rows in ovdi12)
  }
  if (verbose) {
    printf("ovdi12 (%d)\n",n12);
    for (i=1;i<=n12;i++) 
      printf("%d %d %d %d %d\n",ovdi12[i][1],ovdi12[i][2],ovdi12[i][3],ovdi12[i][6],ovdi12[i][7]);
    printf("motif 3 (%d): ",n3);
    for (i=1;i<=n3;i++)
      printf("%d % d %d,",m3_t[i],m3_p[i],m3_s[i]);
    printf("\n");
  }

  /* now call dist_3v_allt_v3 */
  dist_3v_allt_v3(ovdi12,m3_t,m3_p,m3_s,n12,n3,ovdi123,nout123,min_dist,max_dist,n_transcripts,ml,get_positions,order_constraint,max_n_transcripts123,max_n_ovdi123,pos_format,verbose);
  errorcode+=nout123[2];     // errorcode<0 if there are errors
  if (nout123[2]<0) {
    n123=nout123[3];
    if (verbose) 
      printf("ERROR!\tcooc_single_lib.c\tdist_4v_allt_from_m1234_v1\tdist_3v_allt_v3\tnout123[2]=%d\tn123=%d\n",nout123[2],n123);
  } else {
    n123=nout123[0];
  }
  if (verbose) {
    printf("ovdi123 (%d)\n",n123);
    for (i=1;i<=n123;i++) 
      printf("%d %d %d %d %d %d %d\n",ovdi123[i][1],ovdi123[i][2],ovdi123[i][3],ovdi123[i][4],ovdi123[i][8],ovdi123[i][9],ovdi123[i][10]);
    printf("motif 4 (%d): ",n4);
    for (i=1;i<=n4;i++)
      printf("%d % d %d,",m4_t[i],m4_p[i],m4_s[i]);
    printf("\n");
  }

  /* now call dist_4v_allt_v2 */
  dist_4v_allt_v3(ovdi123,m4_t,m4_p,m4_s,n123,n4,ovdi,nout1234,min_dist,max_dist,n_transcripts,ml,get_positions,order_constraint,max_n_transcripts,max_n_ovdi,pos_format);
  n[0]=nout1234[0];
  n[1]=nout1234[1];
  errorcode+=nout1234[2];    // errorcode<0 if there are errors
  n[2]=errorcode;
  n[3]=nout1234[3];
  if (verbose) {
    printf("ovdi1234 (%d)\n",n[0]);
    for (i=1;i<=n[0];i++) 
      printf("%d %d %d %d %d %d %d %d %d\n",ovdi[i][1],ovdi[i][2],ovdi[i][3],ovdi[i][4],ovdi[i][5],ovdi[i][10],ovdi[i][11],ovdi[i][12],ovdi[i][13],ovdi[i][14]);
  }

  free_imatrix(ovdi12,1,max_n_ovdi12,1,7);
  free_imatrix(ovdi123,1,max_n_ovdi123,1,10);  
}

int remove_overlaps_v2(int **m_old,int **m_new,int *ml,int n_entries,int n_motifs)
{
  /* given a matrix with the positions of each motif for multiple occurrences, remove overlaps 
   * (situations with the same motifs in "basically" the same positions)
   * returns the remaining number of entries
   * copies the indices entries in ovdi matrices
   *
   * new in version 2: 01-01-2004: copy also the strand information
   */

  int issimilar;                     /* 1 if two entries are similar and 0 otherwise */
  int i,j,k;                           
  int running_n=1;                   /* running number of entries */
  int tempsimilar;                   /* temporary variable used to compare two entries */
 
  /* the first entry always stays */
  for (j=1;j<=(3*n_motifs);j++) 
    m_new[1][j]=m_old[1][j];  

  for (i=2;i<=n_entries;i++) {
    issimilar=0;                     /* by default, the two entries are not similar */
    k=1;                             /* compare entry number i with all the previous ones starting with k=1 and up to k=running_n or until we find a similarity */
    while ( (issimilar == 0) & (k<=running_n) ) {
      tempsimilar=0;	  
      for (j=1;j<=n_motifs;j++) {	    
	if ( (fabs(m_new[k][j]-m_old[i][j])) <= (ml[j]/2) )
	  tempsimilar++;
      }
      if (tempsimilar == n_motifs) 
	issimilar=1;
      k++;	  
    }	
    if (issimilar==0) {
      //new_n++;
      running_n++;
      for (j=1;j<=(3*n_motifs);j++)                // in remove_overlaps, this was 2*n_motifs
	m_new[running_n][j]=m_old[i][j];
    }
  }

  return running_n;
}

int range_constraint(int **m_old,int **m_new,int range_threshold,int n_entries,int n_motifs)
{
  /* given a matrix with the positions of each motif for multiple occurrences, 
   * remove those entries where the range (maximum position - minimum position) is larger than range_threshold
   * returns the remaining number of entries
   * copies the indices entries in ovdi matrices
   */

  int i,j;                           
  int running_n=0;                   /* running number of entries in m_new */
  int max;                           /* maximum position for current entry */
  int min;                           /* minimum position for current entry */
  int range;                         /* range in current entry             */

  for (i=1;i<=n_entries;i++) {
    max=0;
    min=100000000;
    for (j=1;j<=n_motifs;j++) {
      if (m_old[i][j]<min) 
	min=m_old[i][j];
      if (m_old[i][j]>max) 
	max=m_old[i][j];
    }
    range=max-min;
    if (range<range_threshold) {
      running_n++;
      for (j=1;j<=n_motifs;j++) 
	m_new[running_n][j]=m_old[i][j];
    }
  }

  return running_n;
}

int read_scan_file(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos)
{
  /* here we read a whole scan file in the following format
   * header lines (start with '>' or '%')
   * transcript position position position ...
   * transcript position position position ...
   *
   * returns the total number of entries in transcripts and positions
   * returns also the transcripts in *transcripts and the positions in *positions
   */

  FILE *input_file;
  int t[3];                 // t[0]= current transcript  t[1]=number of positions in current transcript
  int *m_p;                 // positions in current transcript
  int max_occur=10000;      // maximum number of positions per transcript
  int curr_tr;              // current transcript
  int curr_n;               // number of positions in current transcript
  int running_n=0;          // running number of entries in transcripts and positions
  int i;
  int togoon=1;
  float *d;

  input_file=fopen(filename,"r");
  if (!input_file) {
    printf("error! i could not open %s\n",filename);
    return running_n;
  }
  m_p=ivector(1,max_occur);
  read_header(input_file,t,m_p,max_occur);
  curr_tr=t[0];
  curr_n=t[1];
  if (curr_n>0) {
    if (neg2pos==1) {
      for (i=1;i<=curr_n;i++) {
	if (m_p[i]<0) 
	  m_p[i]=-m_p[i];
      }
    }
    if (sortpos>0) {
      /* sort positions of first entry */
      d=vector(1,curr_n);
      quicksorti(curr_n,m_p,d);
      free_vector(d,1,curr_n);
    }
  }

  for (i=1;i<=curr_n;i++) {
    running_n++;
    transcripts[running_n]=curr_tr;
    positions[running_n]=m_p[i];
  }

  while (togoon>=1) {
    read_line(input_file,t,m_p,max_occur,neg2pos,sortpos);
    curr_tr=t[0];
    curr_n=t[1];
    togoon=t[2];
    for (i=1;i<=curr_n;i++) {
      running_n++;
      transcripts[running_n]=curr_tr;
      positions[running_n]=m_p[i];
    }
  }
  fclose(input_file);

  free_ivector(m_p,1,max_occur);

  return running_n;
}

int read_scan_file_v2(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos,int max_overlap)
{
  /* here we read a whole scan file in the following format
   * header lines (start with '>' or '%')
   * transcript position position position ...
   * transcript position position position ...
   *
   * returns the total number of entries in transcripts and positions
   * returns also the transcripts in *transcripts and the positions in *positions
   *
   * 10-17-2003: make it simpler by including the header stuff into the main loop
   * 10-17-2003: ignore overlapping occurrences
   */

  FILE *input_file;
  int t[3];                 // t[0]= current transcript  t[1]=number of positions in current transcript
  int *m_p;                 // positions in current transcript
  int max_occur=10000;      // maximum number of positions per transcript
  int curr_tr;              // current transcript
  int curr_n;               // number of positions in current transcript
  int running_n=0;          // running number of entries in transcripts and positions
  int i;
  int togoon=1;
  float *d;
  int *positions_no_overlap;
  int n_no_overlap;

  input_file=fopen(filename,"r");
  if (!input_file) {
    printf("error! i could not open %s\n",filename);
    return running_n;
  }
  m_p=ivector(1,max_occur);
  positions_no_overlap=ivector(1,max_occur);

  while (togoon>=1) {
    read_line(input_file,t,m_p,max_occur,neg2pos,sortpos);
    curr_tr=t[0];
    curr_n=t[1];
    togoon=t[2];
    if (curr_n>0) {
      if (max_overlap>0) {           // remove overlapping occurrences
	n_no_overlap=vector_overlap(m_p,curr_n,max_overlap,positions_no_overlap);
	for (i=1;i<=n_no_overlap;i++) {
	  running_n++;
	  transcripts[running_n]=curr_tr;
	  positions[running_n]=positions_no_overlap[i];
	}
      } else {
	for (i=1;i<=curr_n;i++) {
	  running_n++;
	  transcripts[running_n]=curr_tr;
	  positions[running_n]=m_p[i];
	}
      }       // close check on max_overlap
    }         // close check on curr_n>0
  }           // close togoon
  fclose(input_file);

  free_ivector(m_p,1,max_occur);
  free_ivector(positions_no_overlap,1,max_occur);

  return running_n;
}

int get_max_n_transcripts(int n_transcripts_cl,int n_cl,double cdf_threshold,int n_transcripts_bck)
{
  // get the maximum number of transcripts that one needs to run to achieve a minimum p value
  // usage: max_n_transcripts=get_max_n_transcripts(k,n,y,n_all);
  // 11-03-2003: make cdf_threshold a double
  // 11-03-2003: make p a double

  double p;                   // probability output for invbinocdf2
  double r=0.001;             // step size
  int max_n_transcripts;      // maximum number of transcripts
  int function_verbose=0;     // if this is 1, then report performace time
  struct timeb ti,tf;         // used to compute performance time
  long tdif;                  // used to compute performance time

  /* double invbinocdf2(double y,int n,int k,double r)
   * returns the probability p in the background so that a binomial with parameters n,k,p will yield a cumulative distribution of y
   * r is the resolution in the search space [default=0.01]
   *
   * e.g. we observed k transcripts in the cluster out of a total of n transcripts in the cluster
   * we wish to know the maximum number of transcripts to process in the background so as to achieve a maximum cdf of y
   */ 

  if (function_verbose==1) {
    printf("get_max_n_transcripts\tn_transcripts_cl=%d\tn_cl=%d\tcdf_threshold=%1.10g\tn_transcripts_bck=%d\n",n_transcripts_cl,n_cl,cdf_threshold,n_transcripts_bck);
    ftime(&ti);
  }
  p=invbinocdf2(cdf_threshold,n_transcripts_cl,n_cl,r); 

  if (p<=0) 
    p=r;

  max_n_transcripts=ceil(p*n_transcripts_bck)+1;

  if (function_verbose==1) {
    ftime(&tf);
    tdif=time_difference(ti,tf);
    printf("p=%1.6g\n",p);
    printf("max_n_transcripts=%d\n",max_n_transcripts);
    printf("\tget_max_n_transcripts:\t%ld ms cpu time\n",tdif);
  }

  return max_n_transcripts;
}

void read_line_v2(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int max_overlap,int seq_ul,int seq_el,int seq_il)
{
  /* given a file handle, read a single line
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   *   with the positions measured with respect to the 3' end of the sequence
   * variable t
   *    t[0] = transcript number
   *    t[1] = number of positions returned in data_p
   *    t[2] = -1 -> end of file   t[2]=0 -> n>max_occur  t[2]=1 -> ok  t[2]=2 -> n=0 (no positions)
   *    t[3] = number of positions before filtering and before removing overlaps
   *    t[4] = number of positions after filtering before removing overlaps
   *
   * new in version 2
   * 10-27-2003: allow to filter the positions by considering only shorter sequences
   *
   * 05-31-2003: read files created with rc=1 (negative positions indicating the opposite strand
   *
   * if neg2pos==1, convert negative positions to positive ones 
   * if sortpos==1, sort the positions
   * if the end of file is reached, returns all t entries as -1, same if it is a comment line with t[2]=2
   */

  char txtline[100000];    // read one line of the input file
  int i,j,k,ll;           // loop counters
  int *v;                 // vector holding the transcripts and positions for running entry
  int start=0;            // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];              // number of entries
  int n;                  // number of entries
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int header;                           // 1 if current line is a header line
  int *temp_positions;                  // temporary positions (contains the raw positions before filtering
  int *filtered_positions;              // filtered positions (contains the positions after filtering)
  int *positions_no_overlap;            // contains the positions after filtering and removing overlaps
  float *indices;                       // dummy variable used in call to sort function
  int curr_transcript;                  // current transcript
  int curr_ul,curr_el,curr_il;          // upstream length, exon length, intron length
  int n_after_filter;                   // number of occurrences after filtering
  int n_no_overlap=0;                   // number of occurrences that do not overlap

  if (feof(input_file)!=0) {
    t[0]=-1;
    t[1]=-1;
    t[2]=-1;
    t[3]=-1;
    t[4]=-1;
  } else {
    fgets(txtline,n_chars,input_file);
    header=is_header(txtline);
    if (header==1) {
      t[0]=-1;
      t[1]=-1;
      t[2]=2;
      t[3]=-1;
      t[4]=-1;
    } else {
      v=ivector(1,max_occur);
      splitline2int(txtline,v,start,np,max_occur);
      n=np[0]-1;
      if (n>0) {
	if (n>max_occur) {
	  n=max_occur;
	  t[2]=0;
	  printf("WARNING! n=%d and max_occur=%d\n",n,max_occur);
	}
 	temp_positions=ivector(1,max_occur);
	filtered_positions=ivector(1,max_occur);
	indices=vector(1,max_occur);
	positions_no_overlap=ivector(1,max_occur);
	curr_transcript=v[1];
	t[0]=curr_transcript;  // transcript number
	curr_ul=transcript_lengths[curr_transcript][1];
	if (curr_ul>seq_ul)
	  curr_ul=seq_ul;
	curr_el=transcript_lengths[curr_transcript][2];
	if (curr_el>seq_el)
	  curr_el=seq_el;
	curr_il=transcript_lengths[curr_transcript][3];
	if (curr_il>seq_il)
	  curr_il=seq_il;

	t[3]=n-1;            // number of positions before filtering and removing overlaps
	t[2]=1;              // ok
	if (neg2pos==1) {
	  for (i=2;i<=n;i++) {
	    j=i-1;
	    if (v[i]>0) {
	      temp_positions[j]=v[i];
	    } else {
	      temp_positions[j]=-v[i];
	    }
	  }
	} else {
	  //for (i=1;i<=n;i++) {                    
	  for (i=2;i<=n;i++) {                 // 01-02-2004 changed the loop to start at i=2 (i=1 is reserved for the transcript number
	    j=i-1;
	    temp_positions[j]=v[i];
	  }
	}
	n_after_filter=filter_positions(temp_positions,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	t[4]=n_after_filter;
	if (n_after_filter>0) {
	  if (sortpos) {
	    quicksorti(n_after_filter,filtered_positions,indices);
	  }
	  n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_overlap,positions_no_overlap);
	  for (k=1;k<=n_no_overlap;k++) {
	    data_p[k]=positions_no_overlap[k];
	  }
	}
	t[1]=n_no_overlap;

	/* free memory */
	free_ivector(temp_positions,1,max_occur);
	free_ivector(filtered_positions,1,max_occur);
	free_ivector(positions_no_overlap,1,max_occur);
	free_vector(indices,1,max_occur);
      } else {
	t[2]=2;
      }   // close check on n>0
      
      free_ivector(v,1,max_occur);
     }  // close if (is_header(txtline)==0)     
  }    // close feof check
}      // end of read_line_v2
   
void read_line_wstrand_v2(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int max_overlap,int seq_ul,int seq_el,int seq_il)
{
  /* given a file positions and a file strands handles, read a single line from each file
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ... for the positions file
   *   transcript\tstrand\tstrand\tstrand ... for the strands file
   *   with the positions measured with respect to the 3' end of the sequence
   * variable t
   *    t[0] = transcript number
   *    t[1] = number of positions returned in data_p
   *    t[2] = -1 -> end of file   t[2]=0 -> n>max_occur  t[2]=1 -> ok  t[2]=2 -> n=0 (no positions)
   *    t[3] = number of positions before filtering and before removing overlaps
   *    t[4] = number of positions after filtering before removing overlaps
   *
   * 01-30-2004 call vector_overlap_v3
   */

  char txtline_pos[100000];             // read one line of the positions input file
  char txtline_strand[10000];           // read one line of the strands input file
  int i,j,k,ll;                         // loop counters
  int *v;                               // vector holding the transcripts and positions for running entry
  int *w;                               // vector holding the transcripts and strands for running entry 
  int start=0;                          // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];                            // number of entries
  int n;                                // number of entries
  int n_chars=1000000;                  // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  float *indices;                       // indices variable used in sorting
  int header;                           // 1 if current line is a header line
  int *temp_positions;                  // temporary positions (contains the raw positions before filtering)
  char *temp_strands;                   // temporary strands (contains the raw strands before filtering)
  int *filtered_positions;              // filtered positions (contains the positions after filtering)
  char *filtered_strands;               // filtered strands (contains the strands after filtering)
  int *positions_no_overlap;            // contains the positions after filtering and removing overlaps
  char *strands_no_overlap;             // contains the strands after filtering and removing overlaps
  int curr_transcript;                  // current transcript
  int curr_ul,curr_el,curr_il;          // upstream length, exon length, intron length
  int n_after_filter;                   // number of occurrences after filtering
  int n_no_overlap=0;                   // number of occurrences that do not overlap

  t[0]=-1;t[1]=-1;t[2]=-1;t[3]=-1;t[4]=-1;
  if (feof(input_file_pos)==0) {
    fgets(txtline_pos,n_chars,input_file_pos);
    fgets(txtline_strand,n_chars,input_file_strand);
    header=is_header(txtline_pos);
    if (header==1) {
      t[2]=2;
    } else {
      v=ivector(1,max_occur);
      w=ivector(1,max_occur);
      splitline2int(txtline_pos,v,start,np,max_occur);
      splitline2int(txtline_strand,w,start,np,max_occur);
      n=np[0]-1;
      if (n>0) {
	if (n>max_occur) {
	  n=max_occur;
	  t[2]=0;
	  printf("WARNING! n=%d and max_occur=%d\n",n,max_occur);
	}
 	temp_positions=ivector(1,max_occur);
	temp_strands=cvector(1,max_occur);
	filtered_positions=ivector(1,max_occur);
	filtered_strands=cvector(1,max_occur);
	indices=vector(1,max_occur);
	positions_no_overlap=ivector(1,max_occur);
	strands_no_overlap=cvector(1,max_occur);
	curr_transcript=v[1];
	t[0]=curr_transcript;  // transcript number
	curr_ul=transcript_lengths[curr_transcript][1];
	if (curr_ul>seq_ul)
	  curr_ul=seq_ul;
	curr_el=transcript_lengths[curr_transcript][2];
	if (curr_el>seq_el)
	  curr_el=seq_el;
	curr_il=transcript_lengths[curr_transcript][3];
	if (curr_il>seq_il)
	  curr_il=seq_il;
	t[3]=n-1;            // number of positions before filtering and removing overlaps
	t[2]=1;              // ok
	if (neg2pos==1) {
	  for (i=2;i<=n;i++) {
	    j=i-1;
	    temp_strands[j]=w[i];
	    if (v[i]>0) {
	      temp_positions[j]=v[i];
	    } else {
	      temp_positions[j]=-v[i];
	    }
	  }
	} else {
	  //for (i=1;i<=n;i++) {          01-02-2004 shouldn't this start at i=2?
	  for (i=2;i<=n;i++) {
	    j=i-1;
	    temp_strands[j]=w[i];
	    temp_positions[j]=v[i];
	  }
	}
	
	//n_after_filter=filter_positions(temp_positions,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	n_after_filter=filter_positions_v3(temp_positions,temp_strands,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
	t[4]=n_after_filter;
	if (n_after_filter>0) {
	  if (sortpos) {
	    for (i=1;i<=n_after_filter;i++)
	      indices[i]=i;
	    quicksorti(n_after_filter,filtered_positions,indices);
	    sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
	  }
	  //n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_overlap,positions_no_overlap);
	  //n_no_overlap=vector_overlap_v2(filtered_positions,n_after_filter,max_overlap,positions_no_overlap,filtered_strands,strands_no_overlap);
	  n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_overlap);
	  for (k=1;k<=n_no_overlap;k++) {
	    data_p[k]=positions_no_overlap[k];
	    data_s[k]=strands_no_overlap[k];
	  }
	}
	t[1]=n_no_overlap;

	/* free memory */
	free_ivector(temp_positions,1,max_occur);
	free_cvector(temp_strands,1,max_occur);
	free_ivector(filtered_positions,1,max_occur);
	free_cvector(filtered_strands,1,max_occur);
	free_ivector(positions_no_overlap,1,max_occur);
	free_cvector(strands_no_overlap,1,max_occur);
	free_vector(indices,1,max_occur);
      } else {
	t[2]=2;
      }   // close check on n>0
      
      free_ivector(v,1,max_occur);
      free_ivector(w,1,max_occur);
     }  // close if (is_header(txtline)==0)     
  }    // close feof check
}      // end of read_line_v2

void read_line_v3(FILE *input_file,int *t,int *data_p,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il)
{
  /* given a file handle, read a single line
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   *   with the positions measured with respect to the 3' end of the sequence
   * variable t
   *    t[0] = transcript number
   *    t[1] = number of positions returned in data_p
   *    t[2] = -1 -> end of file   t[2]=0 -> n>max_occur  t[2]=1 -> ok  t[2]=2 -> n=0 (no positions)  t[2]=3 -> curr_transcript=0,odd
   *    t[3] = number of positions before filtering and before removing overlaps
   *    t[4] = number of positions after filtering before removing overlaps
   *
   * new in version 3
   * 10-30-2003: add sequence ul, el and il so as to report positions with respect to sequence 3' end
   * 10-28-2003: this version is for the scan of all sequences
   * 10-28-2003: use a map of the transcript numbers (because there are many transcripts that are not locus link transcripts)
   *             transcript_map[i][1] = original index
   *             transcript_map[i][2] = new index (locus link file index)
   * 10-28-2003: use different formatting for the scan output (with < 0 values for upstream sequences, etc.)
   *
   * new in version 2
   * 10-27-2003: allow to filter the positions by considering only shorter sequences
   *
   * 05-31-2003: read files created with rc=1 (negative positions indicating the opposite strand
   *
   * if neg2pos==1, convert negative positions to positive ones 
   * if sortpos==1, sort the positions
   * if the end of file is reached, returns all t entries as -1, same if it is a comment line with t[2]=2
   */

  char txtline[100000];    // read one line of the input file
  int i,j,k,ll;           // loop counters
  int *v;                 // vector holding the transcripts and positions for running entry
  int start=0;            // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];              // number of entries
  int n;                  // number of entries
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  //float *d;               // dummy variable used in sorting
  int header;                           // 1 if current line is a header line
  int *temp_positions;                  // temporary positions (contains the raw positions before filtering
  int *filtered_positions;              // filtered positions (contains the positions after filtering)
  int *positions_no_overlap;            // contains the positions after filtering and removing overlaps
  float *dummy;                         // dummy variable used in call to sort function
  int curr_transcript;                  // current transcript
  int curr_ul,curr_el,curr_il;          // upstream length, exon length, intron length
  int n_after_filter;                   // number of occurrences after filtering
  int n_no_overlap=0;                   // number of occurrences that do not overlap
  int curr_transcript_index;

  t[0]=-1;t[1]=-1;t[2]=-1;t[3]=-1;t[4]=-1;

  if (feof(input_file)!=0) {
    /* done with file, exit */
  } else {
    fgets(txtline,n_chars,input_file);
    header=is_header(txtline);
    if (header==1) {
      t[2]=2;
    } else {
      /* debug here       printf("\ntxtline=%s\tmax_occur=%d\n",txtline,max_occur);   */
      v=ivector(1,max_occur);
      splitline2int(txtline,v,start,np,max_occur);
      n=np[0]-1;
      /* debug here      printf("n=%d\n",n);  end debug here */
      if (n>0) {
	if (n>max_occur) {
	  n=max_occur;
	  t[2]=0;
	  printf("WARNING! n=%d and max_occur=%d\n",n,max_occur);
	}
	curr_transcript_index=v[1];
	curr_transcript=transcript_map[curr_transcript_index][2];
	t[0]=curr_transcript;  // transcript number
	if (curr_transcript<=0) {
	  t[2]=3;
	} else {
	  temp_positions=ivector(1,max_occur);
	  filtered_positions=ivector(1,max_occur);
	  dummy=vector(1,max_occur);
	  positions_no_overlap=ivector(1,max_occur);
	  /* printf("curr_transcript=%d\tcurr_transcript_index=%d\n",curr_transcript,curr_transcript_index); debug here */
	  curr_ul=transcript_lengths[curr_transcript][1];
	  if (curr_ul>seq_ul)
	    curr_ul=seq_ul;
	  curr_el=transcript_lengths[curr_transcript][2];
	  if (curr_el>seq_el)
	    curr_el=seq_el;
	  curr_il=transcript_lengths[curr_transcript][3];
	  if (curr_il>seq_il)
	    curr_il=seq_il;
	  /* debug here 	printf("\ncurrent_transcript_index=%d\tcurrent_transcript=%d\tcurr_ul=%d\tcurr_el=%d\tcurr_il=%d\tmax_overlap=%d\n",curr_transcript_index,curr_transcript,curr_ul,curr_el,curr_il,max_overlap); */

	  t[3]=n-1;              // number of positions before filtering and removing overlaps
	  t[2]=1;                // ok
	  if (neg2pos==1) {
	    for (i=2;i<=n;i++) {
	      j=i-1;
	      if (v[i]>0) {
		temp_positions[j]=v[i];
	      } else {
		temp_positions[j]=-v[i];
	      }
	    }
	  } else {
	    // for (i=1;i<=n;i++) {    01-02-2004 changed the loop to start at i=2; i=1 is reserved for the transcript number
	    for (i=2;i<=n;i++) {
	      j=i-1;
	      temp_positions[j]=v[i];
	    }
	  }      // close check on neg2pos==1
	  /* debug here 	printf("n=%d\t",n);	for (i=1;i<=(n-1);i++)	  printf("%d ",temp_positions[i]); end debug here */
	  n_after_filter=filter_positions_v2(temp_positions,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	  /* debug here 	printf("\nn_after_filter=%d\t",n_after_filter);	for (i=1;i<=n_after_filter;i++) 	  printf("%d ",filtered_positions[i]);	printf("\n"); end debug here */
	  t[4]=n_after_filter;
	  if (n_after_filter>0) {
	    if (sortpos)
	      quicksorti(n_after_filter,filtered_positions,dummy);
	    n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_overlap,positions_no_overlap);
	    for (k=1;k<=n_no_overlap;k++)
	      data_p[k]=positions_no_overlap[k];
	  }     // close check on n_after_filter>0
	  t[1]=n_no_overlap;
	  /* debug here 	printf("n_no_overlap=%d\t",n_no_overlap);	for (i=1;i<=n_no_overlap;i++) 	  printf("%d ",positions_no_overlap[i]);	printf("\n"); end debug here */

	  /* free memory */
	  free_ivector(temp_positions,1,max_occur);
	  free_ivector(filtered_positions,1,max_occur);
	  free_ivector(positions_no_overlap,1,max_occur);
	  free_vector(dummy,1,max_occur);
	}     // close check on curr_transcript>0
      } else {
	t[2]=2;
      }   // close check on n>0
      free_ivector(v,1,max_occur);
     }  // close if (is_header(txtline)==0)     
  }    // close feof check
}      // end of read_line_v2

void read_line_wstrand_v3(FILE *input_file_pos,FILE *input_file_strand,int *t,int *data_p,char *data_s,int max_occur,int neg2pos,int sortpos,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int max_overlap,int seq_ul,int seq_el,int seq_il)
{
  /* given a file handle, read a single line
   * the entries have the form:
   *   transcript\tposition\tposition\tposition ...
   *   with the positions measured with respect to the 3' end of the sequence
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
   */

  char txtline_pos[100000];              // read one line of the input positions file
  char txtline_strand[100000];           // read one line of the input strands file 
  int i,j,k,ll;                          // loop counters
  int *w;                                // vector holding the transcripts and strands for running entry
  int *v;                                // vector holding the transcripts and positions for running entry
  int start=0;                           // the 3rd entry is the motif length, the 4th one is the number of sequences, then there is the transcript and position
  int np[1];                             // number of entries
  int n;                                 // number of entries
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
      splitline2int(txtline_strand,w,start,np,max_occur);
      if (n>0) {
	if (n>max_occur) {
	  n=max_occur;
	  t[2]=0;
	  printf("WARNING! n=%d and max_occur=%d\n",n,max_occur);
	}

	curr_transcript_index=v[1];
	curr_transcript=transcript_map[curr_transcript_index][2];
	t[0]=curr_transcript;  // transcript number
	if (curr_transcript<=0) {
	  t[2]=3;
	  printf("WARNING!!! curr_transcript=%d\nthis is odd\n");
	} else {
	  /* get curr_ul,curr_el,curr_el */
	  curr_ul=transcript_lengths[curr_transcript][1];
	  if (curr_ul>seq_ul)
	    curr_ul=seq_ul;
	  curr_el=transcript_lengths[curr_transcript][2];
	  if (curr_el>seq_el)
	    curr_el=seq_el;
	  curr_il=transcript_lengths[curr_transcript][3];
	  if (curr_il>seq_il)
	    curr_il=seq_il;

	  temp_positions=ivector(1,max_occur);               //temp_positions=ivector(1,max_occur);
	  temp_strands=cvector(1,max_occur);                 //temp_strands=cvector(1,max_occur);
	  filtered_positions=ivector(1,max_occur);           //filtered_positions=ivector(1,max_occur);
	  filtered_strands=cvector(1,max_occur);             //filtered_strands=cvector(1,max_occur);
	  indices=vector(1,max_occur);               //indices=vector(1,max_occur);
	  positions_no_overlap=ivector(1,max_occur);         //positions_no_overlap=ivector(1,max_occur);
	  strands_no_overlap=cvector(1,max_occur);           //strands_no_overlap=cvector(1,max_occur);

	  t[3]=n-1;              // number of positions before filtering and removing overlaps
	  t[2]=1;                // ok
	  if (neg2pos==1) {
	    for (i=2;i<=n;i++) {
	      j=i-1;
	      temp_strands[j]=w[i];
	      if (v[i]>0) {
		temp_positions[j]=v[i];
	      } else {
		temp_positions[j]=-v[i];
	      }
	    }
	  } else {
	    // for (i=1;i<=n;i++) {    01-02-2004 changed the loop to start at i=2; i=1 is reserved for the transcript number
	    for (i=2;i<=n;i++) {
	      j=i-1;
	      temp_positions[j]=v[i];
	      temp_strands[j]=w[i];
	    }
	  }      // close check on neg2pos==1
	  //n_after_filter=filter_positions_v2(temp_positions,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions); // filter by position
	  n_after_filter=filter_positions_v4(temp_positions,temp_strands,n-1,max_ul,max_el,max_il,curr_ul,curr_el,curr_il,filtered_positions,filtered_strands);
	  t[4]=n_after_filter;
	  if (n_after_filter>0) {
	    if (sortpos) {
	      for (i=1;i<=n_after_filter;i++)
		indices[i]=i;
	      quicksorti(n_after_filter,filtered_positions,indices);
	      sort_cvector_by_indices(n_after_filter,filtered_strands,indices);
	    }
	    //n_no_overlap=vector_overlap(filtered_positions,n_after_filter,max_overlap,positions_no_overlap);
	    //n_no_overlap=vector_overlap_v2(filtered_positions,n_after_filter,max_overlap,positions_no_overlap,filtered_strands,strands_no_overlap);
	    n_no_overlap=vector_overlap_v3(filtered_positions,filtered_strands,positions_no_overlap,strands_no_overlap,n_after_filter,max_overlap);
	    for (k=1;k<=n_no_overlap;k++) {
	      data_p[k]=positions_no_overlap[k];
	      data_s[k]=strands_no_overlap[k];
	    }
	  }     // close check on n_after_filter>0
	  t[1]=n_no_overlap;
	  /* free memory */
	  free_ivector(temp_positions,1,max_occur);
	  free_cvector(temp_strands,1,max_occur);
	  free_ivector(filtered_positions,1,max_occur);
	  free_cvector(filtered_strands,1,max_occur);
	  free_ivector(positions_no_overlap,1,max_occur);
	  free_cvector(strands_no_overlap,1,max_occur);
	  free_vector(indices,1,max_occur);
	}     // close check on curr_transcript>0
      } else {
	t[2]=2;
      }   // close check on n>0
      free_ivector(v,1,max_occur);
      free_ivector(w,1,max_occur);
     }  // close if (is_header(txtline)==0)     
  }    // close feof check
}      // end of read_line_v2

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
	if (data_format==1) {
	  /* positions measured with respect to the TSS but with the upstream as > 0 */
	  curr_transcript_index=v[1];
	  curr_transcript=transcript_map[curr_transcript_index][2];
	} 
	if (data_format==2) {
	  /* positions measured with respect to the 3' end of the sequence */
	  curr_transcript=v[1];
	}

	t[0]=curr_transcript;  // transcript number
	if (verbose)
	  printf("transcript=%d\n",curr_transcript);
	if (curr_transcript<=0) {
	  t[2]=3;
	  printf("WARNING!!! curr_transcript=%d\nthis is odd\n");
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

int read_scan_file_v3(char *filename,int *transcripts,int *positions,int neg2pos,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il)
{
  /* here we read a whole scan file in the following format
   * header lines (start with '>' or '%')
   * transcript position position position ...
   * transcript position position position ...
   *
   * returns the total number of entries in transcripts and positions
   * returns also the transcripts in *transcripts and the positions in *positions
   *
   * new in v3
   * 10-30-2003: add seq_ul,seq_el,seq_il to allow a default sequence length to measure positions with respect to the 3' end
   * 10-28-2003: add transcript_map, if transcript_map[1][1]>0 then call read_line_v3, otherwise just call read_line_v2 and ignore transcript_map
   * 10-27-2003: call read_line_v2 (therefore allowing filtering of the scan output)
   * 10-27-2003: call to vector_overlap is now performed within read_line_v2
   * new in v2
   * 10-17-2003: make it simpler by including the header stuff into the main loop
   * 10-17-2003: ignore overlapping occurrences
   */

  FILE *input_file;
  int t[5];                 // t[0]= current transcript  t[1]=number of positions in current transcript
  int *m_p;                 // positions in current transcript
  int max_occur=10000;      // maximum number of positions per transcript
  int curr_tr;              // current transcript
  int curr_n;               // number of positions in current transcript
  int running_n=0;          // running number of entries in transcripts and positions
  int i;
  int togoon=1;
  float *d;
  int call2or3=2;           // call read_line_v2 unless there is an entry in the map

  i=1;
  while ( (i<10000) & (call2or3==2) ) {
    if (transcript_map[i][1]>0)
      call2or3=3;
    if (transcript_map[i][2]>0)
      call2or3=3;
    i++;
  }

  input_file=fopen(filename,"r");
  if (!input_file) {
    printf("error! i could not open %s\n",filename);
    return running_n;
  }
  m_p=ivector(1,max_occur);
  //positions_no_overlap=ivector(1,max_occur);

  while (togoon>=0) {
    if (call2or3==3) {
      read_line_v3(input_file,t,m_p,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,transcript_map,max_overlap,seq_ul,seq_el,seq_il);
    } else {
      read_line_v2(input_file,t,m_p,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,max_overlap,seq_ul,seq_el,seq_il);
    }
    curr_tr=t[0];
    curr_n=t[1];
    togoon=t[2];
    /* debug here  printf("read_scan_file_v3: curr_tr=%d\tcurr_n=%d\ttogoon=%d\n",curr_tr,curr_n,togoon); */
    if (curr_n>0) {
      for (i=1;i<=curr_n;i++) {
	running_n++;
	transcripts[running_n]=curr_tr;
	positions[running_n]=m_p[i];
      }       // close i loop
    }         // close check on curr_n>0
  }           // close togoon
  fclose(input_file);

  free_ivector(m_p,1,max_occur);
  //free_ivector(positions_no_overlap,1,max_occur);

  return running_n;
}

int read_scan_file_v4(char *filename_positions,char *filename_strands,int *transcripts,int *positions,char *strands,int neg2pos,int sortpos,int max_overlap,int max_ul,int max_el,int max_il,int **transcript_lengths,int **transcript_map,int seq_ul,int seq_el,int seq_il,int *out)
{
  /* here we read a whole scan file in the following format
   * header lines (start with '>' or '%')
   * transcript position position position ...
   * transcript position position position ...
   *
   * same format with strands instead of positions for the strands file
   *
   * returns the total number of entries in transcripts and positions
   * returns also the transcripts in *transcripts and the positions in *positions
   *
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
  int call2or3=2;           // call read_line_wstrand_v2 unless there is an entry in the map; in that case, call read_line_wstrand_v3
  int n_orig;               // original number of entries
  int n_after_filter;       // number of entries after filtering
  int cumul_filter_out=0;   // cumulative number of entries filtered out
  int cumul_overlap_out=0;  // cumulative number of entries removed (including filter and overlap)

  i=1;
  while ( (i<10000) & (call2or3==2) ) {
    if (transcript_map[i][1]>0)
      call2or3=3;
    if (transcript_map[i][2]>0)
      call2or3=3;
    i++;
  }

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
    if (call2or3==3) {
      //read_line_v3(input_file,t,m_p,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,transcript_map,max_overlap,seq_ul,seq_el,seq_il);
      read_line_wstrand_v3(input_file_pos,input_file_strand,t,m_p,m_s,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,transcript_map,max_overlap,seq_ul,seq_el,seq_il);
    } else {
      //read_line_v2(input_file,t,m_p,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,max_overlap,seq_ul,seq_el,seq_il);
      read_line_wstrand_v2(input_file_pos,input_file_strand,t,m_p,m_s,max_occur,neg2pos,sortpos,max_ul,max_el,max_il,transcript_lengths,max_overlap,seq_ul,seq_el,seq_il);
    }
    curr_tr=t[0];
    curr_n=t[1];
    togoon=t[2];
    n_orig=t[3];
    n_after_filter=t[4];
    cumul_filter_out+=(n_orig-n_after_filter);
    cumul_overlap_out+=(n_orig-curr_n);
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
   *
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
}
