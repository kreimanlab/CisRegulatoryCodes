/* scan_motif_v11.c
 * given a motif weight matrix and a list of upstream sequences, scan for the presence of the motif
 * positions are reported with respect to the 3' end of the sequence
 *
 * gabriel kreiman
 * kreiman@mit.edu
 * last modified: 02-16-2004
 */

/* new in version 11:
 * 02-16-2004: call scan_motif4
 */

/* new in version 10:
 * 02-02-2004: solved a bug in the strand report (needed to assign indices to 1,...,np in each transcript!
 */

/* new in version 9:
 * write the strand information in a separate file
 * 12-28-2003: solved bug, if sortpos=1, then the strands should be printed according to the sorted indices!
 */

/* new in version 8:
 * 08-27-2003: added sortpos
 * 08-17-2003: call functions in regul_search_methods_v7.c
 * 08-17-2003: position and strand information in separate vectors (do not assign negative values to positions in the minus strand
 * 08-17-2003: here we do not use n_seqs4motif
 */

/* new in version 7:
 * 07-17-2003: do not compute logarithm, just add up the corresponding entries of the weight matrix; the log is now done when computing the weight matrix
 * 07-17-2003: use the input weight matrix as is (the fudge correction is done elsewhere
 */

/* new in version 6:
 * allow sequences of arbitrary length
 */

/* new in version 5:
 * allow also scanning of reverse complementary sequence
 * 05-10-2003: retrieve scores for each occurrence
 */

/* new in version 4:
 * scanning of fasta formatted files
 */

/* new in version 3:
 * allow to scan only a segment of the upstream sequences
 */

/* 02-10-2003: change to score>=threshold to accept an occurrence instead of score>threshold. this brought problems, for example, when only perfect matches 
 * were required for a weight matrix with only 1s and 0s.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "my_config.c"
#include "../lib/nr.c"
#include "../lib/my_math_package.c"
#include "../lib/my_file_io.c"
#include "../lib/my_text_package.c"
#include "regul_search_methods_v7.c"

void usage(void);

int verbose=0;
int prog_version=11;
float current_score;
FILE *log_file;
int motif_length;                        // number of nucleotides in the motif
int seq_nr;                              // running number of sequences scanned
char last_modified[200];
int debug_mode=0;

main(int argc, char *argv[])
{
  FILE *output_file;                       // file with position results
  FILE *strand_file;                       // file with strand information
  FILE *upstream_sequence_file;            // file with sequences to scan
  FILE *weight_matrix_file;                // file containing the weigh matrix
  char *num_upstream_sequence;             // upstream sequence in numeric format
  char *strand;
  char info1[1000];
  char info2[1000];
  char infotext[5000];
  char log_filename[200];                  // log filename
  char output_filename[200];               // file name with output position results
  char strand_filename[200];               // file name for the strand of each hit
  char seq[200000];                        // sequence to scan. NOTE: maximum of 200,000 nucleotides in current implementation
  char temptext[5000];
  char upstream_sequence_filename[200];    // file name for sequences
  char weight_matrix_filename[200];        // file name for weight matrix
  float **weight_matrix;                   // weight matrix
  float *indices;                          // used in call to quicksorti
  int curr_index;                          // current index (used when printing the strand information)
  float temp_float;                        // temporary float value
  float threshold;                         // threshold; classify as present if score>=threshold
  int *pos;
  int i,j;
  int itemp;                               // temporary integer value
  int max_hits=10000;                      // maximum number of hits
  int max_upstream_length;                 // maximum sequence length
  int n;
  int n_transcripts_hits=0;                // total number of transcripts with at least one occurrence
  int n_hits=0;                            // total number of occurrences
  int nm;
  int np;
  int rc=0;                                // if rc==1, scan also the reverse complementary sequence
  int sequence_length;                     // sequence length
  int sortpos=1;                           // if sortpos==1, then sort the positions
  time_t    time_now;
  long total_sequence_length=0;            // total sequence length analyzed
  float bp_per_hit;                        // average number of base pairs between hits
  int min_dist;                            // minimum distance between successive occurrences
  int temp_seq_nr;                         // used in call to scan_motif_v4

  sprintf(last_modified,"02_16_2004");
  if (argc<5) 
    usage();
 
  sprintf(upstream_sequence_filename,argv[1]);
  sprintf(weight_matrix_filename,argv[2]);
  motif_length=atoi(argv[3]);
  threshold=atof(argv[4]);
  if (argc>5) {
    max_upstream_length=atoi(argv[5]);
  } else {
    max_upstream_length=50000;
  }
  if (max_upstream_length<=motif_length) {
    printf("ERROR! max_upstream_length=%d\tmotif_length=%d\n",max_upstream_length,motif_length);
    exit(1);
  }
  if (argc>6)
    rc=atoi(argv[6]);
  if (argc>7) 
    sortpos=atoi(argv[7]);
  if (argc>8) 
    verbose=atoi(argv[8]);  

  min_dist=(int)ceil(0.5*motif_length);
  time(&time_now);     /* get time in seconds */
  sprintf(output_filename,"scan_motif_v%d.out.txt",prog_version);
  output_file=fopen(output_filename,"w");
  sprintf(strand_filename,"scan_motif_v%d.strand.txt",prog_version);
  strand_file=fopen(strand_filename,"w");

  sprintf(infotext,"%% scan_motif_v%d.c\n",prog_version);
  sprintf(temptext,"%% %s",asctime(localtime(&time_now)));strcat(infotext,temptext);
  sprintf(temptext,"%% input parameters:\n");strcat(infotext,temptext);
  sprintf(temptext,"%% upstream_sequence_filename=%s\n",upstream_sequence_filename);strcat(infotext,temptext);
  sprintf(temptext,"%% weight_matrix_filename=%s\n",weight_matrix_filename);strcat(infotext,temptext);
  sprintf(temptext,"%% motif_length=%d\n",motif_length);strcat(infotext,temptext);
  sprintf(temptext,"%% threshold=%.4f\n",threshold);strcat(infotext,temptext);
  sprintf(temptext,"%% verbose=%d\n",verbose);strcat(infotext,temptext);
  sprintf(temptext,"%% maximum_upstream_length=%d\n",max_upstream_length);strcat(infotext,temptext);
  sprintf(temptext,"%% rc=%d\n",rc);strcat(infotext,temptext);
  sprintf(temptext,"%% sortpos=%d\n",sortpos);strcat(infotext,temptext);
  sprintf(temptext,"%% min_dist=%d\n",min_dist);strcat(infotext,temptext);
  sprintf(temptext,"%% last modified=%ld",last_modified);strcat(infotext,temptext);
  fprintf(output_file,"%s\n",infotext);
  fprintf(strand_file,"%s\n",infotext);

  if (verbose>=0) {
    sprintf(log_filename,"scan_motif_v%d.log.txt",prog_version);
    log_file=fopen(log_filename,"w");
  } else {
    log_file=0;
  }
  if (verbose>=0)
    printinfo(infotext,verbose,log_file);

  /* memory allocation */
  num_upstream_sequence=cvector(1,max_upstream_length);
  pos=ivector(1,max_hits);
  strand=cvector(1,max_hits);
  if (sortpos==1) {
    indices=vector(1,max_hits);
  }

  /**************************/
  /* load the weight matrix */
  /**************************/
  weight_matrix=matrix(1,4,1,motif_length);
  load_weight_matrix(weight_matrix_filename,motif_length,weight_matrix);
  if (verbose>=0) {
    sprintf(infotext,"%% weight_matrix\n%%");
    for (i=1;i<=4;i++) {
      for (j=1;j<=motif_length;j++) {
	sprintf(temptext,"%.4f\t",weight_matrix[i][j]);strcat(infotext,temptext);
      }
      strcat(infotext,"\n%%");
    }
    printinfo(infotext,verbose,log_file);
  }

  seq_nr=0;
  upstream_sequence_file=fopen(upstream_sequence_filename,"r");
  if (!upstream_sequence_file) {
    printf("ERROR! I could not find the file %s\n",upstream_sequence_filename);
    exit(1);
  }
  //read_fasta_seq2(upstream_sequence_file,info1,seq); // read first information line
  read_fasta_seq3(upstream_sequence_file,info1,seq); // read first information line

  while (feof(upstream_sequence_file)==0) {

    //read_fasta_seq2(upstream_sequence_file,info2,seq);
    read_fasta_seq3(upstream_sequence_file,info2,seq);
    sprintf(info1,"%s",info2);
    sprintf(info2,"");

    seq_nr++;

    sequence_length=strlen(seq);

    if (sequence_length>motif_length) {
      if (sequence_length>max_upstream_length) {
	itemp=sequence_length-max_upstream_length+1;
	total_sequence_length+=max_upstream_length;
      } else {
	itemp=0;
	total_sequence_length+=sequence_length;
      }
      
      seq2num(seq,sequence_length,num_upstream_sequence,itemp);
      
      temp_seq_nr=verbose*seq_nr;
      if (debug_mode)
	printf("processing temp_seq_nr=%d\n",temp_seq_nr);
      np=scan_motif4(weight_matrix,num_upstream_sequence,motif_length,sequence_length,threshold,pos,strand,rc,min_dist,log_file,temp_seq_nr);
      /* if (verbose>=0) {
	 np=scan_motif3(weight_matrix,num_upstream_sequence,motif_length,sequence_length,threshold,pos,strand,rc,log_file,seq_nr);
	 } else {
	 np=scan_motif2(weight_matrix,num_upstream_sequence,motif_length,sequence_length,threshold,pos,strand,rc);
	 }
      */
	
      if (np>0) {
	if (sortpos==1) {
	  for (i=1;i<=np;i++) 
	    indices[i]=(float)i;
	  quicksorti(np,pos,indices);
	}
	
	n_transcripts_hits++;
	n_hits+=np;

	/* print positions */
	sprintf(infotext,"%d\t",seq_nr);
	for (i=1;i<=np;i++) {
	  sprintf(temptext,"%d\t",pos[i]);strcat(infotext,temptext);
	}
	printinfo(infotext,verbose,log_file);
	fprintf(output_file,"%s\n",infotext);
	/* print strand information onlyt to file */
	sprintf(infotext,"%d\t",seq_nr);
	for (i=1;i<=np;i++) {
	  curr_index=(int)indices[i];
	  //sprintf(temptext,"%d\t",strand[i]);strcat(infotext,temptext);
	  sprintf(temptext,"%d\t",strand[curr_index]);strcat(infotext,temptext);
	}
	fprintf(strand_file,"%s\n",infotext);
      }  // close check on np>0
    }    // close check on if sequence_length>0
  }      // another sequence please
  fclose(upstream_sequence_file);
  fclose(output_file);
  fclose(strand_file);
  
  free_matrix(weight_matrix,1,4,1,motif_length);
  free_cvector(num_upstream_sequence,1,max_upstream_length);
  free_ivector(pos,1,max_hits);
  free_cvector(strand,1,max_hits);
  if (sortpos==1)
    free_vector(indices,1,max_hits); // printf("indices\n");

  bp_per_hit=(float)total_sequence_length/(float)n_hits;

  printf("n_transcripts_processed=%d\n",seq_nr);
  printf("n_transcripts_hits=%d\n",n_transcripts_hits);
  printf("n_hits=%d\n",n_hits);
  printf("total_sequence_length=%ld\n",total_sequence_length);
  printf("bp_per_hit=%.2f\n",bp_per_hit);
  printf("system_status=0\n");
  if (verbose>=0) {
    fprintf(log_file,"system_status=0\n");
    fprintf(log_file,"n_transcripts_processed=%d\n",seq_nr);
    fprintf(log_file,"n_transcripts_hits=%d\n",n_transcripts_hits);
    fprintf(log_file,"n_hits=%d\n",n_hits);
    fprintf(log_file,"total_sequence_length=%ld\n",total_sequence_length);
    fprintf(log_file,"bp_per_hits=%.2f\n",bp_per_hit);
    fclose(log_file);
  }
  return 0;
}

void usage(void)
{
  printf("scan_motif_v%d.c\n",prog_version);
  printf("usage\n");
  printf("scan_motif_v%d <sequence> <weight_matrix> <motif_length> <threshold> <max_upstream_length> <rc> <sortpos> <verbose>\n",prog_version);
  printf("\nscans for occurrences of the motif given by the <weight_matrix> file (4 x motif_length) in the <background_sequence> list of sequences\n");
  printf("considers a Stormo-Fields type of score based on the log-odds and uses <threshold> to classify as present/absent\n");
  printf("analyses only the <max_upstream_length> nucleotides to the 3' end (the entire sequence up to 10k if this is ommited)\n");
  printf("<n_seqs4motif> is the number of sequences used to compute the weight matrix\n");
  printf("<max_upstream_length> is the maximum upstream length to consider\n");
  printf("<rc> if rc==1 then scan the reverse complementary sequence as well (default=1)\n");
  printf("<sortpos> if sortpos==1 then sort positions (default=1)\n");
  printf("prints miscellaneous information if verbose=1 (default=0)\n");
  printf("if verbose<2 --> output to log file\n");
  printf("if verbose>0 --> output to screen\n");
  printf("\noutput stored in scan_motif_v%d.out.txt (positions) and scan_motif_v%d.strand.txt (strand)\n\n",prog_version,prog_version);
  gk();
  exit(1);
}
