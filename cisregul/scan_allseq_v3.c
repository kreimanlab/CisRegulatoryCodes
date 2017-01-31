/* scan_allseq.c (version 3)
 * scan one motif across all sequences
 * here we use the files with sequences (upstream, first exon, first intron), separated by chromosome
 */

/* some conventions about the positions:
 *   i) positions in the upstream sequences are > 0
 *  ii) positions downstream of the TSS are < 0
 * iii) if x_i > x_j then x_i is upstream of x_j
 */

/* comes from scan_1m_allseq_v6 */

/* new in version 2:
 * 03-10-2004: call scan_motif6
 */

/* new in version 3:
 * 03-10-2004: allow to provide an input list of motifs to search 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "my_config.c"
#include "../lib/nr.c"
#include "../lib/my_math_package.c"
#include "../lib/my_file_io.c"
#include "regul_search_methods_v7.c"
#include "../lib/my_text_package.c"
//#include "../lib/environment_variables.c"

/* functions */
void usage(void);
int stable_locus_id(char *txtline,char *what2search);
void read_one_genes_file_line(FILE *genes_file);
void parse_line(char *txtline);

/* global variables */
/* char CODE_DIR[200];
   char TEMP_DIR[200];
   char DATA_DIR[200];
*/
char *CODE_DIR;
char *TEMP_DIR;
char *DATA_DIR;
char program_name[200];
int prog_version=3;
char last_modified[200];
char id[1000];               // gene entry id
int current_chr_genes_file;  // current chromosome entry from genes_file
int upstream_length;         // lengths as obtained from the genes_file
int exon_length;
int intron_length;
int debug_mode=0;            // 1 for debug mode

main(int argc, char *argv[])
{
  FILE *exon_file;                           // handle to file with exon information
  FILE *genes_file;                          // handle to file with genes information
  FILE *intron_file;                         // handle to file with intron information
  FILE *log_file;                            // handle to log file
  FILE *out_file;                            // scan positions file
  FILE *strand_file;                         // scan strand file
  FILE *upstream_file;                       // file with sequences
  char chr_name[10];                         // chromosome (text name version)
  char exon_filename[200];                   // file with exon information
  char exon_seq[10000];                      // exon sequence
  char genes_filename[200];                  // list of genes
  char info1_exon[500];                      // exon information line
  char info1_intron[500];                    // intron information line
  char info1_upstream[500];                  // upstream information line
  char info2_exon[500];                      // next sequence exon information
  char info2_intron[500];                    // next sequence intron information
  char info2_upstream[500];                  // next sequence upstream information
  char intron_filename[200];                 // file with intron information
  char intron_seq[10000];                    // intron sequence
  char log_filename[200];                    // log filename
  char strand_filename[200];                 // strand filename
  char motif_filename[200];                  // motif weight matrix filename
  FILE *motif_file;                          // file with the motifs
  char numseq[10001];                        // sequence in numeric format (1234)
  char revnumseq[10001];                     // allocate space for the reverse complement of the sequence (the reverse complement is done in scan_motif_v5)
  char output_filename[200];                 // filename for scan positions
  char infotext[10000];                      // information text (for output)
  char strandtext[10000];                    // used to print out the strand information
  char headertext[10000];                    // header text
  char *seq;                                 // sequence
  char species[50];                          // species (typically a two character variable)
  char upstream_filename[200];               // filename with upstream information
  char upstream_seq[10000];                  // upstream sequence
  float **weight_matrix;                     // weight matrix
  float *mean_weight_matrix;                 // mean for each nucleotide (this is used when there is an "n")
  float *max_weight_matrix;                  // maximum for each nucleotide 
  float motif_threshold;                     // threshold for scanning
  int motif_threshold_column;                // column containing the threshold in the information line 
  char *strand_exon;                         // strand in exon
  char *strand_intron;                       // strand in intron
  char *strand_upstream;                     // strand in upstream
  int n_pos;                                 // number of entries (used in call to remove_vector_overlap)
  int *pos;                                  // generic position used in call to remove_vector_overlap
  int *pos_exon;                             // positions in exon sequences
  int *pos_intron;                           // positions in intron sequences
  int *pos_upstream;                         // positions in upstream sequences
  int actual_exon_length;                    // actual exon length
  int actual_intron_length;                  // actual intron length
  int actual_upstream_length;                // actual upstream length
  int chr;                                   // chromosome number
  int current_exon_length;                   // current exon length (trimmed if longer than maximum allowed)
  int current_intron_length;                 // current intron length (trimmed if longer than maximum allowed)
  int current_pos;                           // current position
  int current_upstream_length;               // current upstream length (trimmed if longer than maximum allowed)
  int i,j;                                   // loop variabels
  int is_locus_id;                           // indicator variable, 1 for locus_id genes
  int locuslink_only=1;                      // consider only entries with a locus link ID if locuslink_only=1
  int max_exon_length=5000;                  // maximum allowed exon length
  int max_hits=50000;                        // maximum number of hits per transcript
  int max_intron_length=5000;                // maximum allowed intron length
  int max_upstream_length=5000;              // maximum allowed upstream length
  int motif_length;                          // motif length
  int n_chromosomes;                         // total number of chromosomes
  int n_exon;                                // number of occurrences in exon sequence in current transcript
  int n_intron;                              // number of occurrences in intron sequence in current transcript
  int n_upstream;                            // number of occurrences in upstream sequence in current transcript
  int rc=1;                                  // consider both strands if rc=1
  int transcript=0;                          // transcript number (corresponding to the entries in the genes_file)
  int chr_transcript;                        // transcript number (within each chromosome; this is only used in the log file...)
  int verbose=2;                             // verbose=0 --> to screen, verbose=2 --> to log file, verbose=1 --> to both (in addition to out file with positions only)
  struct timeb t1,t2;                        // computation time
  long difftime;                             // used to evaluate computation time
  time_t time_now;                           // current time (used in header information)
  char temptext[5000];                       // temporary text
  int sortpos=1;                             // if sortpos==1, then sort positions
  float *indices;                            // used in call to quicksorti
  int curr_index;                            // used when printing the strand information
  int min_dist=0;                            // minimum distance between adjacent motif occurrences
  float min_dist_factor=0.5;                 // min_dist=motif_length*min_dist_factor (if this is <=0 then no overlaps are removed)   
  int progress_report_int=1000;              // report progress every progress_report_int transcripts
  int n_transcripts;                         // total number of transcripts
  int n_transcripts_processed=0;             // total number of transcripts scanned (differs from the above for example if only locusID genes are scanned)
  int n_transcripts_hits=0;                  // number of transcripts with at least one occurrence
  int n_transcripts_hits_up=0;               // number of transcripts wtih at least one occurrence in upstream sequences
  int n_transcripts_hits_ex=0;               // ", exon sequences
  int n_transcripts_hits_in=0;               // ", intron sequences
  //int seq_report=1000;                       // report one every seq_report sequences to check that everything goes ok
  char *search_string;                       // used for sanity check of informatino content
  int n_occurrences=0;                       // cumulative number of occurrences
  int n_occurrences_up=0;                    // cumulative number of occurrences in upstream sequences
  int n_occurrences_ex=0;                    // ", exon
  int n_occurrences_in=0;                    // ", intron
  long cumulative_ul=0;                      // cumulative upstream length
  long cumulative_el=0;                      // cumulative exon length
  long cumulative_il=0;                      // cumulative intron length
  long cumulative_length;                    // cumulative sequence length
  float p_per_transcript_up;                 // proportion of occurrences per transcript, upstream
  float p_per_transcript_ex;                 // ", exon
  float p_per_transcript_in;                 // ", intron
  float p_per_transcript;                    // ", all
  float freq_up;                             // number of base pairs between occurrences, upstream
  float freq_ex;                             // ", exon
  float freq_in;                             // ", intron
  float freq;                                // ", all
  int singlefile=0;                          // if singlefile == 1, then there is a single file with all the transcripts for all chromosomes (one file for upstream sequences, one for first exon, one for first intron)
  int curr_arg;                              // used to retrieve the arguments
  char seq_source[200];                      // gbk, ensembl, etc
  char seq_dir[200];                         // upstream_sequences if source is gbk, ensembl if source is ensembl
  int read_all=0;                            // if read_all == 1, then read all lines, even those that show length = 0
  char envname[200];                         // used to retrieve the environment variables
  char **upstream_sequences;                 // store all upstream sequences
  char **exon_sequences;                     // store all exon sequences
  char **intron_sequences;                   // store all intron sequences
  int *upstream_sequences_in;                // 0 to exclude and actual upstream sequence length to include
  int *exon_sequences_in;                    
  int *intron_sequences_in;
  int max_n_transcripts=20000;               // used for memory allocation
  float memory_req;                          // memory requirement for allocation of a variable 
  float cumul_memory_req=0.0;                // cumulative memory requirement
  int curr_transcript;                       // current transcript
  int *transcript_map;                       // transcript map (because many entries do not have a locus id)
  int *v;                                    // used to split the information line
  int curr_motif;                            // motif number
  char txtline[10000];                       // read one line of the motif input file
  int n_chars=10000;                         // number of chars to read
  int header;                                // 1 if the current line is a header line
  float temp_float;                          // used while reading the motifs file
  int temp_int;                              // used while reading interesting_filename
  char * pch;                                // used by the string tokenizer
  char file_prefix[200];                     // prefix for the output files (including directory information)
  int process_motif;                         // 0 to skip current motif
  int init_motif=1;                          // initial motif
  int finit_motif=1000;                      // final motif
  float *cum_threshold;
  char interesting_filename[200];            // filename with list of motifs to scan
  FILE *interesting_file;                    // used to read interesting_filename;
  int *interesting_list;                     // list of motifs to scan
  int *interesting_indices;                  // used in call to findi
  int n_interesting;                         // whether it is in the interesting list or not
  int interesting=0;                         // 1 if intersesting_filename is an input and 0 otherwise (if 0, scan all)

  sprintf(last_modified,"03_10_2004");
  if (argc<5) 
    usage();

  printf("getting environment variables...\n");
  sprintf(envname,"DATADIR");
  DATA_DIR=getenv(envname);
  sprintf(envname,"CODEDIR");
  CODE_DIR=getenv(envname);
  sprintf(envname,"TEMPDIR");
  TEMP_DIR=getenv(envname);
  if (verbose>=0) {
    printf("DATA_DIR=%s\n",DATA_DIR);
    printf("TEMP_DIR=%s\n",TEMP_DIR);
    printf("CODE_DIR=%s\n",CODE_DIR);
  }
  sprintf(program_name,"scan_allseq_v%d",prog_version);
 
  curr_arg=1;sprintf(motif_filename,argv[curr_arg]);                     // motif weight matrix filename
  curr_arg++;sprintf(species,argv[curr_arg]);                            // species
  curr_arg++;motif_threshold_column=atoi(argv[curr_arg]);                // threshold for scanning
  curr_arg++;sprintf(file_prefix,argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_upstream_length=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_exon_length=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_intron_length=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    locuslink_only=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    rc=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    sortpos=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    min_dist_factor=atof(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    singlefile=atoi(argv[curr_arg]);
  sprintf(seq_source,"gbk");
  curr_arg++;
  if (argc>curr_arg)
    sprintf(seq_source,argv[curr_arg]);
  sprintf(seq_dir,"%s",seq_source);
  if (strcmp(seq_source,"gbk")==0) 
    sprintf(seq_dir,"upstream_sequences");
  if (strcmp(seq_source,"sgd")==0)
    sprintf(seq_dir,"upstream_sequences");
  curr_arg++;
  if (argc>curr_arg) 
    read_all=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    max_n_transcripts=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    init_motif=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg)
    finit_motif=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) 
    verbose=atoi(argv[curr_arg]);
  curr_arg++;
  if (argc>curr_arg) {
    sprintf(interesting_filename,argv[curr_arg]);
    interesting=1;
  }

  printf("read %d arguments\n",argc);
  sprintf(genes_filename,"%s/db/%s/%s/lists/exon_list_v3.txt",DATA_DIR,seq_dir,species);

  if (singlefile==0) {
    if ( (strcmp(species,"hs")==0) )
      n_chromosomes=24;
    if ( (strcmp(species,"mm")==0) )
      n_chromosomes=21;
  } else {
    n_chromosomes=1;             // just go through the loop once
  }
  if (verbose>=0) {
    sprintf(log_filename,"%s.log.txt",program_name);
    log_file=fopen(log_filename,"w");
    if (!log_file) {
      printf("i could not open the file %s for writing\n",log_filename);
      exit(1);
    }
  }

  /* print out parameters */
  time(&time_now);
  sprintf(infotext,"%% scan_allseq_v%d.c [last modified=%s]\n",prog_version,last_modified);
  sprintf(temptext,"%% %s",asctime(localtime(&time_now)));strcat(infotext,temptext);      
  sprintf(temptext,"%% motif_filename=\t%s\n",motif_filename);strcat(infotext,temptext);
  sprintf(temptext,"%% species=%s\n",species);strcat(infotext,temptext);
  //sprintf(temptext,"%% motif_threshold=%.2f\n",motif_threshold);strcat(infotext,temptext);
  //sprintf(temptext,"%% motif_length=%d\n",motif_length);strcat(infotext,temptext);
  sprintf(temptext,"%% max_upstream_length=%d\n",max_upstream_length);strcat(infotext,temptext);
  sprintf(temptext,"%% max_exon_length=%d\n",max_exon_length);strcat(infotext,temptext);
  sprintf(temptext,"%% max_intron_length=%d\n",max_intron_length);strcat(infotext,temptext);
  sprintf(temptext,"%% locuslink_only=%d\n",locuslink_only);strcat(infotext,temptext);
  sprintf(temptext,"%% rc=%d\n",rc);strcat(infotext,temptext);
  sprintf(temptext,"%% sortpos=%d\n",sortpos);strcat(infotext,temptext);
  sprintf(temptext,"%% file_prefix=%s\n",file_prefix);strcat(infotext,temptext);
  sprintf(temptext,"%% min_dist_facotr=%.2f\n",min_dist_factor);strcat(infotext,temptext);
  sprintf(temptext,"%% min_dist=%d\n",min_dist);strcat(infotext,temptext);
  if (singlefile==0) {
    sprintf(temptext,"%% n_chromosomes=%d\n",n_chromosomes);strcat(infotext,temptext);
  } else {
    sprintf(temptext,"%% n_chromosomes=%d (single file with all chromosomes)\n",n_chromosomes);strcat(infotext,temptext);
  }
  sprintf(temptext,"%% singlefiles=%d\n",singlefile);strcat(infotext,temptext);
  sprintf(temptext,"%% seq_source=%s\n",seq_source);strcat(infotext,temptext);
  sprintf(temptext,"%% read_all=%d\n",read_all);strcat(infotext,temptext);
  sprintf(temptext,"%% seq_dir=%s\n",seq_dir);strcat(infotext,temptext);
  sprintf(temptext,"%% max_n_transcripts=%d\n",max_n_transcripts);strcat(infotext,temptext);
  if (interesting>0) {
    sprintf(temptext,"%% interesting_filename=%s\n",interesting_filename);strcat(infotext,temptext);
  }
  sprintf(temptext,"%% verbose=%d\n",verbose);strcat(infotext,temptext);
  sprintf(temptext,"%% debug_mode=%d\n",debug_mode);strcat(infotext,temptext);
  sprintf(temptext,"%%");strcat(infotext,temptext);
  strcpy(headertext,infotext);

  if (verbose>=0)
    printinfo(infotext,verbose,log_file);

  /* memory allocation */
  memory_req=(float)2*max_n_transcripts*max_upstream_length/(float)1000.0;cumul_memory_req+=memory_req;
  printf("upstream_sequences [%d x %d, %.1f kb / %.1f kb]\n",max_n_transcripts,max_upstream_length,memory_req,cumul_memory_req);
  upstream_sequences=cmatrix(1,max_n_transcripts,1,max_upstream_length);
  memory_req=(float)2*max_n_transcripts*max_exon_length/(float)1000.0;cumul_memory_req+=memory_req;
  printf("exon_sequences [%d x %d, %.1f kb / %.1f kb]\n",max_n_transcripts,max_exon_length,memory_req,cumul_memory_req);
  exon_sequences=cmatrix(1,max_n_transcripts,1,max_exon_length);
  memory_req=(float)2*max_n_transcripts*max_intron_length/(float)1000.0;cumul_memory_req+=memory_req;
  printf("intron_sequences [%d x %d, %.1f kb / %.1f kb]\n",max_n_transcripts,max_intron_length,memory_req,cumul_memory_req);
  intron_sequences=cmatrix(1,max_n_transcripts,1,max_intron_length);
  memory_req=(float)4*max_n_transcripts/(float)1000.0;cumul_memory_req+=memory_req;
  printf("upstream_sequences_in,exon_sequences_in,intron_sequences_in [%d x 1, %.1f kb / %.1f kb]\n",max_n_transcripts,memory_req,cumul_memory_req);
  upstream_sequences_in=ivector(1,max_n_transcripts);exon_sequences_in=ivector(1,max_n_transcripts);intron_sequences_in=ivector(1,max_n_transcripts);
  transcript_map=ivector(1,max_n_transcripts);
  for (i=1;i<=max_n_transcripts;i++) {
    upstream_sequences_in[i]=0;
    exon_sequences_in[i]=0;
    intron_sequences_in[i]=0;
    transcript_map[i]=0;
  }
  memory_req=(float)6*4*max_hits/(float)1000.0;cumul_memory_req+=memory_req;
  printf("pos_upstream/strand_upstream [%d x 1, %.1f kb / %.1f kb]\n",max_hits,memory_req,cumul_memory_req);pos_upstream=ivector(1,max_hits);strand_upstream=cvector(1,max_hits);
  printf("pos_exon/strand_exon [%d x 1, %.1f kb / %.1f kb]\n",max_hits,memory_req,cumul_memory_req);pos_exon=ivector(1,max_hits);strand_exon=cvector(1,max_hits);
  printf("pos_intron/strand_intron [%d x 1, %.1f kb / %.1f kb]\n",max_hits,memory_req,cumul_memory_req);pos_intron=ivector(1,max_hits);strand_intron=cvector(1,max_hits);
  memory_req=(float)2*1000000/(float)1000.0;cumul_memory_req+=memory_req;
  printf("seq [char, %d x 1, %.1f kb / %.1f kb]\n",1000000,memory_req,cumul_memory_req);seq=cvector(0,1000000);
  if (sortpos==1) 
    indices=vector(1,max_hits);
  pos=ivector(1,max_hits);

  if (interesting>0) {
    /********************************/
    /* reading interesting_filename */
    /*******************************/
    interesting_file=fopen(interesting_filename,"r");
    if (!interesting_file) {
      printf("error! i could not open %s for reading. exiting...\n",interesting_filename);
      exit(1);
    }
    interesting_list=ivector(1,2000);
    interesting_indices=ivector(1,2000);
    i=0;
    while (feof(interesting_file)==0) {
      fscanf(interesting_file,"%d\n",&temp_int);
      i++;
      interesting_list[i]=temp_int;
    }    
    n_interesting=i;
    fclose(interesting_file);
  }

  /*************************/
  /* reading all sequences */
  /*************************/
  if (verbose>=0)
    fprintf(log_file,"reading all sequences...\n");
  for (chr=1;chr<=n_chromosomes;chr++) {                  // debug: for (chr=1;chr<=1;chr++) {
    chr_transcript=0;
    sprintf(chr_name,"%d",chr);
    if ( (strcmp(species,"hs")==0) & (chr==23) ) 
      sprintf(chr_name,"X");
    if ( (strcmp(species,"hs")==0) & (chr==24) ) 
      sprintf(chr_name,"Y");
    if ( (strcmp(species,"mm")==0) & (chr==20) ) 
      sprintf(chr_name,"X");
    if ( (strcmp(species,"mm")==0) & (chr==21) ) 
      sprintf(chr_name,"Y");
    if (debug_mode==1) {
      printf("chr_name=%s\n",chr_name);
      fprintf(log_file,"chr_name=%s\n",chr_name);
    }
    if (verbose>=0)
      fprintf(log_file,"\treading chr=%d\t(%s)\n",chr,chr_name);
    /*******************************************/
    /* open fasta files for current chromosome */
    /*******************************************/
    if (max_upstream_length>0) {      // open upstream file for current chromosome
      if (singlefile==0) {
	sprintf(upstream_filename,"%s/db/%s/%s/%s_chr%s.%s.upstream.fasta",DATA_DIR,seq_dir,species,species,chr_name,seq_source);
      } else {
	sprintf(upstream_filename,"%s/db/%s/%s/%s.%s.upstream.fasta",DATA_DIR,seq_dir,species,species,seq_source);
      }	
      upstream_file=fopen(upstream_filename,"r");
      if (!upstream_file) {
	printf("error! i could not open the file upstream_filename=%s\n",upstream_filename);
	exit(1);
      }
      read_fasta_seq3(upstream_file,info1_upstream,upstream_seq); // read first information line
      if (debug_mode==1) {
	printf("upstream_filename=%s\n",upstream_filename);fprintf(log_file,"upstream_filename=%s\n",upstream_filename);
      }
    }
    if (max_exon_length>0) {      // open exon file for current chromosome
      if (singlefile==0) {
	sprintf(exon_filename,"%s/db/%s/%s/%s_chr%s.%s.exon1.fasta",DATA_DIR,seq_dir,species,species,chr_name,seq_source);
      } else {
	sprintf(exon_filename,"%s/db/%s/%s/%s.%s.exon1.fasta",DATA_DIR,seq_dir,species,species,seq_source);
      }
      exon_file=fopen(exon_filename,"r");
      read_fasta_seq3(exon_file,info1_exon,exon_seq); // read first information line
      if (debug_mode==1) {
	printf("exon_filename=%s\n",exon_filename);fprintf(log_file,"exon_filename=%s\n",exon_filename);
      }
    }
    if (max_intron_length>0) {      // open intron file for current chromosome
      if (singlefile==0) {
	sprintf(intron_filename,"%s/db/%s/%s/%s_chr%s.%s.intron1.fasta",DATA_DIR,seq_dir,species,species,chr_name,seq_source);
      } else {
	sprintf(intron_filename,"%s/db/%s/%s/%s.%s.intron1.fasta",DATA_DIR,seq_dir,species,species,seq_source);
      }
      intron_file=fopen(intron_filename,"r");
      read_fasta_seq2(intron_file,info1_intron,intron_seq); // read first information line
      if (debug_mode==1) {
	printf("intron_filename=%s\n",intron_filename);fprintf(log_file,"intron_filename=%s\n",intron_filename);
      }
    }

    /**********************/
    /* opening genes_file */
    /**********************/
    current_chr_genes_file=1;
    genes_file=fopen(genes_filename,"r");
    check_file_handle(genes_file,genes_filename);
    if (singlefile==1) 
      current_chr_genes_file=1;
    printf("\tprocesssing file = %s\n",genes_filename);

    while ( (current_chr_genes_file<=chr) & (feof(genes_file)==0) ) {
      read_one_genes_file_line(genes_file);     // results of parsing one line are stored in global variables: current_chr_genes_file, id, upstream/exon/intron_length;
      if (singlefile==1)
	current_chr_genes_file=1;
      if (current_chr_genes_file == chr) {
	transcript++;
	chr_transcript++;
	if ( (transcript % progress_report_int) == 0 )
	  printf("processing chr=%d\ttranscript=%d\n",chr,transcript);             // progress report 
	is_locus_id=stable_locus_id(id,"LocusID");	                           /* search for string LocusID */
	if ( (locuslink_only != 1) | (is_locus_id) ) {                             // if it is a locusID entry or if we do not care about this, process this entry
	  n_transcripts_processed++;
	  transcript_map[n_transcripts_processed]=transcript;
	  if ( (n_transcripts_processed % progress_report_int) == 0 ) 
	    printf("reading sequence = %d\n",n_transcripts_processed);
	  if (debug_mode==1) {
	    sprintf(infotext,"chr=%d\tchr_transcript=%d\ttranscript=%d\tis_locus_id=%d\tid=%s",chr,chr_transcript,transcript,is_locus_id,id);
	    printinfo(infotext,verbose,log_file);
	  }
	  /**********************/
	  /* upstream sequences */
	  /**********************/
	  if (max_upstream_length>0) {	                                         /* read upstream sequence */
	    read_fasta_seq3(upstream_file,info2_upstream,seq);
	    if (!seq) {
	      printf("end of upstream file for chr=%d\n",chr);
	    } else {
	      actual_upstream_length=strlen(seq);
	      if (actual_upstream_length>max_upstream_length) {                  // shorten upstream sequences, extract three prime segment
		substring(seq,upstream_seq,max_upstream_length,-1);              // extract max_upstream_seq nucleotides from the 3' end from the upstream_seq
		current_upstream_length=max_upstream_length;                     // set length to max
		strcpy(seq,upstream_seq);                                        // copy upstream_seq to seq
	      } else {
		current_upstream_length=actual_upstream_length;
	      }
	      seq2num(seq,current_upstream_length,numseq,0);	                 // convert letter sequence to numeric sequence   
	      for (i=1;i<=current_upstream_length;i++) 
		upstream_sequences[n_transcripts_processed][i]=numseq[i];
	      upstream_sequences_in[n_transcripts_processed]=current_upstream_length;
	    }
	    search_string=strstr(info1_upstream,id);   /* sanity check */
	    if (search_string==NULL) {
	      printf("ERROR!\ntranscript=%d\nn_transcripts_processed=%d\n",transcript,n_transcripts_processed);printf("id=%s\n",id);
	      printf("info1_upstream=%s\n",info1_upstream);printf("info2_upstream=%s\n",info2_upstream);printf("seq=%s\n",seq);
	      if ( (strlen(info1_upstream)>1) & (strlen(info2_upstream)>1) )      exit(1);
	    }
	    sprintf(info1_upstream,"%s",info2_upstream);
	    sprintf(info2_upstream,"");
	  }                                                                       // close check on max_upstream_length>0
	  /******************/
	  /* exon sequences */
	  /******************/
	  if ( (max_exon_length>0) & (exon_length>0) ) {	                  /* read exon sequence */
	    read_fasta_seq3(exon_file,info2_exon,seq);
	    if (!seq) {
	      printf("end of exon file for chr=%d\n",chr);
	    } else {
	      actual_exon_length=strlen(seq);
	      if (actual_exon_length>max_exon_length) {	                         // shorten exon sequences, extract five prime segment
		substring(seq,exon_seq,max_exon_length,1);                       // extract max_exon_length nucleotides from the 5' end from the exon_seq
		current_exon_length=max_exon_length;
		strcpy(seq,exon_seq);
	      } else {
		current_exon_length=actual_exon_length;
	      }
	      seq2num(seq,current_exon_length,numseq,0);	                 // convert letter sequence to numeric sequence
	      for (i=1;i<=current_exon_length;i++) 
		exon_sequences[n_transcripts_processed][i]=numseq[i];
	      exon_sequences_in[n_transcripts_processed]=current_exon_length;
	    }
	    search_string=strstr(info1_exon,id);    // sanity check
	    if (search_string==NULL) {
	      printf("ERROR!\ntranscript=%d\nn_transcripts_processed=%d\n",transcript,n_transcripts_processed);printf("id=%s\n",id);
	      printf("info1_exon=%s\n",info1_exon);printf("info2_exon=%s\n",info2_exon);printf("seq=%s\n",seq);printf("exon_length=%d\n",exon_length);
	      if ( (strlen(info1_exon)>1) & (strlen(info2_exon)>1) )
		exit(1);
	    }
	    sprintf(info1_exon,"%s",info2_exon);
	    sprintf(info2_exon,"");
	  } else {                                                               // if max_exon_length <=0 or if exon_length<=0
	    if ( (max_exon_length>0) & (read_all==1) ) {                         // if processing exons and empty sequences are present in the fasta file  
	      read_fasta_seq3(exon_file,info2_exon,seq);
	      sprintf(info1_exon,"%s",info2_exon);
	      sprintf(info2_exon,"");
	    }
	  }
	  /**********************/
	  /* intron sequences   */
	  /**********************/
	  if ( (max_intron_length>0) & (intron_length>0) ) {	                /* read intron sequence */
	    read_fasta_seq3(intron_file,info2_intron,seq);
	    if (!seq) {
	      printf("end of intron file for chr=%d\n",chr);
	    } else {
	      actual_intron_length=strlen(seq);
	      if (actual_intron_length>max_intron_length) {                     // shorten intron sequences, extract five prime segment
		substring(seq,intron_seq,max_intron_length,1);                  // extract max_intron_length nucleotides from the 5' end from the intron_seq
		current_intron_length=max_intron_length;
		strcpy(seq,intron_seq);
	      } else {
		current_intron_length=actual_intron_length;
	      }
	      seq2num(seq,current_intron_length,numseq,0);                       // convert to numeric sequence
	      for (i=1;i<=current_intron_length;i++) 
		intron_sequences[n_transcripts_processed][i]=numseq[i];
	      intron_sequences_in[n_transcripts_processed]=current_intron_length;
	      search_string=strstr(info1_intron,id);   // sanity check
	      if (search_string==NULL) {
		printf("ERROR!\ntranscript=%d\nn_transcripts_processed=%d\n",transcript,n_transcripts_processed);printf("id=%s\n",id);
		printf("info1_intron=%s\n",info1_intron);printf("info2_intron=%s\n",info2_intron);printf("seq=%s\n",seq);printf("intron_length=%d\n",intron_length);
		if ( (strlen(info1_intron)>1) & (strlen(info2_intron)>1) )
		  exit(1);
	      }
	      sprintf(info1_intron,"%s",info2_intron);
	      sprintf(info2_intron,"");
	    }                                                                     // close check on end of intron file	
	  } else {                                                                // close check on max_intron_length>0 and intron_lenght>0
	    if ( (max_intron_length>0) & (read_all==1) ) {                        // if processing introns and empty sequences are present in the fasta file 
	      read_fasta_seq3(intron_file,info2_intron,seq);
	      sprintf(info1_intron,"%s",info2_intron);
	      sprintf(info2_intron,"");
	    }
	  }                                                                       // close check on max_intron_length>0 and intron_length>0
	} else {                                                                  // else to locus link check
	  if (debug_mode==1) {
	    sprintf(infotext,"chr=%d\tchr_transcript=%d\ttranscript=%d\tis_locus_id=%d\tid=%s",chr,chr_transcript,transcript,is_locus_id,id);
	    printinfo(infotext,verbose,log_file);
	  }
	  /* read the files anyway so as to get the right indices matching the right sequences ... */
	  if (max_upstream_length>0) {
	    read_fasta_seq3(upstream_file,info2_upstream,seq);sprintf(info1_upstream,"%s",info2_upstream);sprintf(info2_upstream,"");
	  }
	  if ( (max_exon_length>0) & (exon_length>0) ) {
	    read_fasta_seq3(exon_file,info2_exon,seq);sprintf(info1_exon,"%s",info2_exon);sprintf(info2_exon,"");
	  }
	  if ( (max_intron_length>0) & (intron_length>0) ) {
	    read_fasta_seq3(intron_file,info2_intron,seq);sprintf(info1_intron,"%s",info2_intron);sprintf(info2_intron,"");
	  }
	}                                                                         // close check on locus link
      }                                                                           // check on chromosome from genes_file matching current chr
    }                                                                             // close reading genes_file
    printf("next chromosome please (closing chromosome %s files)\n",chr_name);
    /* close files */
    fclose(genes_file);
    if (max_upstream_length>0) 
      fclose(upstream_file);
    if (max_exon_length>0)
      fclose(exon_file);
    if (max_intron_length>0)
      fclose(intron_file);
  }                                                                               // close chr loop

  /*********************/
  /* sequence scanning */
  /*********************/
  motif_file=fopen(motif_filename,"r");
  if (!motif_file) {
    printf("ERROR! I could not open the file %s for reading; exiting...\n",motif_filename);exit(1);
  }
  curr_motif=0;
  while (feof(motif_file)==0) {
    ftime(&t1);
    curr_motif++;
    process_motif=1;
    if (curr_motif<init_motif)
      process_motif=0;
    if (curr_motif>finit_motif)
      process_motif=0;
    if (process_motif==1) {
      sprintf(strand_filename,"%s_%s_all_%d_%d_%d_%d_%d.strand.txt",file_prefix,species,max_upstream_length,max_exon_length,max_intron_length,curr_motif,motif_threshold_column);
      strand_file=fopen(strand_filename,"r");
      if (!strand_file) {
	printf("%s does not exist\n",strand_filename);
      } else {
	printf("%s already exists; skipping\n",strand_filename);
	process_motif=0;
      }
    }
    printf("processing motif = %d\n",curr_motif);
    if ( (debug_mode==1) & (verbose>=0) )
      fprintf(log_file,"processing motif=%d\n",curr_motif);
    fgets(txtline,n_chars,motif_file);
    header=is_header(txtline);
    if (!header) {
      printf("ERROR! I was expecting to see a header line:\n\tscan_allseq_v%d.c\n\tcurr_motif=%d\n\ttxtline=%s\n",prog_version,curr_motif,txtline);exit(1);
    } else { 
      // get information
      pch = strtok (txtline,"\t");
      i=0;
      while (pch != NULL) {
	if (i==2) {
	  if (pch) 
	    motif_length=atoi(pch);
	}
	if (i==motif_threshold_column) {
	  if (pch)
	    motif_threshold=atof(pch);
	  }
	pch = strtok(NULL,"\t");
	i++;
      }  
      if (min_dist_factor>0)
	min_dist=(int)ceil(motif_length*min_dist_factor);
      // read 4 more lines      
      weight_matrix=matrix(1,4,1,motif_length);
      max_weight_matrix=vector(1,motif_length);
      cum_threshold=vector(1,motif_length);
      mean_weight_matrix=vector(1,motif_length);
      for (j=1;j<=motif_length;j++)
	mean_weight_matrix[j]=0.0;
      for (i=1;i<=4;i++) {
	for (j=1;j<=motif_length;j++) {
	  fscanf(motif_file,"%f\n",&temp_float);
	  weight_matrix[i][j]=temp_float;
	  mean_weight_matrix[j]=mean_weight_matrix[j]+temp_float;
	}
      }
      for (j=1;j<=motif_length;j++)
	mean_weight_matrix[j]=mean_weight_matrix[j]/4.0;
      weight_matrix_max(weight_matrix,motif_length,max_weight_matrix);
      temp_threshold(max_weight_matrix,motif_length,motif_threshold,cum_threshold);
	
      if (process_motif==1) {
	// check whether it is in the interesting_list or not
	process_motif=findi(interesting_list,n_interesting,curr_motif,interesting_indices);
	if (process_motif!=1) {
	  printf("i am sorry but i am currently not interested in motif %d\n",curr_motif);
	} else {
	  if (verbose>=0)  fprintf(log_file,"processing curr_motif=%d\n",curr_motif);
	  sprintf(output_filename,"%s_%s_all_%d_%d_%d_%d_%d.txt",file_prefix,species,max_upstream_length,max_exon_length,max_intron_length,curr_motif,motif_threshold_column);
	  printf("opening filename=%s\n",output_filename);
	  if (verbose>=0)  fprintf(log_file,"opening filename=%s\n",output_filename);
	  out_file=fopen(output_filename,"w");
	  if (!out_file) {
	    printf("i could not open the file %s for writing\n",output_filename);exit(1);
	  }
	  sprintf(strand_filename,"%s_%s_all_%d_%d_%d_%d_%d.strand.txt",file_prefix,species,max_upstream_length,max_exon_length,max_intron_length,curr_motif,motif_threshold_column);
	  printf("opening filename=%s\n",strand_filename);
	  if (verbose>=0)  fprintf(log_file,"opening filename=%s\n",strand_filename);
	  strand_file=fopen(strand_filename,"w");
	  if (!strand_file) {
	    printf("i could not open the file %s for writing\n",strand_filename);
	    exit(1);
	  }
	  fprintf(out_file,"%s\n",headertext);fprintf(strand_file,"%s\n",headertext);
	  sprintf(infotext,"%% motif_length=%d\n",motif_length);
	  sprintf(temptext,"%% motif_threshold=%.4f",motif_threshold);strcat(infotext,temptext);
	  fprintf(out_file,"%s\n",infotext);fprintf(strand_file,"%s\n",infotext);
	  
	  if (verbose>=0) {
	    fprintf(log_file,"%s\n",infotext);
	    sprintf(infotext,"weight_matrix\n");
	    for (i=1;i<=4;i++) {
	      for (j=1;j<=motif_length;j++) {
		sprintf(temptext,"%.4f\t",weight_matrix[i][j]);strcat(infotext,temptext);
	      }
	      strcat(infotext,"\n");
	    }
	    sprintf(temptext,"mean_weight_matrix\n");strcat(infotext,temptext);
	    for (j=1;j<=motif_length;j++) {
	      sprintf(temptext,"%.4f\t",mean_weight_matrix[j]);strcat(infotext,temptext);
	    }
	    printinfo(infotext,verbose,log_file);
	  }
	  
	  printf("sequence scanning [n_transcripts_processed=%d, motif_length=%d, motif_threshold=%.4f...\n",n_transcripts_processed,motif_length,motif_threshold);
	  /********/
	  /* scan */
	  /********/
	  for (curr_transcript=1;curr_transcript<=n_transcripts_processed;curr_transcript++) {
	    n_upstream=0;n_exon=0;n_intron=0;
	    if ( (curr_transcript % progress_report_int) == 0)
	      printf("scanning transcript = %d\n",curr_transcript);
	    current_upstream_length=upstream_sequences_in[curr_transcript];
	    if (current_upstream_length>0) {
	      for (i=1;i<=current_upstream_length;i++)
		numseq[i]=upstream_sequences[curr_transcript][i];
	      //n_upstream=scan_motif5(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_upstream_length,motif_threshold,pos_upstream,strand_upstream,rc,min_dist,log_file,-1);                                                                     // end with -1 means do not print the info to the log_file  
	      n_upstream=scan_motif6(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_upstream_length,motif_threshold,pos_upstream,strand_upstream,rc,min_dist,cum_threshold,log_file,-1);
	      cumulative_ul+=current_upstream_length;
	      if (n_upstream>0) {
		n_transcripts_hits_up++;
		n_occurrences_up+=n_upstream;
	      }
	    }
	    current_exon_length=exon_sequences_in[curr_transcript];
	    if (current_exon_length>0) {
	      for (i=1;i<=current_exon_length;i++)
		numseq[i]=exon_sequences[curr_transcript][i];
	    //n_exon=scan_motif5(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_exon_length,motif_threshold,pos_exon,strand_exon,rc,min_dist,log_file,-1);                                                                                     // end with -1 means do not print the info to the log_file  
	      n_exon=scan_motif6(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_exon_length,motif_threshold,pos_exon,strand_exon,rc,min_dist,cum_threshold,log_file,-1);
	      cumulative_el+=current_exon_length;
	      if (n_exon>0) {
		n_transcripts_hits_ex++;
		n_occurrences_ex+=n_exon;
	      }
	    }
	    current_intron_length=intron_sequences_in[curr_transcript];
	    if (current_intron_length>0) {
	      for (i=1;i<=current_intron_length;i++)
		numseq[i]=intron_sequences[curr_transcript][i];
	      //n_intron=scan_motif5(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_intron_length,motif_threshold,pos_intron,strand_intron,rc,min_dist,log_file,-1);                                                                             // end with -1 means do not print the info to the log_file  
	      n_intron=scan_motif6(weight_matrix,mean_weight_matrix,numseq,revnumseq,motif_length,current_intron_length,motif_threshold,pos_intron,strand_intron,rc,min_dist,cum_threshold,log_file,-1);
	      cumulative_il+=current_intron_length;
	      if (n_intron>0) {
		n_transcripts_hits_in++;
		n_occurrences_in+=n_intron;
	      }
	    }
	    /**************************/
	    /* print out scan results */
	    /**************************/
	    if ( (n_upstream>0) | (n_exon>0) | (n_intron>0) ) {
	      n_transcripts_hits++;
	      n_occurrences=n_occurrences+n_upstream+n_exon+n_intron;
	      sprintf(infotext,"%d\t",transcript_map[curr_transcript]);
	      sprintf(strandtext,"%d\t",transcript_map[curr_transcript]);
	      if (n_intron>0) {                                                               // convert positions (here we assume that all positions are > 0)
		if (sortpos==1) {
		  // for (i=1;i<=max_hits;i++) 
		  for (i=1;i<=n_intron;i++)
		    indices[i]=i; 
		  quicksorti(n_intron,pos_intron,indices);
		}
		for (i=1;i<=n_intron;i++) {
		  pos[i]=current_exon_length+current_intron_length-pos_intron[i];
		  sprintf(temptext,"%d\t",-pos[i]);strcat(infotext,temptext);
		  curr_index=indices[i];
		  sprintf(temptext,"%d\t",strand_intron[curr_index]);strcat(strandtext,temptext);
		}
	      }	    
	      if (n_exon>0) {	                                                            // convert positions (here we assume that all positions are > 0)
		if (sortpos==1) {
		  //for (i=1;i<=max_hits;i++)
		  for (i=1;i<=n_exon;i++)
		    indices[i]=i;
		  quicksorti(n_exon,pos_exon,indices);
		}
		for (i=1;i<=n_exon;i++) {
		  pos[i]=current_exon_length-pos_exon[i];
		  sprintf(temptext,"%d\t",-pos[i]);strcat(infotext,temptext);
		  curr_index=indices[i];
		  sprintf(temptext,"%d\t",strand_exon[curr_index]);strcat(strandtext,temptext);
		}
	      }
	      if (n_upstream>0) {
		if (sortpos==1) {
		  //for (i=1;i<=max_hits;i++)
		  for (i=1;i<=n_upstream;i++)
		    indices[i]=i;
		  quicksorti(n_upstream,pos_upstream,indices);
		}
		for (i=1;i<=n_upstream;i++) {
		  sprintf(temptext,"%d\t",pos_upstream[i]);strcat(infotext,temptext);
		  curr_index=indices[i];
		  sprintf(temptext,"%d\t",strand_upstream[curr_index]);strcat(strandtext,temptext);
		}
	      }
	      if (debug_mode==1) 
		printinfo(infotext,verbose,log_file);
	      fprintf(out_file,"%s\n",infotext);
	      fprintf(strand_file,"%s\n",strandtext);
	    } else {
	      if (debug_mode==1) {
		sprintf(infotext,"transcript=%d\tn_transcripts_processed=%d\tno scan hit\n",curr_transcript,n_transcripts_processed);
		printinfo(infotext,verbose,log_file);
	      }
	    }                                                                                   // close check on whether there was any hit in the scanning
	  }                                                                                     // close curr_transcript loop
	  n_transcripts=transcript;
	  /* print out some statistics */
	  sprintf(infotext,"n_transcripts=%d\n",n_transcripts);
	  sprintf(temptext,"n_transcripts_processed=%d\n",n_transcripts_processed);strcat(infotext,temptext);
	  sprintf(temptext,"n_transcripts_hits=%d\n",n_transcripts_hits);strcat(infotext,temptext);
	  if (verbose>=0)
	    printinfo(infotext,verbose,log_file);
	  printf("%s\n",infotext);                                                              // print to screen anyway because we may need to read these values in perl 
	  cumulative_length=cumulative_el+cumulative_il+cumulative_ul;
	  p_per_transcript_up=(float)n_occurrences_up/(float)n_transcripts_hits_up;
	  p_per_transcript_ex=(float)n_occurrences_ex/(float)n_transcripts_hits_ex;
	  p_per_transcript_in=(float)n_occurrences_in/(float)n_transcripts_hits_in;
	  p_per_transcript=(float)n_occurrences/(float)n_transcripts_hits;
	  freq_up=(float)cumulative_ul/(float)n_occurrences_up;
	  freq_ex=(float)cumulative_el/(float)n_occurrences_ex;
	  freq_in=(float)cumulative_il/(float)n_occurrences_in;
	  freq=(float)cumulative_length/(float)n_occurrences;
	  sprintf(infotext,"\tupstream\texon\tintron\ttotal\n");
	  sprintf(temptext,"n_transcripts_hits\t%d\t%d\t%d\t%d\n",n_transcripts_hits_up,n_transcripts_hits_ex,n_transcripts_hits_in,n_transcripts_hits);strcat(infotext,temptext);
	  sprintf(temptext,"n_occurrences\t%d\t%d\t%d\t%d\n",n_occurrences_up,n_occurrences_ex,n_occurrences_in,n_occurrences);strcat(infotext,temptext);
	  sprintf(temptext,"p_per_transcript\t%.2f\t%.2f\t%.2f\t%.2f\n",p_per_transcript_up,p_per_transcript_ex,p_per_transcript_in,p_per_transcript);strcat(infotext,temptext);
	  sprintf(temptext,"cumulative_len\t%ld\t%ld\t%ld\t%ld\n",cumulative_ul,cumulative_el,cumulative_il,cumulative_length);strcat(infotext,temptext);
	  sprintf(temptext,"frequency\t%.2f\t%.2f\t%.2f\t%.2f\n",freq_up,freq_ex,freq_in,freq);strcat(infotext,temptext);
	  if (verbose>=0) 
	    printinfo(infotext,verbose,log_file);
	  printf("%s\n",infotext);
	  /* cpu time */
	  ftime(&t2);
	  difftime=t2.time-t1.time;
	  printf("cpu time (sec) = %ld\n",difftime);
	  if (verbose>=0) 
	    fprintf(log_file,"cpu time (sec) = %ld\n",difftime);
	  fclose(out_file);
	  fclose(strand_file);
	}                                                                                      // close check on process_motif==1 (whether it is interesting)
      }                                                                                        // close check on process_motif==1
      free_matrix(weight_matrix,1,4,1,motif_length);
      free_vector(mean_weight_matrix,1,motif_length);
      free_vector(max_weight_matrix,1,motif_length);
      free_vector(cum_threshold,1,motif_length);
    }                                                                                          // close check on whether header line while reading the motif_file
  }                                                                                            // close check on feof(motif_file)

  /* free memory */
  free_cmatrix(upstream_sequences,1,max_n_transcripts,1,max_upstream_length);
  free_cmatrix(exon_sequences,1,max_n_transcripts,1,max_exon_length);
  free_cmatrix(intron_sequences,1,max_n_transcripts,1,max_intron_length);
  free_ivector(upstream_sequences_in,1,max_n_transcripts);
  free_ivector(exon_sequences_in,1,max_n_transcripts);
  free_ivector(intron_sequences_in,1,max_n_transcripts);
  free_ivector(transcript_map,1,max_n_transcripts);
  free_ivector(pos_upstream,1,max_hits);  //printf("pos_upstream\n");
  free_ivector(pos_exon,1,max_hits);      //printf("pos_exon\n");
  free_ivector(pos_intron,1,max_hits);    //printf("pos_intron\n");
  free_cvector(seq,0,1000000); //printf("seq\n");
  free_cvector(strand_upstream,1,max_hits); // printf("strand_upstream\n");
  free_cvector(strand_exon,1,max_hits);  //printf("strand_exon\n");
  free_cvector(strand_intron,1,max_hits); // printf("strand_intron\n");
  if (sortpos==1)
    free_vector(indices,1,max_hits); // printf("indices\n");
  free_ivector(pos,1,max_hits);  //printf("pos\n");
  if (interesting>0) {
    free_ivector(interesting_list,1,2000);
    free_ivector(interesting_indices,1,2000);
  }
  /* return system status 0=ok */
  printf("system_status=0\n");
  if (verbose>=0)
    fprintf(log_file,"system_status=0\n");

  /* close files */
  if (verbose>=0)
    fclose(log_file);

  return 0;
}  

void read_one_genes_file_line(FILE *genes_file)
{
  /* read one line from the genes_file
   * ignore header/comments lines
   * store results in global variables
   */

  int header=1;
  char txtline[100000];   
  int n_chars=100000;

  if (feof(genes_file) != 0) {
    current_chr_genes_file=-1;
  } else {
    /* read file header */
    while (header) {
      fgets(txtline,n_chars,genes_file);
      header=is_header(txtline);
    }
    
    parse_line(txtline);
  }
}
	
void parse_line(char *txtline)
{
  /* parse the contents of one entry from the file genes_file 
   * store the results in the global variables:
   * current_chr_genes_file (int)
   * exon_length (int)
   * intron_length (int)
   * upstream_length (int)
   * id   (char)
   */
  int token_counter=1;
  char * pch;
  //int current_chr_genes_file;
  //int exon_length=0;
  //int intron_length=0;
  //int upstream_length=0;

  //printf("txtline=%s\n",txtline);
  pch = strtok (txtline,"\t");
  current_chr_genes_file=atoi(pch);
  
  while (pch != NULL) {
    if (token_counter==6) {
      strcpy(id,pch);
    }
    if (token_counter==16) 
      exon_length=atoi(pch);
    if (token_counter==20)
      intron_length=atoi(pch);
    if (token_counter==24)
      upstream_length=atoi(pch);
    
    token_counter++;
    pch = strtok(NULL,"\t");
  }

  /* printf("id=%s\n",id);
     printf("current_chr_genes_file=%d\texon_length=%d\tintron_length=%d\tupstream_length=%d\n",current_chr_genes_file,exon_length,intron_length,upstream_length);
  */
}

int stable_locus_id(char *txtline,char *what2search) {
  char *search_string;
  int gotit;
  
  search_string=strstr(txtline,what2search);
  if (search_string==NULL) {
    gotit=0;
  } else {
    gotit=1;
  }
  
  return gotit;
}

void usage(void)
{
  printf("scan_allseq_v%d.c\n\n",prog_version);
  printf("usage:\n");
  printf("scan_allseq_v%d.exe <motif_filename> <species> <motif_threshold_column> <file_prefix> (<max_upstream_length> <max_exon_length> <max_intron_length> <locuslink_only> <rc> <sortpos> <min_dist_factor> <singlefile> <seq_source> <read_all> <max_n_transcripts> <init_motif> <finit_motif> <verbose> <interesting_filename>)\n\n",prog_version);
  printf("motif_filename\tfile with motif weight matrix (without any comment lines)\n");
  printf("species\ths or mm\n");
  printf("motif_threshold_column\tcolumn in the information line containing the threshold for motif occurrences (int)\n");
  printf("max_upstream_length\n");
  printf("max_exon_length\n");
  printf("max_intron_length\tmaximum possible intron length, trim 3' if the intron is longer than this value, do not consider intron if this value is < 0\n");
  printf("locuslink_only\tif this value is 1, then consider only entries with LocusID\n");
  printf("sortpos\tsort positions\n");
  printf("min_dist_factor\tmin_dist=min_dist_factor*motif_length  (if min_dist_factor<0, then set remove_overlap to 0 and do not call remove_vector_overlap)\n");
  printf("<singlefile>=if singlefile is 1, then there is a single file with all transcripts for all chromosomes\n");
  printf("<seq_source>=gbk, ensembl, etc.\n");
  printf("<read_all>=1 to read all lines including the ones with 0 length, read_all=0 to skip the zero length lines\n");
  printf("<max_n_transcripts> used only to allocate memory for the sequences [default=20000]\n");
  printf("verbose indicates where and how to report the log output\n");
  printf("<interesting_filename> file with motifs to scan\n");
  printf("\noutput files: %s.txt, %s.strand.txt and %s.log.txt if verbose=1\n",program_name,program_name,program_name);
  printf("\nlast modified=%s\n",last_modified);
  gk();
  exit(1);
}
