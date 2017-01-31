/* compare_cooc_outputs_v2.c
 * given the outputs of cooc_checkresults with different parameters, try to find motif clusters that occur with several parameters 
 *
 * 11-17-2003: allow up to 12 files (previously the maximum was 6)
 * new in v2
 * attempt to solve a bug in call to sort_int_matrix so that we can process 4 motifs
 * 02-01-2004: attempt to solve a bug, file_indices can contain more than n_files elements (multiple occurrences in a file)
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
#include "../lib/my_file_io.c"
#include "../lib/my_text_package.c"
#include "crm_lib.c"
#include "cooc_single_lib.c"

int n_files;                     // number of files to compare 
FILE *log_file;                  // log file
int verbose=0;                   // 0: print to log file, 2: print to screen, 1: print to both
int prog_version=1;              // program version

void usage(void);
int load_file(char *filename,int **data,int filenumber,int max_n_motifs,int pos);
void sort_int_matrix(int **m,int n_rows,int n_cols_sort,int n_cols_write,int **sorted_matrix,int max_entry);
int comp2v(int *v1,int *v2,int vec_length);
int motif_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_list,int n_motifs);
int pair2index(int i1,int i2,int n_motifs);
int motif_pair_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_pairs_list,int n_motifs); 
void index2pair(int index,int n_motifs,int *i12);
void index2triplet(int index,int n_motifs,int *i123);
int triplet2index(int i1,int i2,int i3,int n_motifs);
int motif_triplet_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_triplets_list,int n_motifs);

main(int argc,char *argv[])
{
  int i,j,k,ll,mm;                   // 
  char infotext[10000];              // information for log file
  char temptext[10000];       
  int max_n_motifs;                  // maximum number of motifs (2, 3 or 4)
  char filename[200];
  int **d;
  int **sorted_d;
  int max_n_entries=100000;
  int n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12;
  int *skip;
  int pos=0;
  FILE *log_file;
  char filename1[200];               // scan file for motif 1
  char filename2[200];               // scan file for motif 2
  char filename3[200];               // scan file for motif 3
  char filename4[200];               // scan file for motif 4
  char filename5[200];
  char filename6[200];
  char filename7[200];               // scan file for motif 7
  char filename8[200];               // scan file for motif 8
  char filename9[200];               // scan file for motif 9
  char filename10[200];               // scan file for motif 10
  char filename11[200];
  char filename12[200];
  int n_cols;                        // number of columns in d (max_n_motifs+1)
  int max_entry=1000;                // maximum entry in d (used within the the sort_int_matrix program
  int *file_indices;                 // for each motif cluster, indices of files in which it occurs
  int r;
  int n_current;
  int n_multiple;
  int n_unique;
  int *v1;
  int *v2;
  FILE *unique_file;
  FILE *multiple_file;
  int sameornot;                     // 1 to indicate that two vectors are the same and 0 otherwise
  int n_motifs;                      // total number of motifs analyzed
  int *m1_list;                      // distribution of single motifs
  int nm1;                           // number of single motifs with at least one occurrence
  int *m2_list;                      // distribution of pairs of motifs
  int nm2;                           // number of motif pairs with at least one occurrence
  int *m3_list;                      // distribution of triplets of motifs
  int nm3;                           // number of motif triplets with at least one occurrence
  float *indices1;                   // used to sort m1_list
  float *indices2;                   // used to sort m2_list
  float *indices3;                   // used to sort m3_list
  int cmotif;                        // count of the number of occurrences of one particular motif
  float pmotif;                      // proportion of the number of occurrences of one particular motif
  int curr_index;
  int i12[3];                        // used in call to index2pair to get index1=i12[0] and index2=i12[1]
  int i123[4];                       // used in call to index2triplet to get index1=i123[0], index2=i123[1] and index3=i123[2]
  int i1,i2,i3;                      // index 1, index 2 and index 3
  long max_pairs;                    // maximum number of motif pairs
  long max_triplets;                 // maximum number of motif triplets
  int threshold_report_single=2;     // minimum number of occurrences to report a single motif
  int threshold_report_pair=2;       // minimum number of occurrences to report a pair of motifs
  int threshold_report_triplet=2;    // minimum number of occurrences to report a triplet of motifs
  struct timeb t1,t2;                // computation time
  time_t time_now;
  int curr_arg;

  if (argc<4) 
    usage();
  
  time(&time_now);                                     /* get time in seconds */
  curr_arg=1;n_files=atoi(argv[curr_arg]);             // total number of files
  curr_arg++;max_n_motifs=atoi(argv[2]);
  n_cols=max_n_motifs+1;
  curr_arg++;n_motifs=atoi(argv[3]);
  max_pairs=n_motifs*n_motifs;       // maximum number of motif pairs
  max_triplets=max_pairs*n_motifs;   // maximum number of motif triplets
  curr_arg++;sprintf(filename1,argv[4]);
  curr_arg++;
  if (n_files>1)
    sprintf(filename2,argv[curr_arg]);
  curr_arg++;
  if (n_files>2)
    sprintf(filename3,argv[curr_arg]);
  curr_arg++;
  if (n_files>3)
    sprintf(filename4,argv[curr_arg]);
  curr_arg++;
  if (n_files>4)
    sprintf(filename5,argv[curr_arg]);
  curr_arg++;
  if (n_files>5)
    sprintf(filename6,argv[curr_arg]);
  curr_arg++;
  if (n_files>6)
    sprintf(filename7,argv[curr_arg]);
  curr_arg++;
  if (n_files>7)
    sprintf(filename8,argv[curr_arg]);
  curr_arg++;
  if (n_files>8)
    sprintf(filename9,argv[curr_arg]);
  curr_arg++;
  if (n_files>9)
    sprintf(filename10,argv[curr_arg]);
  curr_arg++;
  if (n_files>10)
    sprintf(filename11,argv[curr_arg]);
  curr_arg++;
  if (n_files>11)
    sprintf(filename12,argv[curr_arg]);

  log_file=fopen("compare_cooc_outputs.log.txt","w");
  if (!log_file) {
    printf("error! i could not open compare_cooc_outputs.log.txt for writing. exiting...\n");
    printf("system_status=-1\n");
    exit(1);
  }

  /* write down list of parameters */
  sprintf(infotext,"%% compare_cooc_outputs_v2.c\n");
  sprintf(temptext,"%% %s",asctime(localtime(&time_now)));strcat(infotext,temptext);      
  sprintf(temptext,"%% n_files = %d\n",n_files);strcat(infotext,temptext);
  sprintf(temptext,"%% max_n_motifs = %d\n",max_n_motifs);strcat(infotext,temptext);
  sprintf(temptext,"%% n_motifs = %d\n",n_motifs);strcat(infotext,temptext);
  sprintf(temptext,"%% max_pairs = %ld\n",max_pairs);strcat(infotext,temptext);
  sprintf(temptext,"%% max_triplets = %ld\n",max_triplets);strcat(infotext,temptext);
  sprintf(temptext,"%% filename1 = %s\n",filename1);strcat(infotext,temptext);
  sprintf(temptext,"%% filename2 = %s\n",filename2);strcat(infotext,temptext);
  sprintf(temptext,"%% filename3 = %s\n",filename3);strcat(infotext,temptext);
  sprintf(temptext,"%% filename4 = %s\n",filename4);strcat(infotext,temptext);
  sprintf(temptext,"%% filename5 = %s\n",filename5);strcat(infotext,temptext);
  sprintf(temptext,"%% filename6 = %s\n",filename6);strcat(infotext,temptext);
  sprintf(temptext,"%% filename7 = %s\n",filename7);strcat(infotext,temptext);
  sprintf(temptext,"%% filename8 = %s\n",filename8);strcat(infotext,temptext);
  sprintf(temptext,"%% filename9 = %s\n",filename9);strcat(infotext,temptext);
  sprintf(temptext,"%% filename10 = %s\n",filename10);strcat(infotext,temptext);
  sprintf(temptext,"%% filename11 = %s\n",filename11);strcat(infotext,temptext);
  sprintf(temptext,"%% filename12 = %s\n",filename12);strcat(infotext,temptext);
  sprintf(temptext,"%% threshold report single = %d\n",threshold_report_single);strcat(infotext,temptext);
  sprintf(temptext,"%% threshold report pair = %d\n",threshold_report_pair);strcat(infotext,temptext);
  sprintf(temptext,"%% threshold report triplet = %d",threshold_report_triplet);strcat(infotext,temptext);
  fprintf(log_file,"%s\n",infotext);

  if (verbose)
    fprintf(log_file,"memory allocation...\n");
  d=imatrix(1,max_n_entries,1,n_cols);
  sorted_d=imatrix(1,max_n_entries,1,n_cols);
  skip=ivector(1,max_n_entries);
  for (i=1;i<=max_n_entries;i++)
    skip[i]=0;
  //file_indices=ivector(1,n_files);
  file_indices=ivector(1,1000);
  v1=ivector(1,max_n_motifs);
  v2=ivector(1,max_n_motifs);
  m1_list=ivector(1,n_motifs);                           // count of occurrences for each possible single motif [memory = 1.4 kb assuming n_motifs=350]
  m2_list=ivector(1,max_pairs);                          // count of occurrences for each possible motif pair   [memory = 490 kb assuming n_motifs=350]
  if (max_n_motifs>2)
    m3_list=ivector(1,max_triplets);                     // count of occurrences for each possible motif triplet [memory = 171.5 Mb assuming n_motifs=350]
  indices1=vector(1,n_motifs);                           // used to sort the values in m1_list [memory = 1.4 kb assuming n_motifs=350]
  indices2=vector(1,max_pairs);                          // used to sort the values in m2_list [memory = 490 kb assuming n_motifs=350]
  if (max_n_motifs>2)
    indices3=vector(1,max_triplets);                     // used to sort the values in m3_list [memory = 171.5 Mb assuming n_motifs=350]

  fprintf(log_file,"%% reading file=%s\t",filename1);
  n1=load_file(filename1,d,1,max_n_motifs,pos);
  fprintf(log_file,"n1=%d\n",n1);
  pos+=n1;

  if (n_files>1) {
    fprintf(log_file,"%% reading file=%s\t",filename2);
    n2=load_file(filename2,d,2,max_n_motifs,pos);
    fprintf(log_file,"n2=%d\n",n2);
    pos+=n2;
    
    if (n_files>2) {
      fprintf(log_file,"%% reading file=%s\t",filename3);
      n3=load_file(filename3,d,3,max_n_motifs,pos);
      fprintf(log_file,"n3=%d\n",n3);
      pos+=n3;
      
      if (n_files>3) {
	fprintf(log_file,"%% reading file=%s\t",filename4);
	n4=load_file(filename4,d,4,max_n_motifs,pos);
	fprintf(log_file,"n4=%d\n",n4);
	pos+=n4;
	
	if (n_files>4) {
	  fprintf(log_file,"%% reading file=%s\t",filename5);
	  n5=load_file(filename5,d,5,max_n_motifs,pos);
	  fprintf(log_file,"n5=%d\n",n5);
	  pos+=n5;
	  
	  if (n_files>5) {
	    fprintf(log_file,"%% reading file=%s\t",filename6);
	    n6=load_file(filename6,d,6,max_n_motifs,pos);
	    fprintf(log_file,"n6=%d\n",n6);
	    pos+=n6;

	    if (n_files>6) {
	      fprintf(log_file,"%% reading file=%s\t",filename7);
	      n7=load_file(filename7,d,7,max_n_motifs,pos);
	      fprintf(log_file,"n7=%d\n",n7);
	      pos+=n7;
	    
	      if (n_files>7) {
		fprintf(log_file,"%% reading file=%s\t",filename8);
		n8=load_file(filename8,d,8,max_n_motifs,pos);
		fprintf(log_file,"n8=%d\n",n8);
		pos+=n8;
		
		if (n_files>8) {
		  fprintf(log_file,"%% reading file=%s\t",filename9);
		  n9=load_file(filename9,d,9,max_n_motifs,pos);
		  fprintf(log_file,"n9=%d\n",n9);
		  pos+=n9;
		
		  if (n_files>9) {
		    fprintf(log_file,"%% reading file=%s\t",filename10);
		    n10=load_file(filename10,d,10,max_n_motifs,pos);
		    fprintf(log_file,"n10=%d\n",n10);
		    pos+=n10;
	      
		    if (n_files>10) {
		      fprintf(log_file,"%% reading file=%s\t",filename11);
		      n11=load_file(filename11,d,11,max_n_motifs,pos);
		      fprintf(log_file,"n11=%d\n",n11);
		      pos+=n11;
	  
		      if (n_files>11) {
			fprintf(log_file,"%% reading file=%s\t",filename12);
			n12=load_file(filename12,d,12,max_n_motifs,pos);
			fprintf(log_file,"n12=%d\n",n12);
			pos+=n12;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  fprintf(log_file,"%% total number of entries (pos) =%d\n",pos);

  if (pos>0) {
    if (verbose) {
      fprintf(log_file,"matrix\n");
      for (i=1;i<=pos;i++) {
	for (j=1;j<=n_cols;j++) 
	  fprintf(log_file,"%d\t",d[i][j]);
	fprintf(log_file,"\n");
      }
    }
    
    fprintf(log_file,"%% computing distribution of single motifs...\n");
    nm1=motif_distribution(d,pos,max_n_motifs,m1_list,n_motifs);
    printf("%% nm1=%d\n",nm1);
    fprintf(log_file,"%% number of single motifs (nm1) = %d\n",nm1);
    /* sort and print the results */
    fprintf(log_file,"%%\tsorting distribution of single motifs...\n");
    for (i=1;i<=n_motifs;i++)
      indices1[i]=(float)i;
    quicksorti(n_motifs,m1_list,indices1);
    cmotif=threshold_report_single+1;
    i=n_motifs+1;
    fprintf(log_file,"%%\tdistribution of single motifs...\n");
    fprintf(log_file,"%%index\toccurrences\tproportion\n");
    while (cmotif>0) {
      i--;
      curr_index=(int)indices1[i];
      cmotif=m1_list[i];
      //if (cmotif>0) {
      if (cmotif>threshold_report_single) {
	pmotif=(float)cmotif/(float)nm1;
	fprintf(log_file,"%d\t%d\t%.6f\n",curr_index,cmotif,pmotif);
      }
      if (i<=0) 
	cmotif=-1;
    }
    
    fprintf(log_file,"%%computing distribution of pairs of motifs...\n");
    nm2=motif_pair_distribution(d,pos,max_n_motifs,m2_list,n_motifs); 
    fprintf(log_file,"%% number of motif pairs (nm2) = %d\n",nm2);
    printf("%% nm2=%d\n",nm2);
    /* sort and print the results */
    fprintf(log_file,"%%\tsorting distribution of motif pairs...\n");
    for (i=1;i<=max_pairs;i++)
      indices2[i]=(float)i;
    quicksorti(max_pairs,m2_list,indices2);
    cmotif=threshold_report_pair+1;
    i=max_pairs+1;
    fprintf(log_file,"%%\tdistribution of motif pairs...\n");
    fprintf(log_file,"%%index\tmotif1\tmotif2\toccurrences\tproportion\n");
    while (cmotif>0) {
      i--;
      curr_index=(int)indices2[i];
      cmotif=m2_list[i];
      index2pair(curr_index,n_motifs,i12);
      i1=i12[0];
      i2=i12[1];
      if (cmotif>threshold_report_pair) {
	pmotif=(float)cmotif/(float)nm2;
	fprintf(log_file,"%d\t%d\t%d\t%d\t%.6f\n",curr_index,i1,i2,cmotif,pmotif);
	if (i<=0)
	  cmotif=-1;
      }
    }

    if (max_n_motifs>2) {
      fprintf(log_file,"%%computing distribution of triplets of motifs...\n");
      nm3=motif_triplet_distribution(d,pos,max_n_motifs,m3_list,n_motifs); 
      fprintf(log_file,"%% number of motif triplets (nm3) = %d\n",nm3);
      printf("%% nm3=%d\n",nm3);
      /* sort and print the results */
      fprintf(log_file,"%%\tsorting distribution of motif triplets...\n");
      for (i=1;i<=max_triplets;i++)
	indices3[i]=(float)i;
      quicksorti(max_triplets,m3_list,indices3);
      cmotif=threshold_report_triplet+1;
      i=max_triplets+1;
      fprintf(log_file,"%%\tdistribution of motif triplets...\n");
      fprintf(log_file,"%%index\tmotif1\tmotif2\tmotif3\toccurrences\tproportion\n");
      while (cmotif>threshold_report_triplet) {
	i--;
	curr_index=(int)indices3[i];
	cmotif=m3_list[i];
	index2triplet(curr_index,n_motifs,i123);
	i1=i123[0];
	i2=i123[1];
	i3=i123[2];
	if (cmotif>0) {
	  pmotif=(float)cmotif/(float)nm3;
	  fprintf(log_file,"%d\t%d\t%d\t%d\t%d\t%.6f\n",curr_index,i1,i2,i3,cmotif,pmotif);
	  if (i<=0)
	    cmotif=-1;
	}
      }
    }
    
    /* sort the matrix */
    //sort_int_matrix(d,pos,n_cols,n_cols,sorted_d,max_entry);
    sort_int_matrix(d,pos,max_n_motifs,n_cols,sorted_d,max_entry);

    unique_file=fopen("compare_cooc_outputs.uni.txt","w");
    if (!unique_file) {
      printf("error! i could not open compare_cooc_outputs.uni.txt for writing. exiting...\n");
      printf("system_status=-1\n");
      exit(1);
    }
    fprintf(unique_file,"%s\n",infotext);
    fprintf(unique_file,"%% motif1\tmotif2\tmotif3\tn_files\tfiles\n");
    multiple_file=fopen("compare_cooc_outputs.mul.txt","w");
    if (!unique_file) {
      printf("error! i could not open compare_cooc_outputs.mul.txt for writing. exiting...\n");
      printf("system_status=-1\n");
      exit(1);
    }
    fprintf(multiple_file,"%s\n",infotext);
    fprintf(multiple_file,"%% motif1\tmotif2\tmotif3\tn_files\tfiles\n");
    r=1;
    n_current=1;                       // number of files with current motif cluster
    n_unique=0;                        // number of motif clusters that occur in only one file
    n_multiple=0;                      // number of motif clusters that occur in >1 file
    for (i=1;i<=max_n_motifs;i++) 
      v1[i]=0;
    while (r<pos) {
      for (i=1;i<=max_n_motifs;i++)
	v2[i]=sorted_d[r][i];
      sameornot=comp2v(v1,v2,max_n_motifs);
      if (sameornot==1) {
	/* the two vectors are the same */
	n_current++;
	file_indices[n_current]=sorted_d[r][n_cols];
      } else {
	if (n_current>1) {
	  /* print as repeated */
	  n_multiple++;
	  for (i=1;i<=max_n_motifs;i++) 
	    fprintf(multiple_file,"%d\t",v1[i]);
	  fprintf(multiple_file,"%d\t",n_current);
	  for (i=1;i<=n_current;i++)
	    fprintf(multiple_file,"%d\t",file_indices[i]);
	  fprintf(multiple_file,"\n");
	} else {
	  if (r>1) {
	    /* print as unique */
	    n_unique++;
	    for (i=1;i<=max_n_motifs;i++) 
	      fprintf(unique_file,"%d\t",v1[i]);
	    fprintf(unique_file,"1\t%d\n",sorted_d[r-1][n_cols]);
	  }
	}
	for (i=1;i<=max_n_motifs;i++)
	  v1[i]=v2[i];
	n_current=1;
	file_indices[n_current]=sorted_d[r][n_cols];
      }
      r++;
      //printf("%d\n",r);
    }
    fclose(unique_file);
    fclose(multiple_file);
  }

  fprintf(log_file,"%% number of motif clusters analyzed = %d\n",pos);
  fprintf(log_file,"%% n_multiple = %d\n",n_multiple);
  fprintf(log_file,"%% n_unique = %d\n",n_unique);

  fprintf(log_file,"%% freeing memory...\n");
  free_imatrix(d,1,max_n_entries,1,n_cols);
  free_imatrix(sorted_d,1,max_n_entries,1,n_cols);
  free_ivector(skip,1,max_n_entries);
  //free_ivector(file_indices,1,n_files);
  free_ivector(file_indices,1,1000);
  free_ivector(v1,1,max_n_motifs);
  free_ivector(v2,1,max_n_motifs);
  free_ivector(m1_list,1,n_motifs);
  free_ivector(m2_list,1,max_pairs);
  if (max_n_motifs>2)
    free_ivector(m3_list,1,max_triplets);
  free_vector(indices1,1,n_motifs);
  free_vector(indices2,1,max_pairs);
  if (max_n_motifs>2)
    free_vector(indices3,1,max_triplets);
  fclose(log_file);

  printf("system_status=0\n");
  printf("%%ready. thanks for flying the friendly skies!\n");
  
  return 0;
}

void usage(void)
{
  printf("compare_cooc_outputs_v2.c\n");
  printf("compare output files from cooc_checkresults\n");
  printf("usage:\n");
  printf("compare_cooc_outputs_v2.exe <n_files> <max_n_motifs> <n_motifs> <filename_1> <filename_2> ...  <filename_n>\n");
  printf("where n<=12\n");

  gk();
  exit(1);
}

int load_file(char *filename,int **data,int filenumber,int max_n_motifs,int start_pos)
{
  FILE *input_file;
  int n=0;
  char txtline[100000];    // read one line of the input file
  int n_chars=1000000;    // number of chars to read (if more than this per line, there is an overflow, here we try to detect it)
  int header;
  int *v;
  int start=0;
  int np[1];
  int max_occur=100;
  float *dummy;
  int filepos;
  int pos;
  int i;
  
  v=ivector(1,max_occur);
  dummy=vector(1,max_n_motifs);
  filepos=max_n_motifs+1;
  pos=start_pos;
  
  input_file=fopen(filename,"r");
  if (!input_file) {
    printf("i could not open %s for reading, exiting...\n",filename);
    n=0;
    return n;
  }
  while (feof(input_file)==0) {
    fgets(txtline,n_chars,input_file);
    header=is_header(txtline);
    if (header!=1) {
      splitline2int(txtline,v,start,np,max_occur);
      quicksorti(max_n_motifs,v,dummy);
      pos++;      
      for (i=1;i<=max_n_motifs;i++)
	data[pos][i]=v[i];
      data[pos][filepos]=filenumber;
    }
  }
  
  n=pos-start_pos;
  return n;
}


void sort_int_matrix(int **m,int n_rows,int n_cols_sort,int n_cols_write,int **sorted_matrix,int max_entry) 
{
  int i,j;
  int **temp_matrix1;
  int **temp_matrix2;
  int *indices;
  int curr_index;
  double g1,g2,g3,g4,g;
  double *vd;

  vd=dvector(1,n_rows);
  indices=ivector(1,n_rows);
  for (i=1;i<=n_rows;i++)
    indices[i]=i;

  if (n_cols_sort==1) {
    for (i=1;i<=n_rows;i++)      
      vd[i]=(double)m[i][1]*(double)max_entry;
  } else {
    if (n_cols_sort==2) {
      for (i=1;i<=n_rows;i++) {
	g1=(double)m[i][1]*(double)max_entry;
	g2=(double)m[i][2];
	g=g1+g2;
	vd[i]=g;
      }
    } else {
      if (n_cols_sort==3) {
	for (i=1;i<=n_rows;i++) {
	  //f1=(float)m[i][1]*(float)max_entry;
	  //f2=(float)m[i][2];
	  //f3=(float)m[i][3]/(float)max_entry;
	  //v[i]=f1+f2+f3;
	  //printf("f1=%.10f\nf2=%.10f\nf3=%.10f\n v=%.10f\n   %.10f\n",f1,f2,f3,v[i],f1+f2+f3);
	  g1=(double)m[i][1]*(double)max_entry;
	  g2=(double)m[i][2];
	  g3=(double)m[i][3]/(double)max_entry;
	  //printf("g1=%.10g\ng2=%.10g\ng3=%.10g\n  =%.10g\n",g1,g2,g3,g1+g2+g3);
	  g=g1+g2+g3;
	  vd[i]=(float)g;
	  //printf("g1=%.10g\ng2=%.10g\ng3=%.10g\n v=%.10g\n",g1,g2,g3,vd[i]);
	  //exit(1);
	  //v[i]=(float)m[i][1]+(float)m[i][2]/(float)max_entry+(float)m[i][3]/(float)(max_entry*max_entry);
	  //v[i]=m[i][1]*max_entry*max_entry+m[i][2]*max_entry+m[i][3];
	  //printf("i=%d\tv=%.10f\tm[1][1]=%d --> %.10f\tm[1][2]=%d --> %.10f \tm[1][3]=%d --> %.10f\n",i,v[i],m[i][1],(float)m[i][1],m[i][2],(float)m[i][2]/(float)max_entry,m[i][3],(float)m[i][3]/(float)(max_entry*max_entry));
	}
      } else {
	if (n_cols_sort==4) {
	  for (i=1;i<=n_rows;i++) {
	    g1=(double)m[i][1]*(double)max_entry;
	    g2=(double)m[i][2];
	    g3=(double)m[i][3]/(double)max_entry;
	    g4=(double)m[i][4]/(double)(max_entry*max_entry);
	    g=g1+g2+g3+g4;
	    vd[i]=g;
	  }
	} else {
	  printf("you shouldn't be here ...\n");
	}
      }
    }
  }
  
  quicksortd(n_rows,vd,indices);
  /* debug here 
     printf("\nv\n");
     for (i=1;i<=n_rows;i++)
     printf("%d %.10g\n",i,vd[i]);
     printf("call quicksortd\n");
     //quicksortl(n_rows,v,indices);
     //quicksort(n_rows,v,indices);
     printf("\nsorted v\n");
     for (i=1;i<=n_rows;i++)
     printf("%d %.10g\n",i,vd[i]);
     printf("\nindices\n");
     for (i=1;i<=n_rows;i++)
     printf("%d %d\n",i,indices[i]);
  */

  /* copy matrix to temp_matrix1 */
  for (i=1;i<=n_rows;i++) {
    curr_index=indices[i];
    for (j=1;j<=n_cols_write;j++)
      sorted_matrix[i][j]=m[curr_index][j];
  }

  //free_lvector(v,1,n_rows);
  //free_vector(v,1,n_rows);
  free_dvector(vd,1,n_rows);
  free_ivector(indices,1,n_rows);
}

int comp2v(int *v1,int *v2,int vec_length)
{
  /* compare two vectors and return 1 if they are the same and 0 otherwise */
  
  int i;
  int sameornot=1;           // they are the same as default

  i=0;
  while (i<vec_length) {
    i++;
    if (v1[i]!=v2[i]) {
      sameornot=0;
      i=vec_length+1;
    }
  }

  return sameornot;
}

int motif_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_list,int n_motifs)
{
  /* given the matrix containing the motif clusters and which file they occur in, count the number of occurrences of each motif 
   * format for matrix
   * motif1 motif2 ... motifn file_index (n_rows x max_n_motifs+1)
   */
  int i,j;
  int motif;
  int n=0;

  for (i=1;i<=n_motifs;i++)
    motifs_list[i]=0;
  /* motifs_list=ivector(1,n_motifs);
     for (i=1;i<=n_motifs;i++)
     motifs_list++;
  */

  /* convert the matrix to a list of motifs */
  for (i=1;i<=n_rows;i++) {
    for (j=1;j<=max_n_motifs;j++) {
      motif=matrix[i][j];
      if (motif>0) {
	motifs_list[motif]=motifs_list[motif]+1;
	n++;
      }
    }
  }

  return n;
}

int motif_pair_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_pairs_list,int n_motifs) 
{
  int index;
  int i,j,k;
  int i1,i2;
  int n=0;
  int max_index;

  max_index=n_motifs*n_motifs;           // maximum possible index (this is the number of rows in motifs_pairs_list)
  for (i=1;i<=max_index;i++) 
    motifs_pairs_list[i]=0;

  /* for (i=1;i<=n_motifs;i++) {
     for (j=1;j<=n_motifs;j++) {
     motifs_pairs_list[i][1]=i;
     motifs_pairs_list[i][2]=j;
     motifs_pairs_list[i][3]=0;
     }
     }
  */

  for (i=1;i<=n_rows;i++) {
    for (j=1;j<max_n_motifs;j++) {
      i1=matrix[i][j];
      if (i1>0) {
	for (k=(j+1);k<=max_n_motifs;k++) {
	  i2=matrix[i][k];
	  if (i2>0) {
	    index=pair2index(i1,i2,n_motifs);
	    //motifs_pairs_list[index][3]=motifs_pairs_list[index][3]+1;
	    motifs_pairs_list[index]=motifs_pairs_list[index]+1;
	    n++;
	  }
	}
      }
    }
  }

  return n;
}

int pair2index(int i1,int i2,int n_motifs)
{
  /* e.g. 1 1            -> 1
   *      1 2            -> 2
   *      ... 
   *      1 n_motifs     -> n_motifs
   *      2 1            -> (1*n_motifs)+1
   *      ...
   *      i j            -> (i-1)*n_motifs+j;
   */

  int index;
  index=(i1-1)*n_motifs+i2;

  return index;
}

void index2pair(int index,int n_motifs,int *i12)
{
  int i1,i2;

  i2=index%n_motifs;
  i1=(index-i2)/n_motifs+1;

  i12[0]=i1;
  i12[1]=i2;
}

int motif_triplet_distribution(int **matrix,int n_rows,int max_n_motifs,int *motifs_triplets_list,int n_motifs) 
{
  int index;
  int i,j,k,ll;
  int i1,i2,i3;
  int n=0;
  int max_index;

  if (max_n_motifs<3) {
    printf("sorry. i cannot compute a distribution of triplets with max_n_motifs=%d\n",max_n_motifs);
    return n;
  }

  /* initialize motifs_triplets_list */
  max_index=n_motifs*n_motifs*n_motifs;           // maximum possible index (this is the number of rows in motifs_triplets_list)
  for (i=1;i<=max_index;i++) 
    motifs_triplets_list[i]=0;

  for (i=1;i<=n_rows;i++) {

    for (j=1;j<max_n_motifs;j++) {
      i1=matrix[i][j];
      if (i1>0) {
	for (k=(j+1);k<=max_n_motifs;k++) {
	  i2=matrix[i][k];
	  if (i2>0) {
	    for (ll=(k+1);ll<=max_n_motifs;ll++) {
	      i3=matrix[i][ll];
	      if (i3>0) {
		index=triplet2index(i1,i2,i3,n_motifs);
		motifs_triplets_list[index]=motifs_triplets_list[index]+1;
		n++;
	      }
	    }
	  }
	}
      }
    }
  }

  return n;
}

int triplet2index(int i1,int i2,int i3,int n_motifs)
{
  /* e.g. 1 1 1           -> 1
   *      1 1 2           -> 2
   *      ... 
   *      1 1 n_motifs    -> n_motifs
   *      1 2 1           -> 0*n_motifs*n_motifs+1*n_motifs+1
   *      ...
   *      1 n_motifs 1    -> 0*n_motifs*n_motifs+(n_motifs-1)*n_motifs+1
   *      2 1 1           -> 1*n_motifs*n_motifs+(0*n_motifs)+1
   *      ...
   *      i j k           -> (i-1)*n_motifs*n_motifs+(j-1)*n_motifs+k;
   */

  int index;

  index=(i1-1)*n_motifs*n_motifs+(i2-1)*n_motifs+i3;

  return index;
}

void index2triplet(int index,int n_motifs,int *i123)
{
  int i1,i2,i3;
  int temp;

  i3=index%n_motifs;
  temp=(index-i3)/n_motifs;
  i2=temp%n_motifs+1;
  i1=(temp+1-i2)/n_motifs+1;

  i123[0]=i1;
  i123[1]=i2;
  i123[2]=i3;
}
