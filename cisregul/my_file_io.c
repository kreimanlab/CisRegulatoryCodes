/* my_file_io.c
 * collection of functions for writing and loading vectors and matrices
 * this is being written as the functions are needed
 *
 * last modified: 09-03-2003
 *
 */

int is_header(char *txtline);
void check_file_handle(FILE *file_handle,char *filename);
void gk(void);
void save_fmatrix(char *output_filename,float **data,int nx,int ny);
void save_fmatrixn(char *output_filename,float **data,int nx,int ny,int n);
void save_fvector(char *output_filename,float *data,int n);
void save_fvector6(char *output_filename,float *data,int n);
void save_fvector8(char *output_filename,float *data,int n);
void save_imatrix(char *output_filename,int **data,int nx,int ny);
void save_ivector(char *output_filename,int *data,int n);
void save_lvector(char *output_filename,long *data,int n);

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void gk(void)
{
  printf("Gabriel Kreiman\n");
  printf("kreiman@mit.edu\n");
  printf("All rights reserved (c)\n");
}

void check_file_handle(FILE *file_handle,char *filename)
{
  if (!file_handle) {
    printf("ERROR! I could not find the file %s",filename);
    printf("exiting...\n");
    exit(1);
  }
}

int is_header(char *txtline)
{
  /* given a line of text, return 1 if it starts with ">" or "%" and 0 otherwise */

  char searchchar;        // search whether it is a comment line
  char searchchar2;       // allow also matlab comments!
  char *gotit;            // used in the search process
  int header;           
  int linelength;         // length of the line (is_header=1 if line is empty)

  searchchar='>';
  searchchar2='\%';
  
  gotit=strchr(txtline,searchchar);
  if (gotit) {
    header=1;
  } else {
    gotit=strchr(txtline,searchchar2);
    if (gotit) {
      header=1;
    } else {
      header=0;
    }  
  }
  if (header==0) {
    linelength=strlen(txtline);
    if (linelength==0)
      header=1;
  }
   
  return header;
}

/* 
 * void splitline2int(char *str,int *v,int start,int *np)
 * this was now moved to my_text_package.c
 {
 
 char * pch;
 int i=1;
 int j=0;
 float f;
 
 pch = strtok (str,"\t");
 f=atof(pch);
 j++;
 v[j]=f;
 while (pch != NULL) {
 pch = strtok(NULL,"\t");
 i++;
 if (i>start) {      
 if (pch) {
 f=atof(pch);
 j++;
 v[j]=f;	
 }
 }
 }  

 np[0]=j;
 
 }
*/

void save_fvector(char *output_filename,float *data,int n)
{
  int i;
  FILE *output_file;
  
  output_file=fopen(output_filename,"w");
  if (!output_file)
    printf("I could not open the file %s for writing\n",output_filename);
  
  for (i = 1; i <= n; i++)	
    fprintf(output_file,"%.4f\n",data[i]);
  
  fclose(output_file);				
}

void save_fvector6(char *output_filename,float *data,int n)
{
  int i;
  FILE *output_file;
  
  output_file=fopen(output_filename,"w");
  if (!output_file)
    printf("I could not open the file %s for writing\n",output_filename);
  
  for (i = 1; i <= n; i++)	
    fprintf(output_file,"%.6f\n",data[i]);
  
  fclose(output_file);				
}

void save_fvector8(char *output_filename,float *data,int n)
{
  int i;
  FILE *output_file;
  
  output_file=fopen(output_filename,"w");
  if (!output_file)
    printf("I could not open the file %s for writing\n",output_filename);
  
  for (i = 1; i <= n; i++)	
    fprintf(output_file,"%.8f\n",data[i]);
  
  fclose(output_file);				
}

void save_ivector(char *output_filename,int *data,int n)
{
	int i;
	FILE *output_file;
	
	output_file=fopen(output_filename,"w");
	if (!output_file)
		printf("I could not open the file %s for writing\n",output_filename);
	
	for (i = 1; i <= n; i++)	
		fprintf(output_file,"%d\n",data[i]);
	
	fclose(output_file);				
}

void save_lvector(char *output_filename,long *data,int n)
{
  int i;
  FILE *output_file;
  
  output_file=fopen(output_filename,"w");
  if (!output_file) {
    printf("I could not open the file %s for writing\n",output_filename);
  } else {	
    for (i = 1; i <= n; i++)	
      fprintf(output_file,"%ld\n",data[i]);
    fclose(output_file);			
  }	
}

void save_fmatrix(char *output_filename,float **data,int nx,int ny)
{
	int i,j;
	FILE *output_file;

	output_file=fopen(output_filename,"w");
	if (!output_file)
		printf("I could not open the file %s for writing\n",output_filename);
	
	for (i = 1; i <= nx; i++)	{
		for (j=1;j<=ny;j++)	{
			fprintf(output_file,"%.4f\t",data[i][j]);
		}
		fprintf(output_file,"\n");
	}
	
	fclose(output_file);				
}

void save_fmatrixn(char *output_filename,float **data,int nx,int ny,int n)
{
  int i,j;
  FILE *output_file;

  output_file=fopen(output_filename,"w");
  if (!output_file)
    printf("I could not open the file %s for writing\n",output_filename);
  
  if (n==6) {
    for (i = 1; i <= nx; i++)	{
      for (j=1;j<=ny;j++)	{
	fprintf(output_file,"%.6f\t",data[i][j]);
      }
      fprintf(output_file,"\n");
    }
  } else {
    if (n==8) {
      for (i = 1; i <= nx; i++)	{
	for (j=1;j<=ny;j++)	{
	  fprintf(output_file,"%.8f\t",data[i][j]);
	}
	fprintf(output_file,"\n");
      }
    } else {
      /* default=4 decimal places */
      for (i = 1; i <= nx; i++)	{
	for (j=1;j<=ny;j++)	{
	  fprintf(output_file,"%.6f\t",data[i][j]);
	}
	fprintf(output_file,"\n");
      }
    }
  }

  fclose(output_file);				
}

void save_imatrix(char *output_filename,int **data,int nx,int ny)
{
	int i,j;
	FILE *output_file;

	output_file=fopen(output_filename,"w");
	if (!output_file)
		printf("I could not open the file %s for writing\n",output_filename);
	
	for (i = 1; i <= nx; i++)	{
		for (j=1;j<=ny;j++)	{
			fprintf(output_file,"%d\t",data[i][j]);
		}
		fprintf(output_file,"\n");
	}
	
	fclose(output_file);				
}
