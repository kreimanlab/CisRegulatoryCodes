/* my_text_package.c 
 * functions used to process text
 */

void splitline2int(char *str,int *v,int start,int *np,int max_j);
void substring(char *string_orig,char *string_new,int n,int from);
void printinfo(char *text,int verbose,FILE *file_handle);

void printinfo(char *text,int verbose,FILE *file_handle)
{
  /* print out information in text variable either to screen and/or to a file handle */

  if (verbose<2)
    printf("%s\n",text);
  if (verbose>0)
    fprintf(file_handle,"%s\n",text);
}

void splitline2int(char *str,int *v,int start,int *np,int max_j)
{
  /* split a tab separated char line into integers starting at position start 
   * maximum of max_j entries allowed
   * np[0] returns the total number of entries
   * np[1] returns -1 if the actual number exceeds max_j 
   */
  char * pch;
  int i=1;
  int j=0;
  float f;
  
  //printf ("Splitting string \"%s\" in tokens:\n",str);
  pch = strtok (str,"\t");
  while (pch != NULL) {
    if (i>start) {      
      if (pch) {
	f=atof(pch);
	j++;
	/* debug here 
	printf("%d %.2f\t",j,f); */
	if (j<=max_j)
	  v[j]=f;	
      }
    }
    pch = strtok(NULL,"\t");
    i++;
  }  
  //printf("\n");

  np[0]=j;
  if (j>max_j) {
    np[1]=-1;
  } else {
    np[1]=0;
  }
}

void substring(char *string_orig,char *string_new,int n,int from)
{
  /* get substring (string_new) from a string (string_orig) by taking n chars either starting from the left if from > 0 or starting from the right if from < 0 */

  int string_orig_length;
  int strncpy_start=0;

  string_orig_length=strlen(string_orig);
  if (from < 0) 
    strncpy_start=string_orig_length-n;

  strncpy(string_new,string_orig+strncpy_start,n);
}
