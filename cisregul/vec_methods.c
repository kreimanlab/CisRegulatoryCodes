/* vec_methods */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "nr.c"
//#include "my_math_package.c"

// 02-28-2003 added max_n_i12 in int_intersect2 to avoid memory overflow

/* function definitions */
void int_intersect2(int *v1,int *v2,int n1,int n2,int **i12,int *f,long max_n_i12);
void intercal2matrix(int *intercalv,int n,int **matrix);
void intercal2vec(int *intercalv,int n,int *col1,int *col2);
int int_intersect(int *v1,int *v2,int n1,int n2);
int int_intersect1(int *v1,int *v2,int n1,int n2);

void int_intersect2(int *v1,int *v2,int n1,int n2,int **i12,int *f,long max_n_i12) 
{
  int i,j;
  int curr_v;
  int n=0;
  
  i=0;
  while ( (i<n1) & (n<max_n_i12) ) {
    //  for (i=1;i<=n1;i++) {
    i++;
    curr_v=v1[i];
    //for (j=1;j<=n2;j++) {
    j=0;
    while ( (j<n2) & (n<max_n_i12) ) {
      j++;
      if (v2[j]==curr_v) {
	n++;
	i12[n][1]=i;
	i12[n][2]=j;
	i12[n][3]=curr_v;
      }      
    }
  }
  if (i<n1)
    printf("WARNING!!!\tvec_methods.c\tint_intersect2\tn=%d\tmax_n_i12=%d\ti=%d\tn1=%d\tOVERFLOW\n",n,max_n_i12,i,n1);
  f[0]=n;
}

int int_intersect1(int *v1,int *v2,int n1,int n2)
{
  /* return the number of intersecting points between v1 and v2 
   * here v1 and v2 are assumed to be sorted in ascending order 
   * this attempts to be faster than int_intersect
   */
  
  int i,j;
  int init_j=1;         // initial value to start searching in vector 2
  int curr_v;
  int n=0;              // number of interesections
  int max2;             // maximum of v2

  max2=maxi(v2,n2);
  i=1;
  while (i<=n1) {
    curr_v=v1[i];
    //printf("i=%d\tcurr_v=%d\n",i,curr_v);
    if (curr_v>max2) {
      i=n1+1;
    } else {
      j=init_j;
      //printf("i=%d\tinit_j=%d\tn=%d\n",i,init_j,n);
      while ( (j<=n2) & (v2[j]<curr_v) ) 
	j++;
      if (v2[j]==curr_v)
	n++;
      init_j=j;
    }
    i++;
  }
  return n;
}

int int_intersect(int *v1,int *v2,int n1,int n2) 
{
  /* return the number of intersecting points between v1 and v2 */
  int i,j;
  int curr_v;
  int n=0;
  
  for (i=1;i<=n1;i++) {
    curr_v=v1[i];
    for (j=1;j<=n2;j++) {
      if (v2[j]==curr_v) 
	n++;
    }    
  }

  return n;
  
}

void intercal2matrix(int *intercalv,int n,int **matrix)
{
  int i;
  int j;
  int k=0;
  
  for (i=1;i<=n;i+=2) {
    j=i+1;
    k++;
    
    matrix[k][1]=intercalv[i];
    matrix[k][2]=intercalv[j];
    
  }  
}

void intercal2vec(int *intercalv,int n,int *col1,int *col2)
{
  int i;
  int j;
  int k=0;
  
  for (i=1;i<=n;i+=2) {
    j=i+1;
    k++;
    
    col1[k]=intercalv[i];
    col2[k]=intercalv[j];
  }
}





