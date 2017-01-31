/****************************************************************/
/*		Numerical Recipes Functions.			*/
/****************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define NR_END 1
#define FREE_ARG char*
#define ITMAX 100             /* maximum allowed number of iterations (used in gser,gcf) */
#define EPS 3.0e-7            /* relative accuracy (user in gser,gcf) */
#define FPMIN 1.0e-30         /* number near the smallest representable floating-point number (used in gcf) */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
char ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
char **cmatrix(long nrl,long nrh,long ncl,long nch);
float **matrix(long nrl,long nrh,long ncl,long nch);
float bico(int n, int k);
float erff(float x);
float erffc(float x);
float factln(int n);
float factrl(int n);
float gammln(float xx);
float gammp(float a,float x);
float gammq(float a,float x);
float select_klargest(unsigned long k, unsigned long n, float arr[]);
int  **imatrix(long nrl,long nrh,long ncl,long nch);
int  *ivector(long nl,long nh);
long **lmatrix(long nrl,long nrh,long ncl,long nch);
long *lvector(long nl,long nh);
unsigned char **ucmatrix(long nrl,long nrh,long ncl,long nch);
unsigned short  **usmatrix(long nrl,long nrh,long ncl,long nch);
unsigned short *usvector(long nl,long nh);
void free_cmatrix(char **m,long nrl,long nrh, long ncl, long nch);
void free_imatrix(int **m,long nrl,long nrh, long ncl, long nch);
void free_ivector(int *v, long nl, long nh);
void free_lmatrix(long **m,long nrl,long nrh, long ncl, long nch);
void free_lvector(long *v, long nl, long nh);
void free_matrix(float **m,long nrl,long nrh, long ncl, long nch);
void free_ucmatrix(unsigned char **m,long nrl,long nrh,long ncl,long nch);
void free_usmatrix(unsigned short **m,long nrl,long nrh, long ncl, long nch);
void free_usvector(unsigned short *v, long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
void free_c3tensor(char ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
void gcf(float *gammcf,float a,float x,float *gln);
void gser(float *gamser,float a,float x,float *gln);
void indexx(int n, float arr[], int indx[]);
double *dvector(long nl,long nh);
void free_dvector(double *v, long nl, long nh);

void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error ...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"... now exiting to system ...\n");
	exit(1);
}

float *vector(long nl,long nh)
{
	float *v;

        v=(float *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

double *dvector(long nl,long nh)
{
  double *v;

  v=(double *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

char *cvector(long nl,long nh)
{
        char *v;

        v=(char *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(char)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

float **matrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));

	if (!m) nrerror("allocation failure 1 in matrix()");

	m +=NR_END;
	m -=nrl;

	m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

void free_vector(float *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_vector; nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_vector; nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(char *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_cvector; nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_matrix; nrh must be >=0");
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

long *lvector(long nl,long nh)
{
	long *v;

	v=(long *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

int *ivector(long nl,long nh)
{
	int *v;

	v=(int *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned short *usvector(long nl,long nh)
{
  unsigned short *v;

  v=(unsigned short *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(unsigned short)));
  if (!v) 
    nrerror("allocation failure in usvector()");
  return v-nl+NR_END;
}

long **lmatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long **m;

	m=(long **)malloc((size_t)((nrow+NR_END)*sizeof(long*)));

	if (!m) nrerror("allocation failure 1 in matrix()");

	m +=NR_END;
	m -=nrl;

	m[nrl]=(long *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

int **imatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));

	if (!m) nrerror("allocation failure 1 in matrix()");

	m +=NR_END;
	m -=nrl;

	m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

unsigned short **usmatrix(long nrl,long nrh,long ncl,long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  unsigned short **m;

  m=(unsigned short **)malloc((size_t)((nrow+NR_END)*sizeof(unsigned short*)));

  if (!m) nrerror("allocation failure 1 in matrix()");

  m +=NR_END;
  m -=nrl;
  
  m[nrl]=(unsigned short *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned short)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

char **cmatrix(long nrl,long nrh,long ncl,long nch)
{
	long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **m;

        m=(char **)malloc((size_t)((nrow+NR_END)*sizeof(char*)));

	if (!m) nrerror("allocation failure 1 in matrix()");

	m +=NR_END;
	m -=nrl;

        m[nrl]=(char *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

unsigned char **ucmatrix(long nrl,long nrh,long ncl,long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  unsigned char **m;
  
  m=(unsigned char **)malloc((size_t)((nrow+NR_END)*sizeof(unsigned char*)));
  
  if (!m) nrerror("allocation failure 1 in matrix()");
  
  m +=NR_END;
  m -=nrl;
  
  m[nrl]=(unsigned char *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned char)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  return m;
}

void free_lvector(long *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_lvector: nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_ivector: nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_usvector(unsigned short *v, long nl, long nh)
{
  if (nh<0)
    nrerror("free_usvector: nh must be >=0");
  free((FREE_ARG) (v+nl-NR_END));
}

void free_lmatrix(long **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_lmatrix: nch must be >=0");
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_imatrix: nch must be >=0");
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_usmatrix(unsigned short **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_usmatrix: nch must be >=0");
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_cmatrix(char **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_cmatrix: nch must be >=0");

  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_ucmatrix(unsigned char **m,long nrl,long nrh, long ncl, long nch)
{
  if ( (nch<0) | (nrh<0) )
    nrerror("free_ucmatrix: nch must be >=0");

  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

float gammln(float xx)
{
  /* return the log of the gamma function for x>o */ 
  /* note: internal arithmetic will be done in double precision, 
   * a nicety that you can omit if five-figure accuracy is good enough
   */

  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) 
    ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

float factrl(int n)
{
  /* returns the value n! as a floating-point number */
  static int ntop=4;
  static float a[33]={1.0,1.0,2.0,6.0,24.0}; /* Fill in table only as required. */
  int j;
  if (n < 0) 
    nrerror("Negative factorial in routine factrl");
  if (n > 32) 
    return exp(gammln(n+1.0));
  /* Larger value than size of table is required. Actually, this big a value is going to over ow
   * on many computers, but no harm in trying.
   */

  while (ntop<n) { 
    /* Fill in table up to desired value. */
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}

float factln(int n)
{
  /* returns ln(n!) */
  //static float a[101]; /* A static array is automatically initialized to zero. */
  static float a[1001];
  if (n < 0) 
    nrerror("Negative factorial in routine factln");
  if (n <= 1) 
    return 0.0;
  //  if (n <= 100) {
  if (n<=1000) {
    /* in range of table */
    return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  } 
  else {
    /* out of range of table */
    return gammln(n+1.0); 
  }
}

double dbico(int n, int k)
{
  float a,b,c,d;
  long double e;
  a=factln(n);
  b=factln(k);
  c=factln(n-k);
  d=a-b-c;
  e=0.5+exp(d);
  e=floor(e);
  //return (floor(0.5+exp(a-b-c)));
  return e;
  /* returns the binomial coefficient nchoosek as a floating number */
  //return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
  /* the floor function cleans up roundoff error for smaller values of n and k. */
}

float bico(int n, int k)
{
  /* returns the binomial coefficient nchoosek as a floating number */
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
  /* the floor function cleans up roundoff error for smaller values of n and k. */
}

float gammp(float a,float x)
{
  /* returns the incomplete gamma function P(a,x) */
  /* uses functions gcf, gser */

  float gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) 
    nrerror("Invalid arguments in routine gammp");
  if (x < (a+1.0)) { 
    /* Use the series representation. */
    gser(&gamser,a,x,&gln);
    return gamser;
  } else { 
    /* Use the continued fraction representation */
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf; /* and take its complement. */
  }
}

float gammq(float a,float x)
{
  /* returns the complement of the incomplete gamma function Q(a,x)=1-P(a,x) */
  /* uses functions gcf, gser */

  float gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) 
    nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) { 
    /* Use the series representation */
    gser(&gamser,a,x,&gln);
    return 1.0-gamser; /* and take its complement. */
  } else { 
    /* Use the continued fraction representation.*/
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}

void gser(float *gamser,float a,float x,float *gln)
{
  /* returns the incomplete gamma function P(a,x) evaluated by its series representation as gamser,
   * also returns the ln [gamma(a)] as gln
   */
  int n;
  float sum,del,ap;
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) 
      nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}

void gcf(float *gammcf,float a,float x,float *gln)
{
  /* returns the incomplete gamma function Q(a,x) evaluated by its continued fraction representation as gammcf,
   * also returns ln [gamma(a)] as gln 
   */

  int i;
  float an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a;   /* Set up for evaluating continued fraction by modified Lentz's method (see NR 5.2), b0=0 */
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) { 
    /* Iterate to convergence. */
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) 
      d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) 
      c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) 
      break;
  }
  if (i > ITMAX) 
    nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h; /* Put factors in front */
}

float erff(float x)
{
  /* returns the error function erf(x) */
  return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

float erffc(float x)
{
  /* returns the complementary error function erfc(x) */
  return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}

float select_klargest(unsigned long k, unsigned long n, float arr[])
{ 
  /* Returns the kth smallest value in the array arr[1..n]. The input array will be rearranged
   * to have this value in location arr[k], with all smaller elements moved to arr[1..k-1] (in
   * arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).
   */

  unsigned long i,ir,j,l,mid;
  float a,temp;
  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {    // Active partition contains 1 or 2 elements.
      if (ir == l+1 && arr[ir] < arr[l]) {    // Case of 2 elements.
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;           // Choose median of left, center, and right el-
      SWAP(arr[mid],arr[l+1]);   // ements as partitioning element a. Also
      if (arr[l] > arr[ir]) {  //  rearrange so that arr[l] >= arr[l+1],
	SWAP(arr[l],arr[ir]);  //  arr[ir]>= arr[l+1].
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1;    // Initialize pointers for partitioning.
      j=ir;
      a=arr[l+1];   //  Partitioning element.
      for (;;) {    // Beginning of innermost loop.
	do i++; 
	while (arr[i] < a);    // Scan up to find element > a.
	do j--; 
	while (arr[j] > a);    // Scan down to find element < a.
	if (j < i) 
	  break;    // Pointers crossed. Partitioning complete.
	SWAP(arr[i],arr[j]);
      }    // End of innermost loop.
      arr[l+1]=arr[j];    // Insert partitioning element.
      arr[j]=a;           // 
      if (j >= k) 
	ir=j-1;    // Keep active the partition that contains the
      if (j <= k) 
	l=i;       // kth element.
    }
  }
}

float median(unsigned long n,float arr[]) {
  /* compute the median of values in arr 
   * note this is my own code calling select_klarges from nr
   * WARNING! The array arr is modified upon calling this function!
   */

  long k;
  float m; /* median value */
  float m1,m2;

  if ((n%2)==0) {
    /* n is even */
    k=n/2;
    select_klargest(k,n,arr);
    m1=arr[k];
    k++;
    select_klargest(k,n,arr);
    m2=arr[k];
    m=(m1+m2)/2;
  } else {
    /* n is odd */
    k=(n+1)/2;
    select_klargest(k,n,arr);
    m=arr[k];
  }

  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  /* return pointer to array of pointers to rows */
  return t;
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
     /* free a float f3tensor allocated by f3tensor() */
{
  if ( (nch<0) | (nrh<0) | (ndh<0) )
    nrerror("free_cmatrix: nch must be >=0");

  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}

char ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a char 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  char ***t;
  /* allocate pointers to pointers to rows */
  t=(char ***) malloc((size_t)((nrow+NR_END)*sizeof(char**)));
  if (!t) nrerror("allocation failure 1 in c3tensor()");
  t += NR_END;
  t -= nrl;
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(char **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char*)));
  if (!t[nrl]) nrerror("allocation failure 2 in c3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(char *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(char)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in c3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  /* return pointer to array of pointers to rows */
  return t;
}

void free_c3tensor(char ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
     /* free a char c3tensor allocated by c3tensor() */
{
  if ( (nch<0) | (nrh<0) | (ndh<0) )
    nrerror("free_cmatrix: nch must be >=0");

  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
