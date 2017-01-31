#include <time.h>
#include <math.h>

/* my_math_package.c
 * elmentary math functions
 */

/* these functions assume that vectors were created according to
 * nr.c standards and start at 1
 */

/* 02-19-2003: changed binopdf_mat and binocdf_mat to require and return double variables, this improves the precision of the algorithm significantly 
 * 09-03-2003: moved gk(void) to my_file_io.c 
 * 10-16-2003: added vector_overlap
 * 01-02-2004: added sort_ivector_by_indices
 * 01-02-2004: added vector_overlap_v2
 * 01-07-2004: added basic_stats_int
 */

void basic_stats_int(int n,int *v,float *out);
char char3min(char x,char y,char z);
double binocdf_mat(int n,int k,double p);
double binopdf_mat(int n,int k,double p);
double gammaln_mat(double x);
float binocdf(int n,int k,float p);
float binopdf(int n,int k,float p);
float compute_percentile_f(float x,float *xv,int n);
float compute_percentile_i(int x,int *xv,int n);
float erf1(float x);
float erfc1(float x);
float erfcore(float x,int jint);
float erfinv(float y);
float fpower(float x,int n);
float maxf(float *x,int n);
float meanf(float *x,int n);
float meani(int *x,int n);
float minf(float *x,int n);
float normf(float *v,int n);
float normi(int *v,int n);
float pe_gauss(float mu0,float sigma0,float mu1,float sigma1,float res);
float range3f(float v1,float v2,float v3);
float rangef(float *x,int n) ;
float spearman(float *data1,float *data2,unsigned long n);
float stdf(float *x,float mean_x,int n);
float stdi(int *x,float mean_x,int n);
int fast_inormalize_vector(int *v,float *nv,float mean_v,float std_v,float sqrtnm1,int n);
int findi(int *vec,int n_vec,int search_value,int *index);
int fnormalize_vector(float *x,int n);
int inormalize_vector(int *x,float *nx,int n);
int max2(int x1,int x2);
int max3(int x1,int x2,int x3);
int max4(int x1,int x2,int x3,int x4);
int maxi(int *x,int n);
int min2(int x1,int x2);
int min3(int x1,int x2,int x3);
int min4(int x1,int x2,int x3,int x4);
int mini(int *x,int n);
int myround(float f);
int vector_overlap(int *original_vector,int n,int max_overlap,int *output_vector);
long maxl(long *x,int n);
long minl(long *x,int n);
long time_difference(struct timeb ti,struct timeb tf);
void crank(unsigned long n,float w[]);
void just_quicksortus(unsigned long n,unsigned short *arr);
void linetrend(float *x,float *y,int n,float *res);
void myrank(long n,float *w);
void quicksort(unsigned long n,float *arr,int *brr);
void quicksortf(unsigned long n,float *arr,float *brr);
void quicksorti(unsigned long n,int *arr,float *brr);
void quicksortl(unsigned long n,long *arr,float *brr);
float ran1(long *idum);
void randperms(int n,int *index,long idum);
void shell_sort(long n,float *arr,int *brr);
void sort_cvector_by_indices(int n,char *v,float *indices);
void sort_ivector_by_indices(int n,int *v,float *indices);
void sort_ivector_by_int_indices(int n,int *v,int *indices);

#define MYONE 1.0                 // used only in other definitions
#define IM 2147483647             // ran1
#define IM1 2147483563            // ran2
#define AM (1.0/IM1)              // ran1,ran2
// #define EPS 1.2e-7                // ran2      already defined in nr.c
#define IA 16807                  // ran1
#define IA1 40014                 // ran2
#define IA2 40692                 // ran2
#define IM2 2147483399            // ran2
#define IMM1 (IM1-1)              // ran2
#define IQ 127773                 // ran1
#define IQ1 53668                 // ran2
#define IQ2 52774                 // ran2
#define IR 2836                   // ran1
#define IR1 12211                 // ran2 
#define IR2 3791                  // ran2
#define M 7                       // quicksort
#define NDIV1 (1+(IM-1)/NTAB)     // ran1
#define NDIV (1+IMM1/NTAB)        // ran2
#define NSTACK 200                // quicksort
#define NTAB 32                   // ran1,ran2
#define PI 3.1415926              // erfcore
#define RNMX (1.0-EPS)            // ran2
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define sqrt2 sqrt(2)

void basic_stats_int(int n,int *v,float *out)
{
  /* return mean, s.d., min, max for a vector of integer values 
   * usage:
   * basic_stats_int(n,v,out);
   * where out[0]=mean, out[1]=s.d., out[2]=min, out[3]=max
   */

  float mean_v;

  mean_v=meani(v,n);
  out[0]=mean_v;
  out[1]=stdi(v,mean_v,n);
  out[2]=mini(v,n);
  out[3]=maxi(v,n);
}

float compute_percentile_f(float x,float *xv,int n)
{
  /* given n values in xv, compute the percentile for x 
   * percentile = percentage of entries in xv that are larger or equal to x 
   * usage: p=compute_percentile_f(x,xv,n);
   */
  
  int i;
  float p=0.0;
  
  for (i=1;i<=n;i++) {
    if (xv[i]>=x) 
      p++;
  }
  p=p/(float)n;
  return p;
}

float compute_percentile_i(int x,int *xv,int n)
{
  /* given n values in xv, compute the percentile for x 
   * percentile = percentage of entries in xv that are larger or equal to x 
   * usage: p=compute_percentile_f(x,xv,n);
   */
  
  int i;
  float p=0.0;
  
  for (i=1;i<=n;i++) {
    if (xv[i]>=x) 
      p++;
  }
  p=p/(float)n;
  return p;
}

void sort_cvector_by_indices(int n,char *v,float *indices)
{
  /* sort a (char) vector v by the indices given in another vector called indices
   * here indices is a float to directly accomodate the output of a program like quicksorti
   * usage:
   * sort_cvector_by_indices(int n,char *v,float *indices);
   */
  
  int i;                 // loop index
  int curr_index;        // current index
  char *temp_v;           // temporary vector to copy the values of v

  temp_v=cvector(1,n);
  for (i=1;i<=n;i++)
    temp_v[i]=v[i];
  for (i=1;i<=n;i++) {
    curr_index=(int)indices[i];
    v[i]=temp_v[curr_index];
  }

  free_cvector(temp_v,1,n);
}

void sort_ivector_by_indices(int n,int *v,float *indices)
{
  /* sort a vector v by the indices given in another vector called indices
   * here indices is a float to directly accomodate the output of a program like quicksorti
   * usage:
   * sort_ivector_by_indices(int n,int *v,float *indices);
   */
  
  int i;                 // loop index
  int curr_index;        // current index
  int *temp_v;           // temporary vector to copy the values of v

  temp_v=ivector(1,n);
  for (i=1;i<=n;i++)
    temp_v[i]=v[i];
  for (i=1;i<=n;i++) {
    curr_index=(int)indices[i];
    v[i]=temp_v[curr_index];
  }

  free_ivector(temp_v,1,n);
}

void sort_ivector_by_int_indices(int n,int *v,int *indices)
{
  /* sort a vector v by the indices given in another vector called indices
   * here indices is a float to directly accomodate the output of a program like quicksorti
   * usage:
   * sort_ivector_by_indices(int n,int *v,float *indices);
   */
  
  int i;                 // loop index
  int curr_index;        // current index
  int *temp_v;           // temporary vector to copy the values of v

  temp_v=ivector(1,n);
  for (i=1;i<=n;i++)
    temp_v[i]=v[i];
  for (i=1;i<=n;i++) {
    curr_index=(int)indices[i];
    v[i]=temp_v[curr_index];
  }

  free_ivector(temp_v,1,n);
}

long time_difference(struct timeb ti,struct timeb tf)
{
  /* compute time difference between two points */
  long tisec;
  int timsec;
  long tfsec;
  int tfmsec;
  long tdif;

  tisec=ti.time;
  timsec=ti.millitm;
  tfsec=tf.time;
  tfmsec=tf.millitm;

  if (timsec > tfmsec) {
    tfsec--;
    tfmsec+=1000;
  }
  tdif=(tfsec-tisec)*1000+(long)tfmsec-(long)timsec;

  return tdif;
}

int vector_overlap(int *original_vector,int n,int max_overlap,int *output_vector)
{
  /* given a vector of n integer values, original_vector
   * keep only entries where the distance between adjacent entries is more than max_overlap
   * this can be used to remove duplicates by setting max_overlap=0
   * original_vector is assumed to be sorted in ascending order
   * returns the number of entries kept and those entries in output_vector
   */

  int i;
  int d;
  int n_kept=1;

  output_vector[n_kept]=original_vector[1];
  for (i=2;i<=n;i++) {
    if ( (original_vector[i]-original_vector[i-1]) > max_overlap ) {
      n_kept++;
      output_vector[n_kept]=original_vector[i];
    }
  }

  return n_kept;
}

int max2(int x1,int x2) 
{
  /* returns the maximum of two integers */
  if (x1>x2) {
    return x1;
  } else {
    return x2;
  }
}

int max3(int x1,int x2,int x3) 
{
  /* returns the maximum of three integers */
  int m;

  m=x1;
  if (x2>m) 
    m=x2;
  if (x3>m) {
    return x3;
  } else {
    return m;
  }
}

int max4(int x1,int x2,int x3,int x4) 
{
  /* returns the maximum of four integers */
  int m;

  m=x1;
  if (x2>m) 
    m=x2;
  if (x3>m) 
    m=x3;
  if (x4>m) {
    return x4;
  } else {
    return m;
  }
}

int min2(int x1,int x2) 
{
  /* returns the minimum of two integers */
  if (x1<x2) {
    return x1;
  } else {
    return x2;
  }
}

int min3(int x1,int x2,int x3) 
{
  /* returns the minimum of three integers */
  int m;

  m=x1;
  if (x2<m) 
    m=x2;
  if (x3<m) {
    return x3;
  } else {
    return m;
  }
}

int min4(int x1,int x2,int x3,int x4) 
{
  /* returns the minimum of four integers */
  int m;

  m=x1;
  if (x2<m) 
    m=x2;
  if (x3<m) 
    m=x3;
  if (x4<m) {
    return x4;
  } else {
    return m;
  }
}

float range3f(float v1,float v2,float v3)
{
  /* range 3
   * computes the range (max-min) among 3 values
   * usage:
   * range=range3(float v1,float v2, float v3);
   */

  float range;
  float max;
  float min;
  
  if (v1>v2) {
    if (v1>v3) {
      max=v1;
      if (v2>v3) {
	min=v3;
      } else {
	min=v2;
      } 
    } else {
      max=v3;
      min=v2;
    }
  } else {
    if (v2>v3) {
      max=v2;
      if (v1>v3) {
	min=v3;
      } else {
	min=v1;
      }
    } else {
      max=v3;
      min=v1;
    }
  }
  
  range=max-min;
  
  return range;
}


float rangef(float *x,int n) 
{
  float minx=1000000;
  float maxx=-1000000;
  int i;
  float range;
  
  for (i=1;i<=n;i++) {
    if (x[i]<minx) 
      minx=x[i];
    if (x[i]>maxx)
      maxx=x[i];
  }
  
  range=maxx-minx;
  
  return range;
}

float meanf(float *x,int n)
{
	/* returns the mean of vector x */
	float m_x;
	int i;
	
	m_x=0;
	for (i=1;i<=n;i++)
	    m_x+=x[i];
	m_x=(float)m_x/(float)n;
	
	return m_x;
}

float meani(int *x,int n)
{
	/* returns the mean of vector x */
	float m_x;
	int i;

	m_x=0;
	for (i=1;i<=n;i++)
	    m_x+=x[i];
	m_x=(float)m_x/(float)n;

	return m_x;
}

float sumi(int *x,int n)
{
	/* returns the mean of vector x */
	float m_x;
	int i;

	m_x=0;
	for (i=1;i<=n;i++)
	    m_x+=x[i];

	return m_x;
}

float stdf(float *x,float m_x,int n)
{
  /* returns the standard deviation of x */
  float s_x;
  int i;
  
  s_x=0;	
  for (i=1;i<=n;i++)
    s_x=s_x+(x[i]-m_x)*(x[i]-m_x);
  s_x=(float)s_x/(float)(n-1);
  s_x=(float)sqrt(s_x);
  
  return s_x;
}

float stdi(int *x,float m_x,int n)
{
	/* returns the standard deviation of x */
	float s_x;
	int i;

	s_x=0;
	for (i=1;i<=n;i++)
	    s_x=s_x+(x[i]-m_x)*(x[i]-m_x);
	s_x=(float)s_x/(float)(n-1);
	s_x=(float)sqrt(s_x);

	return s_x;
}

void shell_sort(long n,float *arr,int *brr)
{
	/* sort an array arr into ascending numerical order by Shell's method
	   arr is replaced by the sorted list
	   brr can be used to return the indices in the original file of the sorted values
	*/
	
	int i,j,inc;
	float v;
	int w;
	inc=1;

	do {
	    inc *= 3;
	    inc++;
	}
	while (inc<=n);
	do {
	  inc/=3;
	  for (i=inc+1;i<=n;i++) {
			v=arr[i];
			w=brr[i];
			j=i;
			while (arr[j-inc]>v) {
		    arr[j]=arr[j-inc];
		    brr[j]=brr[j-inc];
		    j-=inc;
		    if (j<=inc)
		    	break;
			}
			arr[j]=v;
			brr[j]=w;
		}
	}
	while (inc>1);
}

int fnormalize_vector(float *v,int n)
{
  /* normalize input vector v so that v.1=0 and norm(v)=1 */
  float mean_v;
  float std_v;
  float var_v;
  int i;
  int errmsg=0;
  
  /* compute mean value */
  mean_v=meanf(v,n);
  std_v=stdf(v,mean_v,n);
  var_v=std_v*std_v;

  if (var_v==0) {
    errmsg=1;
  }
  else	{
    for (i=1;i<=n;i++) {
      v[i]=(v[i]-mean_v)/sqrt((n-1)*var_v);
    }
  }
  return errmsg;
}

int inormalize_vector(int *v,float *nv,int n)
{
	/* normalize input vector v so that v.1=0 and norm(v)=1 */
	float mean_v;
	float std_v;
	float var_v;
	int i;
	int errmsg=0;
	
	/* compute mean value */
	mean_v=meani(v,n);
	std_v=stdi(v,mean_v,n);
	var_v=std_v*std_v;
	if (var_v==0)	{
		errmsg=1;
	}
	else	{
		for (i=1;i<=n;i++) {
		    nv[i]=(v[i]-mean_v)/sqrt((n-1)*var_v);
		}
	}
	return errmsg;
}

int fast_inormalize_vector(int *v,float *nv,float mean_v,float std_v,float sqrtnm1,int n)
{
	/* normalize input vector v so that v.1=0 and norm(v)=1 */
	int i;
	int errmsg=0;
	float denom;
	
	if (std_v==0)	{
		errmsg=1;
	}
	else	{
		denom=sqrtnm1*std_v;
		for (i=1;i<=n;i++) {
		    nv[i]=(v[i]-mean_v)/denom;
		}
	}
	return errmsg;
}

float normf(float *v,int n)
{
	int i;
	float norm_v=0;
	
	for (i=1;i<=n;i++)
		norm_v+=v[i]*v[i];
	norm_v=sqrt(norm_v);
		
	return norm_v;				
}

float normi(int *v,int n)
{
	int i;
	float norm_v=0;
	
	for (i=1;i<=n;i++)
		norm_v+=v[i]*v[i];
	norm_v=sqrt(norm_v);
		
	return norm_v;				
}

float minf(float *x,int n)
{
	/* returns the min of vector x */
	float m_x;
	int i;
	
	m_x=1000000;
	for (i=1;i<=n;i++)	{
		if (x[i]<m_x)
			m_x=x[i];
	}

	return m_x;
}

int mini(int *x,int n)
{
	/* returns the min of vector x */
	int m_x;
	int i;
	
	m_x=1000000;
	for (i=1;i<=n;i++)	{
		if (x[i]<m_x)
			m_x=x[i];
	}

	return m_x;
}

long minl(long *x,int n)
{
  /* returns the min of vector x */
  long m_x;
  int i;
	
  m_x=10000000;
  for (i=1;i<=n;i++)	{
    if (x[i]<m_x)
      m_x=x[i];
  }
  return m_x;
}

float maxf(float *x,int n)
{
  /* returns the max of vector x */
  float m_x;
  int i;
  
  m_x=-1000000;
  for (i=1;i<=n;i++)	{
    if (x[i]>m_x)
      m_x=x[i];
  }
  
  return m_x;
}

int maxi(int *x,int n)
{
	/* returns the min of vector x */
	int m_x;
	int i;
	
	m_x=-1000000;
	for (i=1;i<=n;i++)	{
		if (x[i]>m_x)
			m_x=x[i];
	}

	return m_x;
}

long maxl(long *x,int n)
{
  /* returns the max of vector x */
  long m_x;
  int i;
	
  m_x=-1000000;
  for (i=1;i<=n;i++)	{
    if (x[i]>m_x)
      m_x=x[i];
  }
  return m_x;
}

float ran1(long *idum)
{
  /* Park and Miller random number generator */
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <=0 || !iy) {
    if (-(*idum)<1) *idum=1;
    else *idum =-(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0) *idum+=IM;
      if (j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum<0) *idum += IM;
  j=iy/NDIV1;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

void randperms(int n,int *index,long idum)
{
  /* generate random permutations of numbers 1...n */
  int i;
  float f;
  /* long idum; */
  /*time_t currtime;*/
  float *r;

  r=vector(1,n);
  /* time(&currtime);
  idum=-currtime; */

  for (i=1;i<=n;i++)
    index[i]=i;

  for (i=1;i<=n;i++){
    f=ran1(&idum);
    r[i]=f;
  }

  shell_sort(n,r,index);
  free_vector(r,1,n);
}

void quicksort(unsigned long n,float *arr,int *brr)
{
  unsigned long i,ir=n,j,k,b,l=1;
  long *istack;
  int jstack=0;
  float a,temp;	
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l]>arr[l+1])	{
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

void quicksortd(unsigned long n,double *arr,int *brr)
{
  unsigned long i,ir=n,j,k,b,l=1;
  long *istack;
  int jstack=0;
  //float a,temp;
  double a,temp;
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l]>arr[l+1])	{
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

void just_quicksortus(unsigned long n,unsigned short *arr)
{
  /* sort the entries in arr in ascending order 
   * here, do not simultaneously rearrange an int as in the quicksort version
   *       also explicitly replace SWAP
   *       also use unsigned short 
   */

  unsigned long i,ir=n,j,k,l=1;
  long *istack;
  int jstack=0;
  float a,temp;	
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      temp=arr[k];arr[k]=arr[l+1];arr[l+1]=temp;
      if (arr[l]>arr[ir]) {
	temp=arr[l];arr[l]=arr[ir];arr[ir]=temp;
      }
      if (arr[l+1]>arr[ir]) {
	temp=arr[l+1];arr[l+1]=arr[ir];arr[ir]=temp;
      }
      if (arr[l]>arr[l+1]) {
	temp=arr[l];arr[l]=arr[l+1];arr[l+1]=temp;
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	temp=arr[i];arr[i]=arr[j];arr[j]=temp;
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

void quicksortf(unsigned long n,float *arr,float *brr)
{
  unsigned long i,ir=n,j,k,l=1;
  float b;
  long *istack;
  int jstack=0;
  float a,temp;	
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l]>arr[l+1])	{
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

void quicksorti(unsigned long n,int *arr,float *brr)
{
  unsigned long i,ir=n,j,k,l=1;
  float b;
  long *istack;
  int jstack=0;
  float a,temp;	
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l]>arr[l+1])	{
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

void quicksortl(unsigned long n,long *arr,float *brr)
{
  unsigned long i,ir=n,j,k,l=1;
  float b;
  long *istack;
  int jstack=0;
  float a,temp;	
  
  istack=lvector(1,NSTACK);
  for (;;)	{
    if (ir-l < M)	{
      for (j=l+1;j<=ir;j++)	{
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--)	{
	  if (arr[i]<=a)
	    break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (jstack==0)
	break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }
    else	{
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l]>arr[l+1])	{
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;)	{
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i)
	  break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack+=2;
      if (jstack>NSTACK)
	nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-1)	{
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }
      else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }	
    }
  }
  free_lvector(istack,1,NSTACK);
}

int myround(float f)
{
	float ceil_f;
	float floor_f;
	int roundf;
	
	ceil_f=ceil(f);
	floor_f=floor(f);
	
	if ( (fabs(f-ceil_f)) < (fabs(f-floor_f)) )
		roundf=(int)ceil_f;
	else
		roundf=(int)floor_f;	

	return roundf;	
}

int factorial (int n)
{
	/* factorial of n */

	int i;
	int temp=1;
	
	for (i=1;i<=n;i++)
		temp*=i;

	return temp;	
}

double binopdf_mat(int n,int k,double p) {
  double res;
  double nk;
  double dn;
  double dk;
  double dnk;
  double lny;
  double cp;
  double nk1,nk2,nk3;

  if (k>n) {
    //printf("error! k must be <= n\n(currently n=%d,k=%d). setting k=n\n",n,k);
    k=n;
  }
  if (p==1.0) {
    res=0.0;
  } else {
    if (p==0.0) {
      res=0.0;
    } else {
      dn=(double)n+1.0;
      dk=(double)k+1.0;
      dnk=dn-dk+1.0;
      nk=gammaln_mat(dn)-gammaln_mat(dk)-gammaln_mat(dnk);
      lny= nk + (double)k*log(p)+((double)n-(double)k)*log(1.0-p);
      res=exp(lny);
    }
  }
  return res;
}

float binopdf(int n,int k,float p) {
  /* returns the probability of obtaining k heads out of n trials if the probability of head is p */
  float q;
  double coef;
  float binop;
  float pk;
  int j;
  float qnk;
  
  if (p==0) {
    binop=0;
  } else {
    j=n-k;
    q=1-p;
    pk=fpower(p,k);  // p^k
    qnk=fpower(q,j); // q^(n-k)
    coef=dbico(n,k);
    
    binop=coef*pk*qnk;
  }
  return binop;
}

int nchoosek(int n,int k)
{
	/* number of combinations of k elements among n */
	int nmk,nf,kf,nmkf;
	int temp;
	
	nmk=n-k;
	nf=factorial(n);
	kf=factorial(k);
	nmkf=factorial(nmk);
	
	temp=(int)nf/(kf*nmkf);
	
	return temp;
}

float erfinv(float y)
{
  
  // erfinv, Inverse error function.
  // x = erfinv(y) is the inverse error function for x.
  // The inverse error functions satisfies y = erf(x), for -1 <= y < 1
  // and -inf <= x <= inf.
  // Translated by G.K. from MATLAB code (The Mathworks)

  float a[4];
  float b[4];
  float c[4];
  float d[2];
  float y0=0.7;
  float z;
  float x=0;

  // coefficients in rational approximations;

  if (y==-1)
    return -1000;
  if (y==1)
    return 1000;
  if (fabs(y)>1)
    nrerror("error in my_math_package/erfinv:\tabs(y) must be <=1\n");

  a[0]=0.886226899;a[1]=-1.645349621;a[2]=0.914624893;a[3]=-0.140543331;
  b[0]=-2.118377725;b[1]=1.442710462;b[2]=-0.329097515;b[3]=0.012229801;
  c[0]=-1.970840454;c[1]=-1.624906493;c[2]=3.429567803;c[3]=1.641345311;
  d[0]=3.543889200;d[1]=1.637067800;

  // central range 
  if (fabs(y)<=y0) {
    z=y*y;
    x=y*(((a[3]*z+a[2])*z+a[1])*z+a[0])/((((b[3]*z+b[2])*z+b[1])*z+b[0])*z+1);
    return x;
  }

  // near end points of range
  if ( (y>y0) & (y<1) ) {
    z=sqrt(-log((1-y)/2));
    x=(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1);
    return x;
  }

  if ( (-y0>y) & (y>-1) ) {
    z = sqrt(-log((1+y)/2));
    x=-(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1);
    return x;
  }

  // Two steps of Newton-Raphson correction to full accuracy.
  // Without these steps, erfinv(y) would be about 3 times
  // faster to compute, but accurate to only about 6 digits.
  // x = x - (erf(x) - y) ./ (2/sqrt(pi) * exp(-x.^2));
  // x = x - (erf(x) - y) ./ (2/sqrt(pi) * exp(-x.^2));

  return x;
}

float erf1(float x)
{
  /* compute error function */
  float y;
  y=erfcore(x,0);
  return y;
}

float erfc1(float x)
{
  float y;
  y=erfcore(x,1);
  return y;
}

float erfcore(float x,int jint)
{ 
  // ERFCORE Core algorithm for error functions.
  //   erf(x) = erfcore(x,0)
  //   erfc(x) = erfcore(x,1)
  //   This is a translation by Gabriel Kreiman of matlab code translating 
  //   itself a FORTRAN program by W. J. Cody, Argonne National Laboratory, NETLIB/SPECFUN, March 19, 1990.
  //   The main computation evaluates near-minimax approximations
  //   from "Rational Chebyshev approximations for the error function"
  //   by W. J. Cody, Math. Comp., 1969, PP. 631-638.

  float xbreak=0.46875;
  float a[5];
  float b[4];
  float c[9];
  float d[8];
  float p[6];
  float q[5];
  float y,z;
  float xnum;
  float xden;
  int i;
  float result;
  float del;

  //    evaluate  erf  for  |x| <= 0.46875
  if (fabs(x)<=xbreak) {
    a[0]=3.16112374387056560e00;a[1]=1.13864154151050156e02;a[2]=3.77485237685302021e02;a[3]=3.20937758913846947e03;a[4]=1.85777706184603153e-1;
    b[0]=2.36012909523441209e01;b[1]=2.44024637934444173e02;b[2]=1.28261652607737228e03;b[3]=2.84423683343917062e03;
    
    y=fabs(x);
    z=y*y;
    xnum=a[4]*z;
    xden=z;
    for (i=0;i<=2;i++) {
      xnum=(xnum+a[i])*z;
      xden=(xden+b[i])*z;
    }
    result=x*(xnum+a[3])/(xden+b[3]);
    if (jint!=0) 
      result=1-result;
    
  }
   
  //	 evaluate  erfc  for 0.46875 <= |x| <= 4.0
  if ( (fabs(x)>xbreak) & (fabs(x)<=4.0) ) {
    c[0]=5.64188496988670089e-1;c[1]=8.88314979438837594e00;c[2]=6.61191906371416295e01;c[3]=2.98635138197400131e02;
    c[4]=8.81952221241769090e02;c[5]=1.71204761263407058e03;c[6]=2.05107837782607147e03;c[7]=1.23033935479799725e03;
    c[8]=2.15311535474403846e-8;
    d[0]=1.57449261107098347e01;d[1]=1.17693950891312499e02;d[2]=5.37181101862009858e02;d[3]=1.62138957456669019e03;
    d[4]=3.29079923573345963e03;d[5]=4.36261909014324716e03;d[6]=3.43936767414372164e03;d[7]=1.23033935480374942e03;
  
    y=fabs(x);
    xnum=c[8]*y;
    xden=y;
    for (i=0;i<=6;i++) {
      xnum=(xnum+c[i])*y;
      xden=(xden+d[i])*y;
    }
    
    result=(xnum+c[7])/(xden+d[7]);
    z=floor(y*16)/16;
    del=(y-z)*(y+z);
    result=exp(-z*z)*exp(-del)*result;
  }     
  
  //  evaluate  erfc  for |x| > 4.0
  if (fabs(x)>4.0) {
    p[0]=3.05326634961232344e-1;p[1]=3.60344899949804439e-1;p[2]=1.25781726111229246e-1;p[3]=1.60837851487422766e-2;
    p[4]=6.58749161529837803e-4;p[5]=1.63153871373020978e-2;
    q[0]=2.56852019228982242e00;q[1]=1.87295284992346047e00;q[2]=5.27905102951428412e-1;q[3]=6.05183413124413191e-2;q[4]=2.33520497626869185e-3;

    y=fabs(x);
    z=1/(y*y);
    xnum=p[5]*z;
    xden=z;
    for (i=0;i<=3;i++) {
      xnum =(xnum+p[i])*z;
      xden =(xden+q[i])*z;
    }
    result=z*(xnum+p[4])/(xden+q[4]);
    result=(1/sqrt(PI)-result)/y;
    z=floor(y*16)/16;
    del = (y-z)*(y+z);
    result=exp(-z*z)* exp(-del)* result;
    
  }
  
  // basic conditions and parameter validation
  if (jint==0) {
    if (x>xbreak)
      result=1-result;
    if (x<-xbreak)
      result=-1+result;
  } else {
    if (x<-xbreak)
      result=2-result;
  }
  
  return result;

}

float roc_gauss(float fa,float mu0,float sigma0,float mu1,float sigma1)
{
  /* computes the receiver operator characteristic assuming that the two random variables
   * are gaussian, according to the formula:
   * y = 1- erf(p1*erfinv(1-fa)-p2)
   * where p1=sigma0/sigma1 and p2=(mu1-mu0)/sigma1
   */

  float p1,p2;
  float phiinv;
  float aux;
  float y;
  float phiy;
  float p;

  p1=sigma0/sigma1;
  p2=( (mu1-mu0)/sigma1 )/sqrt2;

  aux=2*(1-fa)-1;
  phiinv=erfinv(aux);

  y=p1*phiinv-p2;
  phiy=erf1(y);

  p=1-0.5*(1+phiy);
  
  return p;
}

void roc_gauss_vector(float *fa,float *roc,int n_fa,float mu0,float sigma0,float mu1,float sigma1)
{
  /* computes the receiver operator characteristic assuming that the two random variables
   * are gaussian, according to the formula:
   * y = 1- erf(p1*erfinv(1-fa)-p2)
   * where p1=sigma0/sigma1 and p2=(mu1-mu0)/sigma1
   */

  /* this version computes the roc for all values of fa 
   * also note the change from p to roc for the output
   */

  float p1,p2;
  float phiinv;
  float aux;
  float y;
  float phiy;
  float p;
  int i;
  
  p1=sigma0/sigma1;
  p2=( (mu1-mu0)/sigma1 )/sqrt2;

  for (i=1;i<=n_fa;i++) {
    aux=2*(1-fa[i])-1;
    phiinv=erfinv(aux);
    
    y=p1*phiinv-p2;
    phiy=erf1(y);
    
    roc[i]=1-0.5*(1+phiy);
  }
}

float pe_gauss(float mu0,float sigma0,float mu1,float sigma1,float res)
{
  /* computes the probability of misclassification based on ROC analysis 
   * assuming that the two random variables are gaussian, according to the formula:
   * pcd = 1- erf(p1*erfinv(1-pfa)-p2)
   * perror = ( pfa + (1-pcd) )/ 2;
   * pe=min(perror);
   * where p1=sigma0/sigma1 and p2=(mu1-mu0)/sigma1
   */

  float p1,p2;
  float phiinv;
  float aux;
  float y;
  float phiy;
  float perror;
  float pfa;
  float pe=100;
  float pcd;

  p1=sigma0/sigma1;
  p2=( (mu1-mu0)/sigma1 )/sqrt2;

  for (pfa=0;pfa<=1;pfa+=res) {
    aux=2*(1-pfa)-1;
    phiinv=erfinv(aux);
    y=p1*phiinv-p2;
    phiy=erf1(y);
    pcd=1-0.5*(1+phiy);
    perror=(pfa+(1-pcd))/2;
    if (perror<pe)
      pe=perror;
  }
  return pe;
}

char char3min(char x,char y,char z) {
  /* return the minimum of the three char values */

  char minchar;

  if (x<y) {
    if (x<z) {
      minchar=x;
    } else {
      minchar=z;
    }
  } else {
    if (y<z) {
      minchar=y;
    } else {
      minchar=z;
    }
  }
 
  return minchar;
}

int findi(int *vec,int n_vec,int search_value,int *index) 
{
  /* find all occurrences of search_value within vector vec and 
   * return the number of occurrences and the indices in index
   */

  int i;
  int n;

  n=0;
  for (i=1;i<=n_vec;i++) {
    if (vec[i]==search_value) {
      n++;
      index[n]=i;
    }
  }
  return n;
}

void linetrend(float *x,float *y,int n,float *res) 
{
  /* This function returns the linear regression analysis of vectors x,y.
   * The model is y= m * x + b; 
   * The program computes the (m,b) so that (yest-y)^2 is minimized.
   * Input parameters: x,y which are vectors of length n.
   * Output : slope m, intercept b and r^2 value.
   */

  float mean_x,mean_y,std_x,std_y;
  float sum_x,sum_y,sum_x2,sum_y2;
  float cov,sum_xy;
  int i;

  /* compute sum_x, sum_y, sum_x2, sum_y2 */
  sum_x=0;sum_y=0;sum_x2=0;sum_y2=0;
  for (i=1;i<=n;i++) {
    sum_x+=x[i];
    sum_x2+=(x[i]*x[i]);
    sum_y+=y[i];
    sum_y2+=(y[i]*y[i]);
  }
  
  /* compute mean_x,std_x,mean_y,std_y */
  mean_x=sum_x/n;
  mean_y=sum_y/n;
  std_x=sum_x2/n-mean_x*mean_x;
  std_y=sum_y2/n-mean_y*mean_y;

  if ( (std_x<=0) | (std_y<=0) ) {
    res[1]=0;
    res[2]=0;
    res[3]=0;
  } else {
    std_x=sqrt(std_x); 
    std_y=sqrt(std_y);
   /* compute covariance and sum(xy)*/
    cov=0;sum_xy=0;
    for (i=1;i<=n;i++) {
      cov+=(x[i]-mean_x)*(y[i]-mean_y);
      sum_xy+=(x[i]*y[i]);
    }
    cov=cov/(n-1);
    res[1]=cov/(std_x*std_y); // r2, correlation coefficient
    
    res[2]=(sum_xy-sum_x*sum_y/n)/(sum_x2-sum_x*sum_x/n); // m, slope
    res[3]=mean_y-(res[2])*mean_x; // b, intercept
  }
}
  
float spearman(float *data1,float *data2,unsigned long n) 
{
  /* compute spearman correlation coefficient between x and y */
  float *x,*y;
  float mean_x,mean_y,std_x,std_y;
  float sum_x,sum_y,sum_x2,sum_y2;
  float cov,sum_xy;
  unsigned long i;
  float rs;

  x=vector(1,n);
  y=vector(1,n);
  for (i=1;i<=n;i++) {
    x[i]=data1[i];
    y[i]=data2[i];
  }
  myrank(n,x);
  myrank(n,y);
  
  /* compute sum_x, sum_y, sum_x2, sum_y2 */
  sum_x=0;sum_y=0;sum_x2=0;sum_y2=0;
  for (i=1;i<=n;i++) {
    sum_x+=x[i];
    sum_x2+=(x[i]*x[i]);
    sum_y+=y[i];
    sum_y2+=(y[i]*y[i]);
  }
  
  /* compute mean_x,std_x,mean_y,std_y */
  mean_x=sum_x/n;
  mean_y=sum_y/n;
  std_x=sum_x2/n-mean_x*mean_x;
  std_y=sum_y2/n-mean_y*mean_y;

  if ( (std_x <= 0) | (std_y <= 0) ) {
    rs=0;
  } else {
    std_x=sqrt(std_x*n/(n-1));
    std_y=sqrt(std_y*n/(n-1));
    /* compute covariance and sum(xy)*/
    cov=0;sum_xy=0;
    for (i=1;i<=n;i++) {
      cov+=(x[i]-mean_x)*(y[i]-mean_y);
      sum_xy+=(x[i]*y[i]);
    }
    cov=cov/(n-1);
    rs=cov/(std_x*std_y); // r2, correlation coefficient
    
  }

  free_vector(x,1,n);
  free_vector(y,1,n);

  return rs;
}

void crank(unsigned long n,float w[])
{
  /* given a sorted array w[1..n], replaces it by the ranks, including midranking of ties */
  unsigned long j=1,ji,jt;
  float rank;

  while (j<n) {
    //printf("%d\n",j);
    if (w[j+1]!=w[j]) {
      /* not a tie */
      w[j]=j;
      ++j;
    } else {
      for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++)
	w[ji]=rank;
    }
  }
  if (j==n)
    w[n]=n;
}

void myrank(long n,float *w)
{
  /* rank vector w 
   * returns the ranks of each entry in w
   * in the case of ties, returns the mean of the corresponding ranks
   */
  
  int *index;
  float *temp;
  int i1,i2,i,j,jt,ji;
  float rank;

  index=ivector(1,n);
  temp=vector(1,n);
  
  for (i=1;i<=n;i++) {
    temp[i]=w[i];
    index[i]=i;
  }
  quicksort(n,temp,index);

  j=1;
  while (j<n) {
    if (w[index[j]]!=w[index[j+1]]) {
      /* not a tie */
      w[index[j]]=j;
      j++;
    } else {
      for (jt=j+1;jt<=n && w[index[jt]]==w[index[j]];jt++);
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++) 
	w[index[ji]]=rank;
      j=jt;
    }
  }
  if (j==n)
    w[index[j]]=n;

  free_ivector(index,1,n);
  free_vector(temp,1,n);
}

float fpower(float x,int n)
{
  /* returns x^n */
  int i;
  float xn;
  float lx;

  //printf("x=%.2f\tn=%d\n",x,n);
  if (x==0) {
    xn=0;
  } else {
    if (n==0) {
      xn=1;
    } else {
      if (x<0) {
	xn=x;
	for (i=2;i<=n;i++) 
	  xn*=x;
      } else {
	//printf("computing x=%.2f\tn=%d\n",x,n);
	lx=log(x);
	xn=lx*n;
	xn=exp(xn);
	//printf("xn=%.4f\n",xn);
      }
    } 
  }
  if (n<0)
    xn=1/xn;

  return xn;
}

float binocdf(int n,int k,float p) {
  /* returns the cumulative probability distribution that a random variable is <=k following a binomial distribution with parameter p */
  int i;
  float binoc=0;
  float t;
  
  if (k>n) {
    //printf("error! k must be <= n\n(currently n=%d,k=%d). setting k=n\n",n,k);
    k=n;
  }
  for (i=0;i<=k;i++) {
    t=binopdf(n,i,p);
    binoc+=t;
  }
  return binoc;
}

double binocdf_mat(int n,int k,double p) {
  /* returns the cumulative probability distribution that a random variable is <=k following a binomial distribution with parameter p */
  int i;
  double binoc=0.0;
  double t;
  
  if (k>n) {
    //printf("error! k must be <= n\n(currently n=%d,k=%d). setting k=n\n",n,k);
    k=n;
  }
  for (i=0;i<=k;i++) {
    t=binopdf_mat(n,i,p);
    binoc+=t;
  }
  return binoc;
}

double gammaln_mat(double x) 
{
  static double d1 = -5.772156649015328605195174e-1;
  static double p1[9] = {0.0,4.945235359296727046734888e0, 2.018112620856775083915565e2, 2.290838373831346393026739e3, 1.131967205903380828685045e4, 2.855724635671635335736389e4, 3.848496228443793359990269e4, 2.637748787624195437963534e4, 7.225813979700288197698961e3};
  static double q1[9] = {0.0,6.748212550303777196073036e1, 1.113332393857199323513008e3, 7.738757056935398733233834e3, 2.763987074403340708898585e4, 5.499310206226157329794414e4, 6.161122180066002127833352e4, 3.635127591501940507276287e4, 8.785536302431013170870835e3};
  static double d2 = 4.227843350984671393993777e-1;
  static double p2[9] = {0.0,4.974607845568932035012064e0, 5.424138599891070494101986e2, 1.550693864978364947665077e4, 1.847932904445632425417223e5, 1.088204769468828767498470e6, 3.338152967987029735917223e6, 5.106661678927352456275255e6, 3.074109054850539556250927e6};
  static double q2[9] = {0.0,1.830328399370592604055942e2, 7.765049321445005871323047e3, 1.331903827966074194402448e5, 1.136705821321969608938755e6, 5.267964117437946917577538e6, 1.346701454311101692290052e7, 1.782736530353274213975932e7, 9.533095591844353613395747e6};
  static double d4 = 1.791759469228055000094023e0;
  static double p4[9] = {0.0,1.474502166059939948905062e4, 2.426813369486704502836312e6, 1.214755574045093227939592e8, 2.663432449630976949898078e9, 2.940378956634553899906876e10, 1.702665737765398868392998e11, 4.926125793377430887588120e11, 5.606251856223951465078242e11};
  static double q4[9] = {0.0,2.690530175870899333379843e3, 6.393885654300092398984238e5, 4.135599930241388052042842e7, 1.120872109616147941376570e9, 1.488613728678813811542398e10, 1.016803586272438228077304e11, 3.417476345507377132798597e11, 4.463158187419713286462081e11};
  static double c[8] = {0.0,-1.910444077728e-03, 8.4171387781295e-04, -5.952379913043012e-04, 7.93650793500350248e-04, -2.777777777777681622553e-03, 8.333333333333333331554247e-02, 5.7083835261e-03};
  double res;  
  int i;
  double xnum,xden;
  double xm1;
  double xm2;
  double xm4;
  double r;
  double ysq;
  double corr;
  double spi;

  res = x;
  // 0 <= x <= EPS
  if (x<EPS) 
    res = - log(x);

  // EPS < x <= 0.5
  if ( (x>EPS) & (x<=0.5) ) { 
    xden=1;
    xnum=0;
    for (i=1;i<=8;i++) {
      xnum = xnum * x + p1[i];
      xden = xden * x + q1[i];
    }
    res = -log(x) + (x*(d1+x*(xnum/xden)));
  }

  // 0.5 < x <= 0.6796875
  if ( (x>0.5) & ( x < 0.6796875) ) {
    xm1 = x - 0.5 - 0.5;
    xden = 1.0;
    xnum = 0.0;
    for (i=1;i<=8;i++) {
      xnum = xnum * xm1 + p2[i];
      xden = xden * xm1 + q2[i];
    }
    res = -log(x) + xm1 * (d2 + xm1 * (xnum / xden));
  }


  //  0.6796875 < x <= 1.5

  if ( (x > 0.6796875) & (x <= 1.5)) {
    xm1 = x - 0.5 - 0.5;
    xden = 1.0;
    xnum = 0.0;
    for (i = 1;i<=8;i++) {
      xnum = xnum * xm1 + p1[i];
      xden = xden * xm1 + q1[i];
    }
    res = xm1 * (d1 + xm1 * (xnum / xden));
  }

  //  1.5 < x <= 4

  if ((x > 1.5) & (x <= 4)) {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for (i=1;i<=8;i++) {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }
    res = xm2 * (d2 + xm2 * (xnum / xden));
  }

  
  // 4 < x <= 12
  if ((x > 4) & (x <= 12)) {
    xm4 = x - 4.0;
    xden = -1.0;
    xnum = 0.0;
    for (i = 1;i<=8;i++) {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }
    res = d4 + xm4 * (xnum / xden);
  }

  // x > 12
  if (x > 12) {
    r = c[7];
    ysq = x*x;
    for (i=1;i<=6;i++) 
      r = r / ysq + c[i];
    r = r / x;
    corr = log(x);
    spi = 0.9189385332046727417803297;
    res = r + spi - 0.5*corr + x * (corr-1);
  }
  return res;
}

int invbinocdf(double y,int n,double p)
{

//   X = BINOINV(Y,N,P) returns the inverse of the binomial cdf with 
//   parameters N and P. Since the binomial distribution is
//   discrete, BINOINV returns the least integer X such that 
//   the binomial cdf evaluated at X, equals or exceeds Y.
//
//   The size of X is the common size of the input arguments. A scalar input  
//   functions as a constant matrix of the same size as the other inputs.    
//
//   Note that X takes the values 0,1,2,...,N.

  int x=0;
  double cumdist=0;
  double currpdf;
  int count;

  if (y == 1) 
    x=n;
  if ( (n < 0) | (p<0) | (p>1) | (y<0.0) | (y>1.0) )
    x=-1;

  if (x == 0) {
    x=-1;    
    while (cumdist<y) {
      x++;
      currpdf=binopdf_mat(n,x,p);
      cumdist+=currpdf;
    }
  }
  
  return x;
}

double invbinocdf2(double y,int n,int x,double r)
{

//   p = invbinocdf2(Y,N,K) returns the inverse of the binomial cdf
// this is different from invbinocdf:
// invbinocdf: What is the x that will cause a cdf of at least y for a binomial with n,p
// invbinocdf2: What is the p that will cause a cdf of at least y for a binomial with k occurrences out of n possible ones?

  double p=1.0;
  double cumdist=0;
  double currpdf;
  int count;
  
  if ( (y == 1) & (x<n) )
    p=1.0;
  if ( (n < 0) | (p<0) | (p>1) | (y<0.0) | (y>1.0) )
    p=-1.0;

  //printf("n=%d\ty=%1.2g\tx=%d\n",n,y,x);
  if (p == 1.0) {
    p=1.0+r;
    while ( (cumdist<y) & (p>=r) ) {
      p-=r;     
      cumdist=binocdf_mat(n,x,p);
      //cumdist+=currpdf;
      //printf("%1.2g\t%1.2g\n",p,cumdist);
    }
  }
  
  if (p<0)
    p=0;
  
  return p;
}

float ran2(long *idum) {
  /* Long period (> 2 × 10^18 )random number generator of L’Ecuyer with Bays-Durham shuffle
     and added safeguards.Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
     the endpoint values).Call with idum a negative integer to initialize;thereafter,do not alter
     idum between successive deviates in a sequence.RNMX should approximate the largest floating
     value that is less than 1
  */
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) {                  // Initialize.
    if (-(*idum) < 1) 
      *idum=1;                       // Be sure to prevent idum =0
    else 
      *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {        // Load the shuffle table (after 8 warm-ups)
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) 
	*idum += IM1;
      if (j < NTAB) 
	iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                     // Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1;     // Compute idum=(IA1*idum) % IM1 without overflows y Schrage ’s method.
  if (*idum < 0) 
    *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;     // Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) 
    idum2 += IM2;
  j=iy/NDIV;                         // Will be in the range 0..NTAB-1
  iy=iv[j]-idum2;                    // Here idum is shuffled,idum and idum2 are combined to generate output.
  iv[j] = *idum;
  if (iy < 1) 
    iy += IMM1;
  if ((temp=AM*iy) > RNMX) 
    return RNMX;                     // Because users don ’t expect endpoint values.
  else 
    return temp;
}
