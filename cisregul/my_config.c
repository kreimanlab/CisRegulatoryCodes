#include <stdio.h>

void config(void);
void signature(void);

/* global variables */
char gen_dir[100];
char data_dir[100];
int powers_of_four[15];
float sqrt2;

void config()
{
  int i,temp;
  
  // sprintf(gen_dir,"/local/home/gabriel/life");
  sprintf(gen_dir,"c:\\life");
  //sprintf(data_dir,"/data");
  sprintf(data_dir,"c:\\life\\databh");
  temp=1;
  for (i=0;i<15;i++) {
    powers_of_four[i]=temp;
    temp*=4;
  }
  sqrt2=sqrt(2);
}

void signature()
{
  printf("\n\tGabriel Kreiman\n");
  printf("\tgabriel@klab.caltech.edu\n");
  printf("\tcheers\n");
}
