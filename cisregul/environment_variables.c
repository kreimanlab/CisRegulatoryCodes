/* environment_variables.c 
 */

void get_environment_variables(char *DATA_DIR,char *TEMP_DIR,char *CODE_DIR);

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void get_environment_variables(char *DATA_DIR,char *TEMP_DIR,char *CODE_DIR)
{
  sprintf(DATA_DIR,"/data/kreiman");
  sprintf(TEMP_DIR,"/data/kreiman/temp");
  sprintf(CODE_DIR,"/cbcl/cbcl01/kreiman/life");
}

/* void get_environment_variables(char *DATA_DIR,char *TEMP_DIR,char *CODE_DIR)
   {
   char envname[200];
   
   sprintf(envname,"DATADIR");
   DATA_DIR=getenv(envname);
   printf("DATA_DIR=%s\n",DATA_DIR);
   sprintf(envname,"CODEDIR");
   CODE_DIR=getenv(envname);
   sprintf(envname,"TEMPDIR");
   TEMP_DIR=getenv(envname);
   }
*/


