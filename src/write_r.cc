#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

/****************************************************************/
void    write_r(int i1,int i2,char* filename) {
/****************************************************************/

    int i;
    FILE* file;

    if((file = fopen(filename,"w"))==NULL) {
        printf("can't open %s ...\n", filename);
    }

    for(i=i1;i<i2;i++)
        fprintf(file,"%le %le %le %le %le %le\n",
            RS[i].x,RS[i].y,RS[i].z,RSo[i].x,RSo[i].y,RSo[i].z);

}

