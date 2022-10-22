#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/****************************************************************/
void    write_lastr() {
/****************************************************************/

/*  this subroutine is called at the end of each run, in order to
    write the final coordinates of the particles on disk */

    int i,j,k;
    char cord[200];
    char* filename;
    FILE* file;

    cord[0] = 'r';
    cord[1] = '_';
    cord[2] = 'o';
    cord[3] = 'u';
    cord[4] = 't';
    cord[5] = '.';
    cord[6] = '\0';

    filename = strcat(cord,System_name);

    if((file = fopen(cord,"w"))==NULL)
        printf("can't open %s ...\n", cord);
    else {
        for(i=0;i<Nr_Stot;i++)
            fprintf(file,"%le %le %le %le %le %le\n",
                RS[i].x,RS[i].y,RS[i].z,RSo[i].x,RSo[i].y,RSo[i].z);

        fclose(file);
    }
}

