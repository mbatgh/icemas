#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

void        get_r0();


/****************************************************************/
void  get_r0() {
/****************************************************************/

    int k,acc,i,j,l;
    Real dis;
    FILE* file;
    VEKTOR d;

    if((file=fopen(R0_file,"r"))==NULL) {
        printf("can't see no %s\n", R0_file); 
        exit(1); 
    }
    else for(k=0;k<Nr_Stot;k++) {
        fscanf(file,"%le %le %le %le %le %le", &(RS[k].x),&(RS[k].y),
            &(RS[k].z),&(RSo[k].x),&(RSo[k].y),&(RSo[k].z));
    }

    fclose(file);


/* const. ewald-self-energy ...................................... */

    ECself = 0.0;
    acc = 0;

    for(i=0;i<Nr_Stot;i++)
        ECself -= elch[i]*elch[i]*eta/sqrt(PI);

    for(i=0;i<Nr_Spec;i++) {
        for(j=0;j<Nr_Molecules[i];j++) {
                
            for(k=0;k<Nr_Sites[i]-1;k++) {
                for(l=k+1;l<Nr_Sites[i];l++) {
                    d = RS[acc+k]-RS[acc+l];
                    convolute(d,adLh,L);
                    dis = sqrt(norm(d));
                    ECself -= elch[acc+k]*elch[acc+l]*erf(eta*dis)/dis;
                }
            }
            acc += Nr_Sites[i];
        }
    }

    ECself /= (Real)Nr_Stot;

    printf("read coordinates from file %s ...\n", R0_file);
}

