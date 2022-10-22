#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/* =================================================================== */
Real    force() {
/* =================================================================== */

    int     o1,o2,n,a,b,i,j,k,l,spec1,spec2;
    Real    SS2,pot,dljpot,rm_12,rm_6,cpu_time;
    VEKTOR  Rij,forc,dx,dR,dr;
    clock_t t1,t2;

    t1 = clock();

    Virial = 0.0;
    lj_Pot = 0.0;


        for(i=0;i<Nr_Stot-1;i++) {

            spec1 = ljt[i];

            for(j=nextidx[i];j<Nr_Stot;j++) {

                spec2 = ljt[j];

                Rij = RS[i]-RS[j];
                convolute(Rij,adLh,L);
                SS2 = norm(Rij);

                if(SS2<cutoff_2) {

                    rm_6 = 1.0/(SS2*SS2*SS2);
                    rm_12 = rm_6*rm_6;
                    pot = LJP12[spec1][spec2]*rm_12 - LJP6[spec1][spec2]*rm_6;
                    lj_Pot += pot; 
                    dljpot = LJF12[spec1][spec2] * rm_12 
                             - LJF6[spec1][spec2] * rm_6;
                    Virial += dljpot;
                    forc = Rij*(dljpot/SS2);

                    FS[i] = FS[i] + forc;
                    FS[j] = FS[j] - forc;
                }
            }
        }
        
        for(n=0;n<NR_LJIN;n++) {

            o1 = ljnb_i[n];
            o2 = ljnb_j[n];
            spec1 = ljt[o1];
            spec2 = ljt[o2];

            Rij = RS[o1]-RS[o2];
            convolute(Rij,adLh,L);
            SS2 = norm(Rij);

            rm_6 = 1.0/(SS2*SS2*SS2);
            rm_12 = rm_6*rm_6;

            pot = LJP12[spec1][spec2]*rm_12 - LJP6[spec1][spec2]*rm_6;
            lj_Pot += pot;
            dljpot = LJF12[spec1][spec2]*rm_12 - LJF6[spec1][spec2]*rm_6;
            forc = Rij*(dljpot/SS2);

            FS[o1] = FS[o1] + forc;
            FS[o2] = FS[o2] - forc;
        }

    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

