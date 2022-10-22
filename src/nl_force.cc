#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/* =================================================================== */
Real    nl_force() {
/* =================================================================== */

    int     n,a,b,i,j,k,l,o1,o2,firstone,lastone,spec1,spec2,dif;
    Real    SS2,pot,dljpot,rm_12,rm_6,cpu_time;
    VEKTOR  Rij,forc,dx,dR,dr;
    clock_t t1,t2;

    t1 = clock();

    Virial = 0.0;
    lj_Pot = 0.0;


        for(n=0;n<NR_LJN;n++) {

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
            Virial += dljpot;

            forc = Rij*(dljpot/SS2);

            FS[o1] = FS[o1] + forc;
            FS[o2] = FS[o2] - forc;
        }


    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

