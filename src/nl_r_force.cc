#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/* =================================================================== */
Real nl_r_force() {
/* =================================================================== */

/*  this subroutine calculates potential and forces (Lennard-Jones-
    as well as the real part of the electrostatic interactions), using
    a seperate neighborlist for each kind of interaction */

    int n,i,j,l,o1,o2,ww1,ww2,lastone,
        molecule,spec1,spec2,firstone;

    Real    SS2,SS,
            pot,
            dcpot,dljpot,
            qq,pifac,
            rm_12,rm_6,
            cpu_time;

    VEKTOR Rij,forc;

    clock_t  t1,t2;
    t1 = clock();

    c_Potr = 0.0;
    lj_Pot = 0.0;
    Virial = 0.0;

    pifac = 2.0*eta/sqrt(PI);

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

        pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
        dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;

        lj_Pot += pot;
        Virial += dljpot;

        forc = Rij*(dljpot/SS2);
          
        FS[o1] = FS[o1] + forc;
        FS[o2] = FS[o2] - forc;
    }

    for(n=0;n<NR_ELN;n++) {   

        o1 = elnb_i[n];
        o2 = elnb_j[n];

        Rij = RS[o1]-RS[o2];
        convolute(Rij,adLh,L);

        SS2 = norm(Rij);
        SS = sqrt(SS2);

        qq = elch[o1] * elch[o2];
        pot = erfc(eta*SS)*qq/SS;
        dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
        c_Potr += pot;
        Virial += dcpot;

        forc = Rij*(dcpot/SS2);

        FS[o1] = FS[o1] + forc;
        FS[o2] = FS[o2] - forc;
    }

    for(n=0;n<NR_AN;n++) {

        o1 = anb_i[n];
        o2 = anb_j[n];
        spec1 = ljt[o1];
        spec2 = ljt[o2];

        Rij = RS[o1]-RS[o2];
        convolute(Rij,adLh,L);
        SS2 = norm(Rij);
        rm_6 = 1.0/(SS2*SS2*SS2);
        rm_12 = rm_6*rm_6;
        pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
        dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;
        lj_Pot += pot;
        Virial += dljpot;

        SS = sqrt(SS2);
        qq = elch[o1] * elch[o2];
        pot = erfc(eta*SS)*qq/SS;
        dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
        c_Potr += pot;
        Virial += dcpot;
        forc = Rij*((dljpot+dcpot)/SS2);

        FS[o1] = FS[o1] + forc;
        FS[o2] = FS[o2] - forc;
    }

    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

