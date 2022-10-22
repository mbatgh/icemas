#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/* =================================================================== */
Real    r_force() {
/* =================================================================== */

/*  this subroutine calculates potential and forces (Lennard-Jones-
    as well as the real part of the electrostatic interactions), using
    no neighborlist */

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

/*  consider each paticle ... */

    for(i=0;i<Nr_Stot-1;i++) {

        ww1 = ww[i];
        spec1 = ljt[i];

/*  and each possible neighbor (excluding the intramolecular neighbors) */

        for(j=nextidx[i];j<Nr_Stot;j++) {

            ww2 = ww[j];
            spec2 = ljt[j];

            Rij = RS[i]-RS[j];
            convolute(Rij,adLh,L);
            SS2 = norm(Rij);       

            if(SS2<cutoff_2) {    

/*  if the distance between the particles is smaller than the cutoff-radius
    calculate potential and forces depending on the kind of interaction
    (LJ or electrostatic or both) */

                switch(ww1+ww2) {

            case 6:
            case 8:

                SS = sqrt(SS2);

                qq = elch[i] * elch[j];
                pot = erfc(eta*SS)*qq/SS;

                dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
                c_Potr += pot;
                Virial += dcpot;
                forc = Rij*(dcpot/SS2);

                FS[i] = FS[i] + forc;
                FS[j] = FS[j] - forc;

                break;
                
            case 4:

                rm_6 = 1.0/(SS2*SS2*SS2);
                rm_12 = rm_6*rm_6;

                pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
                dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;
                lj_Pot += pot;
                Virial += dljpot;
                SS = sqrt(SS2);

                qq = elch[i] * elch[j];
                pot = erfc(eta*SS)*qq/SS;

                dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
                c_Potr += pot;
                Virial += dcpot;
                forc = Rij*((dljpot+dcpot)/SS2);

                FS[i] = FS[i] + forc;
                FS[j] = FS[j] - forc;

                break;

            case 2:
            case 3:

                rm_6 = 1.0/(SS2*SS2*SS2);
                rm_12 = rm_6*rm_6;

                pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
                dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;

                lj_Pot += pot;
                Virial += dljpot;

                forc = Rij*(dljpot/SS2);
          
                FS[i] = FS[i] + forc;
                FS[j] = FS[j] - forc;
                
                }
            }
        }
    }

/* nonbonded interactions of intramolecular pairs (LJ-IA) ...........  */

    for(n=0;n<NR_LJIN;n++) {

        i = ljnb_i[n];
        j = ljnb_j[n];
        spec1 = ljt[i];
        spec2 = ljt[j];

        Rij = RS[i]-RS[j];
        convolute(Rij,adLh,L);
        SS2 = norm(Rij);       

        rm_6 = 1.0/(SS2*SS2*SS2);
        rm_12 = rm_6*rm_6;
        pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
        dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;

        lj_Pot += pot;
        forc = Rij*(dljpot/SS2);
        FS[i] = FS[i] + forc;
        FS[j] = FS[j] - forc;
    }

/* nonbonded interactions of intramolecular pairs (LJ & electrostatic IA) */

    for(n=0;n<NR_AIN;n++) {

        i = anb_i[n];
        j = anb_j[n];
        spec1 = ljt[i];
        spec2 = ljt[j];
        Rij = RS[i]-RS[j];
        convolute(Rij,adLh,L);
        SS2 = norm(Rij);
        rm_6 = 1.0/(SS2*SS2*SS2);
        rm_12 = rm_6*rm_6;
        pot = LJP12[spec1][spec2] * rm_12 - LJP6[spec1][spec2] * rm_6;
        dljpot = LJF12[spec1][spec2] * rm_12 - LJF6[spec1][spec2] * rm_6;
        lj_Pot += pot;
        Virial += dljpot;

        SS = sqrt(SS2);
        qq = elch[i] * elch[j];
        pot = erfc(eta*SS)*qq/SS;
        dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
        c_Potr += pot;
        Virial += dcpot;
        forc = Rij*((dljpot+dcpot)/SS2);

        FS[i] = FS[i] + forc;
        FS[j] = FS[j] - forc;
    }

/* nonbonded interactions of intramolecular pairs (electrostatic IA) */

    for(n=0;n<NR_ELIN;n++) {
        
        i = elnb_i[n];
        j = elnb_j[n];
        Rij = RS[i]-RS[j];
        convolute(Rij,adLh,L);
        SS2 = norm(Rij);
        SS = sqrt(SS2);
        qq = elch[i] * elch[j];
        pot = erfc(eta*SS)*qq/SS;
        dcpot = pot + qq*pifac*exp(-(eta*eta*SS2));
        c_Potr += pot;
        forc = Rij*(dcpot/SS2);
        FS[i] = FS[i] + forc;
        FS[j] = FS[j] - forc;
    }


    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

