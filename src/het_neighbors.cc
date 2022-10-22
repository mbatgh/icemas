#include "vect_ops.h"
#include "structures.h"
#include "icemas.h" 


/*******************************************************************/
Real    het_neighbors(void) {
/*******************************************************************/

/*  this subroutine calculates the neighborlist, i.e. a pair-list of 
    all the particle-pairs, the distance of which is smaller than
    the cutoff-radius (for systems with both LJ- and electrostatic
    interactions) */

    register int i,na,neigh,IA;
    Real dx2,cpu_time;
    VEKTOR Rij;

    clock_t  t1,t2;
    t1 = clock();

    nnup++;
    GETNGB = 0;
    NR_LJN = NR_LJIN;
    NR_ELN = NR_ELIN;
    NR_AN  = NR_AIN; 

    for(na=0;na<Nr_Stot-1;na++) {
        for(neigh=nextidx[na];neigh<Nr_Stot;neigh++) {

            Rij = RS[na]-RS[neigh];
            convolute(Rij,adLh,L);
            dx2 = norm(Rij);      
            
            if(dx2<far_ctof_2) {
            
                IA = ww[na]+ww[neigh];

                switch(IA) {

                    case 2: 
                    case 3:
                        ljnb_i[NR_LJN] = na;
                        ljnb_j[NR_LJN] = neigh;
                        NR_LJN++;
                        break;
                    case 4: 
                    case 6:
                        anb_i[NR_AN] = na;
                        anb_j[NR_AN] = neigh;
                        NR_AN++;
                        break;
                    case 8:
                        elnb_i[NR_ELN] = na;
                        elnb_j[NR_ELN] = neigh;
                        NR_ELN++;
                }
            }
        }
    }

/*
VdW = IA-1;
VdW = 1/VdW;
Cou = IA-2;
o   = (Real)Cou/((Real)Cou+1.0);
Cou = (int)(o+1.0);
OUT = abs(IA-5)+1;
o   = ((Real)OUT-1.0)/(Real)OUT;
OUT = (int)(o+1.0);
LJOUT = OUT*(IA/8);
ZWEI = IA/2;
DREI  = IA/3;
VIER  = IA/4;
FUENF  = IA/5;
SECHS  = IA/6;
ACHT = IA/8;
*/

/* zero distance accumulators .................................... */

    if(NL_SAI)
        for(i=0;i<Nr_Stot;i++)
            displacement[i] = nulvek;


    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

