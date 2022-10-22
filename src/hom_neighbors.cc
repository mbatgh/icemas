#include "vect_ops.h"
#include "structures.h"
#include "icemas.h" 


/*******************************************************************/
Real    hom_neighbors(void) {
/*******************************************************************/

/*  this subroutine calculates the neighborlist, i.e. a pair-list of 
    all the particle-pairs, the distance of which is smaller than
    the cutoff-radius (just for ensembles of solely LJ-interactions) */

    int i,na,neigh,IN;
    Real dx2,cpu_time;
    VEKTOR Rij;

    clock_t  t1,t2;
    t1 = clock();

    nnup++;
    GETNGB = 0;
    NR_LJN = NR_LJIN;

    for(na=0;na<Nr_Stot-1;na++) {
        for(neigh=nextidx[na];neigh<Nr_Stot;neigh++) {

            Rij = RS[na]-RS[neigh];
            convolute(Rij,adLh,L);
            dx2 = norm(Rij);      
            IN = int(dx2/far_ctof_2) + 1;
            IN = 1/IN;
            ljnb_i[NR_LJN] = na;
            ljnb_j[NR_LJN] = neigh;
            NR_LJN += IN;
        }
    }

/* zero distance accumulators .................................... */

    if(NL_SAI)
        for(i=0;i<Nr_Stot;i++)
            displacement[i] = nulvek;


    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

