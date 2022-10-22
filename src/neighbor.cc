#include "vect_ops.h"
#include "structures.h"
#include "icemas.h" 

Real    neighbor(void);

/*******************************************************************/
Real    neighbor(void) {
/*******************************************************************/

/*  this subroutine calculates the neighborlist, i.e. a pair-list of 
    all the particle-pairs, the distance of which is smaller than
    the cutoff-radius (just for ensembles of solely LJ-interactions) */

    int i,na,neigh,IN,OUT;
    Real dx2,cpu_time,o;
    VEKTOR Rij;

    clock_t  t1,t2;
    t1 = clock();

    nnup++;
    GETNGB = 0;
    NON = Nr_IN;

    for(na=0;na<Nr_Stot-1;na++) {
        for(neigh=nextidx[na];neigh<Nr_Stot;neigh++) {

            Rij = RS[na]-RS[neigh];
            convolute(Rij,adLh,L);
            dx2 = norm(Rij);      
            IN = int(dx2/far_ctof_2) + 1;
            IN = 1/IN;
            OUT = abs(ww[na]+ww[neigh]-5)+1;
            o = (Real(OUT)-1.0)/Real(OUT);
            OUT = int(ceil(o));
            nb_i[NON] = na;
            nb_j[NON] = neigh;
            NON += IN*OUT;
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

