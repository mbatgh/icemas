#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


Real update_dist(void);


/****************************************************************/
Real update_dist(void) {
/****************************************************************/

/*  when using a neighbor-list with self-adjusting interval, 
    this subroutine calculates the largest distance one particle 
    travelled since the last neighborlist-update, in order to decide,
    whether the neighbor-list should be updated again */

    int     i,j;
    Real    maxdis,dis,stop,total_time;
    clock_t   t1,t2;

    t1 = clock();

    maxdis = 0.0;
    stop = far_ctof-cutoff;

    for(i=0;i<Nr_Stot;i++) {
        displacement[i] = displacement[i] + (RS[i]-RSo[i]);
        dis = sqrt(norm(displacement[i]));
        if(dis+maxdis>stop) {
            GETNGB = 1;
            t2 = clock();
            total_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;
            return(total_time);
        }
        if(dis>maxdis) maxdis = dis;
    }

    t2 = clock();
    total_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(total_time);
}

