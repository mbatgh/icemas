#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"




void get_a_dist(Real a) {

    int idx;

    idx = (int)floor(200.0*a/PI);
    if(idx>199) idx = 199;
    a_dist[idx]++;
}


void get_dha_dist(Real dha) {

    int idx;

    idx = 100+(int)floor(100.0*dha/PI);
    dha_dist[idx]++;
}


void get_idha_dist(Real dha) {

    int idx;

    idx = 100+(int)floor(100.0*dha/PI);
    idha_dist[idx]++;
}


