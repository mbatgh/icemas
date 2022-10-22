#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


Real angle_force(void);


Real angle_force(void) {

    clock_t t1,t2;
    int     i,j,nm,ns,na,na0,ns0,aidx,idx1,idx2,idx3;
    Real    cpu_time,
            fc,alf,alf_0, 
            a,b,c,ab,calf,
            pot, dpot, potsum,
            sqrtb2,sinalf;
    VEKTOR d12,d23,d13,f1a1,f1a2,f1a3,
            dum1,dum2;

    t1 = clock();

    a_Pot = 0.0;
    na0 = 0;
    ns0 = 0;

    for(ns=0;ns<Nr_Spec;ns++) {
        for(nm=0;nm<Nr_Molecules[ns];nm++) {
            for(na=0;na<Nr_Angles[ns];na++)  {  

        aidx = na0+na;

        idx1 = ns0 + Angle[aidx].idx1;
        idx2 = ns0 + Angle[aidx].idx2;
        idx3 = ns0 + Angle[aidx].idx3;
        alf_0 = Angle[aidx].a_eqi;
        fc = Angle[aidx].force;

        d12 = RS[idx2]-RS[idx1];
        convolute(d12,adL34,L);
        d23 = RS[idx3]-RS[idx2];
        convolute(d23,adL34,L);
        d13 = RS[idx3]-RS[idx1];
        convolute(d13,adL34,L);

        a = d12*d12;
        b = d23*d23;
        c = d23*d12;

        ab = sqrt(a*b);
        calf = c/ab;

        if(calf>1.0) calf = 1.0;
        if(calf<-1.0) calf = -1.0;
        alf = PI-acos(calf);

        pot = fc*(alf-alf_0)*(alf-alf_0);
        dpot = -2.0*fc*(alf-alf_0);

/* 
    if the equilibrium-bond-angle is near PI i.e. linear,
    an alternative potential may be more approbriate 

        pot = fc*(1.0+cos(alf));
        dpot = fc*sin(alf);
*/
        a_Pot += pot;

        dum1 = d12*(c/a)-d23;
        f1a1 = dum1*(dpot/ab);
        dum1 = d12-d23*(c/b); 
        f1a3 = dum1*(dpot/ab);
        dum1 = f1a1+f1a3;
        f1a2 = nulvek-dum1;

        FS[idx1] = FS[idx1] + f1a1;
        FS[idx2] = FS[idx2] + f1a2;
        FS[idx3] = FS[idx3] + f1a3;

        if((ns==SP_ADIST)&&(na==NR_ADIST))
            if((fmod(it,angdf_step)==0.0)&&(it>begin_avgs))
                get_a_dist(alf);

            }
            ns0 += Nr_Sites[ns];
        }
        na0 += Nr_Angles[ns];
    }

    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

