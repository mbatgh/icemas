#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

void    make_k(Real);

/****************************************************************/
void    make_k() {
/****************************************************************/

    int l,m,n,ksq,nk,nmin,mmin;
    Real rksq,zwapi,rl,rm,rn;

    zwapi = 2.0*PI;

    nk = 0;
    mmin=0;
    nmin=1;

    for(l=0;l<=K_max;l++) {
        rl = zwapi*(Real)l/L;
        for(m=mmin;m<=K_max;m++) {
            rm = zwapi*(Real)m/L; 
            for(n=nmin;n<=K_max;n++) {
                rn = zwapi*(Real)n/L; 

                ksq = l*l+m*m+n*n;
                rksq = rl*rl+rm*rm+rn*rn;

                if(ksq<Ksq_max&&ksq!=0) {
                    nk++;
                    ewfac[nk] = zwapi/(L*L*L)*exp(-0.25*rksq/eta/eta)/rksq;
                }
            }
            nmin = -K_max;
        }
        mmin = -K_max;
    }
    
    if(nk>=KTOT_MAX)
        { printf("make_k needs more memory ...\n"); exit(1); }
    else printf("made %d coefficients for ewald-sum ...\n", nk);
}

