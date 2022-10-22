#include <complex>
#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


Real kforce(void);


/****************************************************************/
Real k_force(void) {
/****************************************************************/

    Real zwapi,pifac,zwapidL,cpu_time,
         xc,yc,zc,
         rl,rm,rn,
         qforce,qpec;
    
    VEKTOR kv;

    int nk,mmin,i,j,nmin,kk,k,l,m,n;

    complex<double> el[NR_SMAX][K_max+1];
    complex<double> em[NR_SMAX][2*K_max+1];
    complex<double> en[NR_SMAX][2*K_max+1];
    complex<double> expikr[NR_SMAX],sum;
    register complex<double> dummy, ii(0.0,1.0), cx_zero(0.0,0.0);

    clock_t  t1,t2;
    t1 = clock();

    pifac = 2.0*eta/sqrt(PI);
    zwapi=2.0*PI;
    zwapidL = 2.0*PI/L;

    for(k=0;k<Nr_Charged;k++) {
        i = el_idx[k];
        el[i][0]   = 1.0;
        em[i][K_max] = 1.0;
        en[i][K_max] = 1.0;
        xc=zwapidL*(RS[i].x);
        yc=zwapidL*(RS[i].y);
        zc=zwapidL*(RS[i].z);
        el[i][1]     = cos(xc)+ii*sin(xc);
        em[i][K_max+1] = cos(yc)+ii*sin(yc);
        en[i][K_max+1] = cos(zc)+ii*sin(zc);
    }

    for(k=0;k<Nr_Charged;k++) {
        i = el_idx[k];
        for(j=2;j<=K_max;j++) {
            el[i][j]     = el[i][j-1]    *el[i][1];
            em[i][K_max+j] = em[i][K_max+j-1]*em[i][K_max+1];
            en[i][K_max+j] = en[i][K_max+j-1]*en[i][K_max+1];
        }
        for(j=-K_max;j<0;j++) {
            em[i][K_max+j] = conj(em[i][K_max-j]);
            en[i][K_max+j] = conj(en[i][K_max-j]);
        }
    }

    qpec=0.0;
    nk=0;
    mmin=0;
    nmin=1;

    for(l=0;l<=K_max;l++) {
        kv.x = zwapidL*(Real)l;
        for(m=mmin;m<=K_max;m++) {
            kv.y = zwapidL*(Real)m;
            for(n=nmin;n<=K_max;n++) {
                kv.z = zwapidL*(Real)n;

                kk=l*l+m*m+n*n;

                if(kk<Ksq_max) {

                    nk++;

                    for(k=0;k<Nr_Charged;k++) {
                        i = el_idx[k];
                        expikr[i] = el[i][l]*em[i][K_max+m]*en[i][K_max+n];
                    }
 
                    sum=cx_zero;

                    for(k=0;k<Nr_Charged;k++) {
                        i = el_idx[k];
                        sum += elch[i]*expikr[i];
                    }

                    dummy = conj(sum)*sum;
                    qpec += ewfac[nk]*dummy.real();

                    for(k=0;k<Nr_Charged;k++) {
                        i = el_idx[k];
                        dummy = sum*conj(expikr[i]);

                        qforce = -4.0*ewfac[nk]*elch[i]*dummy.imag();
                        FS[i] = FS[i] + kv * qforce;
                    }
                }
            }
            nmin = -K_max;
        }
        mmin = -K_max;
    }

    c_Potk = 2.0*qpec;
    Virial += -2.0*qpec;

    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

