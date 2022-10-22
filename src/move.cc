#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

Real    move(void);
void    shake(void);


/* ================================================================= */
Real    move(void) {
/* ================================================================= */

    int         i,j,k;
    Real        NHTfac1,NHTfac2,cpu_time;
    VEKTOR      v;

    KinA = 0.0;
    clock_t t1,t2;
    t1 = clock();

/*_____ Stoermer-Verlet-integration .......................... _____ */

   
    if(it<end_thst) {

        if(NHT) {

/*_____ apply NHT ............................................ _____ */
        
            NHTfac1 = zeta * dt * 0.5;
            NHTfac2 = 1.0 / (1.0 + zeta * dt * 0.5);

            for(j=0;j<Nr_Stot;j++) {

                RSn[j] = RS[j] * 2.0 - RSo[j];
                RSn[j] = RSn[j] + FS[j] * dtqdm[j];
                RSn[j] = (RSn[j] + RSo[j] * NHTfac1) * NHTfac2;
            }    
        }

/*_____ or kinetostat ........................................ _____ */

        else {

            for(j=0;j<Nr_Stot;j++) {

                RSn[j] = RS[j] * 2.0 - RSo[j];
                RSn[j] = RSn[j] + FS[j] * dtqdm[j];
                v = RSn[j]-RSo[j];
                KinA += mass[j]*ad8dtq*norm(v);
            }
            
            scale_new_tmp(KinA,KinAsoll);
            KinA = 0.0;
        }
    }

    else {
    
/*_____ or no thermostat ..................................... _____ */

        for(j=0;j<Nr_Stot;j++) {

            RSn[j] = RS[j] * 2.0 - RSo[j];
            RSn[j] = RSn[j] + FS[j] * dtqdm[j];
        }
    }

/*_____ call SHAKE - routine ................................. _____ */

    shake();

/*_____	calc. kinetic Energies ............................... _____ */

    for(j=0;j<Nr_Stot;j++) {
        v = RSn[j]-RSo[j];
        KinA += mass[j]*ad8dtq*norm(v);
    }

/*_____	apply periodic boundaries ............................ _____ */

    for(j=0;j<Nr_Stot;j++) {
        v = RSn[j];
        convolute(RSn[j],adLh,L);
        RSo[j] = RS[j] +(RSn[j]-v);
        RS[j] = RSn[j];
    }

/*_____	adjust thermostat parameter ........................... _____*/

    if(it<end_thst&&NHT) {
        zetaold = zeta;
        zeta = (KinA/KinAsoll-1.0)*(dt/tausq)+zetaold;
    }


    t2 = clock();
    cpu_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

    return(cpu_time);
}


/* =================================================================== */
/* =                            SHAKE                                = */
/* =================================================================== */

void shake() {

    int         j,z,idx1,idx2,n,MOVING[NR_CMAX],MOVED[NR_CMAX],DONE;
    Real        gamma,nbn,diffsq;
    VEKTOR      bn,db,b[NR_CMAX];

    n = 0;

    for(z=0;z<Nr_Btot;z++) {
        idx1 = zidx1[z];
        idx2 = zidx2[z];
        b[z] = RS[idx2]-RS[idx1];
        convolute(b[z],adL34,L);
        MOVING[idx1] = MOVING[idx2] = 0;
        MOVED[idx1] = MOVED[idx2] = 1;
    }

    do {

        DONE = 1;
        
        for(j=0;j<Nr_Btot;j++) {
            idx1 = zidx1[j];
            idx2 = zidx2[j];
        
            if(MOVED[idx1] || MOVED[idx2]) {

                bn = RSn[idx2]-RSn[idx1];
                convolute(bn,adL34,L);
                nbn = norm(bn);
                diffsq = sq_bnd[j]-nbn;

                if(fabs(diffsq)>sq_bnd[j]*sq_crit) {

                    gamma = diffsq/(((rm[idx1]+rm[idx2])*2.0)*(bn*b[j]));
                    db = b[j]*gamma;
                    RSn[idx2] = RSn[idx2]+db*rm[idx2];
                    RSn[idx1] = RSn[idx1]-db*rm[idx1];
                    DONE = 0;
                    MOVING[idx1] = MOVING[idx2] = 1;
                }
            }    
        }
        
        if(++n>100) {
            printf("stuck in shake ...\n");
            write_means(1.0,1.0,1.0,1.0,1.0,1.0,1.0);
            write_r(0,Nr_Stot,"r_kill");
            exit(1);
        }
        
        for(j=0;j<Nr_Btot;j++) {
            MOVED[j] = MOVING[j];
            MOVING[j] = 0;
        }

    } while(!DONE);

    Nr_SHAKE += n;
}

