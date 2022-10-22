#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


void tdf() {

/* this subroutine is called after each integration-step to update
   the means for the thermodynamic functions, to calculate eventually 
   other physical data and write the data on disk */

    int  a,b,i,j,k,l,tot,o1,o2,idx,root,i0_1,in_1,i0_2,in_2,dp1,dp2,p2,p1; 
    Real EkinA,Epotlj,Epotcr,Epotck,adNA,SS,SS2,
         EAng,EdhAng,EidhAng,T,P,divi,sum,iVir,Vir;
    VEKTOR dr,dum,COM;
    FILE* file;

    tot = 0;
    EAng = 0.0;
    EdhAng = 0.0;
    EidhAng = 0.0;
    iVir = 0.0;
    
/* calculation of the intramolecular virial ................... */

    tot = 0;

    for(i=0;i<Nr_Spec;i++) {
        for(j=0;j<Nr_Molecules[i];j++) {
            COM = cent_of_mass(RS,tot,tot+Nr_Sites[i]-1);
            for(k=0;k<Nr_Sites[i];k++) {
                idx = tot+k;
                dr = RS[idx] - COM;
                convolute(dr,adL34,L);
                iVir -= dr * FS[idx];
            }
            tot += Nr_Sites[i];
        }
    }

/* calculate thermodynamics ................................... */

    T       = 2000.0 * KinA / (Real)Nr_Ftot / Rgas;
    Vir     = iVir+Virial;
    P       = 1.0e-06*((Real)Nr_Mtot*kB*T+1000.0*Vir/3.0/NLo)/L/L/L/1.0e-27;

    divi    = 1.0/(Real)Nr_Mtot;

    Vir     = divi * Vir;
    EkinA   = divi * KinA;
    Epotlj  = divi * lj_Pot;
    Epotcr  = divi * c_Potr;
    Epotck  = divi * c_Potk;

    if(ANGLES)    EAng   = divi * a_Pot;
    if(DHANGLES)  EdhAng = divi * dha_Pot;
    if(IDHANGLES) EidhAng = divi * idha_Pot;

    Etot = Epotlj+Epotcr+Epotck+EkinA+EAng+EdhAng+EidhAng;

/* when turning off the thermostat scale the temperature for a last time, 
   for evaluating any energy-drift remember the total energy ..... */

    if(it==end_thst) scale_tmp(KinA,KinAsoll);
    if(it==end_thst+1) Etot0 = Etot;

/* accumulate energies ........................................... */

    if(it>begin_avgs) {
        S_T += T;
        S_V += Vir;
        S_P += P;
        S_Epotlj += Epotlj;
        S_Epotcr += Epotcr;
        S_Epotck += Epotck;
        S_KinA += EkinA;
        S_APot += EAng;
        S_DHAPot += EdhAng;
        S_IDHAPot += EidhAng;
    }


    if(pcf_step>0&&fmod((double)it,(double)pcf_step)==0.0) {

/* calculate pair-distribution-functions ........................ */

        for(i=0;i<Nr_pcfs;i++) {

            i0_1 = pcf0_1[i];
            in_1 = pcfn_1[i];
            i0_2 = pcf0_2[i];
            in_2 = pcfn_2[i];
            dp1 = pcf_d1[i];
            dp2 = pcf_d2[i];

            if(i0_1!=i0_2) {
                for(p1=i0_1;p1<=in_1;p1+=dp1) {
                    for(p2=i0_2;p2<=in_2;p2+=dp2) {
                        dr = RS[p1]-RS[p2];
                        convolute(dr,adLh,L);
                        SS = betr(dr);
                        if(SS<Lh) {
                            idx = (int)floor((Real)NR_PCF_CHANS*SS/Lh);
                            pcf[i][idx]++;
                        }
                    }
                }
            } else {
                for(p1=i0_1;p1<in_1;p1+=dp1) {
                    for(p2=p1+dp2;p2<=in_2;p2+=dp2) {
                        dr = RS[p1]-RS[p2];
                        convolute(dr,adLh,L);
                        SS = betr(dr);
                        if(SS<Lh) {
                            idx = (int)floor((Real)NR_PCF_CHANS*(SS/Lh));
                            pcf[i][idx]+=2;
                        }
                    }
                }
            }
        }
    }


    if(sep_intv>0&&fmod((double)it,(double)sep_intv)==0.0) {

/* calculate the phase-separation ............................... */

        char psfile[200];
        char* filename;
        Real p,a,m;

        psfile[0] = 'p';
        psfile[1] = 'h';
        psfile[2] = 's';
        psfile[3] = 'p';
        psfile[4] = '.';
        psfile[5] = '\0';

        phase_sep(&p,&a,&m);

        filename = strcat(psfile,System_name);
        if((file = fopen(filename,"a+"))==NULL)
            printf("could not open %s ...\n", filename);
        fprintf(file,"%le %le %le %le\n",(Real)it*dt,p,a,m);  
        fclose(file);
    }


/* write thermodynamic data in file ............................. */

    if(fmod(it,write_therm_step) == 0.0) {

        char data[200];
        char* filename;

        data[0] = 'd';
        data[1] = 'a';
        data[2] = 't';
        data[3] = 'a';
        data[4] = '.';
        data[5] = '\0';

        filename = strcat(data,System_name);

        if((file = fopen(filename,"a+"))==NULL)
            printf("could not open %s ...\n", filename);

        fprintf(file,"%le %le %le %le %le %le %le %le %le %le %le %le\n", 
                (Real)it*dt,T,P,Vir,EkinA,Epotlj,Epotcr,Epotck,
                EAng,EdhAng,EidhAng,Etot);
        fclose(file);
    }
}

