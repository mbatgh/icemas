#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

/*********************************************************************/
Real torsion_force(void) {
/*********************************************************************/

/*  this subroutine calculates the potential of each given
    dihedral angle and the resulting forces */

    int     i,aidx,nd,ns,nm,nd0,ns0,idx1,idx2,idx3,idx4;
    clock_t dt1,dt2;
    Real    cpu_time,cosfi,fi_0,fi,signum,fc,alf,
            a2,b2,cd1,cd2,sqrtb2,period,cosfi_0,sinfi,
            t1,t2,t3,t4,t5,t6,
            c11,c12,c13,c22,c23,c33,
            pot,dpot,dur;
    VEKTOR  fd1,fd2,fd3,fd4,d10,d21,d20,d30,d31,d32,
            dum1,dum2,dum3,sum;

    dt1 = clock();

    dha_Pot = 0.0;
    nd0 = 0;
    ns0 = 0;

/* calculation of dihedral angle potential and forces */

    for(ns=0;ns<Nr_Spec;ns++) {
        for(nm=0;nm<Nr_Molecules[ns];nm++) {
            for(nd=0;nd<Nr_Dieder[ns];nd++) {

                aidx = nd0+nd;
                idx1 = ns0 + Dieder[aidx].idx1;
                idx2 = ns0 + Dieder[aidx].idx2;
                idx3 = ns0 + Dieder[aidx].idx3;
                idx4 = ns0 + Dieder[aidx].idx4;
                fi_0 = Dieder[aidx].d_eqi;
                period = (Real)Dieder[aidx].period;
                fc = Dieder[aidx].force;




/*
RS[idx1] = ((VEKTOR){0.0,1.0,0.0});
RS[idx2] = ((VEKTOR){1.0,0.0,0.0});
RS[idx3] = ((VEKTOR){2.0,0.0,0.0});
RS[idx4] = ((VEKTOR){3.0,1.0,0.0});
FILE* fil;
fil = fopen("anc","w");
alf = 0.0;

for(i=0;i<200;i++) {
RS[idx4].x = 3.0;
RS[idx4].y = cos(alf);
RS[idx4].z = sin(alf);
alf += PI/100.0;
*/

                d10 = RS[idx2]-RS[idx1];
                convolute(d10,adLh,L);
                d21 = RS[idx3]-RS[idx2];
                convolute(d21,adLh,L);
                d32 = RS[idx4]-RS[idx3];
                convolute(d32,adLh,L);
                d31 = RS[idx4]-RS[idx2];
                convolute(d31,adLh,L);
                d20 = RS[idx3]-RS[idx1];
                convolute(d20,adLh,L);
                d30 = RS[idx4]-RS[idx1];
                convolute(d30,adLh,L);

                c11 = d10*d10;
                c12 = d10*d21;
                c13 = d10*d32;
                c22 = d21*d21;
                c23 = d21*d32;
                c33 = d32*d32;

                a2 = c13*c22-c12*c23;
                b2 = (c11*c22-c12*c12)*(c22*c33-c23*c23);

                sqrtb2 = sqrt(b2);
                cd2 = sqrt(c22*c33);

                t1 = c13*c22-c12*c23;
                t2 = c11*c23-c12*c13;
                t3 = c12*c12-c11*c22;
                t4 = c22*c33-c23*c23;
                t5 = c13*c23-c12*c33;
                t6 = -t1;

                cosfi = a2/sqrtb2;
                if(cosfi<-1.0) cosfi = -1.0;
                if(cosfi>1.0) cosfi = 1.0;

                dum1 = d10%d21;
                dum2 = d21%d32;
                dum3 = dum1%dum2;
                signum = ((d21*dum3)>0.0) ? -1.0 : 1.0;

                fi = signum*acos(cosfi);

                if(fi>PI) fi -= 2.0*PI;
                if(fi<-PI) fi += 2.0*PI;


/**************************************************************************
  the following 2 dihedral angle potentials were introduced by toxvaerd 
  and ryckaert resp.
  they yield definitly better results for simple alcanes than the
  harmonic potential does
  
        pot = Rgas*(1.03776+2.42607*cosfi+0.08164*cosfi*cosfi-
                3.12946*cosfi*cosfi*cosfi-0.16328*cosfi*cosfi*cosfi*cosfi-
                0.25273*cosfi*cosfi*cosfi*cosfi*cosfi);

        dpot = -Rgas*(2.42607+0.16328*cosfi-9.38838*cosfi*cosfi-
                0.65312*cosfi*cosfi*cosfi-
                1.26365*cosfi*cosfi*cosfi*cosfi);

printf("tpot:  %le, ",pot);
printf("tdpoz: %le\n",dpot);

***************************************************************************

        pot = Rgas*
(1.116+1.462*cosfi-1.578*cosfi*cosfi-0.368*cosfi*cosfi*cosfi+3.156*cosfi*cosfi*cosfi*cosfi-3.788*cosfi*cosfi*cosfi*cosfi*cosfi);

        dpot = -Rgas*(1.462
                -3.156*cosfi
                -1.104*cosfi*cosfi
                +12.624*cosfi*cosfi*cosfi
                -18.94*cosfi*cosfi*cosfi*cosfi);

printf("rpot:  %le, ",pot);
printf("rdpot: %le\n",dpot);

**************************************************************************/

                sinfi = sin(fi);
                sinfi = (fabs(sinfi)<1.0e-08) ? 1.0e-08 : sinfi;
                pot = fc*(1.0-cos(period*(fi-fi_0)));
/*                pot = fc*(1.0+cos(period*fi-fi_0));*/
                dpot = fc*period*sin(period*fi-fi_0)/sinfi;

                dha_Pot += pot;

                  dum1 = d10*t1;
                  dum2 = dum1+d21*t2;
                  dum3 = dum2+d32*t3;
                fd1 = dum3*(c22/((c11*c22-c12*c12)*sqrtb2));
                  dum1 = d10*t4;
                  dum2 = dum1+d21*t5; 
                  dum3 = dum2+d32*t6;
                fd4 = dum3*(c22/((c22*c33-c23*c23)*sqrtb2));
                  dum1 = fd1*(1.0+c12/c22);
                fd2 = fd4*(c23/c22) - dum1;
                  dum1 = fd1*(c12/c22);
                fd3 = dum1 - (fd4*(1.0+c23/c22));


/*
printf("pot = %le, ", pot);
printf("dpot = %le\n", dpot);
printf("forces: %le %le %le %le\n",betr(fd1*dpot),betr(fd2*dpot),betr(fd3*dpot),betr(fd4*dpot));
printf("fi = %le, fi_0 = %le", fi, fi_0);
presskey(" ...");
FS[idx1] = nulvek;
FS[idx2] = nulvek;
FS[idx3] = nulvek;
FS[idx4] = nulvek;
*/

                FS[idx1] = FS[idx1] + fd1*dpot;
                FS[idx2] = FS[idx2] + fd2*dpot;
                FS[idx3] = FS[idx3] + fd3*dpot;
                FS[idx4] = FS[idx4] + fd4*dpot;




/*
fprintf(fil,"%lf %lf %lf %lf %lf %lf %lf %lf ",alf,fi,pot,dpot,
             betr(FS[idx1]),betr(FS[idx2]),betr(FS[idx3]),betr(FS[idx4]));
dum1 = RS[idx1]-RS[idx2];
convolute(dum1,adLh,L);
dur = FS[idx1]*dum1;
fprintf(fil,"%lf ", dur);

dum1 = RS[idx4]-RS[idx3];
convolute(dum1,adLh,L);
dur = FS[idx4]*dum1;
fprintf(fil,"%lf ", dur);

dum1 = RS[idx3]-RS[idx2];
convolute(dum1,adLh,L);
dur = FS[idx2]*dum1;
fprintf(fil,"%lf ", dur);
dur = -(FS[idx3]*dum1);
fprintf(fil,"%lf\n", dur);

}
fclose(fil);
presskey(" ....");
*/

//printf("tor fi = %le, pot = %le, dpot = %le",fi,pot,dpot);
//presskey(" ...");

                if((ns==SP_DHADIST)&&(nd==NR_DHADIST)&&
                (fmod(it,dhadf_step)==0.0)&&(it>begin_avgs)) 
                    get_dha_dist(fi);
            }
            ns0 += Nr_Sites[ns];
        }
        nd0 += Nr_Dieder[ns];
    }

    dt2 = clock();
    cpu_time = ((Real)(dt2-dt1))/CLOCKS_PER_SEC;

    return(cpu_time);
}
