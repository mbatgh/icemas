#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"



/***************************************************************/
Real impdh_force(void) {
/***************************************************************/

/*  this subroutine calculates the potential of each given
    improper dihedral angle and the resulting forces */

    int     i,j,aidx,nd,ns,nm,nd0,ns0,idx1,idx2,idx3,idx4;
    clock_t dt1,dt2;
    Real    cpu_time,fi,dfi,fi_0,f_c,pot,dpot,
            sinfi,cosfi,signum,
            nof_a,nof_b,nof_c,
            adsinfi,adcosfi,switc;
    VEKTOR  r_ij,r_kl,r_jk,r_a,r_b,r_c,
            dcpa,dcpb,dspa,dspb,dspc,dprjk,
            dum1,dum2,dum3,vf1,vf2,vf3,vf4;

    dt1 = clock();

    idha_Pot = 0.0;
    nd0 = 0;
    ns0 = 0;

    for(ns=0;ns<Nr_Spec;ns++) {
        for(nm=0;nm<Nr_Molecules[ns];nm++) {
            for(nd=0;nd<Nr_iDieder[ns];nd++) {

                aidx = nd0+nd;
                idx1 = ns0 + iDieder[aidx].idx1;
                idx2 = ns0 + iDieder[aidx].idx2;
                idx3 = ns0 + iDieder[aidx].idx3;
                idx4 = ns0 + iDieder[aidx].idx4;
                fi_0 = iDieder[aidx].d_eqi;
                f_c = iDieder[aidx].force;

                r_ij = RS[idx1] - RS[idx2];
                convolute(r_ij,adL34,L);
                r_jk = RS[idx2] - RS[idx3];
                convolute(r_jk,adL34,L);
                r_kl = RS[idx3] - RS[idx4];
                convolute(r_kl,adL34,L);

                r_a = r_ij%r_jk;
                r_a.y = -r_a.y;
                r_b = r_jk%r_kl;
                r_b.y = -r_b.y;
                r_c = r_jk%r_a;
                r_c.y = -r_c.y;

                nof_a = 1.0/sqrt(max(SMALL,norm(r_a)));
                nof_b = 1.0/sqrt(max(SMALL,norm(r_b)));
                nof_c = 1.0/sqrt(max(SMALL,norm(r_c)));

                r_a = r_a*nof_a;
                r_b = r_b*nof_b;
                r_c = r_c*nof_c;

                cosfi = r_a * r_b;
                sinfi = r_c * r_b;

//              fi = copysign(acos(min(1.0,max(-1.0,cosfi))),sinfi);
                fi = PI-copysign(acos(min(1.0,max(-1.0,cosfi))),sinfi);
                if(fi>PI)  fi -= 2.0*PI;
                if(fi<-PI) fi += 2.0*PI;

                dfi = fi - fi_0;
                dfi += 4.0*PI;
                while (dfi > PI) dfi -= 2.0*PI;
                dpot = f_c * dfi;
                pot =  dfi * dpot;
                idha_Pot += pot;
                dpot = 2*dpot;

                switc = (fabs(sinfi)>SMALL) ? -1.0 : 0.0;
                adsinfi = dpot*switc*copysign(1.0,sinfi)/max(fabs(sinfi),SMALL);
                adcosfi = dpot*(switc+1.0)*copysign(1.0,cosfi)/max(fabs(cosfi),SMALL);

                dum1 = r_a * cosfi - r_b;
                dcpa = dum1 * nof_a;
                dum1 = r_b * cosfi - r_a;
                dcpb = dum1 * nof_b;
                dum1 = r_c * sinfi - r_b;
                dspc = dum1 * nof_c;
                dum1 = r_b * sinfi - r_c;
                dspb = dum1 * nof_b;

                dprjk.x = - (    r_jk.y * r_ij.y + r_jk.z * r_ij.z) * dspc.x
	                       + (2 * r_jk.x * r_ij.y - r_ij.x * r_jk.y) * dspc.y
	                       + (2 * r_jk.x * r_ij.z - r_ij.x * r_jk.z) * dspc.z;
                dprjk.y = - (    r_jk.z * r_ij.z + r_jk.x * r_ij.x) * dspc.y
		                   + (2 * r_jk.y * r_ij.z - r_ij.y * r_jk.z) * dspc.z
		                   + (2 * r_jk.y * r_ij.x - r_ij.y * r_jk.x) * dspc.x;
                dprjk.z = - (    r_jk.x * r_ij.x + r_jk.y * r_ij.y) * dspc.z
		                   + (2 * r_jk.z * r_ij.x - r_ij.z * r_jk.x) * dspc.x
                           + (2 * r_jk.z * r_ij.y - r_ij.z * r_jk.y) * dspc.y;

                dum1 = dspb%r_kl;
                dum1.y = -dum1.y;
                dprjk = dprjk-dum1;
                dprjk = dprjk*adcosfi;
                dum1 = dcpa%r_ij;
                dum1.y = -dum1.y;
                dum2 = dcpb%r_kl;
                dum2.y = -dum2.y;
                dum3 = dum1-dum2;
                dprjk = dprjk + dum3 * adsinfi;

		        vf1.x = adsinfi * (r_jk.y * dcpa.z - dcpa.y * r_jk.z) +
      			        adcosfi * ((r_jk.y * r_jk.y + r_jk.z * r_jk.z)*dspc.x - 
				        r_jk.x*r_jk.y*dspc.y - r_jk.x*r_jk.z*dspc.z);
		        vf1.y = adsinfi * (r_jk.z * dcpa.x - dcpa.z * r_jk.x) +
      			        adcosfi * ((r_jk.z * r_jk.z + r_jk.x * r_jk.x)*dspc.y - 
				        r_jk.y*r_jk.z*dspc.z - r_jk.y*r_jk.x*dspc.x);
		        vf1.z = adsinfi * (r_jk.x * dcpa.y - dcpa.x * r_jk.y) +
				        adcosfi * ((r_jk.x * r_jk.x + r_jk.y * r_jk.y)*dspc.z - 
				        r_jk.z*r_jk.x*dspc.x - r_jk.z*r_jk.y*dspc.y);

		        vf2 = dprjk - vf1;

                dum1 = r_jk % dcpb;
                dum2 = r_jk % dspb;
                dum1 = dum1 * adsinfi;
                dum2 = dum2 * adcosfi;
                vf4 = dum1 + dum2;
                vf4.y = -vf4.y;
		        vf3 = (vf4 + dprjk) * -1.0;

                FS[idx1] = FS[idx1] - vf1;
                FS[idx2] = FS[idx2] - vf2;
                FS[idx3] = FS[idx3] - vf3;
                FS[idx4] = FS[idx4] - vf4;

//printf("imp fi = %le, pot = %le, dpot = %le",fi,pot,dpot);
//presskey(" ...");

                if((ns==SP_IDHADIST)&&(nd==NR_IDHADIST)&&
                   (fmod(it,idhdf_step)==0.0)&&it>begin_avgs) 
                    get_idha_dist(fi);
            }
            ns0 += Nr_Sites[ns];
        }
        nd0 += Nr_iDieder[ns];
    }

    dt2 = clock();
    cpu_time = ((Real)(dt2-dt1))/CLOCKS_PER_SEC;

    return(cpu_time);
}

