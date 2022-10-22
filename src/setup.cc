#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

void setup(void);


void setup(void) {

    int         c,acc,i,j,k,l,nljt,ni,ns,nm,nb,m,
                ljcount,elcount,acount,nmtot,IN,IA,
                sp,si,sj,ntot,nstot,nbtot,nctot,sidx,bidx,idxi,idxj,
                a0,b0,d0,i0,mo,idx0;
    Real        dis,eps,sig,e1,e2,s2,s1,rho,maxsig;

    if(NL_CI&&NL_SAI) {
        printf("can't use both fixed and variable neighbor-list-interval\n");
        printf("using NL_CI(10) ...\n");
        NL_SAI = 0;
        nl_intv = 10;
    }

/* parameter ..................................................... */

    dth = 0.5 * dt;
    ad8dtq = 1.0/8.0/dt/dt;
    tausq = tau * tau;
    Lh = 0.5*L;
    adLh = 2.0/L;
    adL34 = 4.0/L/3.0;
    eta = etaL/L;
    elfac = 1.0e+09*NLo*0.001*el_Charge*el_Charge/(4.0*PI*epsilon_0);
    sq_crit = crit*crit;;

/* seed for the random-generator ................................. */

    srandom(r_seed);

/* get the dimension of the system, i.e. the number of sites, ....
   bonds, angles and dihedral angles, allocate memory and ........
   read the data respectivly ..................................... */

//presskey("2.1");

    get_system();

//presskey("2.2");

/* scale the cutoff which was read in units of LJ-sigma .......... */

    maxsig = 0.0;
    for(i=0;i<Nr_LJtypes;i++)
        if(LJtypes[i].sigma>maxsig)
            maxsig = LJtypes[i].sigma;

    cutoff *= maxsig;
    far_ctof *= maxsig;
    cutoff_2 = cutoff*cutoff;
    far_ctof_2 = far_ctof*far_ctof;

    if(far_ctof>Lh) {
        printf("box too small ...\n");
        exit(1);
    }

/* zero .......................................................... */

    for(i=0;i<200;i++) {
        a_dist[i] = 0;
        dha_dist[i] = 0;
        idha_dist[i] = 0;
    }

    nulvek.x = 0.0;
    nulvek.y = 0.0;
    nulvek.z = 0.0;

    zeta = 0.0;
    zetaold = 0.0; 

    S_T = 0.0;
    S_V = 0.0;
    S_Epotlj = 0.0;
    S_Epotcr = 0.0;
    S_Epotck = 0.0;
    S_KinA = 0.0;
    S_APot = 0.0;
    S_DHAPot = 0.0;
    S_IDHAPot = 0.0;

    c_Potr = 0.0;
    c_Potk = 0.0;
    lj_Pot = 0.0;
    Virial = 0.0;
    a_Pot = 0.0;
    dha_Pot = 0.0;
    idha_Pot = 0.0;

    Nr_SHAKE = 0;
    nnup = 0;

    if(Nr_Atot>0) ANGLES = 1;
    else ANGLES = 0;
    if(Nr_Dtot>0) DHANGLES = 1;
    else DHANGLES = 0;
    if(Nr_Itot>0) IDHANGLES = 1;
    else IDHANGLES = 0; 

/* the desired kinetic energy .................................... */

    KinAsoll = 0.001 * (Real)Nr_Ftot * Rgas * Temperature / 2.0;

    ntot = 0;
    nstot = 0;
    c = 0;
    nmtot = 0;
    
/* get several paramerters for the enrgy-calculation*/

    for(sp=0;sp<Nr_Spec;sp++) {
        for(nm=0;nm<Nr_Molecules[sp];nm++) {
        
        MASS[nmtot] = 0.0;
        
            for(ns=0;ns<Nr_Sites[sp];ns++) {

                sidx = nstot+ns;
                
                mass[ntot] = Site[sidx].mass;
                rm[ntot] = 1.0/mass[ntot];
                dtqdm[ntot] = dt*dt/mass[ntot];
                elch[ntot] = Site[sidx].charge * sqrt(elfac);;
                ljt[ntot] = Site[sidx].ljtype;
                MASS[nmtot] += mass[ntot];

                if(ljt[ntot]>0) {
                    ww[ntot] = LJ;
                    if(elch[ntot]!=0.0) {
                        ww[ntot] = LJEL;
                        el_idx[c] = ntot;
                        c++;
                    }
                }
                else if(elch[ntot]!=0.0) {
                    ww[ntot] = EL;
                    el_idx[c] = ntot;
                    c++;
                }

                ntot++;
            }
            nmtot++;
        }
        nstot += Nr_Sites[sp];
    }

    Nr_Charged = c;

/* LJ interaction matrix ............................................... */ 

    for(i=0;i<Nr_LJtypes;i++)
        for(j=0;j<Nr_LJtypes;j++) {
            e1 = LJtypes[i].epsilon;
            e2 = LJtypes[j].epsilon;
            s1 = LJtypes[i].sigma;
            s2 = LJtypes[j].sigma;
            eps = sqrt(e1*e2);
            sig = 0.5*(s1+s2);
            LJP12[i][j] = 4.0*eps*pow(sig,12.0);
            LJP6[i][j] = 4.0*eps*pow(sig,6.0);
            LJF12[i][j] = 48.0*eps*pow(sig,12.0);
            LJF6[i][j] = 24.0*eps*pow(sig,6.0);
        }

/* squared bondlengths and constraint-indices ....................... */    

    i = 0;
    bidx = 0;
    nbtot = 0;
    nstot = 0;

    for(sp=0;sp<Nr_Spec;sp++) {
        for(nm=0;nm<Nr_Molecules[sp];nm++) {
            for(nb=0;nb<Nr_Bonds[sp];nb++) {

                bidx = nbtot+nb;
                zidx1[i] = nstot+Bond[bidx].idx1;
                zidx2[i] = nstot+Bond[bidx].idx2;
                sq_bnd[i] = Bond[bidx].b_eqi*Bond[bidx].b_eqi;
                i++;
            }
            nstot += Nr_Sites[sp];
        }
        nbtot += Nr_Bonds[sp];
    }

/* for looking up the pcf's convert molecular to site-indices */

    for(i=0;i<Nr_pcfs;i++) {

        pcf0_1[i] = tot_index(pcf_sp1[i],0,pcf_si1[i]);
        pcf0_2[i] = tot_index(pcf_sp2[i],0,pcf_si2[i]);
        pcfn_1[i] = tot_index(pcf_sp1[i],Nr_Molecules[pcf_sp1[i]]-1,pcf_si1[i]);
        pcfn_2[i] = tot_index(pcf_sp2[i],Nr_Molecules[pcf_sp2[i]]-1,pcf_si2[i]);
        pcf_d1[i] = Nr_Sites[pcf_sp1[i]];
        pcf_d2[i] = Nr_Sites[pcf_sp2[i]];
    }

/* find the index at which each atom should begin to look for its ...
   neighbors when updating the neighbor-list ( -> in order to exclude
   intramolecular neighbors, the interactions of which are calculated
   seperatly) ....................................................... */

   acc = 0;
   ni = 0;

   for(sp=0;sp<Nr_Spec;sp++) {
       for(nm=0;nm<Nr_Molecules[sp];nm++) {
           ni += Nr_Sites[sp];
           for(ns=0;ns<Nr_Sites[sp];ns++) {
               nextidx[acc] = ni;
               acc++;
           }
       }
   }

/* build the (constant) intramolecular neighborlist, i.e. for each .. */
/* site in each molecular species find the offset-index of all ...... */
/* not excluded neighbors ........................................... */

    ljcount = 0;
    elcount = 0;
    acount = 0;
    idx0 = 0;
    a0 = b0 = d0 = i0 = 0;

    for(sp=0;sp<Nr_Spec;sp++) {
        for(si=0;si<Nr_Sites[sp];si++) {
            for(sj=si+1;sj<Nr_Sites[sp];sj++) {
                IN = 1;
                if(Nr_iDieder[sp]>0)
                for(i=0;i<Nr_iDieder[sp];i++)
                    if(both_in(si,sj,iDieder[i+i0])) { IN = 0; i = Nr_iDieder[sp]; }                  
                if(Nr_Dieder[sp]>0)
                for(i=0;i<Nr_Dieder[sp];i++)
                    if(both_in(si,sj,Dieder[i+d0]))  { IN = 0; i = Nr_Dieder[sp]; }
                if(Nr_Angles[sp]>0)
                for(i=0;i<Nr_Angles[sp];i++)
                    if(both_in(si,sj,Angle[i+a0]))   { IN = 0; i = Nr_Angles[sp]; }
                if(Nr_Bonds[sp]>0)
                for(i=0;i<Nr_Bonds[sp];i++)
                    if(both_in(si,sj,Bond[i+b0]))    { IN = 0; i = Nr_Bonds[sp]; }


                IA = ww[idx0+si]+ww[idx0+sj];
            
                if(IN) {
                    
                    switch(IA) {

                    case 2: 
                    case 3:
                        for(i=0;i<Nr_Molecules[sp];i++) {
                            ljnb_i[ljcount] = idx0+i*Nr_Sites[sp]+si;
                            ljnb_j[ljcount] = idx0+i*Nr_Sites[sp]+sj;
                            ljcount++;
                        }
                        break;
                    case 4: 
                    case 6:
                        for(i=0;i<Nr_Molecules[sp];i++) {
                            anb_i[acount] = idx0+i*Nr_Sites[sp]+si;;
                            anb_j[acount] = idx0+i*Nr_Sites[sp]+sj;;
                            acount++;
                            break;
                        }
                    case 8:
                        for(i=0;i<Nr_Molecules[sp];i++) {
                            elnb_i[elcount] = idx0+i*Nr_Sites[sp]+si;
                            elnb_j[elcount] = idx0+i*Nr_Sites[sp]+sj;
                            elcount++;
                        }
                    }
                }
            }
        }

        idx0 += Nr_Sites[sp]*Nr_Molecules[sp];
        a0 += Nr_Angles[sp];    b0 += Nr_Dieder[sp];
        d0 += Nr_iDieder[sp];   i0 += Nr_Bonds[sp];
    }

    NR_LJIN = ljcount;
    NR_ELIN = elcount;
    NR_AIN = acount;

/* VdW - long range corrections ..................................... */
/* (calc. only if there is just one molecular species) .............. */

    U_lrc = 0.0;
    V_lrc = 0.0;

    if(Nr_Spec==1) {

        Real co;
        if(NL_SAI||NL_CI) co = far_ctof;
        else co = cutoff;
        rho = (Real)Nr_Mtot/L/L/L;

        for(i=0;i<Nr_Sites[0];i++) {
            for(j=0;j<Nr_Sites[0];j++) {

        U_lrc += 8.0*LJtypes[1].epsilon*PI*rho*
              (pow(LJtypes[1].sigma,12.0)*pow(co,-9.0)/9.0
              -pow(LJtypes[1].sigma,6.0)*pow(co,-3.0)/3.0); 
        V_lrc += 16.0*LJtypes[1].epsilon*PI*rho*
              (2.0*pow(LJtypes[1].sigma,12.0)*pow(co,-9.0)/9.0
              -pow(LJtypes[1].sigma,6.0)*pow(co,-3.0)/3.0);
            
            }
        }
    }

    P_lrc = 1000.0*V_lrc/3.0/NLo/L/L/L/1.0e-27;

/* make constant coefficients for the ewald-sum ..................... */

    if(EWALD) make_k();
}

