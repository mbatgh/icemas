/*****************************************************************************/
/*                                                                           */
/*                                  ICEMAS                                   */
/*                                                                           */
/*            Integration of the Classical Equations of motion               */
/*                   of Molecules of Arbitrary Shape                         */
/*                                                                           */
/*****************************************************************************/


#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"
#include <mic_noecho.h>


/* .................................................... 
       global integers 
   .................................................... */
int

/* total number of iterations */       maxit,
/* current number of iterations */     it,
/* seed for the random-generator */    r_seed,
/* iteration-step after which the 
   thermostat is shut down */          end_thst,
/* iteration-step after which summing 
   up for the means starts */          begin_avgs,
/* interval for writing thermodynamic 
   data in file */                     write_therm_step,
/* interval for updating the 
   angle ditributions */               angdf_step,dhadf_step,idhdf_step,
/* parameter for the ewald-sum */      K_max,       
/* parameter for the ewald-sum */      Ksq_max,     
/* the indices of charged sites */     el_idx[NR_SMAX],
/* the number of charged sites */      Nr_Charged,
/* ww indicates the type of 
   interaction(s) (LJ or electrostatic
   or both) for each site */           ww[NR_SMAX], 
/* the Lennard-Jones-type (determines 
   which sigma and epsilon used when 
   calculating the LJ-potential) */    ljt[NR_SMAX],         
/* number of intramolecular 
   neighbors */                        NR_LJIN,NR_ELIN,NR_AIN,
/* total number of neighbors */        NR_LJN,NR_ELN,NR_AN,
/* the lj neighborlist */              ljnb_i[NR_NEIMAX],ljnb_j[NR_NEIMAX],
/* the charged neighborlist */         elnb_i[NR_NEIMAX],elnb_j[NR_NEIMAX],
/* the charged+lj neighborlist */      anb_i[NR_NEIMAX],anb_j[NR_NEIMAX],
/* parameter needed when setting 
   up the neighborlist */              nextidx[NR_SMAX],     
/* total number of 
   neighborlist-updates */             nnup,
                  
/* constraint-indices */               zidx1[NR_CMAX], zidx2[NR_CMAX],
/* constant interval for 
   neighborlist update */              nl_intv,
/* constant interval for calculation 
   of phase separation */              sep_intv,
/* total number of different molecular 
   species' */                         Nr_Spec,
/* total number of molecules */        Nr_Mtot,
/* total number of sites */            Nr_Stot,
/* total number of bonds */            Nr_Btot,     
/* total number of angles */           Nr_Atot,     
/* total number of dihedral angles */  Nr_Dtot,
/* total nr. of impr. dih. angles */   Nr_Itot,
/* total number of 
   Lennard-Jones-types */              Nr_LJtypes,
/* total number degrees of freedom */  Nr_Ftot,
/* total number of SHAKE-loops */      Nr_SHAKE,
/* the above parameters for each 
   species of molecule */              Nr_Molecules[NR_SPECMAX],
                                       Nr_Sites[NR_SPECMAX],
                                       Nr_Bonds[NR_SPECMAX],
                                       Nr_Angles[NR_SPECMAX],
                                       Nr_Dieder[NR_SPECMAX],
                                       Nr_iDieder[NR_SPECMAX],
/* some booleans for the 
   force-calculation */                NL_SAI,
                                       NL_CI,
                                       GETNGB,
                                       ANGLES,
                                       DHANGLES,
                                       IDHANGLES,
                                       EWALD,
                                       NHT,
/* numbers of the molecule and  
  (dihedral) angle, that 
   should be observed */               SP_ADIST,NR_ADIST,
                                       SP_DHADIST,NR_DHADIST,
                                       SP_IDHADIST,NR_IDHADIST,
/* arrays for
   distribution-functions */           a_dist[200],
                                       dha_dist[200],
                                       idha_dist[200],
/* calc pcf every pcf_step its */      pcf_step,
/* indices of the sites for pcf's */   pcf0_1[NR_PCFSMAX],pcf0_2[NR_PCFSMAX],
                                       pcfn_1[NR_PCFSMAX],pcfn_2[NR_PCFSMAX],
                                       pcf_d1[NR_PCFSMAX],pcf_d2[NR_PCFSMAX],
                                       pcf_sp1[NR_PCFSMAX],pcf_si1[NR_PCFSMAX],
                                       pcf_sp2[NR_PCFSMAX],pcf_si2[NR_PCFSMAX],
/* pcf-channels */                     pcf[NR_PCFSMAX][NR_PCF_CHANS],
/* total number of pcfs */             Nr_pcfs;


/* ........................................................ 
       global structures 
   ........................................................ */

/* structure containig information 
   on one site */                      t_site* Site;    
/* structure containig information 
   on one bond */                      t_bond* Bond;    
/* structure containig information 
   on one angle */                     t_angle* Angle;   
/* structure containig information 
   on one dihedral angle */            t_dieder* Dieder;  
/* structure containig information 
   on one improp. dihedral angle */    t_idieder* iDieder;
/* structure containig information 
   on one LJ-type */                   t_ljtype LJtypes[NR_LJTMAX];


/* .........................................................
    global vectors
   ......................................................... */
VEKTOR

/* the current position of a site */   RS[NR_SMAX],
/* the old position of a site */       RSo[NR_SMAX],
/* the new position of a site */       RSn[NR_SMAX],
/* force on a site */                  FS[NR_SMAX],
/* displacement of each site 
   (for the NL) */                     displacement[NR_SMAX],
/* (0,0,0) */                          nulvek;         


/* ......................................................
       global floating points 
   ...................................................... */
Real

/* the desired temperature [K] of the 
   system */                           Temperature,
/* the timestep [ps] and
   related parameters */               dt, ad8dtq, dth,
/* parameter for the 
   SHAKE-algorithm */                  crit, sq_crit,
/* Lennard-Jones-parameters */         epsilon[NR_LJTMAX], 
/* [kJ/mol], [nm] */                   sigma[NR_LJTMAX],
/* the boxlength and 
   related parameters */               L,adLh,adL34,Lh,
/* energies */                         Etot,
                                       Etot0,
                                       lj_Pot, 
                                       c_Potr, 
                                       c_Potk, 
                                       Virial,
                                       ECself, 
                                       KinA,
                                       dha_Pot,
                                       idha_Pot,
                                       a_Pot,
/* the two (squared) cutoffs */        cutoff, cutoff_2, 
                                       far_ctof, far_ctof_2,
/* the mass and related parameters */  mass[NR_SMAX], 
                                       rm[NR_SMAX], 
                                       dtqdm[NR_SMAX],
                                       *MASS,
/* the squared bond lenghts' */        sq_bnd[NR_SMAX],
/* the desired kinetic energy */       KinAsoll,
/* parameters for the
   Nose-Hoover-thermostat */           zeta, zetaold,
                                       tau, tausq,
/* accumulators for the means */       S_KinA,
                                       S_Epotlj,
                                       S_Epotcr,
                                       S_Epotck,
                                       S_T,
                                       S_V,
                                       S_P,
                                       S_APot,
                                       S_DHAPot,
                                       S_IDHAPot,
/* LJ-long-range-corrections */        P_lrc, U_lrc, V_lrc,
/* matrix for the different 
   LJ-interaction-parameters */        LJP12[NR_LJTMAX][NR_LJTMAX],
                                       LJP6[NR_LJTMAX][NR_LJTMAX],
                                       LJF6[NR_LJTMAX][NR_LJTMAX],
                                       LJF12[NR_LJTMAX][NR_LJTMAX],
/* parameter for the ewald-sum */      eta, etaL,
/* scaling-parameter for charges */    elfac,
/* the el. (partial-)charges [e0] */   elch[NR_SMAX],
/* parameter for the ewald-sum */      ewfac[KTOT_MAX];


/* ..........................................................
       global chars
   .......................................................... */

char
/* name of the file with information 
   on the system */                    System_file[100],
/* name of the file with the 
   initial-coordinates */              R0_file[100],
/* ...  */                             System_name[100];

/* a local function */

Real angles(void);



/* =================================================================== */        
/* =                         MAIN                                    = */
/* =================================================================== */


int main(int ARGC, char* ARGV[]) {

    int       n,relstep;
    Real      r_time,k_time,i_time,n_time,a_time,u_time,total_time,
              s,smin;
    clock_t   t1,t2;

/* do run 1 to ARGC-1 as specified on the commandline .............. */

    for(n=1;n<ARGC;n++) {

        a_time      = 0.0;
        r_time      = 0.0;
        k_time      = 0.0;
        i_time      = 0.0;
        n_time      = 0.0;
        u_time      = 0.0;
        total_time  = 0.0;

        GETNGB      = 1;

/* the commandline-argument(s) is (are) the name(s) of the system(s) */

        strcpy(System_name,ARGV[n]);

        get_parameter();
        setup();
        make_data_file();
        get_r0();

        t1 = clock();

        it = 0;

/* the main - loop ............................................. */  


/* the boolean parameters NL_SAI (use a neighboelist with a self adjusting
   interval), NL_CI (use a neighborlist with constant interval) and EWALD 
   (consider electrostatic interactions using the EWALD-method) determine
   the subroutines used in the main-loop */

        if(NL_SAI) {

            if(EWALD) {

                do {

                    it++;

                    if(GETNGB)      n_time += het_neighbors();
                                              zero_forces();
                                    a_time += angles();
                                    r_time += nl_r_force();
                                    k_time += k_force();
                                    i_time += move();
                                    u_time += update_dist();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));
            
            } else {
            
                do {

                    it++;

                    if(GETNGB)      n_time += hom_neighbors();
                                              zero_forces();
                                    a_time += angles();          
                                    r_time += nl_force();
                                    i_time += move();
                                    u_time += update_dist();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));
            }

        } else if(NL_CI) {
            
            if(EWALD) {

                do {

                    if(fmod((double)it,(double)nl_intv)==0.0)
                                    n_time += het_neighbors();
                                              it++;
                                              zero_forces();
                                    a_time += angles();
                                    r_time += nl_r_force();
                                    k_time += k_force();
                                    i_time += move();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));

            } else {

                do {

                    if(fmod((double)it,(double)nl_intv)==0.0)
                                    n_time += hom_neighbors();
                                    it++;
                                              zero_forces();
                                    a_time += angles();       
                                    r_time += nl_force();
                                    i_time += move();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));
            }                       

        } else {
        
            if(EWALD) {

                do {

                    it++;

                                              zero_forces();
                                    a_time += angles();
                                    r_time += r_force();
                                    k_time += k_force();
                                    i_time += move();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));
            
            } else {
            
                do {

                    it++;

                                              zero_forces();
                                    a_time += angles();
                                    r_time += force();
                                    i_time += move();
                                              tdf();

                } while((it<maxit)&&((fopen("the_end","r"))==NULL));
            }
        }        


        t2 = clock();
        total_time = ((Real)(t2-t1))/CLOCKS_PER_SEC;

        write_lastr();
        write_means(r_time,k_time,i_time,n_time,a_time,u_time,total_time);
    }

    exit(1);
}



Real angles(void) {

    Real t;
    t = 0.0;

    if(IDHANGLES)   t += impdh_force();
    if(DHANGLES)    t += torsion_force();
    if(ANGLES)      t += angle_force();

    return(t);
}
