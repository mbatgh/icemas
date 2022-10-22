/**************************************************************
extern declaration of functions and global variables 
***************************************************************/

#include "first.h"


extern void         get_parameter();
extern void         setup(void);
extern void         get_system(void);
extern void         get_r0();
extern void         make_k(void);
extern void         make_data_file();

extern Real         het_neighbors(void);
extern Real         hom_neighbors(void);
extern void         zero_forces(void);
extern Real         angle_force(void);
extern Real         torsion_force(void);
extern Real         impdh_force(void);
extern Real         r_force(void);
extern Real         k_force(void);
extern Real         nl_r_force(void);
extern Real         nl_force(void);
extern Real         force(void);
extern Real         move(void);
extern void         shake(void);
extern void         tdf(void);
extern Real         update_dist(void);

extern void         phase_sep(Real*,Real*,Real*);
extern void         get_a_dist(Real);
extern void         get_dha_dist(Real);
extern void         get_idha_dist(Real);
extern void         write_lastr();
extern void         write_means(Real,Real,Real,Real,Real,Real,Real);
extern void         write_r(int,int,char*);

extern void         presskey(char* message);
extern void         presskey(void);
extern void         scale_tmp(Real,Real);
extern void         scale_new_tmp(Real ist, Real soll);
extern int          both_in(int,int,t_idieder);
extern int          both_in(int,int,t_dieder);
extern int          both_in(int,int,t_angle);
extern int          both_in(int,int,t_bond);
extern VEKTOR       cent_of_mass(VEKTOR*,int,int);
extern int          tot_index(int,int,int);

extern void         readline(FILE*,char*);
extern void         readline(FILE*,float*);
extern void         readline(FILE*,double*);
extern void         readline(FILE*,int*);
extern void         readline(FILE*,short*);
extern void         readline(FILE*,int*,int*,int*,int*,int*,int*);
extern void         readline(FILE*,int*,int*,int*,int*,int*);
extern void         readline(FILE*,int*,int*,int*,int*);
extern void         readline(FILE*,int*,int*,int*,int*,int*,int*,Real*,Real*);
extern void         readline(FILE*,int*,int*,int*,int*,int*,Real*,Real*);
extern void         readline(FILE*,int*,int*,int*,int*,Real*,Real*);
extern void         readline(FILE*,int*,int*,int*,Real*);
extern void         readline(FILE*,int*,int*,Real*,Real*);
extern void         readline(FILE*,int*,Real*,Real*);
extern int          get_data_line(FILE*,char*,char*);
extern int          get_data_line(FILE*,char*,float*);
extern int          get_data_line(FILE*,char*,double*);
extern int          get_data_line(FILE*,char*,int*);
extern int          readmore_4_int(FILE*,char*,int,int*,int*,int*,int*);


extern int          maxit,it,r_seed,end_thst,begin_avgs,
                    write_therm_step,angdf_step,dhadf_step,
                    nnup,nl_intv,idhdf_step,Nr_SHAKE,
                    K_max,Ksq_max,sep_intv,NHT,
                    el_idx[NR_SMAX],Nr_Charged,
                    nextidx[NR_SMAX],ww[NR_SMAX],ljt[NR_SMAX],
                    first[NR_SMAX],last[NR_SMAX],
                    zidx1[NR_CMAX],zidx2[NR_CMAX],
                    Nr_Mtot,Nr_Stot,Nr_Btot,Nr_Atot,Nr_Dtot,
                    Nr_Ftot,Nr_LJtypes,Nr_Spec,Nr_Itot,
                    Nr_Molecules[NR_SPECMAX],
                    Nr_Sites[NR_SPECMAX],
                    Nr_Bonds[NR_SPECMAX],
                    Nr_Angles[NR_SPECMAX],
                    Nr_Dieder[NR_SPECMAX],
                    Nr_iDieder[NR_SPECMAX],
                    a_dist[200],dha_dist[200],idha_dist[200],
                    NL_SAI,NL_CI,GETNGB,
                    NR_LJIN,NR_ELIN,NR_AIN,
                    NR_LJN,NR_ELN,NR_AN,
                    anb_i[NR_NEIMAX],anb_j[NR_NEIMAX],
                    ljnb_i[NR_NEIMAX],ljnb_j[NR_NEIMAX],
                    elnb_i[NR_NEIMAX],elnb_j[NR_NEIMAX],
                    ANGLES,DHANGLES,EWALD,IDHANGLES,
                    SP_ADIST,NR_ADIST,
                    SP_DHADIST,NR_DHADIST,
                    SP_IDHADIST,NR_IDHADIST,
                    pcf0_1[NR_PCFSMAX],pcf0_2[NR_PCFSMAX],
                    pcfn_1[NR_PCFSMAX],pcfn_2[NR_PCFSMAX],
                    pcf_d1[NR_PCFSMAX],pcf_d2[NR_PCFSMAX],
                    pcf[NR_PCFSMAX][NR_PCF_CHANS],
                    pcf_sp1[NR_PCFSMAX],pcf_si1[NR_PCFSMAX],
                    pcf_sp2[NR_PCFSMAX],pcf_si2[NR_PCFSMAX],
                    Nr_pcfs,pcf_step;


extern VEKTOR       RS[NR_SMAX], RSo[NR_SMAX], RSn[NR_SMAX],
                    FS[NR_SMAX], displacement[NR_SMAX], 
                    nulvek, distance[NR_SMAX];

extern Real         Temperature,dt,rho,
                    crit,sq_crit,*MASS,
                    spce_epsilon,spce_sigma,
                    epsilon[NR_LJTMAX],sigma[NR_LJTMAX],
                    adLh,adL34,L,Lh,NoMrzp,ad8dtq,dth,adVm,
                    lj_Pot,c_Potr,c_Potk,ECself,Virial,
                    KinA,dha_Pot,idha_Pot,a_Pot,Etot,Etot0,
                    cutoff_2,cutoff,far_ctof_2,far_ctof,
                    mass[NR_SMAX],rm[NR_SMAX],
                    dtqdm[NR_SMAX],sq_bnd[NR_SMAX],
                    KinAsoll,zeta,zetaold,tausq,tau,
                    S_KinA,S_Epotlj,S_Epotcr,S_Epotck,S_T,S_P,S_V,
                    S_APot,S_DHAPot,S_IDHAPot,P_lrc,U_lrc,V_lrc,
                    LJP12[NR_LJTMAX][NR_LJTMAX],
                    LJP6[NR_LJTMAX][NR_LJTMAX],
                    LJF6[NR_LJTMAX][NR_LJTMAX],
                    LJF12[NR_LJTMAX][NR_LJTMAX],
                    eta,etaL,elfac,elch[NR_SMAX],
                    ewfac[KTOT_MAX];

extern char         System_file[100],
                    System_name[100],
                    R0_file[100];

