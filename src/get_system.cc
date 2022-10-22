#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

void get_system(void);

/**************************************************************/
void get_system(void) {
/**************************************************************/

/*  subroutine for reading the system-data (size of the system, 
    LJ-parameter, etc. ) from an input-file, the name of which 
    is specified in the parameter-file */

    int     mu,sp,ns,nb,na,nd,nid,nm,
            sptot,nstot,nbtot,natot,ndtot,nidtot,
            nlj,so_bond,so_angle,so_dieder,so_idieder,so_site;
    FILE*   sys_data;

    if((sys_data=fopen(System_file,"r"))==NULL) {
        printf("can't see no %s\n", System_file);
        exit(1);
    } else {

    printf("\nreading system data ...\n");

    readline(sys_data,&Nr_Spec);

    if(Nr_Spec>NR_SPECMAX) {
        printf("not more than %d different species of molecules ...\n",NR_SPECMAX);
        exit(1);
    }

    readline(sys_data,&Nr_LJtypes);

    if(Nr_LJtypes>NR_LJTMAX) {
        printf("not more than %d different LJ-types ...\n", NR_LJTMAX);
        exit(1);
    } 

    for(sp=0;sp<Nr_Spec;sp++) {
        readline(sys_data,&(Nr_Molecules[sp]),
                          &(Nr_Sites[sp]),
                          &(Nr_Bonds[sp]),
                          &(Nr_Angles[sp]),
                          &(Nr_Dieder[sp]),
                          &(Nr_iDieder[sp]));
    }

/* allocate the proper amount of memory ........................... */

    so_bond = sizeof(t_bond);
    so_angle = sizeof(t_angle);
    so_dieder = sizeof(t_dieder);
    so_site = sizeof(t_site);
    so_idieder = sizeof(t_idieder);
    
    nstot = 0;
    nbtot = 0;
    natot = 0;
    ndtot = 0;
    nidtot = 0;

    for(sp=0;sp<Nr_Spec;sp++) {

        nstot += Nr_Sites[sp];
        nbtot += Nr_Bonds[sp];
        natot += Nr_Angles[sp];
        ndtot += Nr_Dieder[sp];
        nidtot+= Nr_iDieder[sp];
    }

    for(sp=0;sp<Nr_Spec;sp++) {
        Site = (t_site*)calloc(nstot,so_site);
        Bond = (t_bond*)calloc(nbtot,so_bond);
        Angle = (t_angle*)calloc(natot,so_angle);
        Dieder = (t_dieder*)calloc(ndtot,so_dieder);
        iDieder = (t_idieder*)calloc(nidtot,so_idieder);
    }

    Nr_Stot = 0;
    Nr_Atot = 0;
    Nr_Btot = 0;
    Nr_Dtot = 0;
    Nr_Mtot = 0; 
    Nr_Itot = 0;
    
    for(sp=0;sp<Nr_Spec;sp++) {
        mu = Nr_Molecules[sp];
        Nr_Stot += mu * Nr_Sites[sp];
        Nr_Btot += mu * Nr_Bonds[sp];
        Nr_Atot += mu * Nr_Angles[sp];
        Nr_Dtot += mu * Nr_Dieder[sp];
        Nr_Itot += mu * Nr_iDieder[sp];
        Nr_Mtot += mu;
    }

    Nr_Ftot = 3 * Nr_Stot - Nr_Btot;

    MASS = (Real*)calloc(Nr_Mtot,sizeof(Real));

/* get the data ................................................... */

    nstot = 0;
    nbtot = 0;
    natot = 0;
    ndtot = 0;
    nidtot = 0;

    for(nlj=0;nlj<Nr_LJtypes;nlj++) {
        readline(sys_data,&(LJtypes[nlj].number),
                          &(LJtypes[nlj].epsilon),
                          &(LJtypes[nlj].sigma)); 
    }

    for(sp=0;sp<Nr_Spec;sp++) {

        for(ns=0;ns<Nr_Sites[sp];ns++) {
            readline(sys_data,&(Site[nstot].number),
                              &(Site[nstot].ljtype),
                              &(Site[nstot].mass),
                              &(Site[nstot].charge));
            nstot++;
        }
        
        for(nb=0;nb<Nr_Bonds[sp];nb++) {
            readline(sys_data,&(Bond[nbtot].number),
                              &(Bond[nbtot].idx1),
                              &(Bond[nbtot].idx2),
                              &(Bond[nbtot].b_eqi));
            nbtot++;
        }

        for(na=0;na<Nr_Angles[sp];na++) {
            readline(sys_data,&(Angle[natot].number),
                              &(Angle[natot].idx1),
                              &(Angle[natot].idx2),
                              &(Angle[natot].idx3),
                              &(Angle[natot].a_eqi),
                              &(Angle[natot].force));
            natot++;
        }

        for(nd=0;nd<Nr_Dieder[sp];nd++) {
            readline(sys_data,&(Dieder[ndtot].number),
                              &(Dieder[ndtot].idx1),
                              &(Dieder[ndtot].idx2),
                              &(Dieder[ndtot].idx3),
                              &(Dieder[ndtot].idx4),
                              &(Dieder[ndtot].period),
                              &(Dieder[ndtot].d_eqi),
                              &(Dieder[ndtot].force));
            ndtot++;
        }

        for(nid=0;nid<Nr_iDieder[sp];nid++) {
            readline(sys_data,&(iDieder[nidtot].number),
                              &(iDieder[nidtot].idx1),
                              &(iDieder[nidtot].idx2),
                              &(iDieder[nidtot].idx3),
                              &(iDieder[nidtot].idx4),
                              &(iDieder[nidtot].d_eqi),
                              &(iDieder[nidtot].force));
            nidtot++;
        }
    }
    }
}

