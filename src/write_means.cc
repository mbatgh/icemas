#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/****************************************************************/
void write_means(Real r_t,Real k_t,Real i_t,Real n_t,
                 Real a_t, Real u_t, Real t_t) {
/****************************************************************/

/*  this subroutine is called at the end of each MD-run, to write
    the averages of thermodynamic functions and eventually other data
    on disk */

    int i,j,sum;
    Real nav,etot,r,c,ang,nw,pcfnorm,dr,nn;
    FILE* file;
    char* filename;

    char dataname[200];
    char angdname[200];
    char dhadname[200];
    char idhadname[200];
    char pcfname[200];
    char data[] = "data.\0";
    char angd[] = "angdis.\0";
    char dhad[] = "dhangdis.\0";
    char idhad[] = "idhangdis.\0";
    char pcfd[] = "pcf. .\0";

    strcpy(dataname,data);
    filename = strcat(dataname,System_name);
    if((file = fopen(filename,"a+"))==NULL) {
        printf("can't open %s ...\n", filename);
    }

    r = t_t-i_t-r_t-k_t-n_t-a_t-u_t;
    nav = (Real)(it-begin_avgs);
    if(nav==0.0) nav = 1.0;
    if(nnup==0) nnup = maxit;

    etot = (S_Epotlj+S_Epotcr+S_Epotck+S_KinA+S_APot+S_DHAPot
           +S_IDHAPot)/nav;

    fprintf(file,"# *********************************************\n#\n");
    fprintf(file,"# the means ...\n#\n");
    fprintf(file,"# Temperature:    %lf [K]\n",S_T/nav);
    fprintf(file,"# E.pot.lj:       %lf [kJ/mol]\n",S_Epotlj/nav);
    fprintf(file,"# E.pot.cr:       %lf [kJ/mol]\n",S_Epotcr/nav);
    fprintf(file,"# E.pot.ck:       %lf [kJ/mol]\n",S_Epotck/nav);
    fprintf(file,"# E.pot.angel:    %lf [kJ/mol]\n",S_APot/nav);
    fprintf(file,"# E.pot.dh_angel: %lf [kJ/mol]\n",S_DHAPot/nav);
    fprintf(file,"# E.pot.improp:   %lf [kJ/mol]\n",S_IDHAPot/nav);
    fprintf(file,"# E.kin.A:        %lf [kJ/mol]\n",S_KinA/nav);
    fprintf(file,"# E.pot.lj.lrc:   %lf [kJ/mol]\n",U_lrc);
    fprintf(file,"# E.pot.c.self:   %lf [kJ/mol]\n",ECself);
    fprintf(file,"# E.tot:          %lf [kJ/mol]\n",etot);
    fprintf(file,"# U_lrc:          %lf [kJ/mol]\n",U_lrc);
    fprintf(file,"# EC_self:        %lf [kJ/mol]\n",ECself);
    fprintf(file,"# delta-E:        %le [kJ/mol]\n",Etot0-etot);
    fprintf(file,"# Virial:         %lf [kJ/mol]\n",S_V/nav);
    fprintf(file,"# Vir_lrc:        %lf [kJ/mol]\n",V_lrc);
    fprintf(file,"# P:              %le [MPa]\n",S_P/nav);
    fprintf(file,"# P_lrc:          %le [MPa]\n#\n",P_lrc*1.0e-06);
    fprintf(file,"# total time:     %lf [sec]\n",t_t);
    fprintf(file,"# r time:         %lf [perc]\n",100.0*r_t/t_t);
    fprintf(file,"# k time:         %lf [perc]\n",100.0*k_t/t_t);
    fprintf(file,"# i time:         %lf [perc]\n",100.0*i_t/t_t);
    fprintf(file,"# n time:         %lf [perc]\n",100.0*n_t/t_t);
    fprintf(file,"# a time:         %lf [perc]\n",100.0*a_t/t_t);
    fprintf(file,"# u time:         %lf [perc]\n",100.0*u_t/t_t);
    fprintf(file,"# rest time:      %lf [perc]\n#\n",100.0*r/t_t);
    fprintf(file,"# Nr_SHAKE/maxit: %lf []\n",(Real)Nr_SHAKE/(Real)maxit);
    fprintf(file,"# maxit/nnup:     %lf []\n",(double)it/(double)nnup);

    fclose(file);

/* write angle-distribution-function in file ............... */

    if(angdf_step<maxit) {

        strcpy(angdname,angd);
        filename = strcat(angdname,System_name);

        if((file = fopen(filename,"w"))==NULL) {
            printf("can't open %s ...\n", filename);

        } else {
            
            nw = (Real)(maxit-begin_avgs)/(Real)angdf_step;
            if(nw==0.0) nw = 1.0/(Real)angdf_step;
            ang = 0;

            for(i=0;i<200;i++) {
                ang += PI/200.0;
                c = (Real)a_dist[i]/nw;
                fprintf(file,"%lf %lf\n", ang, c);
            }

            fclose(file);
        }
    }

/* write dihedral-angle-distribution-function in file ........ */

    if(dhadf_step<maxit) {

        strcpy(dhadname,dhad);
        filename = strcat(dhadname,System_name);

        if((file = fopen(filename,"w"))==NULL) {
            printf("can't open %s ...\n", filename);

        } else {

            nw = (Real)(maxit-begin_avgs)/(Real)dhadf_step;
            if(nw==0.0) nw = 1.0/(Real)angdf_step;
            ang = -PI;

            for(i=0;i<200;i++) {
                ang += PI/100.0;
                c = (Real)dha_dist[i]/nw;
                fprintf(file,"%lf %lf\n", ang, c);
            }
            
            fclose(file);
        }
    }

/* write improper dihedral-angle-distribution-function in file ... */

    if(idhdf_step<maxit) {

        strcpy(idhadname,idhad);
        filename = strcat(idhadname,System_name);

        if((file = fopen(filename,"w"))==NULL) {
            printf("can't open %s ...\n", filename);

        } else {

            nw = (Real)(maxit-begin_avgs)/(Real)idhdf_step;
            if(nw==0.0) nw = 1.0/(Real)angdf_step;
            ang = -PI;

            for(i=0;i<200;i++) {
                ang += PI/100.0;
                c = (Real)idha_dist[i]/nw;
                fprintf(file,"%lf %lf\n", ang, c);
            }

            fclose(file);
        }
    }

/* write pair-correlation-functions in file ........ */

    int norm;
    Real rnorm,Vshell,rhoist,rhosoll;

    if((int)Nr_pcfs) {

        char null = '0';
        
        for(i=0;i<Nr_pcfs;i++) {

            pcfd[4] = (null + (char)i);
            strcpy(pcfname,pcfd);
            filename = strcat(pcfname,System_name);

            norm = ((int)(it/pcf_step));
            rnorm = (Real)(norm*Nr_Molecules[pcf_sp1[i]]);
            rhosoll = ((Real)Nr_Molecules[pcf_sp2[i]])/(L*L*L);

            if((file = fopen(filename,"w"))==NULL) {
                printf("can't open %s ...\n", filename);
            } else {
                dr = Lh/(Real)NR_PCF_CHANS;
                for(j=0;j<NR_PCF_CHANS;j++) {
                    Vshell = (pow(((Real)j+1.0)*dr,3.0)-pow((Real)j*dr,3.0))*
                              4.0*PI/3.0;
                    rhoist = ((Real)pcf[i][j])/(rnorm*Vshell);

                    fprintf(file,"%lf %lf\n", (Real)j*dr,rhoist/rhosoll);
                }

                fclose(file);
            }
        }
    }
}
