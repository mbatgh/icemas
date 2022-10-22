#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/****************************************************************/
void    make_data_file() {
/****************************************************************/

    int i;
    char data[200];
    char* filename;
    FILE* file;

    data[0] = 'd';
    data[1] = 'a';
    data[2] = 't';
    data[3] = 'a';
    data[4] = '.';
    data[5] = '\0';

    filename = strcat(data,System_name);

    if((file = fopen(filename,"r"))!=NULL) {
        printf("%s existiert bereits ...\n", filename);
        fclose(file);
        exit(1);
    }
    else {
        file = fopen(data,"w");
        fprintf(file,"#      *** MD-simulation ***\n#\n#\n");
        fprintf(file,"# system_filename:      %s\n", System_file); 
        fprintf(file,"# r0_filename:          %s\n", R0_file);
        for(i=0;i<Nr_Spec;i++)
        fprintf(file,"# Nr_Molecules(%d):      %d\n", i,Nr_Molecules[i]);
        fprintf(file,"# maxit:                %d\n", maxit);
        fprintf(file,"# timestep dt:          %le\n", dt);
        fprintf(file,"# Temperature:          %lf\n", Temperature);
        fprintf(file,"# NHT-parameter tau:    %lf\n", tau);
        fprintf(file,"# L:                    %lf\n", L);
        fprintf(file,"# etaL: (ewald-param.)  %lf\n", etaL);
        fprintf(file,"# crit: (SHAKE-dx)      %le\n", crit);
        fprintf(file,"# EWALD:                %d\n", EWALD);
        fprintf(file,"# KMAX:                 %d\n", K_max);
        fprintf(file,"# KSQMAX:               %d\n", Ksq_max);
        fprintf(file,"# cutoff:                 %lf\n", cutoff);
        fprintf(file,"# far_ctof:             %lf\n", far_ctof);
        fprintf(file,"# random-seed:          %d\n", r_seed);
        fprintf(file,"# write_therm_step:     %d\n", write_therm_step);
        fprintf(file,"# end_thst:             %d\n", end_thst);
        fprintf(file,"# begin_avgs:           %d\n", begin_avgs);
        fprintf(file,"# nl_intv:     %d, ",           nl_intv);
        fprintf(file,"sep_intv:    %d\n",             sep_intv); 
        fprintf(file,"# NL_SAI:      %d\n",             NL_SAI); 
        fprintf(file,"# SP_ADIST:    %d, ",           SP_ADIST);
        fprintf(file,"NR_ADIST:    %d\n",             NR_ADIST);
        fprintf(file,"# SP_DHADIST:  %d, ",           SP_DHADIST);
        fprintf(file,"NR_DHADIST:  %d\n",             NR_DHADIST); 
        fprintf(file,"# SP_IDHADIST: %d, ",           SP_IDHADIST);
        fprintf(file,"NR_IDHADIST: %d\n#\n",             NR_IDHADIST);
        fclose(file);        
    }

    printf("wrote data-file %s\n", data);
}
