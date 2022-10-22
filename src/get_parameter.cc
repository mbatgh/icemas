#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/****************************************************************/
void get_parameter() {
/****************************************************************/

/*  this subroutine reads the parameters for a single MD-run from a 
    file, the name of which is specified on the command-line 
    (for most of the parameters a default is given) */

    int i;
    char c,temp[100];
    FILE* parameter;
    char pars[200];
    char* parameter_filename;

    pars[0] = 'p';
    pars[1] = 'a';
    pars[2] = 'r';
    pars[3] = 's';
    pars[4] = '.';
    pars[5] = '\0';

    parameter_filename = strcat(pars,System_name);

    if((parameter=fopen(parameter_filename,"r"))==NULL) 
        { printf("can't see no %s ...\n", parameter_filename); exit(1); }
    else printf("\ncurrent parameters:\n\n");

/* the defaults ............................................... */

    maxit = 0;
    dt = 0.001;
    Temperature = 300.0;
    tau = 0.02;
    L = 1.9724438;
    etaL = 5.36;
    crit = 1.0e-03;
    K_max = 5;
    Ksq_max = 27;
    cutoff = 2.5;
    far_ctof = 2.7;
    r_seed = 11;
    write_therm_step = 10;
    angdf_step = 100000;
    dhadf_step = 100000;
    idhdf_step = 100000;
    sep_intv = 100000;
    end_thst = 0;
    NHT = 0;
    begin_avgs = 4000;
    EWALD = 0;
    nl_intv = 0;
    sep_intv = 0;
    NL_SAI = 0;
    NL_CI = 0;
    SP_ADIST = 0;
    NR_ADIST = 0;
    SP_DHADIST = 0;
    NR_DHADIST = 0;
    SP_IDHADIST = 0;
    NR_IDHADIST = 0;

/* read parameters from file .................................. */

    if(get_data_line(parameter,"system_filename",System_file))
        printf("System_file:          %s\n",  System_file);
    if(get_data_line(parameter,"r0_filename",R0_file))
        printf("R0_file:              %s\n",  R0_file);
    if(get_data_line(parameter,"maxit",&maxit))
        printf("maxit:                %d\n",  maxit);
    if(get_data_line(parameter,"dt",&dt))
        printf("timestep dt:          %lf\n", dt);
    if(get_data_line(parameter,"temperature",&Temperature))
        printf("Temperature:          %lf\n", Temperature);
    if(get_data_line(parameter,"NHT",&NHT))
        printf("NHT:                  %d\n", NHT);
    if(get_data_line(parameter,"tau",&tau))
        printf("NHT-parameter tau:    %lf\n", tau);
    if(get_data_line(parameter,"boxlength",&L))
        printf("L:                    %lf\n", L);
    if(get_data_line(parameter,"etaL",&etaL))
        printf("etaL: (ewald-param.)  %lf\n", etaL);
    if(get_data_line(parameter,"shake_crit",&crit))
        printf("crit: (SHAKE-dx)      %le\n", crit);
    if(get_data_line(parameter,"k_max",&K_max))
        printf("KMAX:                 %d\n",  K_max);
    if(get_data_line(parameter,"ksq_max",&Ksq_max))
        printf("KSQMAX:               %d\n",  Ksq_max);
    if(get_data_line(parameter,"cutoff",&cutoff))
        printf("cutoff:               %lf\n", cutoff);
    if(get_data_line(parameter,"far_ctf",&far_ctof))
        printf("far_ctof:             %lf\n", far_ctof);
    if(get_data_line(parameter,"random_seed",&r_seed))
        printf("random-seed:          %d\n",  r_seed);
    if(get_data_line(parameter,"write_step",&write_therm_step))
        printf("write_therm_step:     %d\n",  write_therm_step);
    if(get_data_line(parameter,"angdf_step",&angdf_step))
        printf("angdf_step:           %d\n",  angdf_step);
    if(get_data_line(parameter,"dhadf_step",&dhadf_step))
        printf("dhadf_step:           %d\n",  dhadf_step); 
    if(get_data_line(parameter,"idhdf_step",&idhdf_step))
        printf("idhdf_step:           %d\n",  idhdf_step);
    if(get_data_line(parameter,"pcf_step",&pcf_step))
        printf("pcf_step:             %d\n",  pcf_step);
    if(get_data_line(parameter,"end_thermostat",&end_thst))
        printf("end_thst:             %d\n",  end_thst);
    if(get_data_line(parameter,"begin_averages",&begin_avgs))
        printf("begin_avgs:           %d\n",  begin_avgs);
    if(get_data_line(parameter,"EWALD",&EWALD))
        printf("EWALD:                %d\n",  EWALD);
    if(get_data_line(parameter,"NL_SAI",&NL_SAI))
        printf("NL_SAI:               %d\n", NL_SAI);
    if(get_data_line(parameter,"NL_CI",&NL_CI))
        printf("NL_CI:                %d\n", NL_CI); 
    if(get_data_line(parameter,"nl_interval",&nl_intv))
        printf("nl_intv:              %d\n", nl_intv);
    if(get_data_line(parameter,"phase_sep_interval",&sep_intv))
        printf("sep_intv:             %d\n", sep_intv);
    if(get_data_line(parameter,"SP_ADIST",&SP_ADIST))
        printf("SP_ADIST:             %d\n", SP_ADIST);
    if(get_data_line(parameter,"NR_ADIST",&NR_ADIST))
        printf("NR_ADIST:             %d\n", NR_ADIST);
    if(get_data_line(parameter,"SP_DHADIST",&SP_DHADIST))
        printf("SP_DHADIST:           %d\n", SP_DHADIST);
    if(get_data_line(parameter,"NR_DHADIST",&NR_DHADIST))
        printf("NR_DHADIST:           %d\n", NR_DHADIST); 
    if(get_data_line(parameter,"SP_IDHADIST",&SP_IDHADIST))
        printf("SP_IDHADIST:          %d\n", SP_IDHADIST);
    if(get_data_line(parameter,"NR_IDHADIST",&NR_IDHADIST))
        printf("NR_IDHADIST:          %d\n", NR_IDHADIST);
    if((Nr_pcfs = readmore_4_int(parameter,"PCF",NR_PCFSMAX,
                                pcf_sp1,pcf_si1,pcf_sp2,pcf_si2)))
        printf("Nr_pcfs =             %d", Nr_pcfs);

    fclose(parameter);

    printf("\nusing defaults for the remaining parameters ...\n");
}

