#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"



/****************************************************************/
int tot_index(int sp,int mo,int si) {
/****************************************************************/

    int i,j,k,l,tot;
    
    tot = 0;

    for(i=0;i<sp;i++)
        tot += Nr_Molecules[i]*Nr_Sites[i];
    for(i=0;i<mo;i++)
        tot += Nr_Sites[sp];
    tot += si;

    return(tot);
}


/****************************************************************/   
VEKTOR cent_of_mass(VEKTOR* v, int i0, int in) {
/****************************************************************/

    int i,j;
    Real M;
    VEKTOR sum,dpbc;
    VEKTOR bindv[NR_SPMMAX];
    VEKTOR rb[NR_SPMMAX];

    sum = nulvek;
    rb[0] = v[i0];

    for(j=0,i=i0;i<in;j++,i++) {
        bindv[j] = v[i]-v[i+1];
        convolute(bindv[j],adLh,L);
        rb[j+1] = rb[j] + bindv[j];
        M += mass[j];
    }

    for(j=0,i=i0;i<in;j++,i++)
        sum = sum + rb[j] * (mass[j]/M);
    
    convolute(sum,adL34,L);

    return(sum);
}


/****************************************************************/
int both_in(int si,int sj,t_bond B) {
/****************************************************************/

    if(si==B.idx1&&sj==B.idx2||si==B.idx2&&sj==B.idx1)
        return(1);
    else return(0);
}

/****************************************************************/
int both_in(int si,int sj,t_angle A) {
/****************************************************************/

    int i,j;

    if(si==A.idx1||si==A.idx2||si==A.idx3)
        if(sj==A.idx1||sj==A.idx2||sj==A.idx3) return(1);

    return(0);
}

/****************************************************************/
int both_in(int si,int sj,t_dieder D) {
/****************************************************************/

    int i,j;

    if(si==D.idx1||si==D.idx2||si==D.idx3||si==D.idx4)
        if(sj==D.idx1||sj==D.idx2||sj==D.idx3||sj==D.idx4) return(1);

    return(0);
}

/****************************************************************/
int both_in(int si,int sj,t_idieder I) {
/****************************************************************/

    int i,j;

    if(si==I.idx1||si==I.idx2||si==I.idx3||si==I.idx4)
        if(sj==I.idx1||sj==I.idx2||sj==I.idx3||sj==I.idx4) return(1);

    return(0);
}


/****************************************************************/
void        zero_forces() {
/****************************************************************/

    int i,j;

/* Kraefte 0 setzen ........................................... */

    for(i=0;i<Nr_Stot;i++) {
        FS[i] = nulvek;
    }
}


/**************************************************************/
void scale_tmp(Real ist ,Real soll) {
/**************************************************************/

    VEKTOR v;
    Real scal = sqrt(soll/ist);

    for(int j=0;j<Nr_Stot;j++) {
        v = RS[j]-RSo[j];
        RS[j] = RSo[j] + v*scal;
    }
}


/**************************************************************/
void scale_new_tmp(Real ist, Real soll) {
/**************************************************************/

    VEKTOR v;
    Real scal = sqrt(soll/ist);

    for(int j=0;j<Nr_Stot;j++) {
        v = RSn[j]-RS[j];
        RSn[j] = RS[j] + v*scal;
    }
}

