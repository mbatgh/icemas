#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


/* this subroutine calculates a measure for the mean-separation of
polar and apolar groups resp. ......................................... */


/************************************************************************/
void phase_sep(Real* m_poldis,Real* m_apoldis,Real* m_mixdis) {
/************************************************************************/

    int       i,j,wwi,wwj,n_pol,n_mix,n_apol,n_tot,
              l,firstone,lastone;
    Real      distance,a_dist,m_dist,p_dist,meansep;
    VEKTOR    dr;

    n_pol = n_mix = n_apol = 0;
    a_dist = m_dist = p_dist = 0.0;

    for(i=0;i<Nr_Stot-1;i++) {

        wwi = ww[i];

        for(j=nextidx[i];j<Nr_Stot;j++) {

            dr = RS[i]-RS[j];
            convolute(dr,adLh,L);
            distance = betr(dr);

            wwj = ww[j];
            
            switch(wwi+wwj) {
            
                case 2: n_apol++;
                        a_dist += distance; 
                        break;
                case 3: 
                case 5: n_mix++;
                        m_dist += distance; 
                        break;
                case 4:
                case 6:
                case 8: n_pol++;
                        p_dist += distance;
            }
        }
    }

    n_tot = n_pol + n_mix + n_apol;
    if(n_pol==0) n_pol = 1;
    if(n_mix==0) n_mix = 1;
    if(n_apol==0) n_apol = 1;
    
    meansep = (a_dist+m_dist+p_dist)/(Real)n_tot;
    *m_poldis = p_dist/(Real)n_pol/meansep*100.0;
    *m_apoldis = a_dist/(Real)n_apol/meansep*100.0;
    *m_mixdis = m_dist/(Real)n_mix/meansep*100.0;
}
