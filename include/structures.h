/* definition of some structures */


typedef struct str_site     t_site;
typedef struct str_bond     t_bond;
typedef struct str_angle    t_angle;
typedef struct str_dieder   t_dieder;
typedef struct str_idieder  t_idieder;
typedef struct str_ljtype   t_ljtype;

extern              t_site*    Site;
extern              t_bond*    Bond;
extern              t_angle*   Angle;
extern              t_dieder*  Dieder;
extern              t_idieder* iDieder;
extern              t_ljtype   LJtypes[NR_LJTMAX];

struct str_site {
    int number;
    int ljtype;
    Real mass;
    Real charge;
};

struct str_bond {
    int number;
    int idx1;
    int idx2;
    Real b_eqi;
};

struct str_angle {
    int number;
    int idx1;
    int idx2;
    int idx3;
    Real a_eqi;
    Real force;
};

struct str_dieder {
    int number;
    int idx1;
    int idx2;
    int idx3;
    int idx4;
    int period;
    Real d_eqi;
    Real force;
};

struct str_idieder {
    int number;
    int idx1;
    int idx2;
    int idx3;
    int idx4;
    Real d_eqi;
    Real force;
};

struct str_ljtype {
    int number;
    Real epsilon;
    Real sigma;
};
