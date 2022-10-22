
using namespace std;

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <ctype.h>


#define     NLo         6.022045e23
#define     kB          1.380662e-23
#define     Rgas        8.3143
#define     amu         1.6605655e-27
#define     epsilon_0   8.8542e-12
#define     el_Charge   1.6021e-19
#define     PI          3.1415
#define     Real            double

#define     KTOT_MAX        1000
#define     NR_SMAX         2000
#define     NR_SPECMAX      5
#define     NR_NEIMAX       500000
#define     NR_LJTMAX       5
#define     NR_CMAX         2000
#define     NR_SPMMAX       100
#define     NR_PCFSMAX      10
#define     NR_PCF_CHANS    500

#define     LJ          1
#define     LJEL        2
#define     EL          4

#define     max(a,b)   ((a>b) ? a : b)
#define     min(a,b)   ((a<b) ? a : b)
#define     SMALL      (1.0e-10)
#define     BIG        (1.0e+30)
#define     RFACTOR    (1.0e+10)
