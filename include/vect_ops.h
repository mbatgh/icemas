#include "first.h"

typedef struct dob_vektor 	VEKTOR;
typedef struct int_vektor   INTVEK;

struct dob_vektor
{
    double      x,y,z;
};

struct int_vektor
{
    int         x,y,z; 
};

inline VEKTOR operator+(const VEKTOR& a, const VEKTOR& b) {
   return((VEKTOR){a.x+b.x,a.y+b.y,a.z+b.z});
} 
inline VEKTOR operator-(const VEKTOR& a, const VEKTOR& b) {
   return((VEKTOR){(a.x-b.x),(a.y-b.y),(a.z-b.z)});
}
inline VEKTOR operator-(const VEKTOR& a, const INTVEK& b){
   return((VEKTOR){(a.x-(Real)b.x),(a.y-(Real)b.y),(a.z-(Real)b.z)});
}
inline Real operator*(const VEKTOR& a, const VEKTOR& b) {
    return(a.x*b.x+a.y*b.y+a.z*b.z);
}
inline VEKTOR operator*(const Real& s,const VEKTOR& a) {
    return((VEKTOR){a.x*s,a.y*s,a.z*s});
}
inline VEKTOR operator*(const Real& s, const INTVEK& a) {
    return((VEKTOR){(Real)a.x*s,(Real)a.y*s,(Real)a.z*s});
}
inline VEKTOR operator*(const VEKTOR& a,const Real& s) {
    return((VEKTOR){a.x*s,a.y*s,a.z*s});
}
inline VEKTOR operator*(const INTVEK& a,const Real& s) {
    return((VEKTOR){(Real)a.x*s,(Real)a.y*s,(Real)a.z*s});
}
inline VEKTOR operator/(const VEKTOR& a,const Real& s) {
    return((VEKTOR){a.x/s,a.y/s,a.z/s});
}
inline VEKTOR operator%(const VEKTOR& a,const VEKTOR& b) {
    return((VEKTOR){(a.y*b.z-a.z*b.y),(a.x*b.z-a.z*b.x),(a.x*b.y-a.y*b.x)});
}
inline Real betr(const VEKTOR& a) { 
    return(sqrt(a.x*a.x+a.y*a.y+a.z*a.z)); 
}
inline Real norm(const VEKTOR& a) { 
    return(a.x*a.x+a.y*a.y+a.z*a.z); 
}
inline void convolute(VEKTOR& a,const Real& adxh, const Real& x) {
   	a.x = a.x - ((int)(a.x*adxh)) * x;
    a.y = a.y - ((int)(a.y*adxh)) * x;
    a.z = a.z - ((int)(a.z*adxh)) * x;
}
inline void printvek(VEKTOR& v) {
    printf("(%le, %le, %le)",v.x,v.y,v.z);
}

