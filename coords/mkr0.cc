#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#define     NMAX        2000
#define     NSMAX       10

void error_msg(void);


int main(int ARGC, char* ARGV[]) {

    int N,NS,i,j,k,l,COORDS,nb,n,ok;
    double L,v,x,y,z,xe,ye,ze,x0,y0,z0,si,
           dx,dy,dz,sx,sy,sz,NpL,min,dis;
    long t;
    FILE *in, *out;

    double RX[NMAX][NSMAX];
    double RY[NMAX][NSMAX];
    double RZ[NMAX][NSMAX];
    double RSx[NSMAX];
    double RSy[NSMAX];
    double RSz[NSMAX];

    COORDS = 0;
    v = 0.0;    
    NS = 1;
    ok = 0;

    for(i=0;i<NSMAX;i++) {
        RSx[i] = 0.0;
        RSy[i] = 0.0;
        RSz[i] = 0.0;
    }

    srandom(time(&t));

    for(i=1;i<ARGC;i+=2) {
        
        if(ARGV[i][0]!='-') error_msg();
        else switch(ARGV[i][1]) {
        
            case 'n': N = atoi(ARGV[i+1]);
                      ok++;
                      break;
            case 'l': L = atof(ARGV[i+1]);
                      ok++;
                      break;
            case 'v': v = atof(ARGV[i+1]);
                      break;
            case 'f': COORDS = i+1;
                      break;
            default: error_msg();
        }
    }

    if(ok<2) error_msg();

    if((out=fopen("r0","w"))==NULL) {
        printf("can't open r0-file ...\n");
        exit(1);
    }

    if(COORDS) {
    
        if((in=fopen(ARGV[COORDS],"r"))==NULL) {
            printf("can't see no %s\n", ARGV[COORDS]);
            exit(1);
        }
    
        NS = -1;

        do {
            i = fscanf(in,"%lf %lf %lf",&x,&y,&z);
            NS++;
        } while(i==3);

        if(NS>100) error_msg();

        rewind(in);

        for(i=0;i<NS;i++) {
            fscanf(in,"%lf %lf %lf",&x,&y,&z);
            RSx[i] = x;
            RSy[i] = y;
            RSz[i] = z;
        }
    }

    NpL = pow((double)N,(1.0/3.0));

    nb = 0;
    
    xe = L*0.5;
    ye = L*0.5;
    ze = L*0.5;
    x0 = -L*0.5;
    y0 = -L*0.5;
    z0 = -L*0.5;

    dx = L/NpL;
    dy = L/NpL/NpL;
    dz = L/NpL/NpL/NpL;

    x = x0;
    y = y0;
    z = z0; 

    si = 1.0;

    for(i=0;i<N;i++) {

        x += dx;//*(double)(random()/RAND_MAX);
        y += dy;//*(double)(random()/RAND_MAX);
        z += dz;//*(double)(random()/RAND_MAX);

        if(x>L*0.5) x -= L;
        if(y>L*0.5) y -= L;

            for(j=0;j<NS;j++) {
                
                sx = x+RSx[j];
                sy = y+RSy[j];
                sz = z+RSz[j];

                if(sx>L*0.5) sx -= L;
                if(sx<-L*0.5) sx += L;
                if(sy>L*0.5) sy -= L;
                if(sy<-L*0.5) sy += L;
                if(sz>L*0.5) sz -= L;
                if(sz<-L*0.5) sz += L;

                nb += fprintf(out,"%lf %lf %lf %lf %lf %lf\n",
                        sx,sy,sz,sx+v*si,sy+v*si,sz+v*si);

                RX[i][j] = sx;
                RY[i][j] = sy;
                RZ[i][j] = sz;
            }
    
        si = -si;
    }

    min = 1.0e+30;

    for(i=0;i<N-1;i++)
        for(j=i+1;j<N;j++)
                for(k=0;k<NS;k++)
                    for(l=0;l<NS;l++) { 

            dx = RX[i][k]-RX[j][l];
            dy = RY[i][k]-RY[j][l];
            dz = RZ[i][k]-RZ[j][l];

                if(dx>L*0.5) dx -= L;
                if(dx<-L*0.5) dx += L;
                if(dy>L*0.5) dy -= L;
                if(dy<-L*0.5) dy += L;
                if(dz>L*0.5) dz -= L;
                if(dz<-L*0.5) dz += L;

           dis = sqrt(dx*dx+dy*dy+dz*dz);
           if(dis<min) min = dis;
    }

    printf("\n%d bytes written in r0, \n",nb);
    printf("minimum distance = %lf\n\n",min);

    fclose(out);

    if(COORDS) fclose(in);
}



void error_msg(void) {

    printf("\nusage: mkr0 -n number -l length [-v dx][-f file]\n\n");
    printf("* number is the total number of molecules\n");
    printf("  (maximum number of molecules: 2000)\n");
    printf("* length [nm] is the side-length of the cubic box\n");
    printf("* dx [nm] is the displacement of the old coords, which\n");
    printf("  determines the average initial-velocity,\n");
    printf("  (default is zero)\n");
    printf("* file is a file, containing information on the\n");
    printf("  geometry [nm] of a single molecule\n"); 
    printf("  (maximum number of sites per molecule: 10)\n");
    printf("  if no file is given, the molecule is assumed to\n"); 
    printf("  be a single atom\n\n");

    exit(1);
}

