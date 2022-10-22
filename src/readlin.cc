#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"


void    readline(FILE*,char*);
void    readline(FILE*,double*);
void    readline(FILE*,float*);
void    readline(FILE*,int*);

void    readline(FILE*,int*,int*,int*,int*,int*,int*);
void    readline(FILE*,int*,int*,int*,int*);
void    readline(FILE*,int*,int*,int*,int*,int*,int*,double*,double*);
void    readline(FILE*,int*,int*,int*,int*,int*,double*,double*);
void    readline(FILE*,int*,int*,int*,int*,double*,double*);
void    readline(FILE*,int*,int*,int*,double*);
void    readline(FILE*,int*,int*,double*,double*);
void    readline(FILE*,int*,double*,double*);


/****************************************************************/
void readline(FILE *str, int* i1,int* i2,int* i3,int* i4, int* i5, int* i6) {
/****************************************************************/

    const int nn = 6;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k;
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%d",i4);
    sscanf((char*)(text[4]),"%d",i5);
    sscanf((char*)(text[5]),"%d",i6);
}


/****************************************************************/
void readline(FILE *str, int* i1,int* i2,int* i3,int* i4) {
/****************************************************************/

    const int nn = 4;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k;
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%d",i4);
}


/****************************************************************/
void 
readline(FILE *str,int* i1,int* i2,int* i3,int* i4,int* i5,int* i6,
         double* d1,double* d2) {
/****************************************************************/

    const int nn = 8;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k; //before[nw-1] + word[nw-1];
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%d",i4);
    sscanf((char*)(text[4]),"%d",i5);
    sscanf((char*)(text[5]),"%d",i6);
    sscanf((char*)(text[6]),"%lf",d1);
    sscanf((char*)(text[7]),"%lf",d2);
}


/****************************************************************/
void 
readline(FILE *str,int* i1,int* i2,int* i3,int* i4,int* i5,double* d1,double* d2) {
/****************************************************************/

    const int nn = 7;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }

    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k; //before[nw-1] + word[nw-1];
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%d",i4);
    sscanf((char*)(text[4]),"%d",i5);
    sscanf((char*)(text[5]),"%lf",d1);
    sscanf((char*)(text[6]),"%lf",d2);
}


/****************************************************************/
void readline(FILE *str,int* i1,int* i2,int* i3,int* i4,double* d1,double* d2) {
/****************************************************************/

    const int nn = 6;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k; //before[nw-1] + word[nw-1];
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%d",i4);
    sscanf((char*)(text[4]),"%lf",d1);
    sscanf((char*)(text[5]),"%lf",d2);
}


/****************************************************************/
void readline(FILE *str, int* i1,int* i2,int* i3, double* d1) {
/****************************************************************/

    const int nn = 4;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k;
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%d",i3);
    sscanf((char*)(text[3]),"%lf",d1);
}


/****************************************************************/
void readline(FILE *str,int* i1,int* i2,double* d1,double* d2) {
/****************************************************************/

    const int nn = 4;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k; //before[nw-1] + word[nw-1];
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%d",i2);
    sscanf((char*)(text[2]),"%lf",d1);
    sscanf((char*)(text[3]),"%lf",d2);
}


/****************************************************************/
void readline(FILE *str,int* i1,double* d1,double* d2) {
/****************************************************************/

    const int nn = 3;

    char temp[100];
    char text[nn][100];
    int nw,k,j,i,HALT,before[nn],word[nn];

    HALT = 0;
    
    for(j=0;j<100;j++) {
        temp[j] = 32;
        for(i=0;i<nn;i++)
            text[i][j] = 32;
    }
    
    temp[99] = '\0';
    for(i=0;i<nn;i++)
        text[i][99] = '\0';

    for(i=0;i<nn;i++) {
        before[i] = 0;
        word[i] = 0;
    }

    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);

    for(nw=0;nw<nn;nw++) {
        if(nw>0) before[nw] = k;
        while(temp[before[nw]]==32) before[nw]++;
        word[nw] = before[nw];
        while(isgraph(temp[word[nw]])) word[nw]++;
        for(k=before[nw];k<word[nw];k++) text[nw][k-before[nw]] = temp[k];
        text[nw][word[nw]-before[nw]] = '\0';
    }

    sscanf((char*)(text[0]),"%d",i1);
    sscanf((char*)(text[1]),"%lf",d1);
    sscanf((char*)(text[2]),"%lf",d2);
}

/****************************************************************/
void readline(FILE *str, char* text) {
/****************************************************************/

    char temp[100];
    int k,j,i,HALT,before,word;
    
    HALT = 0;
    before = 0;
    word = 0;
    
    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);
    while(temp[before]==32) before++;
    word = before;
    while(isgraph(temp[word])) word++;
    for(k=before;k<word;k++) text[k-before] = temp[k];
    text[word-before] = '\0';
}


/****************************************************************/
void readline(FILE *str, float* number) {
/****************************************************************/

    char temp[100];
    char text[100];
    int k,j,i,HALT,before,word;
    
    HALT = 0;
    before = 0;
    word = 0;
    
    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);
    while(temp[before]==32) before++;
    word = before;
    while(isgraph(temp[word])) word++;
    for(k=before;k<word;k++) text[k-before] = temp[k];
    text[word-before] = '\0';
    sscanf(text,"%f",number);
}


/****************************************************************/
void readline(FILE *str, double* number) {
/****************************************************************/

    char temp[100];
    char text[100];
    int k,j,i,HALT,before,word;
    
    HALT = 0;
    before = 0;
    word = 0;
    
    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);
    while(temp[before]==32) before++;
    word = before;
    while(isgraph(temp[word])) word++;
    for(k=before;k<word;k++) text[k-before] = temp[k];
    text[word-before] = '\0';
    sscanf(text,"%lf",number);
}


/****************************************************************/
void readline(FILE *str, int* number) {
/****************************************************************/

    char temp[100];
    char text[100];
    int k,j,i,HALT,before,word;
    
    HALT = 0;
    before = 0;
    word = 0;
    
    do {
        i = 0;
        do { 
            temp[i] = fgetc(str);
            i++; 
        } while((temp[i-1]!=35)&&(temp[i-1]!=10));

        if(temp[i-1]==35) while(fgetc(str)!=10);

        for(j=0;j<i-1;j++) {
            if(isgraph(temp[j]))
                HALT = 1;  
        }  
    } while(!HALT);
    while(temp[before]==32) before++;
    word = before;
    while(isgraph(temp[word])) word++;
    for(k=before;k<word;k++) text[k-before] = temp[k];
    text[word-before] = '\0';
    sscanf(text,"%d",number);
}

