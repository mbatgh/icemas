#include "vect_ops.h"
#include "structures.h"
#include "icemas.h"

int get_data_line(FILE*,char*,char*);
int get_data_line(FILE*,char*,float*);
int get_data_line(FILE*,char*,double*);
int get_data_line(FILE*,char*,int*);


/*******************************************************************/
int get_data_line(FILE *str, char* word, char* text) {
/*******************************************************************/

    int i,tot_len,GO_ON,FOUND;
    size_t length;
    char tmp[200];
    char* tptr;
    char* cptr;

    length = strlen(word);
    FOUND = 0;
    fseek(str,0,SEEK_SET);

    do {

        GO_ON = 1;
        i = 0;

/* read one line, i.e. terminate when reading line-feed or EOF ... */

        do {
            tmp[i] = fgetc(str);
            i++;
        } while((tmp[i-1]!=10)&&(tmp[i-1]!=EOF));

/* exit when reading EOF ......................................... */

        if(tmp[i-1]==EOF) GO_ON = 0;

/* search the line for word ...................................... */

        tmp[i-1] = '\0';
        tptr = (char*)tmp;

        do {

            if(strncmp(tptr,word,length)==0) {
                cptr = tptr+length+1;
                do {
                    if(isalnum((int)(*cptr))) {
                        sscanf(cptr,"%s",text);
                        FOUND = 1;
                    }
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
            }

            tptr++;

        } while((*(tptr-1) != '\0')&&(!FOUND));

    } while(GO_ON&&(!FOUND));

    return(FOUND);
}


/*******************************************************************/
int get_data_line(FILE *str, char* word, double* number) {
/*******************************************************************/

    int i,tot_len,GO_ON,FOUND;
    size_t length;
    char tmp[200];
    char* tptr;
    char* cptr;
    
    length = strlen(word);
    FOUND = 0;
    fseek(str,0,SEEK_SET);

    do {

        GO_ON = 1;
        i = 0;

/* read one line, i.e. terminate when reading line-feed or EOF ... */

        do {
            tmp[i] = fgetc(str);
            i++;
        } while((tmp[i-1]!=10)&&(tmp[i-1]!=EOF));

/* exit when reading EOF ......................................... */

        if(tmp[i-1]==EOF) GO_ON = 0;

/* search the line for word ...................................... */

        tmp[i-1] = '\0';
        tptr = (char*)tmp;
        
        do {

            if(strncmp(tptr,word,length)==0) {
                cptr = tptr+length+1;
                do {
                    if(isalnum((int)(*cptr))) {
                        sscanf(cptr,"%lf",number);
                        FOUND = 1;
                    }
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
            }

            tptr++;

        } while((*(tptr-1) != '\0')&&(!FOUND));

    } while(GO_ON&&(!FOUND));

    return(FOUND);
}


/*******************************************************************/
int get_data_line(FILE *str, char* word, float* number) {
/*******************************************************************/

    int i,tot_len,GO_ON,FOUND;
    size_t length;
    char tmp[200];
    char* tptr;
    char* cptr;
    
    length = strlen(word);
    FOUND = 0;
    fseek(str,0,SEEK_SET);

    do {

        GO_ON = 1;
        i = 0;

/* read one line, i.e. terminate when reading line-feed or EOF ... */

        do {
            tmp[i] = fgetc(str);
            i++;
        } while((tmp[i-1]!=10)&&(tmp[i-1]!=EOF));

/* exit when reading EOF ......................................... */

        if(tmp[i-1]==EOF) GO_ON = 0;

/* search the line for word ...................................... */

        tmp[i-1] = '\0';
        tptr = (char*)tmp;
        
        do {
    
            if(strncmp(tptr,word,length)==0) {
                cptr = tptr+length+1;
                do {
                    if(isalnum((int)(*cptr))) {
                        sscanf(cptr,"%f",number);
                        FOUND = 1;
                    }
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
            }

            tptr++;

        } while((*(tptr-1) != '\0')&&(!FOUND));

    } while(GO_ON&&(!FOUND));

    return(FOUND);
}


/*******************************************************************/
int get_data_line(FILE *str, char* word, int* number) {
/*******************************************************************/

    int i,tot_len,GO_ON,FOUND;
    size_t length;
    char tmp[200];
    char* tptr;
    char* cptr;
    
    length = strlen(word);
    FOUND = 0;
    fseek(str,0,SEEK_SET);

    do {

        GO_ON = 1;
        i = 0;

/* read one line, i.e. terminate when reading line-feed or EOF ... */

        do {
            tmp[i] = fgetc(str);
            i++;
        } while((tmp[i-1]!=10)&&(tmp[i-1]!=EOF));

/* exit when reading EOF ......................................... */

        if(tmp[i-1]==EOF) GO_ON = 0;

/* search the line for word ...................................... */

        tmp[i-1] = '\0';
        tptr = (char*)tmp;
        
        do {
    
            if(strncmp(tptr,word,length)==0) {
                cptr = tptr+length+1;
                do {
                    if(isalnum((int)(*cptr))) {
                        sscanf(cptr,"%d",number);
                        FOUND = 1;
                    }
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
            }

            tptr++;

        } while((*(tptr-1) != '\0')&&(!FOUND));

    } while(GO_ON&&(!FOUND));

    return(FOUND);
}


/*******************************************************************/
int readmore_4_int(FILE *str, char* word,int max,int* sp_1,int* si_1,int* sp_2,int* si_2) {
/*******************************************************************/

    int i,j,tot_len,GO_ON,LINES,FOUND;
    size_t length;
    char tmp[200];
    char* tptr;
    char* cptr;
    int *sp1 = sp_1;
    int *sp2 = sp_2;
    int *si1 = si_1;
    int *si2 = si_2;

    length = strlen(word);
    fseek(str,0,SEEK_SET);
    LINES = 0;

    do {

        GO_ON = 1;

/* read one line, i.e. terminate when reading line-feed or EOF ... */

        i = 0;

        do {
            tmp[i] = fgetc(str);
            i++;
        } while((tmp[i-1]!=10)&&(tmp[i-1]!=EOF));

/* exit when reading EOF ......................................... */

        if(tmp[i-1]==EOF)
            GO_ON = 0;

/* search the line for word ...................................... */

        tmp[i-1] = '\0';
        tptr = (char*)tmp;

        do {

            if(strncmp(tptr,word,length)==0) {
            
                cptr = tptr+length+1;
                do {
                    FOUND = 0;
                    if(isalnum((int)(*cptr)))
                        FOUND = sscanf(cptr,"%d",sp1);
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
                do {
                    FOUND = 0;
                    if(isalnum((int)(*cptr)))
                        FOUND = sscanf(cptr,"%d",si1);
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
                do {
                    FOUND = 0;
                    if(isalnum((int)(*cptr)))
                        FOUND = sscanf(cptr,"%d",sp2);
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
                do {
                    FOUND = 0;
                    if(isalnum((int)(*cptr)))
                        FOUND = sscanf(cptr,"%d",si2);
                    cptr++;
                } while( (*(cptr-1) != '\0') && (!FOUND) );
   
//printf("LINES = %d, sp1=%d, sp2=%d, si1=%d, si2=%d\n", LINES,
//                     sp_1[LINES],sp_2[LINES],si_1[LINES],si_2[LINES]);
//presskey(" ...");

            LINES++;
            sp1++; sp2++; si1++; si2++;
            }

            tptr++;
            
            if(LINES==max) {
                GO_ON=0;
                break;
            }

        } while(*(tptr-1) != '\0');       

    } while(GO_ON);

    return(LINES);
}

