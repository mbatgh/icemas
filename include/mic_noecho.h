#ifndef STDIOH
#include <stdio.h>
#define STDIOH
#endif
#ifndef TERMIOSH
#include <termios.h>
#define TERMIOSH
#endif
#ifndef STDLIBH
#include <stdlib.h>
#define STDLIBH
#endif
#ifndef UNISTDH
#include <unistd.h>
#define UNISTDH
#endif
#ifndef IOSTREAMH
#include <iostream>
#define IOSTREAMH
#endif


struct termios saved_attributes;

void
reset_input_mode(void) {
    tcsetattr (STDIN_FILENO,TCSANOW,&saved_attributes);
    }

void 
set_input_mode(void) {
    struct termios tattr;
    char* name;

    if(!isatty(STDIN_FILENO)) {
        fprintf(stderr,"\nNot a terminal ...");
        exit(EXIT_FAILURE);
    }

    tcgetattr (STDIN_FILENO, &saved_attributes);
    atexit(reset_input_mode);

    tcgetattr(STDIN_FILENO, &tattr);
    tattr.c_lflag &= ~(ICANON|ECHO);
    tattr.c_cc[VMIN] = 1;
    tattr.c_cc[VTIME] = 0;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &tattr);
}

int ask_stop()
    {
    char c;
    set_input_mode();
    
    while(1)
        {
        cout << "proceed ? (y/n)" << endl;
        read(STDIN_FILENO,&c,1);
        switch(c)
            {
            case 'y': reset_input_mode(); return 0;
            case 'n': reset_input_mode(); return 1;
            default: cout << "y or n !" << endl; break;
            }
        }

    return EXIT_FAILURE;
    }
    
void presskey(void)
    {
    char c;
    set_input_mode();
    cout << "press any key to proceed ..." << endl;
    read(STDIN_FILENO,&c,1);
    reset_input_mode();
    }    
    
void presskey(char* message)
    {
    char c;
    set_input_mode();
    cout << message << endl;
    read(STDIN_FILENO,&c,1);
    reset_input_mode();
    }   
    
char mic_getchar(char* text) 
    {
    char c;
    set_input_mode();
    cout << text << endl;
    read(STDIN_FILENO,&c,1);
    reset_input_mode();
    return(c);
    }    
