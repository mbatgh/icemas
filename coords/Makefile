PROG = 		mkr0
CC = 		g++
OBJS =		mkr0.o

all:$(PROG)

$(PROG):	$(OBJS)
	$(CC) -o $(PROG) $(OBJS)

mkr0.o:			$(HEADER) mkr0.cc
	$(CC) $(CFLAGS) -c mkr0.cc
