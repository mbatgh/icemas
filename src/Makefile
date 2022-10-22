PROG = 		icemas
CC = 		g++
CFLAGS = 	-O2 -finline-functions
HEADER = 	../include/first.h ../include/vect_ops.h \
		../include/structures.h ../include/icemas.h 
OBJS =		force.o angle_force.o main.o setup.o make_data_file.o \
                get_data_line.o make_k.o tdf.o get_parameter.o \
                trash.o get_r0.o update_dist.o phase_sep.o \
                get_system.o move.o write_lastr.o hom_neighbors.o \
                het_neighbors.o write_means.o write_r.o r_force.o k_force.o \
                torsion_force.o distributions.o impdh_force.o \
                readlin.o nl_force.o nl_r_force.o
INCLUDES =	-I../include

all:$(PROG)

$(PROG):	$(OBJS)
	$(CC) -o $(PROG) $(OBJS)

main.o:			$(HEADER) main.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c main.cc

angle_force.o:		$(HEADER) angle_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c angle_force.cc

torsion_force.o:	$(HEADER) torsion_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c torsion_force.cc

k_force.o: 		$(HEADER) k_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c k_force.cc

force.o:		$(HEADER) force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c force.cc

tdf.o:			$(HEADER) tdf.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c tdf.cc

r_force.o:		$(HEADER) r_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c r_force.cc

nl_force.o:		$(HEADER) nl_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c nl_force.cc

nl_r_force.o:		$(HEADER) nl_r_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c nl_r_force.cc

trash.o:		$(HEADER) trash.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c trash.cc

move.o:			$(HEADER) move.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c move.cc

setup.o:		$(HEADER) setup.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c setup.cc

get_parameter.o:	$(HEADER) get_parameter.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c get_parameter.cc

get_data_line.o:	$(HEADER) get_data_line.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c get_data_line.cc

get_r0.o:		$(HEADER) get_r0.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c get_r0.cc

make_data_file.o:	$(HEADER) make_data_file.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c make_data_file.cc

make_k.o:		$(HEADER) make_k.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c make_k.cc

hom_neighbors.o:	$(HEADER) hom_neighbors.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c hom_neighbors.cc

het_neighbors.o:	$(HEADER) het_neighbors.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c het_neighbors.cc

get_system.o:		$(HEADER) get_system.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c get_system.cc

update_dist.o: 		$(HEADER) update_dist.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c update_dist.cc

write_lastr.o:         	$(HEADER) write_lastr.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c write_lastr.cc

write_means.o:         	$(HEADER) write_means.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c write_means.cc

write_r.o:		$(HEADER) write_r.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c write_r.cc

distributions.o:	$(HEADER) distributions.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c distributions.cc

phase_sep.o:		$(HEADER) phase_sep.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c phase_sep.cc

impdh_force.o:		$(HEADER) impdh_force.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c impdh_force.cc

readlin.o:		$(HEADER) readlin.cc
	$(CC) $(INCLUDES) $(CFLAGS) -c readlin.cc

clean:
	rm -f *.o *.bak *~ $(PROG)

install:	$(PROG)
	mv $(PROG) /usr/local/bin

