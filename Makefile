PROG =	minor-gravity

SRCS =	command.for dhms.for get_pid.for minor-gravity.for timestamp2.for

OBJS =	command.o dhms.o get_pid.o minor-gravity.o timestamp2.o

LIBS = -lm	

CC = gcc
CFLAGS = -O
FC = gfortran
#FC = ifort
FFLAGS = -O -mcmodel=medium
F90 = gfortran
#F90 = ifort
F90FLAGS = -O
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.d

.SUFFIXES: $(SUFFIXES) .f .for

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

.for.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

minor-gravity.o: polyhedron.inc polyhedron.inc polyhedron.inc polyhedron.inc \
	polyhedron.inc polyhedron.inc polyhedron.inc polyhedron.inc \
	polyhedron.inc polyhedron.inc polyhedron.inc polyhedron.inc \
	polyhedron.inc polyhedron.inc polyhedron.inc polyhedron.inc \
	polyhedron.inc polyhedron.inc polyhedron.inc polyhedron.inc \
	polyhedron.inc polyhedron.inc
