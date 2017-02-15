.PHONY = default

# Compiling the source code in to an executable is a two step process
default: test
	./test

# The executable is created from the object file. In this process the linker
# must be told where to find the external functions that need to be linked in
# This is done using the -lgsl, -lgslcblas (basic linear algebra sub-routines),
# and -lm (system math subroutines) flags
test: functions.o
	gcc -o test -lgsl -lgslcblas -lm functions.o

# The object file is created simply. In this step only the header files are
# used to get the prototypes for the external functions. The linker links them
# in the next step
functions.o: functions.c
	gcc -c functions.c
