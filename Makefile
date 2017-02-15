.PHONY = default, valgrind

# Define the additional flags used to configure the compiler
# Comment this to speed up both compilation and execution
CFLAGS = -g -O0 -Wall

# -Wall: All warnings have been turned on
# -g: Debugging info is added
# -O0: NO optimization. Compiles and runs slower but allows valgrind to give
# more accurate info


# Compiling the source code in to an executable is a two step process
default: test
	./test

# The executable is created from the object file. In this process the linker
# must be told where to find the external functions that need to be linked in
# This is done using the -lgsl, -lgslcblas (basic linear algebra sub-routines),
# and -lm (system math subroutines) flags
test: functions.o
	gcc $(CFLAGS) -lgsl -lgslcblas -lm functions.o -o test

# The object file is created simply. In this step only the header files are
# used to get the prototypes for the external functions. The linker links them
# in the next step
functions.o: functions.c
	gcc $(CFLAGS) -c functions.c

# Use valgrind to test the program
valgrind:
	valgrind --leak-check=yes ./test
