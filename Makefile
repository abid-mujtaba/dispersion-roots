.PHONY = plot, check, test, data

# Define all object files needed to compile the main test executable
objectfiles = dispersion.o henning.o roots.o hypergeom.o constants.o first.o second.o third.o
headerfiles = dispersion.h henning.h roots.h hypergeom.h constants.h
libraries = -lgsl -lgslcblas -lm -lgmp -lmpfr

# Define the additional flags used to configure the compiler
# Comment this to speed up both compilation and execution
CFLAGS = -g -O0 -Wall -std=gnu11 -pedantic
# CFLAGS = -O3

# -Wall: All warnings have been turned on
# -g: Debugging info is added
# -O0: NO optimization. Compiles and runs slower but allows valgrind to give
# more accurate info
# -O3: Maximum optimization


plot: functions.py
	python3 plots.py

test: test.out
	./test.out

# Compiling the source code in to an executable is a two step process
# The executable is created from the object files.
# In this process the linker must be told where to find the external functions
# (present in the object files) that need to be linked in
# This is done using the -lgsl, -lgslcblas (basic linear algebra sub-routines),
# -lm (system math subroutines), -lgmp (GNU multiple precision), and -lmprf (GNU
# Multiple Precision Floating-point Reliably) flags stored in the $(libraries)
# variable
# Note: In some systems the library MUST be linked at the end after the object
# file
# Note: On some systems the env variable LD_LIBRARY_PATH must be set to the
# library path (On the ultrabook this has been done in the fish config file)
#
# Note: $@ is the automatic variable that translates to the target name, in this
# case test.out
test.out: $(objectfiles) test.o
	gcc $(CFLAGS) $(objectfiles) test.o -o $@ $(libraries)


# Create root data for plotting
data: data.out
	./data.out

data.out: $(objectfiles) data.o
	gcc $(CFLAGS) $(objectfiles) data.o -o $@ $(libraries)


# The object file is created simply. In this step only the header files are
# used to get the prototypes for the external functions. The linker links them
# in the next step
# Note: The use of patterns to declare that any .o files is created as described
# from the corresponding .c file which is referred to using the magic (automatic)
# variable $<
%.o: %.c $(headerfiles)
	gcc $(CFLAGS) -c $<

# test.py imports from libfunctions.so but is NOT created from it. So we declare
# the dependancy but don't provide a rule for building test.py
functions.py: libDroots.so

# The shared object (dynamic library) is created from the relevant c files
libDroots.so: roots.c dispersion.c henning.c hypergeom.c constants.c $(headerfiles)
	gcc $(CFLAGS) -fPIC -shared roots.c dispersion.c henning.c hypergeom.c constants.c -o $@ $(libraries)

# Use valgrind to test the program
check: test.out
	valgrind --leak-check=yes ./test.out

clean:
	rm -f *.o *.so *.gch *.out
