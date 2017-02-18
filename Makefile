.PHONY = default, valgrind

# Define all object files needed to compile the main test executable
objectfiles = functions.o test.o

# Define the additional flags used to configure the compiler
# Comment this to speed up both compilation and execution
CFLAGS = -g -O0 -Wall

# -Wall: All warnings have been turned on
# -g: Debugging info is added
# -O0: NO optimization. Compiles and runs slower but allows valgrind to give
# more accurate info


# the default target is to execute both the C and python test scripts
default: test.out test.py
	./test.out
	python3 test.py

# Compiling the source code in to an executable is a two step process
# The executable is created from the object files and the relevant header files.
# In this process the linker must be told where to find the external functions
# that need to be linked in
# This is done using the -lgsl, -lgslcblas (basic linear algebra sub-routines),
# and -lm (system math subroutines) flags
# Note: In some systems the library MUST be linked at the end after the object
# file
# Note: On some systems the env variable LD_LIBRARY_PATH must be set to the
# library path (On the ultrabook this has been done in the fish config file)
test.out: $(objectfiles)
	gcc $(CFLAGS) $(objectfiles) functions.h -o test.out -lgsl -lgslcblas -lm

# The object file is created simply. In this step only the header files are
# used to get the prototypes for the external functions. The linker links them
# in the next step
# Note: The use of patterns to declare that any .o files is created as described
# from the corresponding .c file which is referred to using the magic variable
# $<
%.o: %.c functions.h
	gcc $(CFLAGS) -c $<

# test.py imports from libfunctions.so but is NOT created from it. So we declare
# the dependancy but don't provide a rule for building test.py
test.py: libfunctions.so

libfunctions.so: functions.c functions.h
	gcc -fPIC -shared functions.c -o libfunctions.so -lgsl -lgslcblas -lm

# Use valgrind to test the program
check:
	valgrind --leak-check=yes ./test.out

clean:
	rm -f *.o *.so *.gch test
