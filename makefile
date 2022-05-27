# Makefile for ios2d
#
# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran

CC=gcc
FC=gfortran

FFLAGS= -fPIC -O3 -funroll-loops -std=legacy -w 
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS)
CPPFLAGS="-I/usr/local/opt/openblas/include" 

CLINK = -lgfortran -lm -ldl

LIBS = -lm

# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1   
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

# For your OS, override the above by placing make variables in make.inc
-include make.inc



# objects to compile
#
# Common objects

DIRN = src/fortran

# Helmholtz objects
OBJS = $(DIRN)/hank103.o $(DIRN)/prini.o $(DIRN)/helm_kernels.o  \
	$(DIRN)/formsysmatbac.o $(DIRN)/kern_mats.o $(DIRN)/lap_kernels.o \
	$(DIRN)/durdec.o $(DIRN)/corrand4.o $(DIRN)/dumb_conres.o \
	$(DIRN)/legeexps.o $(DIRN)/curve_filtering.o $(DIRN)/curve_resampler.o \
	$(DIRN)/dfft.o $(DIRN)/cdjseval2d.o $(DIRN)/h2dcommon.o


.PHONY: usage matlab 

default: usage 

all: examples matlab 

usage:
	@echo "Makefile for ios2d. Specify what to make:"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make mex - generate matlab interfaces (for expert users only, requires mwrap)"
	@echo "  make clean - also remove lib, MEX, and demo executables"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

LIBNAME = libhelmtrans
STATICLIB = $(LIBNAME).a

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

# matlab..
MWRAPFILE = kern_mats
MWRAPFILE2 = helm_kernels
MWRAPFILE3 = lap_kernels
MWRAPFILE4 = curve_resampler
GATEWAY = $(MWRAPFILE)
GATEWAY2 = $(MWRAPFILE2)
GATEWAY3 = $(MWRAPFILE3)
GATEWAY4 = $(MWRAPFILE4)

matlab:	$(STATICLIB) src/matlab/$(GATEWAY).c src/matlab/$(GATEWAY2).c src/matlab/$(GATEWAY3).c src/matlab/$(GATEWAY4).c
	$(MEX) src/matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) -output src/matlab/kern_mats $(MEXLIBS)
	$(MEX) src/matlab/$(GATEWAY2).c lib-static/$(STATICLIB) $(MFLAGS) -output src/matlab/helm_kernels $(MEXLIBS)
	$(MEX) src/matlab/$(GATEWAY3).c lib-static/$(STATICLIB) $(MFLAGS) -output src/matlab/lap_kernels $(MEXLIBS)
	$(MEX) src/matlab/$(GATEWAY4).c lib-static/$(STATICLIB) $(MFLAGS) -output src/matlab/curve_resampler $(MEXLIBS)

mex:  $(STATICLIB)
	cd src; cd matlab;  $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) -v $(GATEWAY).c ../../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) $(MEXLIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw;\
	$(MEX) -v $(GATEWAY2).c ../../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE2) $(MEXLIBS);\
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY3) -mb $(MWRAPFILE3).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY3) -c $(GATEWAY3).c $(MWRAPFILE3).mw;\
	$(MEX) -v $(GATEWAY3).c ../../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE3) $(MEXLIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY4) -mb $(MWRAPFILE4).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY4) -c $(GATEWAY4).c $(MWRAPFILE4).mw;\
	$(MEX) -v $(GATEWAY4).c ../../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE4) $(MEXLIBS);

clean: objclean
	rm -f examples/ext_dir_solver examples/trans_solver

objclean: 
	rm -f $(OBJS) 
	rm -f test/*.o examples/*.o c/*.o
