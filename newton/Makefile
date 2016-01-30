# Makefile for OS X (10.9) and Linux (Ubuntu 13.10) requires Eigen 3.
#
# To install Eigen 3 on OS X:
#   brew install eigen
#
# On Ubuntu:
#   sudo apt-get install libeigen3-dev
#
# Author: Matt Wytock (mwytock@cs.cmu.edu)

UNAME = $(shell uname)
INCLUDES = -I$(EIGEN)
COMMON_CC  = newton_cd_loop.cc newton_cd_loop_mat.cc
COMMON_H   = newton_cd_loop.h  newton_cd_loop_mat.h
COMMON_OBJ = newton_cd_loop.o  newton_cd_loop_mat.o

# For compiling with mex
MEX_CC = $(COMMON_CC) newton_cd_loop_mex.cc
MEX_FLAGS = $(INCLUDES) -largeArrayDims

# For compiling without mex
CFLAGS = $(INCLUDES) -I$(MATLAB)/extern/include -O2 -march=native -g
LDFLAGS = -L$(MATLAB)/bin/$(ARCH) -lmat -lmx

ifeq ($(UNAME), Darwin)
	MATLAB = /Applications/MATLAB_R2013a.app
	EIGEN = /usr/local/include/eigen3
	ARCH = maci64
	SUFFIX = mexmaci64
	CC = clang++
else
	MATLAB = /usr/local/MATLAB/R2013a
	EIGEN = /usr/include/eigen3
	ARCH = glnxa64
	SUFFIX = mexa64
	CC = g++
	CFLAGS += -fopenmp
	LDFLAGS += -Wl,-rpath-link,$(MATLAB)/bin/$(ARCH) -fopenmp
	MEX_FLAGS += CXXFLAGS='$$CXXFLAGS -fopenmp -march=native'
	MEX_FLAGS += LDFLAGS='$$LDFLAGS -fopenmp'
endif

# Compile the MATLAB interface
newton_cd_loop.$(SUFFIX): $(MEX_CC) $(COMMON_H)
	mex -v -output newton_cd_loop $(MEX_FLAGS) $(MEX_CC)

# Compile the standalone test
newton_cd_loop_test: newton_cd_loop_test.o $(COMMON_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

.cc.o: $(COMMON_H)
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -f newton_cd_loop.$(SUFFIX) newton_cd_loop_test *.o
