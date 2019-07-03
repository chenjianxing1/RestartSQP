# Copyright (C) 2019
# All Rights Reserved.
#
# Authors: Xinyi Luo and Andreas Waechter      2019-06-17
#

# C++ compiler
CXX = g++

# Ipopt installation directory
#IPOPTDIR = /Users/xinyiluo/Packages/CoinIpopt/build
#IPOPTDIR = /usr/local/Ipopt
#IPOPTDIR = /home/andreasw/COIN-svn/OBJgcc/Ipopt/stable/3.12
# IT IS BETTER TO JUST SET THIS AS ENVIRONMENT VARIABLE!!

IPOPTLIBDIR = $(IPOPTDIR)/lib 

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
IPOPTINCL = `PKG_CONFIG_PATH=$(IPOPTDIR)/lib/pkgconfig: pkg-config --cflags ipopt`

# Linker flags
IPOPTLIBS = `PKG_CONFIG_PATH=$(IPOPTDIR)/lib/pkgconfig: pkg-config --libs ipopt`

# include qpOASES makefile
QPOASESDIR = ../ThirdParty/qpOASES
#include $(QPOASESDIR)/make.mk

# qpOASES variables
QPOASESINCL = -I$(QPOASESDIR)/include
#QPOASESLIBS = $(QPOASES_LINK) $(LINK_LIBRARIES) #-L$(QPOASESDIR)/bin -lqpOASES
#QPOASESLIBS = -L$(QPOASESDIR)/bin -lqpOASES #-lqpOASES_wrapper
QPOASESLIBS = $(QPOASESDIR)/bin/libqpOASES.a #-lqpOASES_wrapper

# C++ compiler flags
#CXXFLAGS = -g -Wall -std=c++11
CXXFLAGS = -g -Wall -std=c++11 -O

INCLUDES = -I. $(IPOPTINCL) $(QPOASESINCL)
LIBS = $(QPOASESLIBS) $(IPOPTLIBS)

# Library utility command
AR = ar rv

# Set sources and objects
SQPLIB_sources = Algorithm.cpp Matrix.cpp MyNLP.cpp Options.cpp \
	QPhandler.cpp QPsolverInterface.cpp SQPTNLP.cpp \
	Utils.cpp Vector.cpp
SQPLIB_objects = $(SQPLIB_sources:.cpp=.o)

# Set library name
SQPLIB = libSQPhotstart.a

# Executable names
SIMPLE_TEST = simple_test
SIMPLE_TEST_sources = simple_test.cpp
SIMPLE_TEST_objects = simple_test.o

EXES_sources = $(SIMPLE_TEST_sources)

all: $(SIMPLE_TEST)

$(SIMPLE_TEST): $(SIMPLE_TEST_objects) $(SQPLIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(SIMPLE_TEST_objects) $(SQPLIB) $(LIBS)

# Create library
$(SQPLIB): $(SQPLIB_objects)
	$(AR) $@ $?

# Rule for .cpp from .o file
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(SQPLIB_objects) $(SIMPLE_TEST_objects) $(dependlists) *.d.*

# Very clean
veryclean: clean
	rm -f $(SIMPLE_TEST) $(LIB)

# Dependencies (handled automatically using tricks from GNU Make)
dependlists = $(SQPLIB_sources:.cpp=.d) $(EXES_sources:.cpp=.d)
-include $(dependlists)
%.d: %.cpp
	@set -e; \
	rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $(INCLUDES) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
