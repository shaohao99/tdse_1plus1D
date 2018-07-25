# $Id: Makefile 447 2006-08-31 16:38:09Z fhe $
include ../common.mk

.PHONY: all

# make single simple rule for creating each target 
# to have all source files listed in one rule
# to make tar file out of them
#
# here specify the default target:

INCDIR=../$(INCDIRNAME)
FFTWDIR=$(FFTWDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))

all: libnpsf.a

#libnpsf.a: $(libraryfilenames:.cpp=.o)
libnpsf.a: $(objects)
	xiar rs $@ $^
#	ar rs $@ $^

# how to build a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(INCDIR) -I$(FFTWDIR) -I$(dirCoulombwave) -I$(gslInclude) -c $<
#	$(CXX) $(CFLAGS) -I$(INCDIR) -I$(FFTWDIR) -I$(dirCoulombwave) -I$(gslInclude) -I$(omp_include) -c $<
#	$(CXX) $(CFLAGS) -I$(INCDIR) -I$(FFTWDIR) -I$(gslInclude) -I$(omp_include) -c $<

# auto-generate the dependency tree
#include $(subst .cpp,.d,$(libraryfilenames))

#%.d: %.cpp
#	$(CXX) -M -I$(INCDIR) -I$(FFTWDIR) -I$(dirCoulombwave) -I$(gslInclude) -I$(omp_include) $(CPPFLAGS) $< > $@.$$$$; \
#	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
#	rm -f $@.$$$$

# auto-depend done

#doc:
#	$(DOXYGEN)

clean: 
	@rm -f *.o *.a *.d
 
.PHONY: all clean
