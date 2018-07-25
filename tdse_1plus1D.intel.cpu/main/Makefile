# $Id: Makefile $

npsflib_path=/project/scv/shaohao/tdse_1plus1D.intel.cpu/npsflib-parallel

include $(npsflib_path)/common.mk

INCDIR=$(npsflib_path)/$(INCDIRNAME)
LIBDIR=$(npsflib_path)/$(LIBDIRNAME)
LIBSRCDIR=$(npsflib_path)/$(LIBSRCDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))
libraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfilenames))

# tell make where libnpsf.a is lying
libnpsf.a: $(LIBDIR)/libnpsf.a
$(LIBDIR)/libnpsf.a:
	$(MAKE) -C $(LIBDIR)

tdse2d: main_volkov.cpp libnpsf.a
	$(CXX) -I$(INCDIR) -L$(LIBDIR) -I$(FFTWDIRNAME) ${OPTFLAGS} main_volkov.cpp -lnpsf $(FFTWLIBDIR)/libfftw3.a -o $@
#	$(CXX) -I$(INCDIR) -L$(LIBDIR) -I$(FFTWDIRNAME) -O3 -qopenmp -xMIC-AVX512 main_volkov.cpp -lnpsf $(FFTWLIBDIR)/libfftw3.a -o $@
#	$(CXX) $(CFLAGS) -I$(INCDIR) main_volkov.cpp -o $@ -L$(LIBDIR) -lnpsf -I$(FFTWDIRNAME) $(FFTWLIBDIR)/libfftw3.a -I$(omp_include)
	$(call gettar)
	$(call lib_msg, $@)

clean: 	
	rm -f *.o
#	rm -f $(npsflib_path)/lib/*.a
#	rm -f $(npsflib_path)/libsrc/*.a
#	rm -f $(npsflib_path)/libsrc/*.o
#	rm -f $(npsflib_path)/libsrc/*.d


.PHONY: all clean

