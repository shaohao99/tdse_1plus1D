SHELL=/bin/bash

CXX=icpc
CC=icc

npsflibdir=/project/scv/shaohao/tdse_1plus1D.intel/npsflib-parallel

headerfilenames=constant.h utils.h Cartesian2D.h wavefunction.h DefaultParameters.h ParameterMap.h ConfigFileParser.h Observable.h Laser.h InputStuff.h 
libraryfilenames=Cartesian2D.cpp  Parameter.cpp ConfigFileParser.cpp InputStuff.cpp Laser.cpp 
objects=Cartesian2D.o  Parameter.o ConfigFileParser.o InputStuff.o Laser.o
commonfilenames=$(headerfilenames) $(libraryfilenames)
docfilenames=README INSTALL

#OPTFLAGS=-O2 -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -mfpmath=sse -funroll-loops -fopenmp -nostartfiles
#OPTFLAGS=-fast -openmp # work, but slow
#OPTFLAGS=-O3 -qopenmp -xMIC-AVX512 -qopt-report=2
OPTFLAGS=-O3 -qopenmp -xmic-avx512 -qopt-report=2
#OPTFLAGS=-O3 -xMIC-AVX512 -qopenmp -align -fma -ftz -finline-functions
#OPTFLAGS=-qopenmp-simd -xMIC-AVX512
#OPTFLAGS=-O2 -qopenmp-simd -xcore-avx2
#OPTFLAGS=-qopenmp -xMIC-AVX512
#OPTFLAGS=-openmp -xMIC-AVX512
#OPTFLAGS=-qopenmp -xMIC-AVX512 -qopt-report=2
#OPTFLAGS=-O3 -qopenmp
#OPTFLAGS=-openmp 
#OPTFLAGS=-O3 -fast -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -funroll-loops -openmp
#OPTFLAGS=-O3 -fast -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -funroll-loops -openmp -nostartfiles
#OPTFLAGS=-O3 -fast -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -mfpmath=sse -funroll-loops -openmp -nostartfiles
#OPTFLAGS=-O4 -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -mfpmath=sse -funroll-loops 
#CFLAGS= -fno-common $(OPTFLAGS)
CFLAGS= $(OPTFLAGS)
#CFLAGS= $(OPTFLAGS)

LIBSRCDIRNAME=libsrc
LIBDIRNAME=lib
INCDIRNAME=include

FFTWDIRNAME=/share/pkg/fftw/3.3.4/install/include   # for fftw-3.2
FFTWLIBDIR=/share/pkg/fftw/3.3.4/install/lib

gslInclude=/share/pkg/gsl/1.16/install/include
gsllibdir=/share/pkg/gsl/1.16/install/lib

omp_include=/usr/local/compilers/Intel/cluster_studio_xe_2013.1.046/composer_xe_2013_sp1.2.144/compiler/include

dirCoulombwave=$(npsflibdir)/Coulomb_wave  # for Coulomb wave

#DOXYGEN=/data/npsf/Tools/doxygen-1.4.7/bin/doxygen

define lib_msg
	@echo "--------------------------------------------------------------"
	@echo "    To run $(1), please set (just one time)"
	@echo "    For bash:"
	@echo 
	@echo '  source /data/npsf/share/setenv.sh'
	@echo ""
	@echo "    For csh:"
	@echo 
	@echo '  source /data/npsf/share/setenv.csh'
	@echo
	@echo "(you may also consider putting that line into either .bash_profile or .login)"
endef

define gettar
	@echo "(should be) Creating tarball..."
	@echo "(but it's actually not done at this time.. fix it in common.mk)"
endef

