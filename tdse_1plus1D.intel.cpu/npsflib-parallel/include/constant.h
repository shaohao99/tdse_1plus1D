//Constants
#ifndef CONSTANT_H
#define CONSTANT_H
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>	
#include <stdlib.h>  // for using exit
#include <complex>
#define complex complex<double>
#define DIM 3           // the number of dimensions
const std::complex I=std::complex(0.,1.); 
const double lightC_au = 137.036;
const double one_by_lightC_au = 1./lightC_au;
const double pi = 3.141592653589793238463383;
const double charge_e1_au = -1.;
const double charge_e2_au = -1.;
const double mass_e1_au = 1.;
const double mass_e2_au = 1.;
const double charge_n1_au=1.;
const double charge_n2_au=1.;
const double mass_n1_au=1836.;
const double mass_n2_au=1836.;
const double eta=1.0472; //pi/2.1;
const int    N_X1_left=20;
const int    N_X1_right=20;
const int    N_X2_left=20;
const int    N_X2_right=20;
const int    N_X3_left=20;
const int    N_X3_right=20;
// shaohao 2008.04
const double tran_au_eV=27.2113845;
const double tran_au_fs=0.024188843265;

const double e_Euler = 2.71828247254;

enum gauge_t { lengthgauge, velocitygauge, othergauge };

// shaohao2009.04
enum orient { rho, z };

enum dump_type { no_dump, exp_dump, poly_dump };
enum ft_type { odd, even, forward, backward };

//void* __dso_handle;

#endif	/* CONSTANT_H */
