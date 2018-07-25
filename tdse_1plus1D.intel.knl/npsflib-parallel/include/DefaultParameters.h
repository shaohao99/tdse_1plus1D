//Parameters
#ifndef DEFAULTPARAMETERS_H
#define DEFAULTPARAMETERS_H
#include "ParameterMap.h"

double soft_param_ee_a = 0.1; //softening parameter in the ee coulomb-potential
double soft_param_en1_b= 0.2; //softening parameter in the en1 coulomb-potential
double soft_param_en2_b= 0.2; //softening parameter in the en2 coulomb-potential

//laser parameters.

double intensity_wcm2 = 5.e+14;
double cep_rad =pi/2.; //Carrier envelope phase
double number_of_cycles = 4.;
double time_after_pulse_in_cycles = 2.;
double omega_au=5.;
double time_step_au = 0.01;
string shape_pulse="gauss";

//Grid parameters

Parameter<vector<int> > n_points( "n_points", vector<int>(3,50), "number of points for (X1,X2,X3)" );

Parameter<vector<double> > spatial_step( "spatial_step", vector<double>(3,0.3), "spatial_steps for (X1,X2,X3)" );

// Gaussian Initial guess

Parameter<vector<double> > sigma("sigma", vector<double>(3,1.), "sigma");
Parameter<vector<double> > gauss_center(
    "gauss_center", vector<double>(3,0.), "Center of Gaussian");
Parameter<vector<double> > initial_momentum(
    "initial_momentum", vector<double>(3,0.), "Initial momentum of wavepacket");

//ECS Complex Scaling parameters.
double angle_eta_X1;
double angle_eta_X2;
double angle_eta_X3;

double fraction_of_grid_X1=0.2;
double fraction_of_grid_left_X2=0.2;
double fraction_of_grid_left_X3=0.2;
double fraction_of_grid_right_X2=0.2;
double fraction_of_grid_right_X3=0.2;

// The Name of the job
Parameter<string> jobname("jobname", "Job", "The Name of the job");
Parameter<string> outputFolder("outputFolder", "out", "The basename of the output directory");  // just a default, is set in Input.cpp from jobname()
Parameter<string> wavefunctionFolder("wavefunctionFolder", "/data/npsf/lib/states/tmp", "The wavefunction archive");  // just a default, is set in Input.cpp from jobname()



// Paramters for the Hamiltonian:
Parameter<double> A_x1("A_x1",-1., "Parameter for the second derivative of X1");
Parameter<double> A_x2("A_x2",-1., "Parameter for the second derivative of X2");
Parameter<double> A_x3("A_x3",-1./4., "Parameter for the second derivative of X3");

Parameter<double> B_x1("B_x1",-1., "Parameter for the first derivative of X1");
Parameter<double> B_x2("B_x2",0., "Parameter for the first derivative of X2");
Parameter<double> B_x3("B_x3",0., "Parameter for the first derivative of X3");




Parameter<string>potential("potential","He","The type of the potential");
Parameter<vector<double> > potentialParam("potentialParam", vector<double>(3,1.), "Parameters");



#endif /* DEFAULTPARAMETERS_H */
