// $Id:$
#ifndef CARTESIAN2D_H
#define CARTESIAN2D_H
#include "constant.h"
#include "wavefunction.h"
#include "potentials.h"

using namespace std; // to save the std:: in front of vector and cout

namespace Cartesian_2D{

/**
 * Declare everything what is needed for the Crank-Nicholson-Scheme.  
 * independent of geometry namespace
 */

struct temp_data_t {
/**
 * temporary storage for Tridag
 */
  vector<complex> gam;

/**
 * the 3 diagonals
 */
  vector<complex> tridag_upp;
  vector<complex> tridag_mid;
  vector<complex> tridag_low;

/**
 * the potential with respect to the coordinate
 */
  vector<complex> v_1D;
  vector<complex> wf_1D;
  vector<complex> wf_1D_rightside;
  vector<complex> wf_1D_solution;
  
};
  
/**
 *  The namespace Cartesian deals only Hamiltonians H=Dxx + V.
 */

  
/*--------------------------  Hamiltonian class  ----------------------------*/

class Hamiltonian {
  public:

  /**
   * The crucial grid parameters
   * (introduced here to check if wf dimensions match grid dimensions)
   */
  double dx[DIM+1];/**< Spatial Step x[j]( i ) */  
  int n[DIM+1];/**< Number of points in the N[j] dimension */  

  temp_data_t temp_data[DIM+1]; /* provide workspace for DIM dimensions */

  /**
   * Constructor just Initializes the workspace
   */
  Hamiltonian( const wavefunction &wf );
  
  /**
   * Propagation Operators
   */

  void X1( const complex time_step, wavefunction &wf, const ABVparam &p );
  void X2( const complex time_step, wavefunction &wf, const ABVparam &p );
  

  /*
   * Propagation of one time step
   */
  void operator()( const complex time_step, wavefunction &wf, const ABVparam &p)
	  {
		  X1( time_step*.25 , wf , p );
		  X2( time_step*.5 , wf , p );
		  X1( time_step*.25 , wf , p );
	  }

  void X1_Laser( const complex time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X2_Laser( const complex time_step, wavefunction &wf, const ABVparam &p,
		 const double field, const gauge_t gauge );
  void X1_Laser_OMP(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge );
  void X2_Laser_OMP( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge);
  void X1_Laser_dx4_OMP(const complex dt, wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge );
  void X2_Laser_dx4_OMP(const complex dt, wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge );
  

 void X1_ECS_AP_OMP(complex time_step , wavefunction &wf, const ABVparam &p, double laser_vector, double frac_n1, double eta0);
 void X2_ECS_AP_OMP(complex time_step , wavefunction &wf, const ABVparam &p, double laser_vector, double frac_n2, double eta0);
  /*void ECS( const complex time_step, wavefunction &wf, const ABVparam &p );
  void X1_ECS( const complex time_step, wavefunction &wf, const ABVparam &p );  
  void X2_ECS( const complex time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS( const complex time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS_E_R( const complex time_step, wavefunction &wf, const ABVparam &p,
     const double electric_field ); 
  void X3_ECS_A_P( const complex time_step, wavefunction &wf, const ABVparam &p,
      const double electric_field ); */


 /**
   * Adaptative mesh
   */

  //void X1_AM( const complex time_step , wavefunction &wf , const ABVparam &p );
  //void X2_AM(const complex time_step , wavefunction &wf  , const ABVparam &p );
  
  
  //void X2_Laser_AM( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge);
  //void X3_Laser_AM( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge );

}; // end of class Hamiltonian
 
// Shaohao 2012.05
  void Open_File( ofstream &filename, string name );
//  void cft_1D (const vector<complex> &input, const vector<double> &input_x, const int n_inp, const double inpstep, vector<complex> &output, const vector<double> &out_x, const int n_out, const double outstep, const ft_type type, const dump_type dump, const double damp_f);

 
  /**
   * Calculation of Observables:
   */

  double Obs_Norm( wavefunction &wf );
  double Obs_Energy( wavefunction &wf, ABVparam &p );
  double Obs_Energy_X2( wavefunction &wf, ABVparam &p);
  double Obs_Expectation_Value_X1( wavefunction &wf );
  double Obs_Expectation_Value_X2( wavefunction &wf );
  //double Obs_dipole( wavefunction &wf );
  void Obs_dipole( wavefunction wf, complex &dip);
  void Obs_dipole_gex( wavefunction wf, wavefunction small_wf, complex &dip);
 
  double Obs_Expectation_Value_Width_X1( wavefunction &wf );
  double Obs_Expectation_Value_Width_X2( wavefunction &wf );
 
  double Obs_Expectation_Value_Inv_X1( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X2( wavefunction &wf );
 
  double Obs_Expectation_Value_r( wavefunction &wf );
  double Obs_Expectation_Value_Inv_r( wavefunction &wf );
  // vector<double>  H2p_nuclei( Potential_H2_plus_e3D_n2D &pot );

  vector<complex> Obs_Serialize( wavefunction &wf );
  void Mask_Function( wavefunction &wf , double frac_n1_right , double frac_n2_right, double frac_n1_left, double frac_n2_left  , double exponent);
  void Mask_Function_OMP(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n1_left, double frac_n2_left ,double exponent );
  void Gauss_partition(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n1_left, double frac_n2_left, double v0, int n_SI, vector< vector<complex> > &psi_1r, vector< vector<complex> > &psi_1l);
  void Psi_edge( wavefunction wf, double frac_n1_right, double frac_n1_left, double x_SI, double x_edge, vector<complex> &psi_right, vector<complex> &psi_left);
  void TSURFF( wavefunction wf, double frac_n1_right, double frac_n1_left, double x_SI, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left);
   void TSURFF_PRJ( wavefunction wf, wavefunction wf_ion, double frac_n1_right, double frac_n1_left, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left);
   void TSURFF_1D( wavefunction wf, double frac_n1, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left);

  /**
   * Calculation of Derivatives:
   */
  complex First_Derivative_X1( wavefunction &wf , int j , int i );
  complex First_Derivative_X2( wavefunction &wf , int j , int i );
 
  complex Second_Derivative_X1( wavefunction &wf , int j , int i );
  complex Second_Derivative_X2( wavefunction &wf , int j , int i );
 

/*------------------  Stuff needed for FFT and Masks  -----------------------*/


  void Initialize_Momentum( wavefunction &wf, wavefunction &wf_mom);
  void FFT( wavefunction &wf , wavefunction &wf_transformed);
  void Reduce_Wavefunction_DI( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_DI_Mask( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_SI_Mask( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_Bound_Mask( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_DI_both_left( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_DI_both_right( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_DI_opposite( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_SI( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_SI_one_electron( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );
  void Reduce_Wavefunction_Bound( wavefunction &wf , wavefunction &wf_reduced , vector<double> &parameters );

  vector<double> Distribution_X1( wavefunction &wf );
  vector<double> Distribution_X2( wavefunction &wf );
  vector<double> Distribution_X3( wavefunction &wf );
  void Write_2DGraph(vector<double> &x, vector<double> &y, string name );


/*------------------  Wavefunction related functions  -----------------------*/

  // Attention !!!!!
     
  
  /**
   * A 3D initializer (void).  
   * Simply allocate the arrays and instanciate the number of points and the spatial step.   
   * Allocate the wavefunction (wave) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the potential (v) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the x1 (x1) as an array of (n1+2).   
   */

  void Initialize( wavefunction &wf, vector<int> N_points , vector<double> spatial_steps, int symmetry_x1 = 0, int symmetry_x2 = 0 );
  //void Initialize_part( wavefunction &wf, vector<int> N_points , vector<double> spatial_steps, int symmetry_x1 = 0, int symmetry_x2 = 0, int N_left );
  void Initialize_part( wavefunction &wf, vector<int> N_points , vector<double> spatial_steps, int symmetry_x1 = 0, int symmetry_x2 = 0);
  //void Runge_Kutta( wavefunction &wf, Potential_H2_plus_e3D_n2D &pot, double dt, const double Ex, const double Ey );
  //vector<double> H2p_Force(wavefunction &wf , Potential_H2_plus_e3D_n2D &pot, const double Ex, const double Ey);

  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta );
  void Even_Parity_Function( wavefunction &wf , vector<double> sigmas );
  void Interact_Init_Wave( wavefunction &wf, wavefunction wf0);
  void Normalize( wavefunction &wf );

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf ); 
  complex Project( wavefunction &wf , wavefunction &wf2 );
  double Overlap( wavefunction &wf , wavefunction &wf2);
  complex Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );

 /**
   * Adaptative mesh
   */
  void SizeGrid_X1(wavefunction &wf, double criteria, int grow);
  void SizeGrid_X2(wavefunction &wf, double criteria, int grow);
  void SizeGrid_X3(wavefunction &wf, double criteria, int grow);
  void SizeGrid(wavefunction &wf, double criteria, int grow);



/*----------------  Classical Nuclear Motion functions  ---------------------*/

  /* double Calc_Distance( double velocity ); */
/*   double Calc_Velocity( const wavefunction &wf , double distance , Potential_H2_planar_parallel_deriv &pot ); */
/*   void Newton_Equation_Moving_Nuclei( const wavefunction &wf , vector<double> &initial_values , double time_step , Potential_H2_planar_parallel_deriv &pot ); */
  
  void Initial_Eng( wavefunction &wf );


  void X1_AM( const complex time_step, wavefunction &wf, const ABVparam &p );
  void X2_AM( const complex time_step, wavefunction &wf, const ABVparam &p );
  void X3_AM( const complex time_step, wavefunction &wf, const ABVparam &p );

  //void X2_Laser_AM( const complex time_step, wavefunction &wf, const ABVparam &p,
  //    const double field, const gauge_t gauge );
  //void X3_Laser_AM( const complex time_step, wavefunction &wf, const ABVparam &p,
  //    const double field, const gauge_t gauge );





  
 }

#endif	/* CARTESIAN3D_H */
