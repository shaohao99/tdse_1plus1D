// $Id: Cartesian3D.cpp 234 2006-06-23 11:25:42Z silvio $ 
#include "Cartesian2D.h"
#include "utils.h"
#include "potentials.h"
#include <math.h>
#include <fftw3.h>
#include <complex>
#define complex complex<double>
#include "ParameterMap.h" 
#include "wavefunction.h"
#include "constant.h"
#include <omp.h>

namespace Cartesian_2D{

/*--------------------------  Hamiltonian class  ----------------------------*/

Hamiltonian::Hamiltonian(const wavefunction &wf) {
  /**
   *  Note the grid dimensions that were used to initialize the
   * ABV and temp data vectors
   */
  for (int coord_index=1; coord_index<=2; ++coord_index) {
    n[coord_index]  = wf.n[coord_index];
    dx[coord_index] = wf.x[coord_index].delta;
  }

  /**
   * Allocate the temporary storage for Tridag
   */
  for ( int coord_index=1; coord_index <= 2; ++coord_index ) {
    int num = wf.n[coord_index] + 2;
    temp_data_t &d=temp_data[coord_index];

    d.gam.resize(num, complex(0., 0.) );
    d.tridag_low.resize(num, complex(0., 0.) );
    d.tridag_mid.resize(num, complex(0., 0.) );
    d.tridag_upp.resize(num, complex(0., 0.) );
    d.v_1D.resize(num, complex(0., 0.) );
    d.wf_1D.resize(num, complex(0., 0.) );
    d.wf_1D_rightside.resize(num, complex(0., 0.) );
    d.wf_1D_solution.resize(num, complex(0., 0.) );
  }
};

  /**
   * One time step in Real time (NO LASER), the default value for a in a*(d/dx^2) is a=1/2
   * IMPORTANT, the time step used comes from the Cranck-Nichols equation
   * The auxiliary vector are allocated and initialized inside the function, it uses the tridag (www.nr.com) function implemented in utils.h 
   */  

    /**
     * Propagation of one time step in X1 direction
     * Look at doc/abv.tex Eq. (6) in doc/ABV/abv.tex for details.
     */

  void Hamiltonian::X1( const complex time_step , wavefunction &wf , const ABVparam &p )
  {
    temp_data_t &d=temp_data[1]; 
    if(n[1] != wf.n1 || dx[1] != wf.dx1) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    
    const complex idt=I*time_step;
    const complex arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x1;
    const complex arg_B = complex(0.,0.); //( idt*.25*wf.one_by_dx1   )*p.ABV_B_x1;
    complex arg_V;
    complex tridag_low_Fast =arg_A-arg_B;
    complex tridag_upp_Fast =arg_A+arg_B;

    for( int j = 1 ; j <= wf.n2 ; j++ )  //interation for x2
    {
       for( int i = 0 ; i < wf.n1+2 ; i++ )//Copy the wavefunction
       {
		    /**
		     * Creates a 1D potential and a 1D wavefunction out of the 2D ones
		     */
	    const int index=wf.in2(j,i);
	    d.v_1D[ i ] = p.ABV_V[ index ];
	    d.wf_1D[ i ] = wf.wave[  index ];
       }
	    
       for( int i = 1 ; i <= wf.n1 ; i++ ) //interation for x1
       {
		    /**
		     * Definition of the arguments in Eq. (7-9) in doc/ABV/abv.tex
		     */
	    arg_V = ( idt*.5 )*d.v_1D[ i ];
		    
		    /**
		     * Define the 3 diagonals X of the left side of Eq. (6)  in doc/ABV/abv.tex => look at Eqs. (7-9) in doc/ABV/abv.tex
		     */
	    d.tridag_mid[ i ] =1.-2.*arg_A+arg_V;
		    
	    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ]+( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//general case, symmetry==1
		    
		    /**
		     * for (anti-)symmetric wavefunctions one has to take care of the right boundary conditions (phi(0)=+/-phi(1). That leads to modified values.
		     */
	    if (i==1) {
		    if (wf.symmetry_x1==1)
		    {
			    d.tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
			    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
		    }
			    
		    if (wf.symmetry_x1==-1)
		    {
			    d.tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
			    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
		    }
	    }  // end of the first point
		    
      }  // end of x1
	    
	    /**
	     * Solve the system of linear equation of Eq. (6) in doc/ABV/abv.tex
	     */
      Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
	    
	    /**
	     * Copy back the 1D wavefunction to the 3D one
	     */
      for( int i = 1 ; i <= wf.n1 ; i++ )
      {
	   wf.wave[ wf.in2(j,i) ] = d.wf_1D_solution[ i ];
      }

    } // end of x2
      
    
  } // end of X1_3D
	

	
    /**
     * Propagation of one time step in X2 direction
     * Look at X1_3D for more details
     */
 void Hamiltonian::X2(const complex time_step , wavefunction &wf , const ABVparam &p )
  {
	  temp_data_t &d=temp_data[2]; 
	  if(n[2] != wf.n2 || dx[2] != wf.dx2) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
	  
	  const complex idt=I*time_step;
	  const complex arg_A = idt*.5*wf.one_by_dx2sqr  *p.ABV_A_x2;
	  const complex arg_B = complex(0.,0.); //( idt*.25*wf.one_by_dx2 )*p.ABV_B_x2;
	  complex arg_V;
          complex tridag_low_Fast =arg_A-arg_B;
          complex tridag_upp_Fast =arg_A+arg_B;
	 

	  for( int i = 1 ; i <= wf.n1 ; i++ ) // iteration of x1
	  {
		  
	     for( int j = 0 ; j < wf.n2+2 ; j++ )
	     {
		  const int index=wf.in2(j,i);
		  d.v_1D[ j ]  = p.ABV_V[ index ];
		  d.wf_1D[ j ] = wf.wave[ index ];
	     }
		  
	     for( int j = 1 ; j <= wf.n2 ; j++ )
	     {
			  
		  arg_V = ( idt*.5  )*d.v_1D[ j ];
		  d.tridag_mid[ j ] = 1.-2.*arg_A+arg_V;
			  
		  d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] + ( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ] - tridag_upp_Fast*d.wf_1D[ j+1 ];//sym==1
			  
		  if (j==1) {
			  if (wf.symmetry_x2==1)
			  {
				  d.tridag_mid[ j ] = 1.-arg_A-arg_B+arg_V;
				  d.wf_1D_rightside[ j ]=-tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ j ]-tridag_upp_Fast*d.wf_1D[ j+1 ];
			  }
			  if (wf.symmetry_x2==-1)
			  {
				  d.tridag_mid[ j ] = 1.-3.*arg_A+arg_B+arg_V;
				  d.wf_1D_rightside[ j ]=-tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ j ]-tridag_upp_Fast*d.wf_1D[ j+1 ];
			  }
		  }
			  
	  }  // end of x2
		  
	  Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );

	  for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
		  wf.wave[ wf.in2(j,i) ] = d.wf_1D_solution[ j ];
	  }

     }  // end of x1
  
  } // end of X2_2D
  

  
  void Hamiltonian::X1_Laser(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
	  
	  const complex idt=I*time_step;
	  const complex arg_A = idt*.5*wf.one_by_dx1sqr *p.ABV_A_x1;
	  const complex f_arg_B=time_step*.25*wf.one_by_dx1; 
	  complex arg_B;
	  complex arg_V;
	  complex tridag_upp_Fast;
	  complex tridag_low_Fast;
	  
	  for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
		  for( int i = 0 ; i < wf.n1+2 ; i++ )
		  {
			  const int index=wf.in2( j , i );
			  d.v_1D[ i ]  = p.ABV_V[ index ];
			  d.wf_1D[ i ] = wf.wave[ index ];
		  }
		  
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  
			  if ( gauge == lengthgauge )  // - E dot r
			  {
				  arg_B = 0.; //f_arg_B       *p.ABV_B_x2;  // no first order derivative
				  //arg_V = idt*.5*( p.ABV_B_x1/p.ABV_A_x1*0.5*wf.x1[ i ]*field + d.v_1D[ i ] ); 
				  arg_V = idt*.5*( -wf.x1[ i ]*field + d.v_1D[ i ] ); 
			  }
			  else if ( gauge == velocitygauge )  // - A dot p
			  {
				  //arg_B = f_arg_B*( -p.ABV_B_x1*field*one_by_lightC_au );
				  arg_B = f_arg_B*( -field*one_by_lightC_au );  // using p=-i*d/dx
				  arg_V = idt*.5*d.v_1D[ i ];
			  }
			  
			  tridag_low_Fast    = arg_A-arg_B;
			  d.tridag_mid[ i ] = 1.-2.*arg_A+arg_V;
			  tridag_upp_Fast    = arg_A+arg_B;	     
			  
			  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//sym==1
			  
			  if (i==1) 
			  {
				  if (wf.symmetry_x1==1)
				  {
					  d.tridag_mid[ i ] = 1.-arg_A-arg_B+arg_V;
					  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
				  }
				  
				  if (wf.symmetry_x1==-1)
				  {
					  d.tridag_mid[ i ] = 1.-3.*arg_A+arg_B+arg_V;
					  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
				  }
			  }
			  
		  }
		  
		  Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
		  
		  
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  wf.wave[ wf.in2( j , i ) ] = d.wf_1D_solution[ i ];
		  }
		  
	  }
  } // end of X1_2D_Laser
  


  void Hamiltonian::X2_Laser( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge)
  {
    temp_data_t &d=temp_data[2]; 
    if(n[2] != wf.n2 || dx[2] != wf.dx2) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx2sqr *p.ABV_A_x2;
    const complex f_arg_B=time_step*.25*wf.one_by_dx2; 
    complex arg_B;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;

    for( int i = 1 ; i <= wf.n1 ; i++ )
    {
	    for( int j = 0 ; j < wf.n2+2 ; j++ )
	    {
		    const int index=wf.in2( j , i );
		    d.v_1D[ j ]  = p.ABV_V[ index ];
		    d.wf_1D[ j ] = wf.wave[ index ];
	    }
	    
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	    {
		    
		    if ( gauge == lengthgauge )
		    {
			    arg_B = 0.;  //f_arg_B       *p.ABV_B_x2;
			    //arg_V = idt*.5*( p.ABV_B_x2/p.ABV_A_x2*0.5*wf.x2[ j ]*field + d.v_1D[ j ] ); 
	                    arg_V = idt*.5*( -wf.x2[ j ]*field + d.v_1D[ j ] ); 
		    }
		    else if ( gauge == velocitygauge )
		    {
			    //arg_B = f_arg_B*( -p.ABV_B_x2*field*one_by_lightC_au );
			    arg_B = f_arg_B*( -field*one_by_lightC_au );  // using p=-i*d/dx
			    arg_V = idt*.5*d.v_1D[ j ];
		    }
		    
		    tridag_low_Fast    = arg_A-arg_B;	     
		    d.tridag_mid[ j ] = 1.-2.*arg_A+arg_V;
		    tridag_upp_Fast    = arg_A+arg_B;	     
		    
		    d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];//sym==1
		    
		    if (j==1) {
			    if (wf.symmetry_x2==1)
			    {
				    d.tridag_mid[ j ] = 1.-arg_A-arg_B+arg_V;
				    d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];
			    }
			    
			    if (wf.symmetry_x2==-1)
			    {
				    d.tridag_mid[ j ] = 1.-3.*arg_A+arg_B+arg_V;
				    d.wf_1D_rightside[ j ]= -tridag_low_Fast*d.wf_1D[ j-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ j ] -tridag_upp_Fast*d.wf_1D[ j+1 ];
			    }
		    }
		    
	    }
	    
	    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
	    
	    
	    for( int j = 1 ; j <= wf.n2 ; j++ )
	    {
		    wf.wave[ wf.in2(  j , i ) ] = d.wf_1D_solution[ j ];
	    }
	  
    }
  } // end of X2_3D_Laser

// openmp version
  void Hamiltonian::X1_Laser_OMP(const complex time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge )
  {
	  //temp_data_t &d=temp_data[1]; 
	  const complex idt=I*time_step;
	  const complex arg_A = idt*.5*wf.one_by_dx1sqr *p.ABV_A_x1;
	  const complex f_arg_B=time_step*.25*wf.one_by_dx1; 
	  complex arg_B;
	  complex arg_V;
	  complex tridag_upp_Fast;
	  complex tridag_low_Fast;

          int j=1, i1=0, i2=1, i3=1;
	  
//#pragma omp parallel shared(idt, arg_A, f_arg_B)
#pragma omp parallel
{
    #pragma omp for private(j, i1, i2, i3, arg_V, arg_B, tridag_low_Fast, tridag_upp_Fast)

//#pragma omp simd private(j, i1, i2, i3, arg_V, arg_B, tridag_low_Fast, tridag_upp_Fast)
	  for( j = 1 ; j <= wf.n2 ; j++ )
	  {
                  //int ip=omp_get_thread_num();
                  //cout<<"--X1-- Thread "<<ip<<": j="<<j<<endl;
            
		  vector<complex> v_1D(wf.n1+2), wf_1D(wf.n1+2);    // arrays as private data
		  for( i1 = 0 ; i1 < wf.n1+2 ; i1++ )
		  {
			  const int index=wf.in2( j , i1 );
			  v_1D[ i1 ]  = p.ABV_V[ index ];
			  wf_1D[ i1 ] = wf.wave[ index ];
		  }
		  
                  vector<complex> tridag_mid(wf.n1+2),  wf_1D_rightside(wf.n1+2);    // arrays as private data
		  for( i2 = 1 ; i2 <= wf.n1 ; i2++ )
		  {
			  
			  if ( gauge == lengthgauge )  // - E dot r
			  {
				  arg_B = 0.; //f_arg_B       *p.ABV_B_x2;  // no first order derivative
				  //arg_V = idt*.5*( p.ABV_B_x1/p.ABV_A_x1*0.5*wf.x1[ i ]*field + d.v_1D[ i ] ); 
				  arg_V = idt*.5*( -wf.x1[ i2 ]*field + v_1D[ i2 ] ); 
			  }
			  else if ( gauge == velocitygauge )  // - A dot p
			  {
				  //arg_B = f_arg_B*( -p.ABV_B_x1*field*one_by_lightC_au );
				  arg_B = f_arg_B*( -field*one_by_lightC_au );  // calc -A(t)/c*p using p=-i*d/dx and i^2=-1
				  arg_V = idt*.5*v_1D[ i2 ];
			  }
			  
			  tridag_low_Fast    = arg_A-arg_B;
			  tridag_mid[ i2 ] = 1.-2.*arg_A+arg_V;
			  tridag_upp_Fast    = arg_A+arg_B;	     
			  
			  wf_1D_rightside[ i2 ]= -tridag_low_Fast*wf_1D[ i2-1 ] +( 1.+2.*arg_A-arg_V )*wf_1D[ i2 ] -tridag_upp_Fast*wf_1D[ i2+1 ];//sym==1
			  
			  if (i2==1) 
			  {
				  if (wf.symmetry_x1==1)
				  {
					  tridag_mid[ i2 ] = 1.-arg_A-arg_B+arg_V;
					  wf_1D_rightside[ i2 ]= -tridag_low_Fast*wf_1D[ i2-1 ] +( 1.+arg_A+arg_B-arg_V )*wf_1D[ i2 ] -tridag_upp_Fast*wf_1D[ i2+1 ];
				  }
				  
				  if (wf.symmetry_x1==-1)
				  {
					  tridag_mid[ i2 ] = 1.-3.*arg_A+arg_B+arg_V;
					  wf_1D_rightside[ i2 ]= -tridag_low_Fast*wf_1D[ i2-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*wf_1D[ i2 ] -tridag_upp_Fast*wf_1D[ i2+1 ];
				  }
			  }
			  
		  }
		  
            	  vector<complex> wf_1D_solution(wf.n1+2), gam(wf.n1+2);    // arrays as private data
		  Tridag_Fast( tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_1D_rightside, wf_1D_solution, gam );
		  
		  
		  for( i3 = 1 ; i3 <= wf.n1 ; i3++ )
		  {
			  wf.wave[ wf.in2( j , i3 ) ] = wf_1D_solution[ i3 ];
		  }
		  
	  } // end iteration of j

}  // end omp

  } // end  X1_2D_Laser


  void Hamiltonian::X2_Laser_OMP( const complex time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge)
  {
    //temp_data_t &d=temp_data[2]; 
    const complex idt=I*time_step;
    const complex arg_A = idt*.5*wf.one_by_dx2sqr *p.ABV_A_x2;
    const complex f_arg_B=time_step*.25*wf.one_by_dx2; 
    complex arg_B;
    complex arg_V;
    complex tridag_upp_Fast;
    complex tridag_low_Fast;

    int i=1, j1=0, j2=1, j3=1;

//#pragma omp parallel shared(idt, arg_A, f_arg_B)
#pragma omp parallel 
{
    #pragma omp for private(i, j1, j2, j3, arg_V, arg_B, tridag_low_Fast, tridag_upp_Fast)
//#pragma omp simd private(i, j1, j2, j3, arg_V, arg_B, tridag_low_Fast, tridag_upp_Fast)

    for( i = 1 ; i <= wf.n1 ; i++ )
    {
            //int ip=omp_get_thread_num();
            //cout<<"--X2-- Thread "<<ip<<": i="<<i<<endl;
		  
	    vector<complex> v_1D(wf.n2+2), wf_1D(wf.n2+2);    // arrays as private data
	    for( j1 = 0 ; j1 < wf.n2+2 ; j1++ )
	    {
		    const int index=wf.in2( j1 , i );
		    v_1D[ j1 ]  = p.ABV_V[ index ];
		    wf_1D[ j1 ] = wf.wave[ index ];
	    }
	    
            vector<complex> tridag_mid(wf.n2+2),  wf_1D_rightside(wf.n2+2);    // arrays as private data
	    for( j2 = 1 ; j2 <= wf.n2 ; j2++ )
	    {
		    
		    if ( gauge == lengthgauge )
		    {
			    arg_B = 0.;  //f_arg_B       *p.ABV_B_x2;
			    //arg_V = idt*.5*( p.ABV_B_x2/p.ABV_A_x2*0.5*wf.x2[ j ]*field + d.v_1D[ j ] ); 
	                    arg_V = idt*.5*( -wf.x2[ j2 ]*field + v_1D[ j2 ] ); 
		    }
		    else if ( gauge == velocitygauge )
		    {
			    //arg_B = f_arg_B*( -p.ABV_B_x2*field*one_by_lightC_au );
			    arg_B = f_arg_B*( -field*one_by_lightC_au );  // using p=-i*d/dx
			    arg_V = idt*.5*v_1D[ j2 ];
		    }
		    
		    tridag_low_Fast    = arg_A-arg_B;	     
		    tridag_mid[ j2 ] = 1.-2.*arg_A+arg_V;
		    tridag_upp_Fast    = arg_A+arg_B;	     
		    
		    wf_1D_rightside[ j2 ]= -tridag_low_Fast*wf_1D[ j2-1 ] +( 1.+2.*arg_A-arg_V )*wf_1D[ j2 ] -tridag_upp_Fast*wf_1D[ j2+1 ];//sym==1
		    
		    if (j2==1) {
			    if (wf.symmetry_x2==1)
			    {
				    tridag_mid[ j2 ] = 1.-arg_A-arg_B+arg_V;
				    wf_1D_rightside[ j2 ]= -tridag_low_Fast*wf_1D[ j2-1 ] +( 1.+arg_A+arg_B-arg_V )*wf_1D[ j2 ] -tridag_upp_Fast*wf_1D[ j2+1 ];
			    }
			    
			    if (wf.symmetry_x2==-1)
			    {
				    tridag_mid[ j2 ] = 1.-3.*arg_A+arg_B+arg_V;
				    wf_1D_rightside[ j2 ]= -tridag_low_Fast*wf_1D[ j2-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*wf_1D[ j2 ] -tridag_upp_Fast*wf_1D[ j2+1 ];
			    }
		    }
		    
	    }
	    
            vector<complex> wf_1D_solution(wf.n2+2), gam(wf.n2+2);    // arrays as private data
	    Tridag_Fast(tridag_low_Fast , tridag_mid, tridag_upp_Fast , wf_1D_rightside, wf_1D_solution, gam );
	    
	    
	    for( j3 = 1 ; j3 <= wf.n2 ; j3++ )
	    {
		    wf.wave[ wf.in2(  j3 , i ) ] = wf_1D_solution[ j3 ];
	    }
	  
    }  // end iteration of i
}  // end omp
  } // end of X2_3D_Laser


// propagation with order dx^4, keeping tridiagonal matrice feature
  void Hamiltonian::X1_Laser_dx4_OMP(const complex dt, wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge )
  {
	  //temp_data_t &d=temp_data[1]; 
          double a=-5./3., b=-1./6.; 
          complex idt2=I*dt/2.;
          complex idt2_dxsquare = idt2/(wf.dx1*wf.dx1); 
          complex c1=b+idt2_dxsquare;
          complex c2=b-idt2_dxsquare;
          complex c3=a+2.*idt2_dxsquare;
          complex c4=a-2.*idt2_dxsquare;
          complex d1=a*idt2;
          complex d2=b*idt2;
	  complex Mr_upp, Mr_low, Mr_mid, Ml_upp, Ml_low;
          int j=1;

          vector<double> vt(wf.n1+2);
          int ii;
	  //if(gauge == lengthgauge ){  // - E dot r
//#pragma omp parallel
//{
//    #pragma omp for private(ii)
          for(ii=1; ii<wf.n1+1 ; ii++ ) vt[ii] = -wf.x1[ii]*field;
//}  // end omp
          vt[0]=-(wf.x1[1]-wf.dx1)*field;
          vt[wf.n1+1]=-(wf.x1[wf.n1]+wf.dx1)*field;
          //}	
	  
//#pragma omp parallel shared(c1, c2 ,c3 ,c4, d1, d2, vt)
//{
//    #pragma omp for private(j, Mr_upp, Mr_low, Mr_mid, Ml_upp, Ml_low)

	  for( j = 1 ; j <= wf.n2 ; j++ )
	  {
            
		  vector<complex> v_1D(wf.n1+2), wf_1D(wf.n1+2);    // arrays as private data
		  for( int i = 0 ; i < wf.n1+2 ; i++ )
		  {
			  const int index=wf.in2( j , i );
			  if ( gauge == lengthgauge ) v_1D[ i ]  = vt[i] + p.ABV_V[ index ];
			  wf_1D[ i ] = wf.wave[ index ];
		  }
		  
                  vector<complex> Ml_mid(wf.n1+2),  wf_1D_rightside(wf.n1+2);    // arrays as private data
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  
			  if ( gauge == lengthgauge )  // - E dot r
			  {
			    Ml_upp = c1 + d2*v_1D[i+1];
			    Ml_mid[i] = c4 + d1*v_1D[i];
			    Ml_low = c1 + d2*v_1D[i-1];
			    Mr_upp = c2 - d2*v_1D[i+1];
			    Mr_mid = c3 - d1*v_1D[i];
			    Mr_low = c2 - d2*v_1D[i-1];
			  }
			  else if ( gauge == velocitygauge )  // - A dot p
			  {
			  }
			  
			  wf_1D_rightside[ i ]= Mr_low*wf_1D[ i-1 ] + Mr_mid*wf_1D[ i ] + Mr_upp*wf_1D[ i+1 ];
			  
		  } // end of i
		  
            	  vector<complex> wf_1D_solution(wf.n1+2), gam(wf.n1+2);    // arrays as private data
		  Tridag_Fast( Ml_low, Ml_mid, Ml_upp, wf_1D_rightside, wf_1D_solution, gam );
		  
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  wf.wave[ wf.in2( j , i ) ] = wf_1D_solution[ i ];
		  }
		  
	  } // end iteration of j

//}  // end omp

  } // end  X1_2D_Laser


// propagation with order dx^4, keeping tridiagonal matrice feature
  void Hamiltonian::X2_Laser_dx4_OMP(const complex dt, wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge )
  {
	  //temp_data_t &d=temp_data[1]; 
          double a=-5./3., b=-1./6.; 
          complex idt2=I*dt/2.;
          complex idt2_dxsquare = idt2/(wf.dx2*wf.dx2); 
          complex c1=b+idt2_dxsquare;
          complex c2=b-idt2_dxsquare;
          complex c3=a+2.*idt2_dxsquare;
          complex c4=a-2.*idt2_dxsquare;
          complex d1=a*idt2;
          complex d2=b*idt2;
	  complex Mr_upp, Mr_low, Mr_mid, Ml_upp, Ml_low;
          int i=1;

          vector<double> vt(wf.n2+2);
          int jj;
	  //if(gauge == lengthgauge ){  // - E dot r
//#pragma omp parallel
//{
//    #pragma omp for private(jj)
            for(jj=1; jj<wf.n2+1 ; jj++ ) vt[jj] = -wf.x2[jj]*field;
//} // end omp	
            vt[0]=-(wf.x2[1]-wf.dx2)*field;
            vt[wf.n2+1]=-(wf.x2[wf.n2]+wf.dx2)*field;
//          }	
	  
//#pragma omp parallel shared(c1, c2 ,c3 ,c4, d1, d2, vt)
//{
//    #pragma omp for private(i, Mr_upp, Mr_low, Mr_mid, Ml_upp, Ml_low)

	  for( i = 1 ; i <= wf.n1 ; i++ )
	  {
            
		  vector<complex> v_1D(wf.n2+2), wf_1D(wf.n2+2);    // arrays as private data
		  for( int j = 0 ; j < wf.n2+2 ; j++ )
		  {
			  const int index=wf.in2( j , i );
			  if ( gauge == lengthgauge ) v_1D[ j ]  = vt[j] + p.ABV_V[ index ];
			  wf_1D[ j ] = wf.wave[ index ];
		  }
		  
                  vector<complex> Ml_mid(wf.n2+2),  wf_1D_rightside(wf.n2+2);    // arrays as private data
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  
			  if ( gauge == lengthgauge )  // - E dot r
			  {
			    Ml_upp = c1 + d2*v_1D[j+1];
			    Ml_mid[j] = c4 +d1*v_1D[j];
			    Ml_low = c1 + d2*v_1D[j-1];
			    Mr_upp = c2 - d2*v_1D[j+1];
			    Mr_mid = c3 - d1*v_1D[j];
			    Mr_low = c2 - d2*v_1D[j-1];
			  }
			  else if ( gauge == velocitygauge )  // - A dot p
			  {
			  }
			  
			  wf_1D_rightside[ j ]= Mr_low*wf_1D[ j-1 ] + Mr_mid*wf_1D[ j ] + Mr_upp*wf_1D[ j+1 ];
			  
		  } // end of j
		  
            	  vector<complex> wf_1D_solution(wf.n2+2), gam(wf.n2+2);    // arrays as private data
		  Tridag_Fast( Ml_low, Ml_mid, Ml_upp, wf_1D_rightside, wf_1D_solution, gam );
		  
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  wf.wave[ wf.in2( j , i ) ] = wf_1D_solution[ j ];
		  }
		  
	  } // end iteration of i

//}  // end omp

  } // end  X2_2D_Laser



  /********************************************************************************/
  /* ECS */
  /********************************************************************************/

 void Hamiltonian::X1_ECS_AP_OMP(complex time_step , wavefunction &wf, const ABVparam &p, double laser_vector, double frac_n1, double eta0)
  {

    int n1_right=int(wf.n1*(1.-frac_n1));
    int n1_left=int(wf.n1*frac_n1);
    complex arg_A = I*time_step/( 2.*wf.dx1*wf.dx1 )*p.ABV_A_x1;
    complex arg_B = I*time_step/( 4.*wf.dx1 )       *(-I*laser_vector/lightC_au );
    //complex tridag_low, tridag_up;
    int j;

//#pragma omp parallel shared(n1_left,n1_right,arg_A,arg_B,eta0,time_step,I)
//#pragma omp parallel shared(n1_left,n1_right,arg_A,arg_B,eta0,time_step)
//{
//    #pragma omp for private(j)
    
    for( j = 1 ; j <= wf.n2 ; j++ )
      {	

	vector<complex> v_1D(wf.n1+2), wf_1D(wf.n1+2);    // arrays as private data
	for( int i = 0 ; i < wf.n1+2 ; i++ )
	  {
	    v_1D[ i ]  = p.ABV_V[ wf.in2( j , i ) ];
	    wf_1D[ i ] = wf.wave[  wf.in2(  j , i ) ];
	  }
	
        vector<complex> tridag_mid(wf.n1+2), tridag_low(wf.n1+2), tridag_up(wf.n1+2), wf_1D_rightside(wf.n1+2);
	for( int i = 1 ; i <= wf.n1 ; i++ )
	  {
	    complex arg_V = I*time_step/2. *v_1D[ i ];
	    
	    if (i>n1_left && i<n1_right) // center of X1
	      {
		tridag_low[ i ]     = arg_A-arg_B;
		tridag_mid[ i ]     = 1.-2.*arg_A+arg_V;
		tridag_up[ i ]     = arg_A+arg_B;
		wf_1D_rightside[ i ]=
		  -tridag_low[ i ]*wf_1D[ i-1 ]
		  +( 2.- tridag_mid[ i ] )*wf_1D[ i ]
		  -tridag_up[ i ]*wf_1D[ i+1 ];
	      }
	    
	    if (i<n1_left | i>n1_right )  // left of x1
	      {
		tridag_low[ i ]     = arg_A*exp(-2.*I*eta0)-arg_B;
		tridag_mid[ i ]  = 1.-2.*arg_A*exp(-2.*I*eta0)+arg_V;
		tridag_up[ i ]     = arg_A*exp(-2.*I*eta0)+arg_B;
		
		wf_1D_rightside[ i ]=
		  -tridag_low[ i ]*wf_1D[ i-1 ]
		  +( 2.-tridag_mid[ i ] )*wf_1D[ i ]
		  -tridag_up[ i ]*wf_1D[ i+1 ];
	      }
	    
	    if (i==n1_left) // left discontinuity point
	      {     
		tridag_low[ i ] = arg_A*2.*exp(-I*eta0)/(1.+exp(I*eta0))+2.*arg_B*(-exp(-I*eta0)/(1.+exp(I*eta0)));
		tridag_mid[ i ] = 1.+arg_A*(-2.*exp(-I*eta0))+arg_V-2.*arg_B*(exp(I*eta0)-1.)/exp(I*eta0);
		tridag_up[ i ] = arg_A*2./(1.+exp(I*eta0))+2.*arg_B*(exp(I*eta0)/(1.+exp(I*eta0)));
		
		wf_1D_rightside[ i ]=
		  -tridag_low[ i ]*wf_1D[ i-1 ]
		  +( 2. - tridag_mid[ i ] )*wf_1D[ i ]
		  -tridag_up[ i ]*wf_1D[ i+1 ];	
		
	      }
	    
	    if (i==n1_right) // right discontinuity point
	      {
		tridag_low[ i ]  = arg_A*2./(1.+exp(I*eta0))+2.*arg_B*(-exp(I*eta0)/(1.+exp(I*eta0)));
		tridag_mid[ i ]  = 1.+arg_A*(-2.*exp(-I*eta0))+arg_V+2.*arg_B*(exp(I*eta0)-1.)/exp(I*eta0);
		tridag_up[ i ]  = arg_A*(2.*exp(-I*eta0))/(1.+exp(I*eta0))+2.*arg_B*(exp(-I*eta0)/(1.+exp(I*eta0)));
		
		wf_1D_rightside[ i ]=
		  -tridag_low[ i ]*wf_1D[ i-1 ]
		  +( 2.-tridag_mid[ i ] )*wf_1D[ i ]
		  -tridag_up[ i ]*wf_1D[ i+1 ];
	      }
	  }

        vector<complex> wf_1D_solution(wf.n1+2), gam(wf.n1+2);    // arrays as private data
	//Tridag_Fast( tridag_low , tridag_mid , tridag_up , wf_1D_rightside , wf_1D_solution , gam );
	Tridag( tridag_low , tridag_mid , tridag_up , wf_1D_rightside , wf_1D_solution , gam );

	for( int i = 1 ; i <= wf.n1 ; i++ )  wf.wave[ wf.in2( j , i ) ] = wf_1D_solution[ i ];
	
      } // end i
//}  // end omp

  } // end of Hamiltonian_X2_ECS


 void Hamiltonian::X2_ECS_AP_OMP(complex time_step , wavefunction &wf, const ABVparam &p, double laser_vector, double frac_n2, double eta0)
  {

    int n2_right=int(wf.n2*(1.-frac_n2));
    int n2_left=int(wf.n2*frac_n2);
    complex arg_A = I*time_step/( 2.*wf.dx2*wf.dx2 )*p.ABV_A_x2;
    complex arg_B = I*time_step/( 4.*wf.dx2 )       *(-I*laser_vector/lightC_au );
    //complex tridag_low, tridag_up;
    int i;

//#pragma omp parallel shared(n2_left,n2_right,arg_A,arg_B,eta0,time_step,I)
//#pragma omp parallel shared(n2_left,n2_right,arg_A,arg_B,eta0,time_step)
//{
//    #pragma omp for private(i)
    
    for( i = 1 ; i <= wf.n1 ; i++ )
      {	

	vector<complex> v_1D(wf.n2+2), wf_1D(wf.n2+2);    // arrays as private data
	for( int j = 0 ; j < wf.n2+2 ; j++ )
	  {
	    v_1D[ j ]  = p.ABV_V[ wf.in2( j , i ) ];
	    wf_1D[ j ] = wf.wave[  wf.in2(  j , i ) ];
	  }
	
        vector<complex> tridag_mid(wf.n2+2), tridag_low(wf.n2+2), tridag_up(wf.n2+2), wf_1D_rightside(wf.n2+2);
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    complex arg_V = I*time_step/2. *v_1D[ j ];
	    
	    if (j>n2_left && j<n2_right) // center of X2
	      {
		tridag_low[ j ]     = arg_A-arg_B;
		tridag_mid[ j ]     = 1.-2.*arg_A+arg_V;
		tridag_up [ j ]    = arg_A+arg_B;
		wf_1D_rightside[ j ]=
		  -tridag_low[ j ]*wf_1D[ j-1 ]
		  +( 2.- tridag_mid[ j ] )*wf_1D[ j ]
		  -tridag_up[ j ]*wf_1D[ j+1 ];
	      }
	    
	    if (j<n2_left | j>n2_right )  // left of x2
	      {
		tridag_low[ j ]     = arg_A*exp(-2.*I*eta0)-arg_B;
		//d.tridag_low[ j ]     = arg_A*exp(-2.*I*eta)-arg_B*exp(-I*eta);
		tridag_mid[ j ]  = 1.-2.*arg_A*exp(-2.*I*eta0)+arg_V;
		tridag_up[ j ]     = arg_A*exp(-2.*I*eta0)+arg_B;
		//d.tridag_up[ j ]     = arg_A*exp(-2.*I*eta)+arg_B*exp(-I*eta);
		
		wf_1D_rightside[ j ]=
		  -tridag_low[ j ]*wf_1D[ j-1 ]
		  +( 2.-tridag_mid[ j ] )*wf_1D[ j ]
		  -tridag_up[ j ]*wf_1D[ j+1 ];
	      }
	    
	    if (j==n2_left) // left discontinuity point
	      {     
		tridag_low[ j ] = arg_A*2.*exp(-I*eta0)/(1.+exp(I*eta0))+2.*arg_B*(-exp(-I*eta0)/(1.+exp(I*eta0)));
		tridag_mid[ j ] = 1.+arg_A*(-2.*exp(-I*eta0))+arg_V-2.*arg_B*(exp(I*eta0)-1.)/exp(I*eta0);
		tridag_up [ j ]= arg_A*2./(1.+exp(I*eta0))+2.*arg_B*(exp(I*eta0)/(1.+exp(I*eta0)));
		
		wf_1D_rightside[ j ]=
		  -tridag_low[ j ]*wf_1D[ j-1 ]
		  +( 2. - tridag_mid[ j ] )*wf_1D[ j ]
		  -tridag_up[ j ]*wf_1D[ j+1 ];	
		
	      }
	    
	    if (j==n2_right) // right discontinuity point
	      {
		tridag_low[ j ]  = arg_A*2./(1.+exp(I*eta0))+2.*arg_B*(-exp(I*eta0)/(1.+exp(I*eta0)));
		tridag_mid[ j ]  = 1.+arg_A*(-2.*exp(-I*eta0))+arg_V+2.*arg_B*(exp(I*eta0)-1.)/exp(I*eta0);
		tridag_up[ j ]  = arg_A*(2.*exp(-I*eta0))/(1.+exp(I*eta0))+2.*arg_B*(exp(-I*eta0)/(1.+exp(I*eta0)));
		
		wf_1D_rightside[ j ]=
		  -tridag_low[ j ]*wf_1D[ j-1 ]
		  +( 2.-tridag_mid[ j ] )*wf_1D[ j ]
		  -tridag_up[ j ]*wf_1D[ j+1 ];
	      }
	  }

        vector<complex> wf_1D_solution(wf.n2+2), gam(wf.n2+2);    // arrays as private data
	//Tridag_Fast( tridag_low , tridag_mid , tridag_up , wf_1D_rightside , wf_1D_solution , gam );
	Tridag( tridag_low , tridag_mid , tridag_up , wf_1D_rightside , wf_1D_solution , gam );

	for( int j = 1 ; j <= wf.n2 ; j++ )  wf.wave[ wf.in2( j , i ) ] = wf_1D_solution[ j ];
	
      } // end i
//}  // end omp

  } // end of Hamiltonian_X2_ECS




/*------------------  Wavefunction related functions  -----------------------*/

  /**
   * Initialize wavefunction: partial grids
   */
	void Initialize_part( wavefunction &wf , vector<int> n_points , vector<double> spatial_steps , int symmetry_x1, int symmetry_x2)
//	void Initialize_part( wavefunction &wf , vector<int> n_points , vector<double> spatial_steps , int symmetry_x1, int symmetry_x2, int N_left)
	{
		wf.verbose=false;
		/**
		 * Initialize symmetries of the wavefunction.  
		 */
		wf.symmetry_x1=symmetry_x1;
		wf.symmetry_x2=symmetry_x2;
	
                wf.n1 = n_points[ 0 ];
                wf.n2 = n_points[ 1 ];
		
		/**
		 * Allocate the wavefunction and the potential
		 */
		wf.wave.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		
		
		/**
		 * Allocate the grid
		 */
                wf.x1.type=PartAxis;
                wf.x2.type=FullAxis;
		//wf.x1.Init_partial("x1", n_points[0], spatial_steps[0], wf.x1.type, N_left);
		wf.x1.Init("x1", n_points[0], spatial_steps[0], wf.x1.type);
		wf.x2.Init("x2", n_points[1], spatial_steps[1], wf.x2.type);
		wf.dp1 = 2*pi/(wf.dx1*wf.n1);
		wf.dp2 = 2*pi/(wf.dx2*wf.n2);
		wf.p1.Init("p_x1", n_points[0], wf.dp1, FFTWAxis);
		wf.p2.Init("p_y2", n_points[1], wf.dp2, FFTWAxis);
    
		/**
		 * Initialize the spatial grid helper constants.
		 */
		wf.one_by_dx1=1. / wf.dx1;     // wf.dx1=x[1].delta (which itself is set in Axis.Init)
		wf.one_by_dx2=1. / wf.dx2;
		wf.one_by_dx1sqr=1. / ( wf.dx1*wf.dx1 );
		wf.one_by_dx2sqr=1. / ( wf.dx2*wf.dx2 );

        } //end


	void Initialize( wavefunction &wf , vector<int> n_points , vector<double> spatial_steps , int symmetry_x1, int symmetry_x2)
	{
		wf.verbose=false;
		/**
		 * Initialize symmetries of the wavefunction.  
		 */
		wf.symmetry_x1=symmetry_x1;
		wf.symmetry_x2=symmetry_x2;
	
                wf.n1 = n_points[ 0 ];
                wf.n2 = n_points[ 1 ];
		
		
		/**
		 * Allocate the wavefunction and the potential
		 */
		wf.wave.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		
		
		/**
		 * Allocate the grid
		 */
		if (abs(symmetry_x1) == 1) wf.x1.type=HalfAxis;
		else wf.x1.type=FullAxis;
		if (abs(symmetry_x2) == 1) wf.x2.type=HalfAxis;
		else wf.x2.type=FullAxis;
		wf.x1.Init("x1", n_points[0], spatial_steps[0], wf.x1.type);
		wf.x2.Init("x2", n_points[1], spatial_steps[1], wf.x2.type);
		wf.dp1 = 2*pi/(wf.dx1*wf.n1);
		wf.dp2 = 2*pi/(wf.dx2*wf.n2);
		wf.p1.Init("p_x1", n_points[0], wf.dp1, FFTWAxis);
		wf.p2.Init("p_y2", n_points[1], wf.dp2, FFTWAxis);
		
    
		/**
		 * Initialize the spatial grid helper constants.
		 */
		wf.one_by_dx1=1. / wf.dx1;     // wf.dx1=x[1].delta (which itself is set in Axis.Init)
		wf.one_by_dx2=1. / wf.dx2;
		wf.one_by_dx1sqr=1. / ( wf.dx1*wf.dx1 );
		wf.one_by_dx2sqr=1. / ( wf.dx2*wf.dx2 );
 
	
	} // Initialize_wavefunction
	
	void Initialize_Momentum( wavefunction &wf , wavefunction &wf_mom)
	{
		
		wf_mom.verbose=false;
		
		/**
		 * Initialize symmetries of the wavefunction.  
		 */
		wf_mom.symmetry_x1=wf.symmetry_x1;
		wf_mom.symmetry_x2=wf.symmetry_x2;
	
		//symmetry=0: X starts from -X_max, no symmetry or antisymmetry for wavefunction necessary
		//symmetry=1: X starts from 0, wavefunction is symmetric
		//symmetry=-1: X starts from 0, wavefunction is antisymmetric
		
		/**
		 * Allocate the grid
		 */
		// !!!!!!!!!!! only symmetry =0 !!!!!!!!!!!!
		wf_mom.x1.Init("p_x", wf_mom.n1, wf_mom.dx1, FullAxis);
		wf_mom.x2.Init("p_y", wf_mom.n2, wf_mom.dx2, FullAxis);
		wf_mom.x3.Init("p_z", wf_mom.n3, wf_mom.dx3, FullAxis);


		/**
		 * Allocate the wavefunction and the potential
		 */
		wf_mom.wave.resize(( wf_mom.n1+2 )*( wf_mom.n2+2 ), 0.);
		
		
		/**
		 * Allocate the grid
		 */
		wf_mom.x1.resize(wf_mom.n1+2, 0.);
		wf_mom.x2.resize(wf_mom.n2+2, 0.);
		wf_mom.x3.resize(wf_mom.n3+2, 0.);
    
		/**
		 * Define the real grid coordinates
		 */
		
		//    if(wf.symmetry_x1==1 || wf.symmetry_x1==-1)
		/**
		 * symmetric(=1) or antisymmetric(=-1) wavefunction; only half the grid is needed
		 */
		



    // !!!!!!!!!!! only symmetry =0 !!!!!!!!!!!!

		for( int i = 1 ; i <= int(wf_mom.n1/2) +1 ; i++ )
		{	
			wf_mom.x1[ i ] =   (i-1)*wf_mom.dx1;
		}  
		for( int i = int(wf_mom.n1/2)+2 ; i <= wf_mom.n1 ; i++ )
		{	
			wf_mom.x1[ i ] =   (i-1-wf_mom.n1)*wf_mom.dx1;
		}  
		for( int j = 1 ; j <= int(wf_mom.n2/2) +1 ; j++ )
		{	
			wf_mom.x2[ j ] =   (j-1)*wf_mom.dx2;
		}  
		for( int j = int(wf_mom.n2/2)+2 ; j <= wf_mom.n2 ; j++ )
		{	
			wf_mom.x2[ j ] =   (j-1-wf_mom.n2)*wf_mom.dx2;
		}  
		
		
	} // Initialize_momentum
	
	
	/**
	 * Initialize a Gaussian function as a initial wavepacket for imaginary time propagation
	 */
   void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta ){
		
        cout<<"---dx1,dx2"<<'\t'<<wf.dx1<<'\t'<<wf.dx2<<endl;
        cout<<"---n1,n2"<<'\t'<<wf.n1<<'\t'<<wf.n2<<endl;
      
	for( int j = 1 ; j <= wf.n2 ; j++ )
	{
		for( int i = 1 ; i <= wf.n1 ; i++ )
		{
			double arg_x1 = ( wf.x1[ i ]-gauss_centers[ 0 ] )/( 2.*sigmas[ 0 ] );
			double arg_x2 = ( wf.x2[ j ]-gauss_centers[ 1 ] )/( 2.*sigmas[ 1 ] );
				
			complex arg_p1 = I*initial_momenta[ 0 ]*wf.x1[ i ];
			complex arg_p2 = I*initial_momenta[ 1 ]*wf.x2[ j ];
				
			wf.wave[ wf.in2(j,i) ] = exp( -arg_x1*arg_x1-arg_x2*arg_x2 )*exp( arg_p1+arg_p2 );

		}
	}
		
		
   } // end of Guess_Function_Gauss
  
  
   void Even_Parity_Function( wavefunction &wf , vector<double> sigmas ){
      
     for( int j = 1 ; j <= wf.n2 ; j++ )
     { 
	for( int i = 1 ; i <= wf.n1 ; i++ )
	{
	   wf.wave[wf.in2(j,i)] = sin(wf.x1[i])*sin(wf.x2[j])*exp( -wf.x1[i]*wf.x1[i]/sigmas[0] -wf.x2[j]*wf.x2[j]/sigmas[1] );
	}
     }
		
   } // end of Even_Parity_Function

   
   void Interact_Init_Wave( wavefunction &wf, wavefunction wf0  ){
     
     for( int j = 0 ; j <= wf.n2+1 ; j++ )  // set 0 on large grids
     { 
	for( int i = 0 ; i <= wf.n1+1 ; i++ )
	{
	   wf.wave[wf.in2(j,i)] = 0.; 
	}
     }

	int shifter_n1=0;
	int shifter_n2=0;
		
	if (wf.symmetry_x1==0)
	{
		shifter_n1=(wf.n1-wf0.n1)/2;
	}    
	if (wf.symmetry_x2==0)
	{
		shifter_n2=(wf.n2-wf0.n2)/2;
	}

        int N_left=240; //280;
        if (wf.symmetry_x1==2)  //partial grids
	{
		shifter_n1=int(N_left-wf0.n1/2);  //even number
	}    


     for( int j = 1 ; j <= wf0.n2 ; j++ )  // set values on small grids
     { 
	for( int i = 1 ; i <= wf0.n1 ; i++ )  // psi=-(x1+x2)*psi0
	{
	   wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)] = -(wf0.x1[i]+wf0.x2[j])*wf0.wave[wf0.in2(j,i)];
	   //wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)] = -wf0.x1[i]*wf0.wave[wf.in2(j,i)];
	}
     }

   }  //end of Interact_Init_Wave


    /**
     * Normalize of the wavefunction phi' = phi/sqrt( <phi|phi> )
     */
	void Normalize( wavefunction &wf )
	{
		double norm = Obs_Norm( wf );
		//	if(wf.verbose==true)
		//	cout << "Normalizing: Input Norm = " << norm << "\n";
		
		for( int index = 1 ; index < ( wf.n1+2 )*( wf.n2+2 ) ; index++ ) 
		{
			wf.wave[ index ] /= norm; 
		}
		
		norm = Obs_Norm( wf );
		
		//	if(wf.verbose==true)
		//	cout << "Normalizing: Output Norm = " << norm << "\n";
		
	} // end of Normalize3D
	
	
    /**
     * Calculate the norm ( = sqrt( <phi|phi> ) = sqrt( int( x1*dx1*dx2* ( conj( phi )*phi ) ) ) )
     */
	double Obs_Norm( wavefunction &wf )
	{
		double obs=0.;
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += norm( wf.wave[ wf.in2(j,i) ] );
			}
		}
		
		obs *= wf.dx1*wf.dx2;
		obs= sqrt( obs );
		return obs;
		
	} // end of Obs_Norm
	

	
	/**
	 * Calculate the energy ( = <phi|H|phi> with the Hamiltonian in Eq. (4) in doc/ABV/abv.tex
	 */
	double Obs_Energy( wavefunction &wf, ABVparam &p)
	{
		double obs=0.;
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += real( conj( wf.wave[ wf.in2(j,i) ] )*(
						     wf.wave[ wf.in2(j,i) ]*( p.ABV_V[ wf.in2(j,i) ] * 2. )
						     +Second_Derivative_X1( wf , j , i )*p.ABV_A_x1
						     +Second_Derivative_X2( wf , j , i )*p.ABV_A_x2
						     
						     ) );
			
			}
		}
		obs *= wf.dx1*wf.dx2;
		return obs;
		
	} // end of Obs_Energy
	

	 // Calculate the energy ( = <phi|H|phi> with the Hamiltonian in Eq. (4) in doc/ABV/abv.tex
	double Obs_Energy_X2( wavefunction &wf, ABVparam &p)
	{
		double obs=0.;
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += real( conj( wf.wave[ wf.in2(j,i) ] )*(
						   wf.wave[ wf.in2(j,i) ] * p.ABV_V[ wf.in2(j,i) ] 
						   + Second_Derivative_X2( wf , j , i )*p.ABV_A_x2
						  ) );
			
			}
		}
		obs *= wf.dx1*wf.dx2;
		return obs;
		
	} // end of Obs_Energy


	/**
	 * Calculate the expectation value <x1> = <phi|x1|phi>
	 */
	double Obs_Expectation_Value_X1( wavefunction &wf )
	{
		double obs = 0.;
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += wf.x1[ i ]
					*norm( wf.wave[ wf.in2(j,i) ] );
			}
		}
	
		obs *= wf.dx1*wf.dx2;
		return obs;
		
		
	} // end of Obs_Expectation_Value_X1



	/**
	 * Calculate the expectation value <x2> = <phi|x2|phi>
	 */
	double Obs_Expectation_Value_X2( wavefunction &wf )
	{
		double obs = 0.;
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += wf.x2[ j ]
					*norm( wf.wave[ wf.in2(j,i) ] );
			}
		}
		
		obs *= wf.dx1*wf.dx2;
		return obs;
		
		
	} // end of Obs_Expectation_Value_X2

	


    /**
     * Calculate dx1 = <phi|sqrt( <x1^2> -<x1>^2 )|phi>
     */
  double Obs_Expectation_Value_Width_X1( wavefunction &wf )
  {
    double obs;
    double exp_x1 = Obs_Expectation_Value_X1( wf );
    double exp_x1_2 = 0.;
    
 
	for( int j = 1 ; j <= wf.n2 ; j++)
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		exp_x1_2 += wf.x1[ i ]*wf.x1[ i ]
				  *norm( wf.wave[ wf.in2(j,i) ] );
	      }
	  }
      
    exp_x1_2 *=wf.dx1*wf.dx2;
    obs = sqrt( exp_x1_2-exp_x1*exp_x1 );
    return obs;
    
  } // end of Obs_Expectation_Value_Width_X1

  
  
    /**
     * Calculate dx2 = <phi|sqrt( <x2^2> -<x2>^2 )|phi>
     */
  double Obs_Expectation_Value_Width_X2( wavefunction &wf )
  {
    double obs;
    double exp_x2 = Obs_Expectation_Value_X2( wf );
    double exp_x2_2 = 0.;
    
  
	for( int j = 1 ; j <= wf.n2 ; j++)
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		exp_x2_2 += wf.x2[ j ]*wf.x2[ j ]
				  *norm( wf.wave[ wf.in2(j,i) ] );
	      }
	  }
    
    exp_x2_2 *=wf.dx1*wf.dx2;
    obs = sqrt( exp_x2_2-exp_x2*exp_x2 );
    return obs;
    
  } // end of Obs_Expectation_Value_Width_X2

  



    /**
     * Dump the wavefunction into a serial vector of complex values
     */
  vector<complex> Obs_Serialize( wavefunction &wf )
  {
    vector<complex> result;
    result.reserve((wf.n2+2)*(wf.n1+2));  // pre-allocate memory for all the elements
    
   
	for( int j = 0 ; j <= wf.n2+1 ; j++)
	  {
	    for( int i = 0 ; i <= wf.n1+1 ; i++ )
	      {
		result.push_back(wf.wave[ wf.in2(j,i) ]);
	      }
	  }
      
    return result;
  } // end of Obs_Serialize



  /**
     * Project over the Coordinate X2
     * No volume element was used in this case. The size of the oputcoming vector is (n3+2)*(n1+2)
     */
  vector<double> Obs_Projection_X2 (wavefunction &wf)
  {
    vector<double> result;
    result.reserve((wf.n1+2));  // pre-allocate memory for all the elements

   
    for( int i = 0 ; i <= wf.n1+1 ; i++)
    {
	    double sum = 0.;
	    for( int j = 0 ; j <= wf.n2+1 ; j++ )
		    sum += abs( wf.wave[ wf.in2(j,i) ] )*wf.dx2;
	    result[i]=sum;
    }
    
    return result;
    
  } // end of Obs_Projection_X2





	
    /**
     * Calculate first derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
	complex First_Derivative_X1( wavefunction &wf , int j , int i )
	{
		complex deriv;
		
		if (wf.symmetry_x1==1 && i ==1 )//symmetry
		{
			
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx1 );
			
		}
		if (wf.symmetry_x1==-1 && i ==1 )//antisymmetry
		{
			
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]+wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx1 );
			
		}
		else
		{
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(  j , i-1 ) ] ) / ( 2.0*wf.dx1 );
		}//sym==1
		return deriv;
		
	} // end of First_Derivative_X1
	
	
	
	/**
	 * Calculate first derivative in x2-direction - look at Eq. (5) in doc/ABV/abv.tex
	 */
  complex First_Derivative_X2( wavefunction &wf  , int j , int i )
  {
    complex deriv;

    if (wf.symmetry_x2==1 && j==1)//symmetry
      {
	deriv = ( wf.wave[ wf.in2(  j+1 , i ) ]-wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx2 );
	
      }
    if (wf.symmetry_x2==-1 && j==1)//antisymmetry
      {
	deriv = ( wf.wave[ wf.in2(  j+1 , i ) ]+wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx2 );
	
      }
    else
      {
	deriv = ( wf.wave[ wf.in2(  j+1 , i ) ]-wf.wave[ wf.in2(  j-1 , i ) ] ) / ( 2.0*wf.dx2 );
      }

    return deriv;
  } // end of First_Derivative_X2
  


  


    /**
     * Calculate second derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
  complex Second_Derivative_X1( wavefunction &wf , int j , int i )
  {

    if (i==1) {
    if (wf.symmetry_x1==1)//symmetry
      {
	return ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx1sqr;
      }
    if (wf.symmetry_x1==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in2( j , i+1 ) ]-3.*wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx1sqr;
      }
    //if (wf.symmetry_x1==1)
    } else {

      return ( wf.wave[ wf.in2(  j , i+1 ) ]-2.*wf.wave[ wf.in2(j,i) ] + wf.wave[ wf.in2(  j , i-1 ) ] ) * wf.one_by_dx1sqr;//sym==1

    }
      
  } // end of Second_Derivative_X1
  


    /**
     * Calculate second derivative in x2-direction - look at Eq. (5) in doc/ABV/abv.tex
     */
  complex Second_Derivative_X2( wavefunction &wf , int j , int i )
  {

    if (j==1) {
    if (wf.symmetry_x2==1)//symmetry
      {
	return ( wf.wave[ wf.in2(  j+1 , i ) ]-wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx2sqr;
	
      }
    if (wf.symmetry_x2==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in2(  j+1 , i ) ]-3.*wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx2sqr;
	
      }
    // if (wf.symmetry_x2==1)
    } else {

      return ( wf.wave[ wf.in2(  j+1 , i ) ]-2.*wf.wave[ wf.in2(j,i) ] + wf.wave[ wf.in2(  j-1 , i ) ] ) * wf.one_by_dx2sqr;  //sym==1

    }

  } // end of Second_Derivative_X2
  

// gaussian absorber
	void Gauss_partition(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n1_left, double frac_n2_left, double v0, int n_SI, vector< vector<complex> > &psi_1r, vector< vector<complex> > &psi_1l)
{

     int start_n1_right=int((wf.n1)*(1.-frac_n1_right));
     int start_n2_right=int((wf.n2)*(1.-frac_n2_right));
     int start_n1_left=int((wf.n1)*frac_n1_left);
     int start_n2_left=int((wf.n2)*frac_n2_left);
     double start_x1_right=wf.x1[start_n1_right];
     double start_x2_right=wf.x2[start_n2_right];
     double start_x1_left=wf.x1[start_n1_left+1];  // i=0 is not counted for x1
     double start_x2_left=wf.x2[start_n2_left+1];

     /*int n1_right=int((wf.n1)*frac_n1_right)+2;
     int n2_right=int((wf.n2)*frac_n2_right)+2;
     int n1_left=int((wf.n1)*frac_n1_left)+2;
     int n1_right=int((wf.n2)*frac_n2_left)+2;*/
      
     //int n1_SI_low = int(wf.n1/2) - n_SI;
     //int n1_SI_up = int(wf.n1/2) + n_SI +1;
     int n2_SI_low = int(wf.n2/2) - n_SI;
     int n2_SI_up = int(wf.n2/2) + n_SI +1;

     psi_1r.resize(wf.n1-start_n1_right+2); //wf in exchanged regoin in SI reoin
     for(int i=0; i<wf.n1-start_n1_right+2; i++) psi_1r[i].resize(n2_SI_up-n2_SI_low+2);
     psi_1l.resize(start_n1_left+2);
     for(int i=0; i<start_n1_left+2; i++) psi_1l[i].resize(n2_SI_up-n2_SI_low+2);

     for(int j=1;j<=wf.n2;j++)
     for(int i=1;i<=wf.n1;i++)
     {

	if( i > start_n1_right ){
          double xx = (wf.x1[i]-start_x1_right) / (wf.x1[wf.n1]-start_x1_right); 
          double vdt = v0*xx*xx;
          if( j>=n2_SI_low && j<=n2_SI_up ){  
              psi_1r[i-start_n1_right-1][j-n2_SI_low] = (1.-exp(-vdt))*wf.wave[wf.in2(j,i)];  //set up outer wf for SI
          }
          wf.wave[wf.in2(j,i)]*= exp(-vdt); // absorber for inner wf for both SI and DI regoins
        }

	if( i < start_n1_left ){
          double xx = (start_x1_left-wf.x1[i]) / (start_x1_left-wf.x1[1]); 
          double vdt = v0*xx*xx;
          if( j>=n2_SI_low && j<=n2_SI_up ){   
              psi_1l[i-1][j-n2_SI_low] = (1.-exp(-vdt))*wf.wave[wf.in2(j,i)];
          }
          wf.wave[wf.in2(j,i)]*= exp(-vdt);
        }

	if( j > start_n2_right ){  // only absorb in x2
          double xx = (wf.x2[j]-start_x2_right) / (wf.x2[wf.n2]-start_x2_right); 
          double vdt = v0*xx*xx;
          wf.wave[wf.in2(j,i)]*= exp(-vdt);
        }

	if( j < start_n2_left ){
          double xx = (start_x2_left-wf.x2[j]) / (start_x2_left-wf.x2[1]); 
          double vdt = v0*xx*xx;
          wf.wave[wf.in2(j,i)]*= exp(-vdt);
        }

     }  // end of grids

} // end Gauss


  
// mask function	
	void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n1_left, double frac_n2_left ,double exponent )
{
	if ( frac_n1_right>1. || frac_n1_right<0. )
	{
		cout<< "frac_n1_right is wrong";
	 	exit(1);
	}
	if ( frac_n2_right>1. || frac_n2_right<0. )
	{
		cout<< "frac_n2_right is wrong";
		exit(1);
	}
	if ( frac_n1_left>1. || frac_n1_left<0. )
	{
		cout<< "frac_n1_left is wrong";
		exit(1);
	}
	if ( frac_n2_left>1. || frac_n2_left<0. )
	{
		cout<< "frac_n2_left is wrong";
		exit(1);
	}
		
     double mask_start_x1_right=wf.x1[int((wf.n1)*(1.-frac_n1_right))];
     double mask_start_x2_right=wf.x2[int((wf.n2)*(1.-frac_n2_right))];
     double mask_start_x1_left=wf.x1[int((wf.n1)*frac_n1_left)+1];
     double mask_start_x2_left=wf.x2[int((wf.n2)*frac_n2_left)+1];
     double argument_x1_right;
     double argument_x2_right;
     double argument_x1_left;
     double argument_x2_left;
     double mask_x1_right;
     double mask_x2_right;
     double mask_x1_left;
     double mask_x2_left;
		
     for(int j=1;j<=wf.n2;j++)
     for(int i=1;i<=wf.n1;i++)
     {
				
	if(wf.symmetry_x1==0 || wf.symmetry_x1==2)
	{
	  if(i< int(wf.n1*frac_n1_left)){
	    argument_x1_left=(pi/2.)*(wf.x1[i]-mask_start_x1_left)/(wf.x1[1]-mask_start_x1_left+1.e-20);
	    mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
	    wf.wave[wf.in2(j,i)]*=mask_x1_left;
          }
	} //if sym==1 or sym==-1, we needn't this part.
				
	if(i> int(wf.n1*(1.-frac_n1_right))){
	  argument_x1_right=(pi/2.)*(wf.x1[i]-mask_start_x1_right)/(wf.x1[wf.n1]-mask_start_x1_right+1.e-20);
	  mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
	  wf.wave[wf.in2(j,i)]*=mask_x1_right;
        }
				
	if (wf.symmetry_x2==0 || wf.symmetry_x2==2)
	{
	  if(j< int(wf.n2*frac_n2_left)){
	    argument_x2_left=(pi/2.)*(wf.x2[j]-mask_start_x2_left)/(wf.x2[1]-mask_start_x2_left+1.e-20);
	    mask_x2_left=pow(fabs(cos(argument_x2_left)),exponent);
	    wf.wave[wf.in2(j,i)]*=mask_x2_left;
          }
	}//if sym==1 or sym==-1, we needn't this part.
			
	if (j> int(wf.n2*(1.-frac_n2_right))){
	   argument_x2_right=(pi/2.)*(wf.x2[j]-mask_start_x2_right)/(wf.x2[wf.n2]-mask_start_x2_right+1.e-20);
	   mask_x2_right=pow(fabs(cos(argument_x2_right)),exponent);
	   wf.wave[wf.in2(j,i)]*=mask_x2_right;
        }
				
     }  // end of i,j

  }  // end of mask


	void Mask_Function_OMP(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n1_left, double frac_n2_left ,double exponent )
{
		
     double mask_start_x1_right=wf.x1[int((wf.n1)*(1.-frac_n1_right))];
     double mask_start_x2_right=wf.x2[int((wf.n2)*(1.-frac_n2_right))];
     double mask_start_x1_left=wf.x1[int((wf.n1)*frac_n1_left)+1];
     double mask_start_x2_left=wf.x2[int((wf.n2)*frac_n2_left)+1];
     double argument_x1_right;
     double argument_x2_right;
     double argument_x1_left;
     double argument_x2_left;
     double mask_x1_right;
     double mask_x2_right;
     double mask_x1_left;
     double mask_x2_left;
     int j=1;
		
//#pragma omp parallel shared(mask_start_x1_right, mask_start_x2_right, mask_start_x1_left, mask_start_x2_left)
//{
//    #pragma omp for private(j, argument_x1_right, argument_x2_right, argument_x1_left, argument_x2_left, mask_x1_right, mask_x2_right, mask_x1_left, mask_x2_left)

     for(j=1;j<=wf.n2;j++){
      //int ip=omp_get_thread_num();
      //cout<<"--Mask-- Thread "<<ip<<": j="<<j<<endl;

      for(int i=1;i<=wf.n1;i++)
      {
				
	if(wf.symmetry_x1==0 || wf.symmetry_x1==2)
	{
	  if(i< int(wf.n1*frac_n1_left)){
	    argument_x1_left=(pi/2.)*(wf.x1[i]-mask_start_x1_left)/(wf.x1[1]-mask_start_x1_left+1.e-20);
	    mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
	    wf.wave[wf.in2(j,i)]*=mask_x1_left;
          }
	} //if sym==1 or sym==-1, we needn't this part.
				
	if(i> int(wf.n1*(1.-frac_n1_right))){
	  argument_x1_right=(pi/2.)*(wf.x1[i]-mask_start_x1_right)/(wf.x1[wf.n1]-mask_start_x1_right+1.e-20);
	  mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
	  wf.wave[wf.in2(j,i)]*=mask_x1_right;
        }
				
	if (wf.symmetry_x2==0 || wf.symmetry_x2==2)
	{
	  if(j< int(wf.n2*frac_n2_left)){
	    argument_x2_left=(pi/2.)*(wf.x2[j]-mask_start_x2_left)/(wf.x2[1]-mask_start_x2_left+1.e-20);
	    mask_x2_left=pow(fabs(cos(argument_x2_left)),exponent);
	    wf.wave[wf.in2(j,i)]*=mask_x2_left;
          }
	}//if sym==1 or sym==-1, we needn't this part.
			
	if (j> int(wf.n2*(1.-frac_n2_right))){
	   argument_x2_right=(pi/2.)*(wf.x2[j]-mask_start_x2_right)/(wf.x2[wf.n2]-mask_start_x2_right+1.e-20);
	   mask_x2_right=pow(fabs(cos(argument_x2_right)),exponent);
	   wf.wave[wf.in2(j,i)]*=mask_x2_right;
        }
		
      }		
     }  // end of i,j

//}  // end omp

  }  // end of mask

	void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf)
	{
		
		//int shifter_n1=(wf.n1-small_wf.n1)/2;
		for (int i=1; i<=2; i++) {
			if(small_wf.n[i] > wf.n[i]) {
				cerr<<"cannot load bigger wavfunction into small grid"<<endl;
				exit(1);
			}
		}
		
		int shifter_n1=0;
		int shifter_n2=0;
		
		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;  //even number
		}    
		if (wf.symmetry_x2==0)
		{
			shifter_n2=(wf.n2-small_wf.n2)/2;
		}
		
                int N_left=240; //280;
		if (wf.symmetry_x1==2 && small_wf.symmetry_x1==0)  //partial grids of wf and sym of small wf
		{
			shifter_n1=int(N_left-small_wf.n1/2);  //even number
		}    
		if (wf.symmetry_x1==2 && small_wf.symmetry_x1==2)  //partial grids for two wfs
		{
			shifter_n1=0;  //no shift
		}    
		
		for(int j=0;j<=wf.n2+1;j++)  // set 0
			for(int i=0;i<=wf.n1+1;i++)
			{
				wf.wave[wf.in2(j,i)]=complex(0.,0.);
			}
		
		for(int j=1;j<=small_wf.n2;j++)
			for(int i=1;i<=small_wf.n1;i++)
			{
				wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)]=small_wf.wave[small_wf.in2(j,i)];
			}
		
	}


	complex Project( wavefunction &wf , wavefunction &wf2)
	{
		
		complex projection=complex(0.,0.);
		
		for(int j=1;j<=wf.n2;j++)
			for(int i=1;i<=wf.n1;i++)
			{
				projection+= conj(wf2.wave[wf2.in2(j,i)])*wf.wave[wf.in2(j,i)];
			}
		projection*=wf.dx1*wf.dx2;
		return projection;
	}
	
	double Overlap( wavefunction &wf , wavefunction &wf2)
	{
		double overlap_;
		complex proj = Project ( wf , wf2 );
		overlap_ = sqrt( real( proj*conj( proj ) ) );
		return overlap_;
	}
	
	complex Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf)
	{
		int shifter_n1=0;
		int shifter_n2=0;
		int shifter_n3=0;//if sym==0 or 2, shift_n==0;

		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;
		}    
		if (wf.symmetry_x2==0)
		{
			shifter_n2=(wf.n2-small_wf.n2)/2;
		}
		
                int N_left=240; //280;
		if (wf.symmetry_x1==2)  //partial grids
		{
			shifter_n1=int(N_left-small_wf.n1/2);  //even number
		}    
		
		complex projection=complex(0.,0.);
		
		
		for(int j=1;j<=small_wf.n2;j++)
			for(int i=1;i<=small_wf.n1;i++)
			{
				projection+= conj(small_wf.wave[small_wf.in2(j,i)])
					*wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)];
			}
		projection*=small_wf.dx1*small_wf.dx2;
		return projection;
	}
	
	
	void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf)
	{
		int shifter_n1=wf.n1-small_wf.n1;
		int shifter_n2=wf.n2-small_wf.n2;
	
		
		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;
		}    
		if (wf.symmetry_x2==0)
		{
			shifter_n2=(wf.n2-small_wf.n2)/2;
		}

                int N_left=240; //280;
		if (wf.symmetry_x1==2)  //partial grids
		{
			shifter_n1=int(N_left-small_wf.n1/2);  //even number
		}    
		
		complex projection=complex(0.,0.);   
		
		
		for(int j=1;j<=small_wf.n2;j++)
			for(int i=1;i<=small_wf.n1;i++)
			{
				projection+= conj(small_wf.wave[small_wf.in2(j,i)])
					*wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)];
			}
		projection*=small_wf.dx1*small_wf.dx2;
		
		
		for(int j=1;j<=small_wf.n2;j++)
			for(int i=1;i<=small_wf.n1;i++)
			{
				wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)]-=projection*small_wf.wave[small_wf.in2(j,i)];
			}
		
		//    return projection;
	}
	


vector<double> Distribution_X1( wavefunction &wf )
{
  vector<double> dist(wf.n1+2,0.);
  double sum;
  for( int i = 1 ; i <= wf.n1 ; i++) 
    {
      sum = 0.;
      
      for( int j = 1 ; j <= wf.n2 ; j++ ) 
      { 
	      sum += norm(wf.wave[ wf.in2(  j , i ) ]);
      }
	
      dist[i] = wf.dx2*sum;
    }
  return dist;
  
}


vector<double> Distribution_X2( wavefunction &wf )
{
  vector<double> dist(wf.n2+2,0.);
  double sum;
  for( int j = 1 ; j <= wf.n2 ; j++) 
    {
      sum = 0.;
      for( int i = 1 ; i <= wf.n1 ; i++ ) 
      { 
	      sum += norm(wf.wave[ wf.in2( j , i ) ]);
      }
      
      dist[j] = wf.dx1*sum;
    }
  return dist;
  
}




/*void Write_2DGraph(vector<double> &x, vector<double> &y, string name )
{

  std::ofstream txt_ofs;

  if(x.size()!= y.size()) { cerr<<"x.size()!= y.size() in write_xy !"<<endl; exit(1);}
  string txtname=name;
  txt_ofs.open(txtname.c_str());
  if(txt_ofs.fail()) {
    cerr<<"Could not open output file "<< txtname<<endl;
    exit(1);
  }
  txt_ofs.setf(std::ios::showpoint|std::ios::scientific);
  ParameterMap pMap;
  txt_ofs << pMap;              // write all parameters in the header.
  //  txt_ofs << "#  time " << label << endl;
  for(int i=1; i<x.size(); ++i) {
    txt_ofs << x[i] << "\t" << y[i] << endl;
  }
  txt_ofs.close();
}*/




double Obs_Expectation_Value_r( wavefunction &wf )
  {
    double obs;
    double exp_r = 0.;
    
    
    for( int j = 1 ; j <= wf.n2 ; j++)
    {
	    for( int i = 1 ; i <= wf.n1 ; i++)
	    {
		    
		    exp_r += sqrt(wf.x1[i]*wf.x1[i]+wf.x2[j]*wf.x2[j])*wf.x1[i]*norm( wf.wave[ wf.in2(  j , i ) ] );
	    }
    }
  
    exp_r *=wf.dx1*wf.dx2;
    obs = exp_r;
    return obs;
    
  } // end of Obs_Expectation_Value_r



double Obs_Expectation_Value_Inv_r( wavefunction &wf )
  {
    double obs;
    double exp_inv_r = 0.;
    
    for( int j = 1 ; j <= wf.n2 ; j++)
    {
	    for( int i = 1 ; i <= wf.n1 ; i++)
	    {
		    //Assuming that sqrt(wf.x1[i]*wf.x1[i]+wf.x2[j]*wf.x2[j])=zero is not a point on the grid
		    exp_inv_r += 1./sqrt(wf.x1[i]*wf.x1[i]+wf.x2[j]*wf.x2[j])*norm( wf.wave[ wf.in2(  j , i ) ] );
	    }
    }
  
    exp_inv_r *=wf.dx1*wf.dx2;
    obs = exp_inv_r;
    return obs;
    
  } // end of Obs_Expectation_Value_Inv_r



// shaohao 2012.05
// open file
   void Open_File( ofstream &filename, string name )
   {
     filename.open(name.c_str());
     filename.setf(ios::showpoint|ios::scientific);
   }

// for energy spec
   void Psi_edge( wavefunction wf, double frac_n1_right, double frac_n1_left, double x_SI, double x_edge, vector<complex> &psi_right, vector<complex> &psi_left){

     int i_right = int(wf.n1*(1.-frac_n1_right))-int(x_edge/wf.dx1);  // edge point
     int i_left = int(wf.n1*frac_n1_left)+int(x_edge/wf.dx1);  // edge point
     int index=0;

     for(int j=1; j<wf.n2; j++){  // interation in x2
       if(abs(wf.x2[j])<x_SI){
          psi_right[index] = wf.wave[wf.in2(j,i_right)];
          psi_left[index] = wf.wave[wf.in2(j,i_left)];
          index += 1;
       } // end if 
     }  //end j

   }

// tsurff treatment at a time instant
   void TSURFF( wavefunction wf, double frac_n1_right, double frac_n1_left, double x_SI, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left){

     // retrive 3 points of wf before the boundary in SI regoin    
     int i_right = int(wf.n1*(1.-frac_n1_right))-int(x_edge/wf.dx1);  // edge point
     int i_left = int(wf.n1*frac_n1_left)+int(x_edge/wf.dx1);  // edge point
     int n_SI=2*int(x_SI/wf.dx2);
     int n2_SI_low = int(wf.n2/2) - n_SI/2 + 1;
     int n2_SI_up = int(wf.n2/2) + n_SI/2;
     vector<vector<complex> > psi_right, psi_left; 
     psi_right.resize(n_SI);
     for(int i=0; i<n_SI; i++) psi_right[i].resize(3);
     psi_left.resize(n_SI);
     for(int i=0; i<n_SI; i++) psi_left[i].resize(3);

     for(int j=n2_SI_low; j<=n2_SI_up; j++){  // in SI regoin
         psi_right[j-n2_SI_low][0] = wf.wave[wf.in2(j,i_right-1)];
         psi_left[j-n2_SI_low][0] = wf.wave[wf.in2(j,i_left-1)];
         psi_right[j-n2_SI_low][1] = wf.wave[wf.in2(j,i_right)];
         psi_left[j-n2_SI_low][1] = wf.wave[wf.in2(j,i_left)];
         psi_right[j-n2_SI_low][2] = wf.wave[wf.in2(j,i_right+1)];
         psi_left[j-n2_SI_low][2] = wf.wave[wf.in2(j,i_left+1)];
     }  //end j

     double fac=pow(2.*pi, -1.5);     

/*     for(int k=0; k<np; k++){ 
        amp_right[k]=complex(0.,0.);
        amp_left[k]=complex(0.,0.);
     }

     for(int isi=0; isi<n_SI; isi++ ){  // in SI regoin
        complex right_2nd = (psi_right[isi][2] - 2.*psi_right[isi][1] + psi_right[isi][0])/(wf.dx1*wf.dx1);
        complex left_2nd = (psi_left[isi][2] - 2.*psi_left[isi][1] + psi_left[isi][0])/(wf.dx1*wf.dx1);
        complex right_1st = (psi_right[isi][2] - psi_right[isi][0])/(2.*wf.dx1);
        complex left_1st = (psi_left[isi][2] - psi_left[isi][0])/(2.*wf.dx1);
        int k;
#pragma omp parallel shared(fac,I,right_2nd,right_1st,left_2nd,left_1st,i_right,i_left)
{
    #pragma omp for private(k)
        for(k=0; k<np; k++){  // coherent sum over SI regoin
          amp_right[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_right])*(-0.5*right_2nd+I*Afield*right_1st);
          amp_left[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_left])*(-0.5*left_2nd+I*Afield*left_1st);
        } // end of momentum
}  // end omp 
     } // end SI */

     int isi, k;
     vector<complex> right_2nd(n_SI), left_2nd(n_SI), right_1st(n_SI), left_1st(n_SI);
//#pragma omp parallel shared(n_SI, right_2nd, left_2nd, right_1st, left_1st, psi_right, psi_left)
//{
//    #pragma omp for private(isi)
     for(isi=0; isi<n_SI; isi++ ){  // in SI regoin
        right_2nd[isi] = (psi_right[isi][2] - 2.*psi_right[isi][1] + psi_right[isi][0])/(wf.dx1*wf.dx1);
        left_2nd[isi] = (psi_left[isi][2] - 2.*psi_left[isi][1] + psi_left[isi][0])/(wf.dx1*wf.dx1);
        right_1st[isi] = (psi_right[isi][2] - psi_right[isi][0])/(2.*wf.dx1);
        left_1st[isi] = (psi_left[isi][2] - psi_left[isi][0])/(2.*wf.dx1);
     } // end of SI
//}  // end omp

//#pragma omp parallel shared(n_SI,I,fac,right_2nd,left_2nd,right_1st,left_1st,i_right,i_left)
//#pragma omp parallel shared(n_SI,fac,right_2nd,left_2nd,right_1st,left_1st,i_right,i_left)
//{
//    #pragma omp for private(k)
     for(k=0; k<np; k++){  // coherent sum over SI regoin
        amp_right[k]=complex(0.,0.);
        amp_left[k]=complex(0.,0.);
        for(int i=0; i<n_SI; i++ ){  // in SI regoin
     amp_right[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_right])*(-0.5*right_2nd[i]+I*Afield*right_1st[i]);
     amp_left[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_left])*(-0.5*left_2nd[i]+I*Afield*left_1st[i]);
        } // end SI
     } // end momentum
//}  // end omp 

   }  // end tsurff


// tsurff treatment at a time instant
   void TSURFF_PRJ( wavefunction wf, wavefunction wf_ion, double frac_n1_right, double frac_n1_left, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left){

     // retrive 3 points of wf before the boundary
     int i_right = int(wf.n1*(1.-frac_n1_right))-int(x_edge/wf.dx1);  // edge point
     int i_left = int(wf.n1*frac_n1_left)+int(x_edge/wf.dx1);  // edge point
     vector<complex> psi_right(3), psi_left(3); 
     int j_shift = int((wf.n1-wf_ion.n1)/2);  // shift for diff size in X1

     for(int i=0; i<3; i++){  // only for 3 points near surface
       psi_right[i]=complex(0.,0.);
       psi_left[i]=complex(0.,0.);
       for(int j=wf_ion.n2/4; j<=3*wf_ion.n2/4; j++){  // project to ion ground state in X2
          psi_right[i] += conj(wf_ion.wave[wf_ion.in2(j,1)]) * wf.wave[wf.in2(j+j_shift,i_right+i-1)];
          psi_left[i] += conj(wf_ion.wave[wf_ion.in2(j,1)]) * wf.wave[wf.in2(j+j_shift,i_left+i-1)];
       }  //end j
       psi_right[i] *= wf_ion.dx2;
       psi_left[i] *= wf_ion.dx2;
     }  //end i

     // initial for tsurff
     double fac=pow(2.*pi, -1.5);     
     int k;

     // calc FD
     complex right_2nd = (psi_right[2] - 2.*psi_right[1] + psi_right[0])/(wf.dx1*wf.dx1);
     complex left_2nd = (psi_left[2] - 2.*psi_left[1] + psi_left[0])/(wf.dx1*wf.dx1);
     complex right_1st = (psi_right[2] - psi_right[0])/(2.*wf.dx1);
     complex left_1st = (psi_left[2] - psi_left[0])/(2.*wf.dx1);
     complex temp3r=-0.5*right_2nd+I*Afield*right_1st;
     complex temp3l=-0.5*left_2nd+I*Afield*left_1st;
     double epsilon=1e-15;
     if(norm(temp3r)>epsilon && norm(temp3l)>epsilon) cout<<"### output test ###"<<endl;

//#pragma omp parallel shared(fac,right_2nd,left_2nd,right_1st,left_1st,i_right,i_left)
//{
//    #pragma omp for private(k)
    for(k=0; k<np; k++){  // momentum
      complex temp1=exp(I*volkov_phase[k]);
      complex temp2r=exp(-I*p[k]*wf.x1[i_right]);
      complex temp2l=exp(-I*p[k]*wf.x1[i_left]);
      amp_right[k] = fac*temp1*temp2r*temp3r;
      amp_left[k] = fac*temp1*temp2l*temp3l;
      if(norm(temp3r)>epsilon && norm(temp3l)>epsilon && k%100==0) cout<<p[k]<<' '<<temp1<<' '<<temp2r<<' '<<temp3r<<' '<<temp2l<<' '<<temp3l<<endl;
      //amp_right[k] = fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_right])*(-0.5*right_2nd+I*Afield*right_1st);
      //amp_left[k]  = fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_left])*(-0.5*left_2nd+I*Afield*left_1st);
    } // end of momentum
//}  // end omp

   }  // end tsurff prj


// tsurff treatment at a time instant
   void TSURFF_1D( wavefunction wf, double frac_n1, double x_edge, vector<double> p, int np, vector<double> volkov_phase,  double Afield, vector<complex> &amp_right, vector<complex> &amp_left){

     // retrive 3 points of wf before the boundary in SI regoin    
     int i_right = int(wf.n1*(1.-frac_n1))-int(x_edge/wf.dx1);  // edge point
     int i_left = int(wf.n1*frac_n1)+int(x_edge/wf.dx1);  // edge point
     vector<complex> psi_right, psi_left; 
     psi_right.resize(3);
     psi_left.resize(3);

     psi_right[0] = wf.wave[wf.in2(1,i_right-1)];
     psi_left[0] = wf.wave[wf.in2(1,i_left-1)];
     psi_right[1] = wf.wave[wf.in2(1,i_right)];
     psi_left[1] = wf.wave[wf.in2(1,i_left)];
     psi_right[2] = wf.wave[wf.in2(1,i_right+1)];
     psi_left[2] = wf.wave[wf.in2(1,i_left+1)];

     // initial for tsurff
     double fac=pow(2.*pi, -1.5);     
     for(int k=0; k<np; k++){
        amp_right[k]=complex(0.,0.);
        amp_left[k]=complex(0.,0.);
     }

     // tsurff treatment
    complex right_2nd = (psi_right[2] - 2.*psi_right[1] + psi_right[0])/(wf.dx1*wf.dx1);
    complex left_2nd = (psi_left[2] - 2.*psi_left[1] + psi_left[0])/(wf.dx1*wf.dx1);
    complex right_1st = (psi_right[2] - psi_right[0])/(2.*wf.dx1);
    complex left_1st = (psi_left[2] - psi_left[0])/(2.*wf.dx1);
    for(int k=0; k<np; k++){  // coherent sum over SI regoin
       amp_right[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_right])*(-0.5*right_2nd+I*Afield*right_1st);
       amp_left[k] += fac*exp(I*volkov_phase[k])*exp(-I*p[k]*wf.x1[i_left])*(-0.5*left_2nd+I*Afield*left_1st);
    } // end of momentum

   }  // end tsurff_1D


// calc dipole
	void Obs_dipole( wavefunction wf, complex &dip)
	{
                dip = complex(0.,0.);
		
		for( int j = 1 ; j <= wf.n2 ; j++ )
		{
			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				dip -= conj(wf.wave[ wf.in2(j,i) ]) * ( wf.x1[ i ] + wf.x2[ i ] )
					* wf.wave[ wf.in2(j,i) ];
			}
		}
	
		dip *= wf.dx1*wf.dx2;
		
	} // end of Obs_Expectation_Value_X1

// calc dipole between ground and excited states
        void Obs_dipole_gex( wavefunction wf, wavefunction small_wf, complex &dip)
        {
                dip = complex(0.,0.);
		int shifter_n1=wf.n1-small_wf.n1;
		int shifter_n2=wf.n2-small_wf.n2;

                ProjectOUT_Diff_Sizes(wf, small_wf);  // wf2 grids <= wf1 grids
                
                for( int j = 1 ; j <= small_wf.n2 ; j++ )
                {
                        for( int i = 1 ; i <= small_wf.n1 ; i++ )
                        {
                                dip -= ( small_wf.x1[ i ] + small_wf.x2[ i ] )
                                        * conj(wf.wave[wf.in2(shifter_n2+j,shifter_n1+i)]) * small_wf.wave[ small_wf.in2(j,i) ];
                        }
                }
        
                dip *= wf.dx1*wf.dx2;
                
        } // end of Obs_Expectation_Value_X1 


// transfer wavefuntion from length gauge to velocity gauge
   void length_to_velocity (wavefunction ){


   }

} // end of namespace Cartesian_2D
