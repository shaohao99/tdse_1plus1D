// $Id: Laser.h 315 2006-07-18 21:41:15Z arvid $
#ifndef POTENTIALS_H
#define POTENTIALS_H

using namespace std;


namespace Cartesian_2D {

  struct ABVparam {
  
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x1;
    double ABV_B_x1;
  
    double ABV_A_x2;
    double ABV_B_x2;
  
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;
    vector<double> ABV_V_deriv;

  };


  /***********************************************************************************************/
  /***********************************************************************************************/
  
  //Potential for Hydorgen like  2D

  struct Potential_H_2D : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  double mass_nucleus;
	  double softening_en;
	  
	  /* Constructor to initialize the default values: */
	  Potential_H_2D(const wavefunction &wf) {
		    
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;
		  
		  charge_nucleus=1;
		  mass_nucleus=1837;
		  softening_en=1;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;
		  
		  updateABV(wf);
	  };

	  /* Constructor to input values: */
	  Potential_H_2D(const wavefunction &wf, const double charge_ion, const double soft_en) {
		
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;

		  /* These parameters can be changed */
		  charge_nucleus=charge_ion;		
		  softening_en=1;
		  
		  /* The rest stays the same */
		  mass_nucleus=1837;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;
		  
		  updateABV(wf);
	  };


	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ] +softening_en );	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };//End Potential Hydrogen 2D

  
  /***********************************************************************************************/
  /***********************************************************************************************/
// shaohao2009.04
  //Potential for Helium in 1D in Jacobi Coordinates: x1 is Z=(z1+z2)/2., x2 is z=z1-z2
  struct Helium_2D_Jacobi : ABVparam
  {
	  double charge_1,charge_2;
	  double charge_nucleus;
	  double mass_1;
	  double mass_2;
	  double mass_nucleus;
	  double softening_en;
	  double softening_ee;
	  
//          cout << "---sh1---" ;
	  /* Constructor to initialize the default values: */
	  Helium_2D_Jacobi(const wavefunction &wf) {
		  
		  charge_1=charge_2=-1.;
		  charge_nucleus=2;
		  mass_1=mass_2=1;
		 
		  softening_en=1.;
		  softening_ee=1.;
		  
		  double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
		  
		  double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */
		  
		  ABV_A_x1 =-0.5/center_of_mass;
		  ABV_A_x2 =-0.5/reduced_mass;
		  
		  ABV_B_x1 = center_of_mass_charge/center_of_mass;
		  ABV_B_x2 = reduced_charge/reduced_mass;	 
		  
		  updateABV(wf);
	  };

//          cout << "---sh2---" ;
	  /* Constructor to initialize the default values: */
	  Helium_2D_Jacobi(const wavefunction &wf, const double q1, const double q2, const double m1, const double m2, 
			   const double soft_en, const double soft_ee, const double Q ) {
		  
		  charge_1=q1;
		  charge_2=q2;
		  charge_nucleus=Q;
		  mass_1=m1;
		  mass_2=m2;
	
		  softening_en=soft_en;
		  softening_ee=soft_ee;
		  
		  double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
		  
		  double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */
		  
		  ABV_A_x1 =-0.5/center_of_mass;
		  ABV_A_x2 =-0.5/reduced_mass;
		  
		  ABV_B_x1 = center_of_mass_charge/center_of_mass;
		  ABV_B_x2 = reduced_charge/reduced_mass;	 
		  
		  updateABV(wf);
	  };

//          cout << "---sh3---" ;
	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] = 
					  charge_1*charge_2/sqrt( wf.x2[ j ]*wf.x2[ j ]+softening_ee )
					  +charge_1*charge_nucleus/sqrt( ( wf.x1[ i ]+wf.x2[ j ]/2. )*( wf.x1[ i ]+wf.x2[ j ]/2. )+softening_en )
					  +charge_2*charge_nucleus/sqrt( ( wf.x1[ i ]-wf.x2[ j ]/2. )*( wf.x1[ i ]-wf.x2[ j ]/2. )+softening_en );
				  
				  ABV_V[ wf.in2( j , i ) ] /= 2.;		
			  }
		  }
	  }
  };
  
// shaohao2012.5

  //short range Coulomb potential
  struct He_Short_Coulomb : ABVparam
  {
    double c_en,c_ee,R0;
    double mass_1, mass_2;
	  
	  /* Constructor to initialize the default values: */
    He_Short_Coulomb(const wavefunction &wf) {
          c_en=0.4786;
          c_ee=0.3673;
          R0=75.;  // a.u.
          //v00=1./wf.dx1+1.;  //1./sqrt(0.399);
	  mass_1=mass_2=1.;
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
	  updateABV(wf);
    };

	  /* Constructor for input of parameters: */
//    He_Short_Coulomb(const wavefunction &wf, const double m1, const double m2, const double c_en_, const double c_ee_, const double v00_, const double R0_) {
    He_Short_Coulomb(const wavefunction &wf, const double m1, const double m2, const double c_en_, const double c_ee_, const double R0_) {
	  c_en=c_en_;
	  c_ee=c_ee_;
          R0=R0_; 
          //v00=v00_;
	  mass_1=m1;
	  mass_2=m2;
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
	  updateABV(wf);
   };

   void updateABV(const wavefunction &wf) {
		  
	  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
          
          double R_en=sqrt(R0*R0+c_en);
          double R_ee=sqrt(R0*R0+c_ee);

	  for( int j=1 ; j <= wf.n2 ; j++ ) {
	     for( int i=1 ; i <= wf.n1 ; i++ ) {
               double v_en1=0., v_en2=0., v_ee=0.;
               /*if(abs(wf.x1[i])<=R0) v_en1=c_en*( -2./abs(wf.x1[i])-wf.x1[i]*wf.x1[i]/(R0*R0*R0)+3./R0 );
               if(abs(wf.x2[i])<=R0) v_en2=c_en*( -2./abs(wf.x2[j])-wf.x2[j]*wf.x2[j]/(R0*R0*R0)+3./R0 );
               if(abs(wf.x1[i]-wf.x2[j])<=R0 && (wf.x1[i]-wf.x2[j])!=0.) v_ee=c_ee*( 1./abs(wf.x1[i]-wf.x2[j])+(wf.x1[i]-wf.x2[j])*(wf.x1[i]-wf.x2[j])/(2.*R0*R0*R0)-3./(2.*R0) );
               if((wf.x1[i]-wf.x2[j])==0.) v_ee=v00;*/
               if(abs(wf.x1[i])<=R0){
                double x=sqrt(wf.x1[i]*wf.x1[i]+c_en);
                v_en1= -2./x - x*x/(R_en*R_en*R_en) + 3./R_en;
               }
               if(abs(wf.x2[j])<=R0){
                double x=sqrt(wf.x2[j]*wf.x2[j]+c_en);
                v_en2= -2./x - x*x/(R_en*R_en*R_en) + 3./R_en;
               }
               if(abs(wf.x1[i]-wf.x2[j])<=R0){
                double x=sqrt((wf.x1[i]-wf.x2[j])*(wf.x1[i]-wf.x2[j])+c_ee);
                v_ee= 1./x + x*x/(2.*R_ee*R_ee*R_ee) + 3./(2.*R_ee);
               }
	       ABV_V[ wf.in2(j,i) ] = v_en1 + v_en2 + v_ee;
	       ABV_V[ wf.in2(j,i) ] /= 2.;		
	     }
	  }
   } // end of updateABV

  };  // end short Coulomb


  //short range Coulomb potential for He+
  struct Heplus_Short_Coulomb : ABVparam
  {
    double c_en,R0;
    double mass_1, mass_2;
	  
	  /* Constructor to initialize the default values: */
    Heplus_Short_Coulomb(const wavefunction &wf) {
          c_en=0.4786;
          R0=75.;  // a.u.
	  mass_1=mass_2=1.;
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
	  updateABV(wf);
    };

	  /* Constructor for input of parameters: */
    Heplus_Short_Coulomb(const wavefunction &wf, const double m1, const double m2, const double c_en_, const double R0_) {
	  c_en=c_en_;
          R0=R0_; 
	  mass_1=m1;
	  mass_2=m2;
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
	  updateABV(wf);
   };

   void updateABV(const wavefunction &wf) {
		  
	  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);

          double R_en=sqrt(R0*R0+c_en);

	  for( int j=1 ; j <= wf.n2 ; j++ ) {
	     for( int i=1 ; i <= wf.n1 ; i++ ) {
               double v_en1=0., v_en2=0.;
               if(abs(wf.x1[i])<=R0){
                double x=sqrt(wf.x1[i]*wf.x1[i]+c_en);
                v_en1= -2./x - x*x/(R_en*R_en*R_en) + 3./R_en;
               }
               if(abs(wf.x2[j])<=R0){
                double x=sqrt(wf.x2[j]*wf.x2[j]+c_en);
                v_en2= -2./x - x*x/(R_en*R_en*R_en) + 3./R_en;
               }
	       ABV_V[ wf.in2(j,i) ] = v_en1 + v_en2;
	       ABV_V[ wf.in2(j,i) ] /= 2.;		
	     }
	  }
   } // end of updateABV

  };  // end short Coulomb for He+


  //Potential for Helium in 1D in Cartesian Coordinates: x1 is z1, x2 is z2
  struct Helium_1plus1D : ABVparam
  {
    double charge_1,charge_2;
    double charge_nucleus;
    double charge_repulsive;
    double mass_1;
    double mass_2;
    double softening_en;
    double softening_ee;
	  
	  /* Constructor to initialize the default values: */
    Helium_1plus1D(const wavefunction &wf) {
		  
	  charge_1=charge_2=-1.;
	  charge_nucleus=2;
          charge_repulsive=1.;
	  mass_1=mass_2=1.;
	  softening_en=0.5;
	  softening_ee=0.339;
		  
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
		  
	  updateABV(wf);
    };

	  /* Constructor for input of parameters: */
    Helium_1plus1D(const wavefunction &wf, const double q1, const double q2, const double m1, const double m2, const double Q, const double C, const double soft_en, const double soft_ee ) {
		  
	  charge_1=q1;
	  charge_2=q2;
	  charge_nucleus=Q;
          charge_repulsive=C;
	  mass_1=m1;
	  mass_2=m2;
	  softening_en=soft_en;
	  softening_ee=soft_ee;
		  
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;

	  updateABV(wf);
   };

   void updateABV(const wavefunction &wf) {
		  
	  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);

	  for( int j=1 ; j <= wf.n2 ; j++ ) {
	     for( int i=1 ; i <= wf.n1 ; i++ ) {
	       ABV_V[ wf.in2(j,i) ] = 
	         charge_repulsive*charge_1*charge_2/sqrt( (wf.x1[i]-wf.x2[j])*(wf.x1[i]-wf.x2[j]) + softening_ee)
	         + charge_1*charge_nucleus/sqrt( wf.x1[i]*wf.x1[i] + softening_en )
	         + charge_2*charge_nucleus/sqrt( wf.x2[j]*wf.x2[j] + softening_en );
	       ABV_V[ wf.in2( j , i ) ] /= 2.;		
	     }
	  }
   } // end of updateABV

  };  // end of He_1plus1D


// shaohao2012.5
  struct He_ion_1plus1D : ABVparam
  {
    double charge_1,charge_2;
    double charge_nucleus;
    double mass_1;
    double mass_2;
    double softening_en;
	  
	  /* Constructor to initialize the default values: */
    He_ion_1plus1D(const wavefunction &wf) {
		  
	  charge_1=charge_2=-1.;
	  charge_nucleus=2;
	  mass_1=mass_2=1.;
	  softening_en=0.5;
		  
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
		  
	  updateABV(wf);
    };

	  /* Constructor for input of parameters: */
    He_ion_1plus1D(const wavefunction &wf, const double q1, const double q2, const double m1, const double m2, const double Q, const double soft_en ) {
		  
	  charge_1=q1;
	  charge_2=q2;
	  charge_nucleus=Q;
	  mass_1=m1;
	  mass_2=m2;
	  softening_en=soft_en;
		  
	  ABV_A_x1 =-0.5/mass_1;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5/mass_2;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;

	  updateABV(wf);
   };

   void updateABV(const wavefunction &wf) {
		  
	  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);

	  for( int j=1 ; j <= wf.n2 ; j++ ) {
	     for( int i=1 ; i <= wf.n1 ; i++ ) {
	       ABV_V[ wf.in2(j,i) ] = charge_1*charge_nucleus/sqrt( wf.x1[i]*wf.x1[i] + softening_en )
	         + charge_2*charge_nucleus/sqrt( wf.x2[j]*wf.x2[j] + softening_en );
	       ABV_V[ wf.in2( j , i ) ] /= 2.;		
	     }
	  }
   } // end of updateABV

  };  // end of He_ion_1plus1D
  

// shaohao2012.10
  struct He_ion_X2: ABVparam
  {
    double charge;
    double softening_en;
	  
	  /* Constructor to initialize the default values: */
    He_ion_X2(const wavefunction &wf) {
		  
	  charge=2.;
	  softening_en=0.339;
		  
	  ABV_A_x1 =0.; //-0.5;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;
		  
	  updateABV(wf);
    };

	  /* Constructor for input of parameters: */
    He_ion_X2(const wavefunction &wf, const double Q, const double soft_en ) {
		  
	  charge=Q;
	  softening_en=soft_en;
		  
	  ABV_A_x1 =-0.;  // coefficient for d^2/dx^2
	  ABV_A_x2 =-0.5;
	  ABV_B_x1 = 0.;         // coefficient for d/dx, without field
	  ABV_B_x2 = 0.;

	  updateABV(wf);
   };

   void updateABV(const wavefunction &wf) {
		  
	  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);

	  for( int j=1 ; j <= wf.n2 ; j++ ) {
	     for( int i=1 ; i <= wf.n1 ; i++ ) {  // i=1 only
	       ABV_V[ wf.in2(j,i) ] = -charge/sqrt( wf.x2[j]*wf.x2[j] + softening_en );
	     }
	  }
   } // end of updateABV

  };  // end of He_ion_X2


  /***********************************************************************************************/
  /***********************************************************************************************/


  
} // end of Namespace Cartesian2D




#endif	/* POTENTIALS_H */


