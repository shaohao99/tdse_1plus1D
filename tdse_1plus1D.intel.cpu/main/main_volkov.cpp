#include <fftw3.h>
#include <string>
#include <vector>
#include "constant.h"
#include "Cartesian2D.h" 
#include "ParameterMap.h"
#include "InputStuff.h"
#include "Laser.h"
#include <fstream>
#include "Observable.h"
#include "utils.h"
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
//#include <gsl/gsl_sf_coulomb.h>

using namespace std;
using namespace Cartesian_2D;
using namespace Laser_fields;

Parameter<vector<int> > n_points( "n_points", vector<int>(2,100), "number of points for (X1,X2,X3)" );
Parameter<vector<int> > n_points_g( "n_points_g", vector<int>(2,100), "number of points for (X1,X2,X3)" );
Parameter<vector<double> > spatial_step( "spatial_step", vector<double>(2,0.3), "spatial_steps for (X1,X2,X3)" );

Parameter<double> dt_i("dt_i",0.02, "Time step for imaginary time propagation");
Parameter<double> dt_p("dt_p",0.02, "Time step for pulse 1");

Parameter<vector<double> > sigma("sigma", vector<double>(2,1.), "width of the pulse");
Parameter<vector<double> > gauss_center("gauss_center", vector<double>(2,0.), "Center of Gaussian");
Parameter<vector<double> > initial_momentum("initial_momentum", vector<double>(2,0.), "Initial momentum of wavepacket");

Parameter<double> soft_en("soft_en",0.5, "e-n soft coulomb parameter");
Parameter<double> soft_ee("soft_ee",0.339, "e-e soft coulomb parameter");
Parameter<double> Cr("Cr",1., "repulsive charge, 1 for He atom");

Parameter<double> GroundStateEnergyThreshold("GroundStateEnergyThreshold",1.e-9,"how long should we evolve the ground state"); 

Parameter<double> lambda1("lambda1",69.0, "Wavelength of the Laser (nm)"); 
Parameter<double> i1("i1",1.0e13, "Intensity of the Laser (W/cm2)");
Parameter<double> cep1("cep1",0., "fraction of Pi for the initial phase");
Parameter<double> n_c1("n_c1",15., "Number of cycles of the as pulse");
Parameter<double> i2("i2",1.0e13, "Intensity of the Laser (W/cm2)");
Parameter<double> n_c2("n_c2",3., "Number of cycles of the as pulse");
Parameter<double> cep2("cep2",0., "fraction of Pi for the initial phase");
Parameter<int> nenv("nenv",2, "Sin^n envelope");
Parameter<int> nhhg("nhhg",11, "Number of HHG");
Parameter<double> ohhg0("ohhg0",75.75, "Starting number of HHG");
Parameter<double> delay("delay",55., "Delay between pulses");  // unit: cycle
Parameter<double> t_wait("t_wait",10., "Time for free propogation");

Parameter<double> x_SI("x_SI",5., "Single ionization regoin.");
Parameter<double> fact("fact",2., "Multiple of time for FFT.");
Parameter<double> facout("facout",0.75, "Fraction of time for FFT output.");
Parameter<int> fft_sign("fft_sign",-1, "forward or backward fft");

Parameter<int> facinter("facinter",150, "How many time steps absorber is appled?");
Parameter<double> epsilon("epsilon",1.e-15, "Absorbing potential.");
Parameter<double> v0("v0",1., "Absorbing potential.");
Parameter<double> tf("tf",1000., "Final time for Volkov propagation."); // a.u.
Parameter<double> facp("facp",2., "Multiple of momentum for FFT.");
Parameter<double> Emin("Emin",1., "Min of energy spectrum");
Parameter<double> Emax("Emax",1., "Max of energy spectrum");
Parameter<double> dE("dE",0.02, "Max of energy spectrum");

Parameter<double> a_percent1("absorb_percent1",0.02,"on how many percents of the grids are the wf absorbed?");
Parameter<double> a_percent2("absorb_percent2",0.02,"on how many percents of the grids are the wf absorbed?");
Parameter<double> mask_exp("mask_exponent",0.167,"the exponent for cos function");

Parameter<string> jobname("jobname", "He_apt", "The Name of the job");
Parameter<string> outputFolder("outputFolder", "He_Cartesian2D/", "The basename of the output directory");  // just a default, is set in Input.cpp from jobname()
Parameter<string> wavefunctionFolder("wavefunctionFolder",".", "path to the directory of the waveFunction to load"); 



// FFT functions based on fftw3 libs
void fft_real_1D (const int n, vector<double> &inf, vector<double> &outf)
{
  fftw_plan p_real_1d = fftw_plan_r2r_1d(n, &inf[0], &outf[0], FFTW_RODFT11, FFTW_ESTIMATE);

  fftw_execute(p_real_1d);

  fftw_destroy_plan(p_real_1d);

}

// fft for complex 1D
void fft_complex_1D (const int &n, vector<complex> &input_f, vector<complex> &output_f, double dt, int fft_sign)
{

  fftw_complex *in;
  fftw_complex *out;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

  fftw_plan p_complex_1d = fftw_plan_dft_1d(n, &in[0], &out[0], FFTW_FORWARD, FFTW_ESTIMATE);
  //fftw_plan p_complex_1d = fftw_plan_dft_1d(n, &in[0], &out[0], FFTW_BACKWARD, FFTW_ESTIMATE);
  //fftw_plan p_complex_1d = fftw_plan_dft_1d(n, &in[0], &out[0], fft_sign, FFTW_ESTIMATE); 
                     //-1 forward, +1 backward

  for(int i=0; i<n; i++){
    in[i][0]=real(input_f[i]);
    in[i][1]=imag(input_f[i]);
  }

  fftw_execute(p_complex_1d);

  double dw=2.*pi/(double(n-1)*dt);
  double fac=dw; //dw/n;  // normalized factor
  
  /*for(int i=0; i<n; i++){ 
    w[i]=( i - int((n+1)/2) )*dw;
  }*/

  for(int i=0; i<n; i++){
    output_f[i] = out[i][0] + I*out[i][1];
  }

  fftw_destroy_plan(p_complex_1d);
  free(in);  free(out);

}


// main program
int main( int argc , char *argv[] )
{
	
// measure time
  struct timeval start_time, stop_time, elapsed_time;

  /*int ithread=omp_get_thread_num();
  cout<<"thread number ="<<ithread<<endl;*/

  cout << outputFolder() << endl;
  parseInput(argc, argv);           // input data
  cout << outputFolder() << endl;

// ----- wavefunction and Hamiltonian
  wavefunction psi_ground;
  Initialize( psi_ground, n_points_g(), spatial_step(), 0, 0);  //0: full axis; 1: half axis
  Helium_1plus1D  pot_ground(psi_ground, -1., -1., 1., 1., 2., Cr(), soft_en(), soft_ee());  
  Hamiltonian H_free(psi_ground);

  Guess_Function_Gauss( psi_ground , sigma() , gauss_center() , initial_momentum() );
  Normalize( psi_ground );
  cout << "Norm(guess)=" << Obs_Norm(psi_ground) << endl;

// ----- virtual time propagation to obtain ground state
  Observable<double, wavefunction, ABVparam>  energy_ground( "Energy_ground" , 1 , Obs_Energy , psi_ground, pot_ground );
  double energy_old = 0.;
  double energy_new = 1000.;
  int counter = 0;
  double error=0.;
  do{
      H_free.X2( -I*0.5*dt_i(), psi_ground, pot_ground ); 
      H_free.X1( -I*dt_i(), psi_ground, pot_ground );
      H_free.X2( -I*0.5*dt_i(), psi_ground, pot_ground ); 
      Normalize( psi_ground);
       
      energy_old = energy_new;
      energy_new = energy_ground.measure();  // get energy of ground state.  shaohao 2008.09
      error=abs(energy_new-energy_old);

       if ( counter%50 == 0 ){
           cout <<counter*dt_i()<<"\t"<< "Norm(ground)="<<Obs_Norm(psi_ground)<<"\t"<<"E(ground)="<<energy_new<<"\t"<<"error="<<error<<endl;
       }
       counter += 1;
   } while ( error > GroundStateEnergyThreshold() );

   //psi_ground.save("wf_ground");
   /*ofstream output_grnd;
   Open_File(output_grnd, outputFolder()+"wf_ground.dat");
   for(int i=1; i<=psi_ground.n1; i++)
     for(int j=1; j<=psi_ground.n2; j++){
     output_grnd<<psi_ground.x1[i]<<"\t"<<psi_ground.x2[j]<<"\t"<<
                  norm(psi_ground.wave[psi_ground.in2(j,i)]) <<endl;
   }*/


// ----- Laser
   vector<double> ohhg(nhhg()), ahhg(nhhg()), cephhg_pi(nhhg());
   for(int j=0; j<nhhg(); j++){
     ohhg[j]=ohhg0()+2.*double(j); // odd
     ahhg[j]=1./sqrt(double(nhhg()));
     //ahhg[j]=double((nhhg()-j)/nhhg()); //1./sqrt(double(nhhg()));  // 1/nhhg of intensity
     cephhg_pi[j]=cep2();
   }
   vector<double> time, E1, E2, Efield, Afield;
   int iend_pulses;
   IR_APT ( lambda1(), i1(), n_c1(), cep1(), nenv(), i2(), n_c2(), nhhg(), ohhg, ahhg, cephhg_pi, delay(), dt_p(),t_wait(), time, E1, E2, Efield, Afield, iend_pulses); 

   iend_pulses=iend_pulses+2;
   int nt=time.size();  // points of time

   /*vector<double> Amid_test(nt);
   for( int i = 0 ; i < nt-1 ; i++ )
       Amid_test[i]=(Afield[i+1]+Afield[i])/2.;

   ofstream output_field; 
   Open_File(output_field, outputFolder()+"/field_origin.dat");
   for(int i=0; i<nt; i++)
       //output_field<<time[i]<<'\t'<<E1[i]<<'\t'<<E2[i]<<"\t"<<Efield[i]<<"\t"<<Afield[i]<<endl;
       output_field<<time[i]<<' '<<Efield[i]<<' '<<Afield[i]<<' '<<Amid_test[i]<<endl;
   //exit(1);*/

/*
   // fft for freq domain
   int nfft=int(nt*fact());   // number of points for fft
   if(nfft%2!=0) nfft+=1; // must be even
   int nout=int(nfft*facout());
   int nmid=int(nfft/2);
   double dw=2.*pi/(double(nfft-1)*dt_p());
   vector<complex> E1_t(nfft), E1_spec(nfft), E2_t(nfft), E2_spec(nfft), E_t(nfft), E_spec(nfft);

   vector<double> freq(nfft);
   for(int i=0; i<=nmid; i++){  // n must be even
      freq[i]=i*dw;
   }
   for(int i=nmid+1; i<=nfft-1; i++){  // n must be even, shifted to negative w
      freq[i]=(i-nfft)*dw;
   }

   for(int i=0; i<nt; i++){
      E1_t[i]=E1[i];  // real part
      E2_t[i]=E2[i];  // real part
      E_t[i]=Efield[i];
   }
   for(int i=nt; i<nfft; i++){  // fill 0
      E1_t[i]=0.;
      E2_t[i]=0.;
      E_t[i]=0.;
   }
*/

   /*fft_complex_1D(nfft, E1_t, E1_spec, dt_p(), -1); //-1 forward, 1 backward
   fft_complex_1D(nfft, E2_t, E2_spec, dt_p(), -1);
   fft_complex_1D(nfft, E_t, E_spec, dt_p(), -1);

   ofstream output_freq;  
   Open_File(output_freq, outputFolder()+"/field_freq.dat");
   //for(int i=nmid; i<nout; i++){ 
   for(int i=nfft-nout; i<nfft; i++){ 
      output_freq<<freq[i]<<" "<<norm(E1_spec[i])<<" "<<norm(E2_spec[i])<<" "<<norm(E_spec[i])<<endl;
   }*/

   //exit(1);

// ----- set initial values for parameters used in propagation
   wavefunction psi; // to store the time-dependent wf
   Initialize( psi, n_points() , spatial_step(), 0, 0 );  
   //Initialize_part( psi, n_points() , spatial_step(), 2, 0);  //  //0: full; 1: half; 2: partial
   PlaceWaveFunction( psi , psi_ground);
   Helium_1plus1D pot(psi, -1., -1., 1., 1., 2., Cr(), soft_en(), soft_ee());
   Hamiltonian H(psi);

   /*ofstream output_x1end;    // output coordinates
   Open_File(output_x1end, outputFolder()+"/x1.dat");
   for(int i=1; i<=psi.n1; i++) { output_x1end<<psi.x1[i]<<endl; }
   ofstream output_x2end;  
   Open_File(output_x2end, outputFolder()+"/x2.dat");
   for(int j=1; j<=psi.n2; j++) { output_x2end<<psi.x2[j]<<endl; }*/
   ofstream output_fs;   
   Open_File(output_fs, outputFolder()+"/field_select.dat");
   double tmid = 0., Emid=0., Amid=0.;

/*
// for partition
   double tfinal=time[nt-1]+tf();  // final time for Volkov propagation
   int nfinal=int(tfinal/dt_p());
   double deltat=facinter()*dt_p();   // time instant for partition
   double v_dt=v0()*deltat;
   vector< vector<complex> > psi_1r, psi_1l;
   int n_SI=int(x_SI()/psi.dx2);
   int np = psi.n1*facp();  // number of fft for momentum space
   if(np%2!=0) np+=1; // must be even
   vector<complex> psi_x1r(np), phi_p1r(np), psi_x1l(np), phi_p1l(np);
   vector<double> p(np), E(np);
   double dp = 2.*pi /(np*psi.dx1);
   int npmid=int(np/2);
   for(int k=0; k<=npmid; k++){  // n must be even
      p[k]=k*dp;
      E[k]=p[k]*p[k]/2.*27.211;  // energy in eV
   }
   for(int k=npmid+1; k<=np-1; k++){  // n must be even, shifted to negative w
      p[k]=(k-np)*dp;
      E[k]=p[k]*p[k]/2.*27.211;  // energy in eV
   }

   int nselect=0;
   for(int k=0; k<np; k++){  // set selected energy regoin
      if(E[k]>Emin() && E[k]<Emax() ){
       nselect++;
      }
   } 
   vector<int> ip(nselect);
   int iselect=0;
   for(int k=0; k<np; k++){  // set selected energy regoin
      if(E[k]>Emin() && E[k]<Emax() ){
       ip[iselect]=k;
       iselect++;
      }
   } 
   vector<complex> phi_p1l_sum(nselect), phi_p1r_sum(nselect), Uvolkov(nselect);
   for(int k=0; k<nselect; k++){  // set 0 for sum
      phi_p1r_sum[k]=(0.,0.);
      //phi_p1l_sum[k]=(0.,0.);
   }
*/

   /*ofstream output_test1;  
   Open_File(output_test1, outputFolder()+"/psi_exchange.dat");
   //ofstream output_test2;  
   //Open_File(output_test2, outputFolder()+"/psi_mom.dat");
   ofstream output_test3;  
   Open_File(output_test3, outputFolder()+"/volkov_phase.dat");*/
   int icount=0;

   gettimeofday(&start_time,NULL);  // start measure time

// ----- start propagation with field
   for( int i = 0 ; i < nt-1 ; i++ )
     {

       tmid=(time[i+1]+time[i])/2.;
       Emid=(Efield[i+1]+Efield[i])/2.;
       Amid=(Afield[i+1]+Afield[i])/2.;

       /*H.X2_Laser_OMP( 0.5*dt_p(), psi, pot, Emid, lengthgauge );   
       H.X1_Laser_OMP( dt_p(), psi, pot, Emid, lengthgauge );     
       H.X2_Laser_OMP( 0.5*dt_p(), psi, pot, Emid, lengthgauge );  */
       H.X2_Laser_OMP( 0.5*dt_p(), psi, pot, Amid, velocitygauge ); 
       H.X1_Laser_OMP( dt_p(), psi, pot, Amid, velocitygauge ); 
       H.X2_Laser_OMP( 0.5*dt_p(), psi, pot, Amid, velocitygauge ); 
     
    /*
       // partition for Volkov, and absorber
       if( i%facinter() == 0){  // every facinter() steps, noted as t_s

         icount++;
         cout<<"---"<<i<<' '<<icount<<endl;
        
         // absorber for inner wf in all regoin, and set up outer wf in SI regoin
         Gauss_partition(psi, a_percent1(),a_percent2(),a_percent1(),a_percent2(), v_dt, n_SI, psi_1r, psi_1l);

         int nx=psi_1r.size();  // number of exchange grids
         int nvolkov=int((tfinal-time[i])/dt_p());  // steps of volkov propagation from t_s to t_f
         
         // clac volkov propagator from t_s to t_f: U(t_s,t_f)*phi_0 = U(t_s,t_s+dt)*...*U(t_f-dt,t_f)*phi_0
         for(int k=0; k<nselect; k++){  // iteration of p for selected regoin
            int ii = ip[k];
            double sum1=0.;
            for(int it=0; it<nvolkov; it++){  // iteration of time from t_s to t_f, with time step = deltat
              double Delta=0.;  // 0 for absence of field
              if(i+it<nt-1) Delta=Afield[i+it]*one_by_lightC_au; // A(t)/c after partition time during the presence of field
              //sum1 += (p[ii]-Delta)*(p[ii]-Delta); // (p-A)^2, gauge
              sum1 += p[ii]*p[ii]-2.*Delta*p[ii]; // p^2-2pA
              //sum1 += p[ii]*p[ii]; // p^2, plane wave
            }
            Uvolkov[k] = exp(-I/2.*dt_p()*sum1);
            //if(icount==1) output_test3<< p[ii] <<' '<< arg(Uvolkov[k]) <<endl;  // output volkov phase at the 2nd partition
         } 

         //start SI
         int isi;
//#pragma omp parallel shared(n_SI, psi_1r, psi_1l, psi, icount, np, nselect, Uvolkov)  // parallel in OpenMP
//{
//    #pragma omp for private(isi, psi_x1r, psi_x1l, phi_p1r, phi_p1l)
         for(isi=0; isi<2*n_SI; isi++){

           for(int k=0; k<nx; k++){
             psi_x1r[k]=psi_1r[k][isi];   // extract 1D wf at right side: need wf in exchanged regoin
             //psi_x1l[k]=psi_1l[k][isi];   // extract 1D wf at left side
             //if(isi==n_SI-1) output_test1<< norm(psi_x1r[k])<<' '<<norm(psi_x1l[k])<<endl;
           }
           for(int k=nx; k<np; k++){
             psi_x1r[k]=0.;  // fill 0 for FFT resolusion
             //psi_x1l[k]=0.;  // fill 0 for FFT resolusion
           }

           fft_complex_1D(np, psi_x1r, phi_p1r, psi.dx1, -1);  // wf in momuntum space, normalization?
           //fft_complex_1D(np, psi_x1l, phi_p1l, psi.dx1, -1);  // wf in momuntum space, normalization?

           for(int k=0; k<nselect; k++){  // iteration of p
             int ii = ip[k];
             phi_p1r_sum[k] += Uvolkov[k]*phi_p1r[ii];  //coherent sum in SI regoin and for each t_s
             //phi_p1l_sum[k] += Uvolkov[k]*phi_p1l[ii];  //coherent sum in SI regoin and for each t_s
           }

         } //end  SI
//} // end parallel

       }  // end partition

       //Mask_Function_OMP(psi, a_percent1(), a_percent2(), a_percent1(), a_percent2(), mask_exp);
       if(i % 20 == 0)  output_fs << tmid <<" "<< Emid <<" "<< Amid << endl;
     */
       if(i % 100 == 0)  cout<<tmid<<"\t"<<"Norm="<< Obs_Norm(psi)<<endl;

    }  // end of time evolution

    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);


/*
// output espec
   ofstream output_espec;  
   Open_File(output_espec, outputFolder()+"/espec.dat");
   for(int k=1; k<=nselect; k++){
      int ii=ip[k]; 
      //output_espec<<p[ii]<<' '<<E[ii]<<' '<<norm(phi_p1r_sum[k])<<' '<<norm(phi_p1l_sum[k])<<endl; 
      output_espec<<p[ii]<<' '<<E[ii]<<' '<<norm(phi_p1r_sum[k])<<endl; 
   }

// ----- save final excited wavefunction
   //psi.save("wf_end");
   ProjectOUT_Diff_Sizes(psi, psi_ground);
   cout <<"Norm(end)-Norm(ground)="<< Obs_Norm(psi)<<endl;
*/
/*
   ofstream output_psiend;  
   Open_File(output_psiend, outputFolder()+"/psi_end.dat");
   for(int i=1; i<=psi.n1; i++)
     for(int j=1; j<=psi.n2; j++){
     output_psiend<<norm(psi.wave[psi.in2(j,i)]) <<endl;
   }
*/

}  // end of program

