// $Id: Laser.h 536 2007-05-02 14:30:33Z arvid $
#ifndef LASER_H
#define LASER_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>      // for Field::Field(filename)
#include "utils.h"      // for skip_comments(std::istream & is)
#include "constant.h"   // for value of pi

using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::endl;


namespace Laser_fields {


// void Two_pulses ( vector<double> wavelength, vector<double> intensity, vector<double> n_cycle, double dt, vector<double> cep, double delay,  double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period1);
 void Two_pulses ( vector<double> wavelength, vector<double> intensity, vector<double> n_cycle, double dt, vector<double> cep, double delay, double t_wait, vector<double> &time, vector<double> &Efield1, vector<double> &Efield2, vector<double> &Efield, vector<double> &Afield, double &period1, double &tstart2);
// void Two_pulses ( vector<double> wavelength, vector<double> intensity, vector<double> n_cycle, double dt, vector<double> cep, double delay, double t_wait, vector<double> &time, vector<double> &Efield1, vector<double> &Efield2, vector<double> &Efield, vector<double> &Afield, double &period1);

 void IR_APT ( double lambda1, double i1, double n_c1, double cep1_pi, int nenv, double i2, double n_c2, int nhhg,vector<double> ohhg, vector<double> ahhg, vector<double> cephhg_pi, double delay, double dt, double t_wait, vector<double> &time, vector<double> &E1, vector<double> &E2, vector<double> &Efield, vector<double> &Afield, int &iend_pulses);

 void One_pulse ( double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period);

 void SINn_SINm_pulse ( double n_front, double n_back, double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period);

 void Static_plus_1pulse ( double static_E, double wavelength, double intensity, double n_cycle, vector<double> dt, double cep, double t_wait1,  double t_wait2, vector<double> &time, vector<double> &Efield, vector<double> &Afield, double &period);

}



#endif	/* LASER_H */
