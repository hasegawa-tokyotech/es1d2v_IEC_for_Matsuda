//
//  main.cpp
//  create_input_txt
//
//  Created by Jun on 2019/07/19.
//  Copyright © 2019 Jun. All rights reserved.
//
//	This code was developed to create an input file for the 1D-PIC simulation
//	code (es1des) specialized for IEC plasma potential analysis.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "phys_const.h"

int main()
{
	using namespace std;
	
	double n0;
	cout << "number density of initial plasma (m^-3) = ";
	cin >> n0;
	
	double kte = 1.0;
	//cout << "electron temperature (eV) = ";
	//cin >> kte;
	
	double mi = 4.0;	// D2 (deuterium molecule)
	//cout << "ion mass (amu) = ";
	//cin >> mi;
	
	int qi = 1;
	//cout << "ion charge = ";
	//cin >> qi;
	
	double vd = 0.0;
	//cout << "drifting energy (eV) = ";
	//cin >> vd;
	
	double lx = 0.2;
	//cout << "gap distance (cm) = ";
	//cin >> lx;
	//lx *= 1e-2;
	
	double x0 = 0.0;
	//cout << "left boundary position of initial plasma (cm) = ";
	//cin >> x0;
	//x0 *= 1e-2;
	
	double x1 = 0.0;
	//cout << "right boundary position of initial plasma (cm) = ";
	//cin >> x1;
	//x1 *= 1e-2;
	
	int nx = 1000;
	//cout << "number of cells = ";
	//cin >> nx;
	
	double phi;
	cout << "discharge voltage (-kV) = ";
	cin >> phi;
	phi *= -1e3;    // -kV -> V

	double cur;
	cout << "discharge current (mA) = ";
	cin >> cur;
	cur *= 1e-3;    // mA -> A
	
	double rp;	// plasma radius (m)
    cout << "plasma radius (cm) = ";
    cin >> rp;
    rp *= 1e-2; // cm -> m
	
	// normalization factors
	double wp = sqrt(n0/si::me/si::ep0)*si::e;	// plasma frequency (rad/s)
	double x_norm = si::c/wp;	// distance
	double t_norm = 1.0/wp;	// time
	
	double dx_ = lx/nx/x_norm;	// cell size (m)
	int nc_ = 1000;	// number of particles per cell
	double n0_ = nc_/dx_;	// number density of super particles

	double phi_norm = si::me*si::c*si::c/si::e;
	double ef_norm = si::me*si::c*wp/si::e;
	double bf_norm = si::me*wp/si::e;
	double n_norm = n0/n0_;
	double num_norm = n_norm*x_norm;
	double j_norm = n_norm*si::e*si::c;
	//double rho_norm = n_norm*si::e;
    double ene_norm = si::me*si::c*si::c;
	double ene_norm_ev = ene_norm/si::e;
	
	// output normalized parameters to the input file
	double me_ = 1.0;	// electron mass
	int qe_ = -1;	// electron charge
	double x0_ = x0/x_norm;	// left boundary position (initial plasma)
	double x1_ = x1/x_norm;	// right boundary position (initial plasma)
	double dt_ = 0.1;	// time step
	long np_ = n0_*(x1_-x0_);	// initial number of plasma particles
	double kte_ = kte/ene_norm_ev;	//	electron temperature
	double vde_ = sqrt(2*vd*si::e/si::me);	// drift velocity
	double jd = cur/si::pi/(rp*rp);	// discharge current density
	long np0_ = jd/j_norm*dt_;	// number of particles generated per step
	
	cout << "number of particles generated per step = " << np0_ << endl;
	if ( np0_ < 1 ) cout << "too small np0, give larger initial density";
	
	ofstream fout( "input.txt", ios::out );
	fout << "#pic_1des\n\n";

	fout << scientific;
	fout << "$particles\n";	// electron
	fout << "np=\t" << np_ << '\n'
		 << "m=\t" << me_ << '\n'
		 << "q=\t" << qe_ << '\n'
		 << "x0=\t" << x0_ << '\n'
		 << "x1=\t" << x1_ << '\n'
		 << "kt=\t" << kte_ << '\n'
		 << "vd=\t" << vde_ << '\n'
		 << "np0=\t" << np0_ << '\n'
		 << "kt0=\t" << kte_ << '\n'
		 << "vd0=\t" << vde_ << '\n';
	fout << '\n';

	double mi_ = 2*si::amu/si::me;	// electron mass
	int qi_ = 1;	// electron charge
	double kti_ = 300*si::kb/ene_norm;	//	room temperature
	double vdi_ = vd;
	
	fout << "$particles\n";	// ion
	fout << "np=\t" << np_ << '\n'
		 << "m=\t" << mi_ << '\n'
		 << "q=\t" << qi_ << '\n'
		 << "x0=\t" << x0_ << '\n'
		 << "x1=\t" << x1_ << '\n'
		 << "kt=\t" << kti_ << '\n'
		 << "vd=\t" << vdi_ << '\n'
		 << "np0=\t" << np0_ << '\n'
		 << "kt0=\t" << kti_ << '\n'
		 << "vd0=\t" << vdi_ << '\n';
	fout << '\n';

	fout << "$field\n";	// ion
	fout << "nx=\t" << nx << '\n'
		 << "lx=\t" << lx/x_norm << '\n'
		 << "ns=\t" << 2 << '\n'
		 << "n0=\t" << n0_ << '\n'
		 << "phi_n=\t" << phi/phi_norm << '\n'
		 << "bz=\t" << 0.0 << '\n'
		 << "bc0=\t" << 0 << '\n'
         << "bc1=\t" << 0 << '\n';
	fout << '\n';

	fout << "$world\n";	// ion
	fout << "tend=\t" << 10000 << '\n'
		 << "tout=\t" << 100 << '\n'
		 << "dt=\t" << dt_ << '\n';
	fout << '\n';
	
	fout.close();
	fout.open( "norm_params.txt", ios::out );
	fout << "x_norm=\t" << x_norm << '\n'
		 << "t_norm=\t" << t_norm << '\n'
		 << "phi_norm=\t" << phi_norm << '\n'
		 << "wp=\t" << wp << '\n'
         << "n_norm=\t" << n_norm << '\n'
         << "num_norm=\t" << num_norm << '\n'
         << "ef_norm=\t" << ef_norm << '\n'
		 << "bf_norm=\t" << bf_norm << '\n'
         << "ene_norm=\t" << ene_norm << '\n'
         << "ene_norm_ev=\t" << ene_norm_ev << '\n'
         << "j_norm=\t" << j_norm << '\n';
	fout.close();
	
	return 0;
}
