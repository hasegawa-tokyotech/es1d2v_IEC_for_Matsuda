#ifndef PHYS_CONST_H
#define PHYS_CONST_H

namespace si{
	const double c   = 2.9979246e8;		// light speed [m/s]
	const double e   = 1.6021773e-19;	// charge of electron [C]
	const double mu0 = 1.2566371e-6;	// magnetic permeability for vacuum [H/m]
	const double ep0 = 8.8541878e-12;	// dielectric constant for vacuum [F/m]
	const double h   = 6.6260755e-34;	// Planck constant [Js]
	const double h2  = 1.05457266e-34;	// h/2pi [Js]
	const double me  = 9.1093897e-31;	// mass of electron [kg]
	const double mp  = 1.6726231e-27;	// mass of proton [kg]
	const double mn  = 1.6749286e-27;	// mass of neutron [kg] 
	const double amu = 1.6605402e-27;	// atomic mass unit [kg]
	const double na  = 6.0221367e23;	// Avogadro number [1/mol]
	const double kb  = 1.380658e-23;	// Boltzmann constant [J/K]
	const double sg  = 5.67051e-8;		// Stefan-Boltzmann coefficient [W/m^2/K^4] 
	const double a0  = 5.2917725e-11;	// Bohr radius [m] 
	const double vb  = 2.1876914e6;		// Bohr velocity [m/s] 
	const double pi  = 3.1415926535;	// circumference
	const double af  = 7.29735308e-3;	// fine structure constant
	
	const double km  = 1e3;				// kilometer
	const double m   = 1;				// meter
	const double cm  = 0.01;			// centimeter
	const double mm  = 1e-3;			// millimeter
	const double um  = 1e-6;			// micrometer
	const double nm  = 1e-9;			// nanometer
	const double fm  = 1e-12;			// femtometer
	const double ang = 1e-10;			// angstrom
	const double kg  = 1;				// kilogram
	const double g   = 1e-3;			// gram
	const double mg  = 1e-6;			// miligram
	const double ug  = 1e-9;			// microgram
	const double ng  = 1e-12;			// nanogram
	const double sec = 1;				// second
	const double ms  = 1e-3;			// millisecond
	const double us  = 1e-6;			// microsecond
	const double ns  = 1e-9;			// nanosecond
	const double fs  = 1e-12;			// femtosecond
	const double rad = 1;				// radian
	const double deg = pi/180;			// degree
	const double eV  = e;				// electron volt [J/eV]
	const double keV = e*1e3;
	const double MeV = e*1e6;
	const double GeV = e*1e9;
}

namespace cgs{
	const double c   = 2.9979246e10;	// light speed [cm/s]
	const double e   = 4.80325e-10;		// charge of electron [esu]
	const double mu0 = 1;				// magnetic permeability for vacuum [H/m]
	const double ep0 = 1;				// dielectric constant for vacuum [F/m]
	const double h   = 6.6260755e-27;	// Planck constant [erg s]
	const double h2  = 1.05457266e-27;	// h/2pi [erg s]
	const double me  = 9.1093897e-28;	// mass of electron [g]
	const double mp  = 1.6726231e-24;	// mass of proton [g]
	const double mn  = 1.6749286e-24;	// mass of neutron [g] 
	const double amu = 1.6605402e-24;	// atomic mass unit [g]
	const double na  = 6.0221367e23;	// Avogadro number [1/mol]
	const double kb  = 1.380658e-16;	// Boltzmann constant [erg/K]
	const double sg  = 5.67051e-5;		// Stefan-Boltzmann coefficient [erg/s/cm^2/K^4] 
	const double a0  = 5.2917725e-9;	// Bohr radius [cm] 
	const double vb  = 2.1876914e8;		// Bohr velocity [cm/s] 
	const double pi  = 3.1415926535;	// circumference
	const double af  = 7.29735308e-3;	// fine structure constant

	const double km  = 1e5;				// kilometer
	const double m   = 1e2;				// meter
	const double cm  = 1;				// centimeter
	const double mm  = 1e-1;			// millimeter
	const double um  = 1e-4;			// micrometer
	const double nm  = 1e-7;			// nanometer
	const double fm  = 1e-10;			// femtometer
	const double ang = 1e-8;			// angstrom
	const double kg  = 1e3;				// kilogram
	const double g   = 1;				// gram
	const double mg  = 1e-3;			// miligram
	const double ug  = 1e-6;			// microgram
	const double ng  = 1e-9;			// nanogram
	const double hr  = 3600;			// hour
	const double min = 60;				// minite
	const double s	 = 1;				// second
	const double ms  = 1e-3;			// millisecond
	const double us  = 1e-6;			// microsecond
	const double ns  = 1e-9;			// nanosecond
	const double fs  = 1e-12;			// femtosecond
	const double rad = 1;				// radian
	const double deg = pi/180;			// degree
	const double eV  = e*1e7;			// electron volt [erg/eV]
	const double keV = e*1e10;
	const double MeV = e*1e13;
	const double GeV = e*1e16;
}

#endif

