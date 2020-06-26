//
//  HydrogenCrossSection.hpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  20 Jan 2018 - Initial Revision (Jun Hasegawa)
//
//	References
//	[1] T. Tabata and T. Shirai, "Analytic Cross Sections for Collisions of H+, H2+, H3+, H, H2, and H-
//  with Hydrogen Molecules", Atomic Data and Nuclear Data Tables 76, 1-25 (2000).
//

#ifndef HYDROGENCROSSSECTION_HPP
#define HYDROGENCROSSSECTION_HPP

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

const double s0 = 1e-16;	// sigma_0 = 1e-16 [cm^2]
const double er = 1.361e-2;	// Rydberg constant [keV]

class HydrogenCrossSection {
public:
    // constructor
    //HydrogenCrossSection( double emin, double emax, int eq, double eth, std::vector<double>& a, std::string proc = "" );
    // emin: minimum energy of the recommended data [keV]
    // emax: maximum energy of the recommended data [keV]
    // eq: the identifying number of the equation to be used
    // eth: Threshold energy of the reaction [keV]
    // a: coefficients, a_[1...n] = a1, a2, ..., an
    
    // constructor
    HydrogenCrossSection( const char* fn, int no );
    // fn: data file name
	// no: the identifying number of a reaction process in Table A in ref.[1]
    
    // return an interpolated cross section value [m^-2]
    double sigma( double ek ) const;
    // ek : kinetic energy of projectile [eV]
    
    // output interpolated cross secion data [m^-2] to an output stream.
    void output( std::ostream& out, int n ) const;
    // n: number of energy grids
    
    // return process name
    std::string proc( void ) { return proc_; }
    
private:
    double sig( double ) const;
    
    double f1( double x, double c1, double c2 ) const { return s0*c1*pow(x/er,c2); }
    double f2( double x, double c1, double c2, double c3, double c4 ) const { return f1(x,c1,c2)/(1+pow(x/c3,c2+c4)); }
    double f3( double x, double c1, double c2, double c3, double c4, double c5, double c6 ) const { return f1(x,c1,c2)/(1+pow(x/c3,c2+c4)+pow(x/c5,c2+c6)); }
    double f4( double x, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8 ) const { return f1(x,c1,c2)*(1+pow(x/c3,c4-c2))/(1+pow(x/c5,c4+c6)+pow(x/c7,c4+c8)); }
    
    double s1( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]); }
    double s2( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + a_[5]*f2(e1/a_[6],a_[1],a_[2],a_[3],a_[4]); }
    double s3( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f2(e1,a_[5],a_[6],a_[7],a_[8]); }
    double s4( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f2(e1,a_[5],a_[6],a_[7],a_[8]) + a_[9]*f2(e1/a_[10],a_[5],a_[6],a_[7],a_[8]); }
    double s5( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f2(e1,a_[5],a_[6],a_[7],a_[8]) + f2(e1,a_[9],a_[10],a_[11],a_[12]); }
    double s6( double e1 ) const { return f3(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]); }
    double s7( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f3(e1,a_[5],a_[2],a_[6],a_[7],a_[8],a_[4]); }
    double s8( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f3(e1,a_[5],a_[2],a_[6],a_[7],a_[8],a_[9]); }
    double s9( double e1 ) const { return f2(e1,a_[1],a_[2],a_[3],a_[4]) + f3(e1,a_[5],a_[6],a_[7],a_[8],a_[9],a_[4]); }
    double s10( double e1 ) const { return f3(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]) + a_[7]*f3(e1/a_[8],a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]); }
    double s11( double e1 ) const { return f3(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]) + f2(e1,a_[7],a_[8],a_[9],a_[10]); }
    double s12( double e1 ) const { return f3(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]) + f2(e1,a_[7],a_[8],a_[9],a_[10]) + a_[11]*f2(e1/a_[12],a_[7],a_[8],a_[9],a_[10]); }
    double s13( double e1 ) const { return f3(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6]) + f3(e1,a_[7],a_[8],a_[9],a_[10],a_[11],a_[12]); }
    double s14( double e1 ) const { return f4(e1,a_[1],a_[2],a_[3],a_[4],a_[5],a_[6],a_[7],a_[8]); }

    int no_;	// process number
    std::string proc_;	// process
    double emin_;	// minimum energy [keV]
    double emax_;	// maximum energy [keV]
    double delta_rms_;
    double delta_max_;
    double e_delta_max_;
    int eq_;	// equation number
    int n_;		// number of coefficients
    std::vector<double> a_;	// coefficients: a[0] = Eth, a[1] = a1, a[2] = a2, ...
};

#endif /* HYDROGENCROSSSECTION_HPP */
