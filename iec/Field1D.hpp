//
//  Field1D.hpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  19 Jan 2018 - Initial Revision (Jun Hasegawa)
//

#ifndef IEC_FIELD_HPP
#define IEC_FIELD_HPP

#include "Ip1d.h"

#include <iostream>
#include <vector>
#include <memory>

//#define IEC_USE_COMSOL

class Field1D {
    
public:

//#ifdef IEC_USE_COMSOL   // use COMSOL field data
    
    // constructor for COMSOL field
    Field1D( const char* fn, double phi0, double phi );
    // fn: file name of COMSOL field data
    // phi0: cathode voltage in COMSOL [V]
    // phi: actual cathode voltage [V]
    
    // constructor for PIC field
    Field1D( const char* fn );
    // fn: file name of field data created by PIC code

	Field1D( double lx0, double lx1, double vc, int nc );
	// lx0: cathode half length [m]
	// lx1: gap length [m]
	// vc: cathode voltage [V]
	// nc:  number of cells

	double x( int i ) const { return x_.at(i); }
    double phi( int i ) const { return phi_.at(i); }
    double ex( int i ) const { return ex_.at(i); }
    double size() const { return x_.size(); }
    
    double xmin() const { return x_.front(); }
	// return minimum position in x [m]
    double xmax() const { return x_.back(); }
	// return maximum position in x [m]

	virtual double phi( double x ) const;
	virtual double ex( double x ) const;
	
//#else // use field in vacuum
    

    
    //double lx0() const { return lx0_; }
    //double lx1() const { return lx1_; }

//#endif /* IEC_USE_COMSOL */
    


    void output( std::ostream& out ) const;
    void output( std::ostream& out, long n ) const;

private:

//#ifdef IEC_USE_COMSOL

    std::vector<double> x_;
    std::vector<double> phi_;    // electrostatic potential along center axis [V]
    std::vector<double> ex_;    // electric field along center axis [V/m]
    
    std::unique_ptr<Ip1d> p_;   // phi
    std::unique_ptr<Ip1d> e_;   // ex

//#else

    //double lx0_; // cathode half length [m]
    //double lx1_; // gap length [m]
    //double ex_;	// constant electric field [V/m]

//#endif /* IEC_USE_COMSOL */

};

#endif /* IEC_FIELD_HPP */
