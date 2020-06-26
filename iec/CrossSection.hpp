//
//  CrossSection.hpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  20 Jan 2018 - Initial Revision (Jun Hasegawa)
//

#ifndef CROSSSECTION_HPP
#define CROSSSECTION_HPP

#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "Ip1d.h"

const int nc_max = 10;  // maximum number of data file columns

class CrossSection {
    
public:
    // constructor
    CrossSection( const char* fn, int ce = 0, int cc = 1 );
    // fn: file name of cross sction data
    // ce: column number of energy values
    // cc: column number of cross section values
    
    // return an interpolated cross section value [m^-2]
    double sigma( double ek ) const;
    // ek : kinetic energy of particle [eV]
    
    // output all cross secion data [m^-2] to output stream.
    void output( std::ostream& out ) const;
    
    // output interpolated cross secion data [m^-2] to output stream.
    void output( std::ostream& out, int n ) const;
    // n: number of energy grids
    
    // output interpolated cross secion data [m^-2] to output stream.
    void output( std::ostream& out, double emin, double emax, int n ) const;
    // emin : minimum energy of range of output cross section [eV]
    // emax : maximum energy of range of output cross section [eV]
    // n: number of energy grids

private:
    double emin_;
    double emax_;
    std::vector<double> ek_;
    std::vector<double> sigma_;
    std::string title_[2];
    
    std::unique_ptr<Ip1d> s_;   // for interpolation
};

#endif /* CROSSSECTION_HPP */
