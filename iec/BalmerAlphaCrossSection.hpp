//
//  BalmerAlphaCrossSection.hpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  20 Jan 2018 - Initial Revision (Jun Hasegawa)
//

#ifndef BALMERALPHACROSSSECTION_HPP
#define BALMERALPHACROSSSECTION_HPP

#include <vector>
#include <memory>

#include "Ip1d.h"

class BalmerAlphaCrossSection {
public:
    BalmerAlphaCrossSection( const char* fn );
    
    double sigma3s( double ek ) const { return s3s_->value(ek); }
    double sigma3d( double ek ) const { return s3s_->value(ek); }
	// These functions return an interpolated cross section value [m^-2]
    // ek : kinetic energy of particle [eV]
    
private:
    std::vector<double> ek_;
    std::vector<double> sigma3s_;
    std::vector<double> sigma3d_;
    
    std::unique_ptr<Ip1d> s3s_;
    std::unique_ptr<Ip1d> s3d_;
};

#endif /* BALMERALPHACROSSSECTION_HPP */
