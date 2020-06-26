//
//  BalmerAlphaCrossSection.cpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  20 Jan 2018 - Initial Revision (Jun Hasegawa)
//

#include "BalmerAlphaCrossSection.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdexcept>

BalmerAlphaCrossSection::BalmerAlphaCrossSection( const char* fn )
{
    using namespace std;
    
    ifstream fin( fn, ios::in );
    if ( !fin ) throw runtime_error( "can't open file in CrossSection()" );
    
    string dum;
    char buf[255];
    
    // check file format
    fin >> dum;
    if ( dum != "Energy(eV)" )
        throw runtime_error( "illegal file format in CrossSection()" );
    
    // skip the rest of the first line
    fin.getline( buf, 255 );
    
    int count = 0;
    for (;;) {
        double ek, sigma_t, sigma3s, sigma3d;
        fin >> ek >> sigma_t >> sigma3s >> sigma3d;
        if ( fin.eof() ) break;
        
        ek_.push_back(ek);
        sigma3s_.push_back(sigma3s);
        sigma3d_.push_back(sigma3d);
        count++;
    }
    
    cout << fn << ": ";
    cout << count << " cross section data loaded." << endl;

    // prepare a member for interpolation
    s3s_ = make_unique<Spline>( ek_, sigma3s_ );
    s3d_ = make_unique<Spline>( ek_, sigma3d_ );

}

