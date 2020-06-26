//
//  HydrogenCrossSection.cpp
//  iec
//
//  Created by Jun on 2017/12/01.
//  Copyright © 2017年 Tokyo Tech. All rights reserved.
//

#include "HydrogenCrossSection.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>

/*
HydrogenCrossSection::HydrogenCrossSection( double emin, double emax, int eq, double eth, std::vector<double>& a, std::string proc )
: emin_(emin), emax_(emax), eq_(eq), a_(a), proc_(proc)
{
    ;
}
*/

HydrogenCrossSection::HydrogenCrossSection( const char* fn, int no )
{
    using namespace std;
    ifstream fin( fn, ios::in );
    if ( !fin ) throw runtime_error( "can't open data file in HydrogenCrossSection()" );
    
    char buf[255];
    string dum;
    
    fin >> dum;
    if ( dum != "No." ) throw runtime_error( "wrong data file format" );
    fin.getline( buf, 255 );	// skip first line
    
    for (;;) {
        fin >> no_;
        if ( no_ == no ) {	// load data
            fin >> proc_ >> emin_ >> emax_ >> delta_rms_ >> delta_max_ >> e_delta_max_ >> eq_ >> n_;
            a_.resize(n_+1);
            for ( int i = 0; i <= n_; i++ ) {
                fin >> a_.at(i);
            }
            break;
        }
        fin.getline( buf, 255 );    // skip rest of line
    }
    //cout << "-- cross section data loaded:" << proc_ << endl;
}

double HydrogenCrossSection::sigma( double ek ) const
{
    ek *= 1e-3; // eV -> keV
    if ( ek < emin_ || ek > emax_ ) {
        //throw std::out_of_range( "HydrogenCrossSection::sigma()" );
        return 0;
    }
    return sig( ek ) * 1e-4;  // cm^2 -> m^2
}

double HydrogenCrossSection::sig( double ek ) const
{
    double e1 = ek-a_[0];
	switch( eq_ ) {
        case 1:
            return s1( e1 );
            break;
        case 2:
            return s2( e1 );
            break;
        case 3:
            return s3( e1 );
            break;
        case 4:
            return s4( e1 );
            break;
        case 5:
            return s5( e1 );
            break;
        case 6:
            return s6( e1 );
            break;
        case 7:
            return s7( e1 );
            break;
        case 8:
            return s8( e1 );
            break;
        case 9:
            return s9( e1 );
            break;
        case 10:
            return s10( e1 );
            break;
        case 11:
            return s11( e1 );
            break;
        case 12:
            return s12( e1 );
            break;
        case 13:
            return s13( e1 );
            break;
        case 14:
            return s14( e1 );
            break;
        default:
            throw std::runtime_error( "illegal equation number in HydrogenCrossSection::sig()" );
            break;
    }
}

void HydrogenCrossSection::output( std::ostream& out, int n ) const
{
    double de = (log(emax_)-log(emin_))/n;  // keV
    
    double ek;
    for ( int i = 0; i < n; i++ ) {
        ek = (exp(log(emin_)+i*de));
        out << ek*1e3 << '\t' << sig( ek ) << '\n';
    }
    ek = emax_;
    out << ek*1e3 << '\t' << sig( ek );
}
