//
//  CrossSection.cpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  20 Jan 2018 - Initial Revision (Jun Hasegawa)
//

#include "CrossSection.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

CrossSection::CrossSection( const char* fn, int ce, int cc )
{
    using namespace std;
    
    ifstream fin( fn, ios::in );
    if ( !fin ) throw runtime_error( "can't open file in CrossSection()" );
    
    char buf[255];
    fin.getline( buf, 255 );    // load first line
    stringstream strin( buf );
    
    // count number of columns
    int nc = 0;
    while ( !strin.eof() ) {
        char c;
        strin.get(c);
        if ( c == '\t' ) nc++;
    }
    nc++;
    
    if ( nc > nc_max ) throw runtime_error( "too many columns in CrossSection()" );
    
    // load title
    string title;
    stringstream strin2( buf );
    for ( int i = 0; i < nc; ++i ) {
        strin >> title;
        if ( i == ce ) title_[0] = title;
        else if ( i == cc ) title_[1] = title;
    }
    
    // load data
    int count = 0;
    for (;;) {
        double val;
        double ek, sigma;
        for ( int i = 0; i < nc; ++i ) {
            fin >> val;
            if ( i == ce ) ek = val;
            else if ( i == cc ) sigma = val;
        }
        if ( fin.eof() ) break; // end of file?
        
        ek_.push_back(ek);
        sigma_.push_back(sigma);
        
        count++;
    }
    
    cout << fn << ": ";
    cout << count << " cross section data loaded." << endl;
    
    emin_ = ek_.at(0);
    emax_ = ek_.at(ek_.size()-1);

    // prepare a member for interpolation
    s_ = make_unique<Spline>( ek_, sigma_ );

}

double CrossSection::sigma( double ek ) const
{
    if ( ek < emin_ || ek > emax_ ) return 0;
    return s_->value(ek);
}


void CrossSection::output( std::ostream& out ) const
{
    out << title_[0] << '\t' << title_[1] << '\n';
    auto n = ek_.size();
    for ( int i = 0; i < n; ++i ) {
        out << ek_.at(i) << '\t' << sigma_.at(i) << '\n';
    }
}

void CrossSection::output( std::ostream& out, int n ) const
{
    out << title_[0] << '\t' << title_[1] << '\n';

    double de = (log(emax_) - log(emin_))/(n-1);
 
    for ( int i = 0; i < n; ++i ) {
        double ek = emin_*exp(i*de);
        out << ek << '\t' << s_->value( ek ) << '\n';
    }
}

void CrossSection::output( std::ostream& out, double emin, double emax, int n ) const
{
    out << title_[0] << '\t' << title_[1] << '\n';
    
    double de = (log(emax) - log(emin))/(n-1);
    
    for ( int i = 0; i < n; ++i ) {
        double ek = emin*exp(i*de);
        out << ek << '\t' << s_->value( ek ) << '\n';
    }
}
