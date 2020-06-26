//
//  Field1D.cpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  19 Jan 2018 - Initial Revision (Jun Hasegawa)
//	08 Aug 2019 - Updated (Jun Hasegawa)
//

#include "Field1D.hpp"

#include <stdexcept>
#include <string>
#include <fstream>

Field1D::Field1D( const char* fn, double vc0, double vc )
{
    std::ifstream fin( fn, std::ios::in );
    if ( !fin ) throw std::runtime_error( "can't open file in Field1D()" );
    
    std::string dum;
    char buf[255];
    
    // check file format
    fin >> dum;
    if ( dum != "%" ) {
        throw std::runtime_error( "illegal file format in Field1D()" );
    }
    else {
        fin >> dum;
        if ( dum != "Model:" ) {
            throw std::runtime_error( "illegal file format in Field1D()" );
        }
        else {
            for ( int i = 0; i < 9; i++ )
                fin.getline( buf, 255 );    // skip header
        }
    }
    
    double f = vc/vc0;    // correction factor
    
    // load data
    for (;;) {
        double x, phi, ex;
        fin >> buf >> x >> buf >> phi >> buf >> ex >> buf;
        if ( fin.eof() ) break;
        
        x_.push_back( x );
        phi_.push_back( f*phi );
        ex_.push_back( f*ex );
    }
    
    // prepare a member for interpolation
    p_ = std::make_unique<Linear>(x_,phi_);
    e_ = std::make_unique<Linear>(x_, ex_);
    
}

Field1D::Field1D( const char* fn )
{
    std::ifstream fin( fn, std::ios::in );
    if ( !fin ) throw std::runtime_error( "can't open file in Field1D()" );
    
    std::string dum;
    char buf[255];
    
    fin.getline( buf, 255 );    // skip first line
    
    // load data
    for (;;) {
        double x, phi, ex;
        fin >> x >> phi >> ex;
        if ( fin.eof() ) break;
        
        x_.push_back( x );
        phi_.push_back( phi );
        ex_.push_back( ex );
    }
    
    // prepare a member for interpolation
    p_ = std::make_unique<Linear>(x_, phi_);
    e_ = std::make_unique<Linear>(x_, ex_);
    
}

Field1D::Field1D( double lx0, double lx1, double vc, int nc )
{
	double lx = lx0 + lx1;
	double dx = 2*lx/nc;
	
	for ( int i = 0; i <= nc; i++ ) {
		double x = -lx + i*dx;
		double phi;
		if ( -lx <= x && x < -lx0 ) phi = vc*(x+lx)/lx1;
		else if ( -lx0 <= x && x <= lx0 ) phi = vc;
		else phi = vc*(lx-x)/lx1;

		x_.push_back( x );
		phi_.push_back( phi );
		
		//std::cout << x << '\t' << phi << '\n';	// for debug
	}
	
	double ex = vc/lx1;
	ex_.push_back( ex );	// left boundary
	for ( int i = 1; i < nc; i++ ) {
		ex = -(phi_[i+1]-phi_[i-1])/(2*dx);
		ex_.push_back( ex );
	}
	ex_.push_back( ex );	// right boundary

	// prepare a member for interpolation
	p_ = std::make_unique<Linear>(x_, phi_);
	e_ = std::make_unique<Linear>(x_, ex_);
}

void Field1D::output( std::ostream& out ) const
{
    out << "x[m]\tphi[V]\tex[V/m]\n";
    auto n = x_.size();
    for ( long i = 0; i < n; ++i )
        out << x_.at(i) << '\t' << phi_.at(i) << '\t' << ex_.at(i) << '\n';
}

void Field1D::output( std::ostream& out, long n ) const
{
    out << "x[m]\tphi[V]\tex[V/m]\n";
    
    double x0 = x_.front();
    double x1 = x_.back();
    double dx = (x1-x0)/n;
    
    for ( int i = 0; i < n; ++i ) {
        double x = x0 + i*dx;
        out << x << '\t' << phi(x) << '\t' << ex(x) << '\n';
    }

}


double Field1D::phi( double x ) const
{
	try {
		return p_->value( x );
	}
	catch( std::out_of_range ) {
		return 0;
	}
}

double Field1D::ex( double x ) const
{
	try {
		return e_->value( x );
	}
	catch( std::out_of_range ) {
		return 0;
	}
}
