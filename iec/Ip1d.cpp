/*
 *  Ip1d.cpp
 *
 *	ver. 181210
 *
 */

#include "Ip1d.h"

#include <stdexcept>
using namespace std;

void Ip1d::locate( const vector<double>& xx, double x, long& j ) const
{
	long ju,jm,jl;
	bool ascnd;

	long n=xx.size();
	jl = -1;
	ju = n;
	ascnd = (xx[n-1] >= xx[0]);
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl = jm;
		else
			ju = jm;
	}
	if (x == xx[0]) j = 0;
	else if (x == xx[n-1]) j = n-2;
	else j=jl;
}

//------------------ linear ---------------------

double Linear::value( double x ) const
{
	long j;
	locate( x_, x, j );

	if ( j == -1 || j == x_.size()-1 )
		throw out_of_range( "Linear::value()" );

	// linear interpolation
	double m = (x-x_[j])/(x_[j+1]-x_[j]);
	
	return (1-m)*y_[j] + m*y_[j+1];
}

double Linear::slope( double x ) const
{
	long j;
	locate( x_, x, j );
    
	if ( j == -1 || j == x_.size()-1 )
		throw out_of_range( "Linear::slope()" );
    
	// calculate slope
	double a = (y_[j+1]-y_[j])/(x_[j+1]-x_[j]);
	
	return a;
}

//------------------- spline ----------------------

Spline::Spline( const vector<double>& x, const vector<double>& y )
: Ip1d( x, y )
{
	y2_.resize( y.size() );
	spline( x, y, 1e30, 1e30, y2_ );
}

double Spline::value( double x ) const
{
	double y;
	splintv( x_, y_, y2_, x, y );
	return y;
}

double Spline::slope( double x ) const
{
	double dydx;
	splints( x_, y_, y2_, x, dydx );
	return dydx;
}

void Spline::spline( const vector<double>& x, const vector<double>& y, double yp1, double ypn, vector<double>& y2 ) const
{
	long n = y2.size();
	vector<double> u;
	u.resize( y2.size()-1 );
	
	if ( yp1 > 0.99e30 )
		y2[0] = u[0] = 0.0;
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for ( int i = 1; i < n-1; ++i ) {
		double sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		double p = sig*y2[i-1] + 2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
	}
	double qn, un;
	if ( ypn > 0.99e30 )
		qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for ( long k = n-2; k >= 0; --k )
		y2[k] = y2[k]*y2[k+1]+u[k];
}

void Spline::splintv( const vector<double>& xa, const vector<double>& ya, const vector<double>& y2a, double x, double& y ) const
{
	long n = xa.size();
	long klo = 0;
	long khi = n-1;
	while ( khi-klo > 1 ) {
		long k = (khi+klo) >> 1;
		if ( xa[k] > x ) khi = k;
		else klo = k;
	}
	double h = xa[khi] - xa[klo];
	if ( h == 0.0 ) runtime_error( "Bad xa input to routine splint" );
	double a = (xa[khi]-x)/h;
	double b = (x-xa[klo])/h;
	y = a*ya[klo]+b*ya[khi] + ((a*a*a-a)*y2a[klo]
		+ (b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void Spline::splints( const vector<double>& xa, const vector<double>& ya, const vector<double>& y2a, double x, double& dydx ) const
{
	long n = xa.size();
	long klo = 0;
	long khi = n-1;
	while ( khi-klo > 1 ) {
		long k = (khi+klo) >> 1;
		if ( xa[k] > x ) khi = k;
		else klo = k;
	}
	double h = xa[khi] - xa[klo];
	if ( h == 0.0 ) runtime_error( "Bad xa input to routine splint" );
	double a = (xa[khi]-x)/h;
	double b = (x-xa[klo])/h;
    
	dydx = (ya[khi]-ya[klo])/h - (3*a*a-1)*h*y2a[klo]/6.0 + (3*b*b-1)*h*y2a[khi]/6.0;   // cf. NRC, p.105, Eq.(3.3.5)
}

