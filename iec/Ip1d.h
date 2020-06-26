/*
 *  Ip1d.h
 *
 *	ver. 181210
 *
 *  One-dimensional interpolation class developed based on
 *  C functions given in Nurmerical Recipes in C, Section 3.3.
 *
 *  20 Jan 2018 - Minor revision with new features of C++11 (Jun Hasegawa)
 *  10 Dec 2018 - added virtual destructor  (Jun Hasegawa)
 *
 */

#ifndef IP1D_H
#define IP1D_H

#include <vector>

class Ip1d {
public:
	Ip1d( const std::vector<double>& x, const std::vector<double>& y )
	: x_( x ), y_( y ) {}
    virtual ~Ip1d() {}
	virtual double value( double x ) const = 0;
    virtual double slope( double x ) const = 0;
protected:
	const std::vector<double> x_;
	const std::vector<double> y_;
	void locate( const std::vector<double>& xx, double x, long& j ) const;
};

class Linear : public Ip1d {    // linear interpolation
public:
	Linear( const std::vector<double>& x, const std::vector<double>& y ) : Ip1d( x, y ) {}
	virtual double value( double x ) const; // return interpolated value at given x
    virtual double slope( double x ) const; // return slope value at given x 
};

class Spline : public Ip1d {    // Cubic spline interpolation
public:
	Spline( const std::vector<double>& x, const std::vector<double>& y );
	virtual double value( double x ) const; // return interpolated value at given x
    virtual double slope( double x ) const; // return slope value at given x
protected:
	std::vector<double> y2_;

	void splintv( const std::vector<double>& xa, const std::vector<double>& ya,
                 const std::vector<double>& y2a, double x, double& y ) const;
	void splints( const std::vector<double>& xa, const std::vector<double>& ya,
                 const std::vector<double>& y2a, double x, double& dydx ) const;
	void spline( const std::vector<double>& x, const std::vector<double>& y,
                double yp1, double ypn, std::vector<double>& y2 ) const;
};

#endif
