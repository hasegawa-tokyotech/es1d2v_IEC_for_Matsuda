//
//  Particle1D.hpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  19 Jan 2018 - Initial Revision (Jun Hasegawa)
//	12 Sep 2019 - Integration with boost odeint library was introduced
//

#ifndef IEC_PARTICLE1D_HPP
#define IEC_PARTICLE1D_HPP

#include <memory>
#include <vector>
#include <array>
#include <iostream>

#include <boost/numeric/odeint.hpp>

#include "Matrix.h"
#include "phys_const.h"

#include "Field1D.hpp"

using namespace boost::numeric::odeint;

class Particle1D {
public:
	// constructor
	Particle1D( int m, int q, double x, double vx, std::shared_ptr<Field1D> f )
	: m_(m), q_(q), s_({x,vx}), f_(f), dist_(0) {}
	// m : mass [amu]
	// q : charge state
	// x : position in x [m]
	// vx: velocity in x [m/s]
	// f : shared pointer to Field1D instance

    // constructor (with stream output)
	//Particle1D( int m, int q, double x, double vx, std::shared_ptr<Field1D> f, std::ostream& out )
	//: m_(m), q_(q), s_({x,vx}), dist_(0), f_(f), out_(out) {}
	// m : mass [amu]
    // q : charge state
    // x : position in x [m]
    // vx: velocity in x [m/s]
	// f : shared pointer to Field1D instance
	// out : output stream for monitoring intermediate state
    
    virtual void move( double dt );
    // dt : time step [s]
    
    virtual void accelerate( double dt );
    // dt : time step [s]
	
	virtual void advance( double dt, double t = 0 );
	// dt : time step [s]
	// t  : initial time [s] (default = 0)
	
    void set_m( int m ) { m_ = m; }
    // m : mass [amu]
    
    void set_q( int q ) { q_ = q; }
    // q : charge state [1]
    
    // member functions to access private members
    int m() const { return m_; }        // [amu]
    int q() const { return q_; }	    //
    double x() const { return s_[0]; }     // [m]
    double vx() const { return s_[1]; }	// [m/s]
    double dist() const { return dist_; }   // [m]
    
    void set_event( int event ) { event_.push_back( event ); }
    // event : kind of event
    
    // record particle history
    void record() {
        hm_.push_back(m_);
        hq_.push_back(q_);
        hx_.push_back(s_[0]);
        hvx_.push_back(s_[1]);
    }

    // access to an element of record arrays
    int event( long i ) const { return event_.at(i); }
    int hm( long i ) const { return hm_.at(i); }
    int hq( long i ) const { return hq_.at(i); }
    double hx( long i ) const { return hx_.at(i); }
    double hvx( long i ) const { return hvx_.at(i); }
    
    // return number of events
    long ne() const { return event_.size(); }
    
    // return kinetic energy (classical)
    double ek() const { return 0.5*m_*si::amu*s_[1]*s_[1]/si::e; }	// [eV]

	// for boost odeint
	using state = std::array<double, 2>;	// definition of type of state
	void operator()( const state& x, state& dx, double t ) {
		dx[0] = x[1];
		dx[1] = q_*si::e/(m_*si::amu)*f_->ex(x[0]);
	}

private:
    
    int m_;  // mass number [amu]
    int q_;  // charge state
	state s_;	// particle state s_[0] = x, s_[1] = vx
	double dist_;   // total flight distance
	//std::ostream& out_;	// output stream for monitoring intermediate state (for debug)
	
	std::shared_ptr<Field1D> f_;	// pointer to an Field1D instance
    
    // containers for recording particle history
    std::vector<int> event_;
    std::vector<int> hm_;
    std::vector<int> hq_;
    std::vector<double> hx_;
    std::vector<double> hvx_;
	
};

#endif /* IEC_PARTICLE_HPP */
