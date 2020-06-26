//
//  Particle1D.cpp
//
//  (C) Copyright Jun Hasegawa 2017.
//
//  19 Jan 2018 - Initial Revision (Jun Hasegawa)
//	12 Sep 2019 - Integration with boost odeint library was introduced
//



#include "Particle1D.hpp"

#include <iostream>
#include <fstream>
#include <cmath>


void Particle1D::advance( double dt, double t )
{
	using namespace boost::numeric::odeint;

	// for debug
	//integrate_adaptive( make_controlled( 1.0e-10 , 1.0e-6 , runge_kutta_dopri5<state>() ),
	//				   *this, s_, t, dt, dt, [this]( const state& x , double t ) { out_ << t << '\t' << x[0] << '\t' << x[1] << '\n'; } );

	integrate_adaptive( make_controlled( 1.0e-10 , 1.0e-6 , runge_kutta_dopri5<state>() ), *this, s_, t, dt, dt );
}


void Particle1D::move( double dt )
{
    double dx = s_[1]*dt;	// dx = vx*dt
    s_[0] += dx;	// x = x+dx
    
    dist_ += fabs( dx );    // record total flight distance
}

void Particle1D::accelerate( double dt )
{    
    double qm = q_*si::e/(m_*si::amu);
	double ex = f_->ex(s_[0]);
	
    // acceleration
    double dv = qm*ex*dt;
    s_[1] += dv;
}

