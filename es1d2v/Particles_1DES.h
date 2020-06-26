//
//  Particles_1DES.h
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/05.
//  Last Updated by Jun Hasegawa on 2019/07/12.
//
//  Copyright (c) 2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#ifndef Particles_1DES_h
#define Particles_1DES_h

#include <iostream>
#include <fstream>
#include <vector>

#include "Matrix.h"

#include "Particle_1DES.h"
class Field_1DES;

class Particles_1DES {
public:
	Particles_1DES( long np, double m, double q, double x0, double x1, double kt, double vd, long np0 = 0, double kt0 = 0, double vd0 = 0 );
	// constructor
    // np: number of particles initially generated
    // m:  mass
    // q:  charge
    // x0: left boundary of region where particles are initially generated
    // x1: right boundary of region where particles are initially generated
    // kt: temperature of particles initially generated
    // vd: drift velocity of particles initially generated
	// np0:	number of particles per step generated uniformly
	// kt0: temperature of particles generated uniformly
	// vd0: drift velocity of particles generated uniformly

	~Particles_1DES();  // destructor
    
    virtual void advance( double dt );
	// advance particles
	// dt: time step
	
    virtual void accelerate( double dt );
	// accelerate particles (used only for initial push back)
	// dt: time step
	
	virtual void weight( int s ) const;
	// weight particle charges and currents to grids
	// s: index of particle specie.
	
	//virtual void inject( double dt );
	// inject particles
	// dt: time step
	
	//virtual void generate_uniform( long np, double x0, double x1, double dt );
	// generate plasma particles uniformly
	// np: number of ion and electron pairs generated per step
	// x0: left boundary of a region where particles are generated
	// x1: right boundary of a region where particles are generated
	// dt: time step

	virtual void generate_uni( double dt );
	// generate plasma particles uniformly
	// dt: time step
	
	virtual void generate_exp( double dt, double alpha );
	// generate plasma particles based on exponential distribution
	// dt: time step
	
    virtual void add_particle( double m, double q, double x, double vx, double vy );
	// add a particle
	// m: mass
	// q: charge
	// x: position in x
	// vx: velocity in x
	// vy: velocity in y
	
	virtual void add_field( Field_1DES* field );
	// add reference to an instance of Field_1DES class

    virtual void snapshot( std::ostream& out, bool title = true ) const;
	// output data of particles to stream

    // member functions to access private members
    virtual long np() const { return ps_.size(); }
	virtual double j0( int i ) const { return j0_.at(i); }
	virtual double j1( int i ) const { return j1_.at(i); }
	virtual double m() { return m_; }
	virtual double q() { return q_; }
	virtual long np0() const { return np0_; }
	
private:
	virtual int get_region( double x ) const;
	// returns an index of the region containing the position x
	
	double pdist_exp( double x, double a );
	// probability distribution function for particle generation with exponential profile

    long pid_;  // particle ID (initial value)
    int nr_;	// number of region
    double xmin_;	// left boundary
	double xmax_;	// right boundary
	
	double m_;	// mass
	double q_;	// charge
	
	long np0_;	// number of particles per step generated uniformly
	double kt0_;	// temperature of particles generated uniformly
	double vd0_;	// drift velocity of particles generated uniformly
	
	std::vector<double> boundary_;	// boundary positions
	std::vector<Field_1DES*> field_;
    std::vector<Particle_1DES*> ps_;
	std::vector<double> j0_, j1_;	// conduction current at both boundaries
};

#endif
