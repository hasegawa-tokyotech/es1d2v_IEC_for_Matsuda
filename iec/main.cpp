//
//  main.cpp
//  iec
//
//  Created by Jun on 2017/10/28.
//  Copyright © 2017年 Tokyo Tech. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <random>
#include <stdexcept>

#include "phys_const.h"
#include "Particle1D.hpp"
#include "Field1D.hpp"
#include "CrossSection.hpp"
#include "HydrogenCrossSection.hpp"

#define IEC_DEUTERIUM // comment out in hydrogen case
//#define IEC_SCATTERING  // comment out if ignore scattering
//#define IEC_USE_COMSOL_FIELD //comment out if use PIC field
#define IEC_ADAPTIVE_TIME_STEP  // enable adaptive time step algorithm

int main()
{
    using namespace std;
    
    // prepare a random generator using C++ standard library
    mt19937 mt(12345);
    uniform_real_distribution<> rand_uniform(0,1);	// prepare a random generator for uniform distribution

    // IEC operation parameters
#ifdef IEC_USE_COMSOL_FIELD
    double phi;	// cathode potential [V]
    cout << "cathode voltage [-kV] = ";
    cin >> phi;
    phi *= -1e3;    // -kV -> V
    auto f = make_shared<Field1D>( "phi_3d.txt", -30e3, phi );	// generate an instance of Field class
#else
    auto f = make_shared<Field1D>( "pic_field.txt" );
    double phi = f->phi(0.0);  // potential value at x = 0
#endif
	
    double ek_max = fabs(phi);    // maximum kinetic energy in eV
	
	// get maximum electric field
	double ex_max = 0;
    for ( int i = 0; i < f->size(); i++ ) {
        double ex = f->ex(f->x(i));
        if ( ex > ex_max ) ex_max = ex;
    }
    
    double p;	// background H2 gas pressure [Pa]
    cout << "background gas pressure [Pa] = ";
    cin >> p;

    double tg = 353;    // gas temperature (<- glass tube temperature ~ 80ºC)
    double ng = p/si::kb/tg;	// gas molecular number density [m^-3]
    double kt = tg*si::kb;    // gas temperature [J]
#ifdef IEC_DEUTERIUM
    double mg = 4*si::amu;    // mass of gas molecular (D2)
#else
    double mg = 2*si::amu;    // mass of gas molecular (H2)
#endif
    double sig_v = sqrt(kt/mg);    // standard deviation of velocity distribution of gas molecule
    normal_distribution<> rand_norm(0,sig_v);	// prepare a random generator for Gaussian (normal) distribution


    // cross sections for H+
    HydrogenCrossSection cs06( "cs_tabata.txt", 6 );  // 1: H+ + H2 -> fast H
    HydrogenCrossSection cs08( "cs_tabata.txt", 8 );  // 2: H+ + H2 -> H_alpha
    HydrogenCrossSection cs10( "cs_tabata.txt", 10 );  // 3: H+ + H2 -> momentum transfer
    
    // cross sections for H2+
    HydrogenCrossSection cs11( "cs_tabata.txt", 11 );  // 21: H2+ + H2 -> H3+ + H
    HydrogenCrossSection cs12( "cs_tabata.txt", 12 );  // 4: H2+ + H2 -> slow H2+
    HydrogenCrossSection cs14( "cs_tabata.txt", 14 );  // 5: H2+ + H2 -> fast H+
    HydrogenCrossSection cs16( "cs_tabata.txt", 16 );  // 6: H2+ + H2 -> H_alpha

    // cross sections for H3+
    HydrogenCrossSection cs18( "cs_tabata.txt", 18 );  // 7: H3+ + H2 -> fast H+
    HydrogenCrossSection cs19( "cs_tabata.txt", 19 );  // 8: H3+ + H2 -> fast H2+
    HydrogenCrossSection cs20( "cs_tabata.txt", 20 );  // 9: H3+ + H2 -> fast H
    HydrogenCrossSection cs21( "cs_tabata.txt", 21 );  // 10: H3+ + H2 -> fast H2
    HydrogenCrossSection cs23( "cs_tabata.txt", 23 );  // 11: H3+ + H2 -> H_alpha
    HydrogenCrossSection cs25( "cs_tabata.txt", 25 );  // 12: H3+ + H2 -> momentum transfer
    
    // cross sections for H
    HydrogenCrossSection cs30( "cs_tabata.txt", 30 );  // 13: H + H2 -> slow H2+
    HydrogenCrossSection cs31( "cs_tabata.txt", 31 );  // 14: H + H2 -> fast H+
    HydrogenCrossSection cs33( "cs_tabata.txt", 33 );  // 15: H + H2 -> H_alpha
    HydrogenCrossSection cs36( "cs_tabata.txt", 36 );  // 16: H + H2 -> momentum transfer

    // cross sections for H2
    HydrogenCrossSection cs38( "cs_tabata.txt", 38 );  // 17: H2 + H2 -> fast H2+
    HydrogenCrossSection cs43( "cs_tabata.txt", 43 );  // 18: H2 + H2 -> fast H+
    HydrogenCrossSection cs44( "cs_tabata.txt", 44 );  // 19: H2 + H2 -> H_alpha
    HydrogenCrossSection cs46( "cs_tabata.txt", 46 );  // 20: H2 + H2 -> momentum transfer

    // cross section for D-D fusion
    CrossSection csf( "dd_fusion.txt" );	// D + D -> 3He + n / T + p

#ifdef IEC_DEUTERIUM
    double cd = 0.5;	// energy correction for deuterium
    int md = 2;	// mass correction for deuterium
    double cf = 1e8;	// artificial enhancement factor of fusion reaction
#else
    double cd = 1;
    int md = 1;
    double cf = 0;
#endif

#ifdef IEC_SCATTERING
    int cs = 1;
#else
    int cs = 0;
#endif

#ifdef IEC_ADAPTIVE_TIME_STEP
    const double de_max = fabs(phi)*1e-3;
    //const double dx_max = de_max/ex_max;
	//const double dx_min = dx_max/10.0;	// minimum distance per step [m]
	const double dx_max = 1e-4;	// maximum flight distance per step [m]
	const double dx_min = 1e-5;	// minimum flight distance per step [m]
#endif
    const long ni = 100;	// <- change this number for changing total number of iterations
    const long np = 1000;

    ofstream fout1( "out.txt", ios::out );  // for all data
    ofstream fout2( "ha.txt", ios::out );   // for Ha emission
    ofstream fout3( "fusion.txt", ios::out );   // for fusion reaction
    ofstream fout4( "error.txt", ios::out );   // for error event

    // region for calculation
#ifdef IEC_USE_COMSOL_FIELD
    double xmin = -0.2; // m
    double xmax = 0.2;  // m
#else
	double xmin = f->xmin();
	double xmax = f->xmax();
#endif
	double lx = xmax - xmin;	// region length
    
    // containers for recording Ha emission
    vector<int> pa;
    vector<double> xa;
    vector<double> va;
    // containers for recording fusion reaction
    vector<int> pf;
    vector<double> xf;
    vector<double> vf;
	// containers for recording forced termination
	vector<int> pt;
	vector<double> xt;
	vector<double> vt;
    
    // containers for counting events
    int num_of_events = 27;
    vector<int> event(num_of_events,0);
    
    double eth;
    cout << "theshold ion energy for initial distribution [eV] = ";
    cin >> eth;
    
    
    for ( long i = 0; i < ni; i++ ) {
        cout << '[' << i << ']';
        cout.flush();
        
        unique_ptr<Particle1D> p[np];	// pointers to Particle instances
        
        // prepare initial set of particles
        for ( long j = 0; j < np; j++ ) {
            
            int m;
			//double gen_ratio = 20;	// generation ratio of H2+ to H+ (= N_H2+/N_H+)
            //double r = rand_uniform(mt);
            //if ( r < 1/(gen_ratio+1) ) m = md;	// generate H+/D+
            //else m = 2*md;	// generate H2+/D2+
			m = 2*md;	// always generate H2+/D2+
			
            double xi;
            for (;;) {
	            xi = xmin + lx*rand_uniform(mt);	// initial position (uniform)
                if ( f->phi(xi)-phi > eth ) break;	// threshold of initial potential energy
            }
            double vi = rand_norm(mt);	// initial velocity
            p[j] = make_unique<Particle1D>( m, 1, xi, vi, f ); // create an instance of Particle
            
			
			double dt = 1e-9;	// initial time step [s]
			p[j]->accelerate( -0.5*dt );    // push back a half step for the leap frog method
			p[j]->set_event( 0 );   // initial state
			p[j]->record();
			event[0]++;

            for ( long k = 0; ; k++ ) {
				
				// adjust time step so that particle can move minimum distance by one step
#ifdef IEC_ADAPTIVE_TIME_STEP
				double vx = p[j]->vx();
				if ( fabs(vx*dt) > dx_max ) {
					for (;;) {
						dt /= 2;
						if ( fabs(vx*dt) <= dx_max ) break;
					}
				}
				else if ( fabs(vx*dt) < dx_min ) {
					for (;;) {
						dt *= 2;
						if ( fabs(vx*dt) >= dx_min ) break;
					}
				}
#endif
                double ngvdt = ng*fabs(p[j]->vx())*dt;

            	// check collision
                double r1 = rand_uniform(mt);
                
                if ( p[j]->m() == md && p[j]->q() == 1 ) {   // H+
                    double ek = p[j]->ek();
                    double s1 = cs06.sigma( ek*cd );
                    double s2 = cs08.sigma( ek*cd );
                    double s3 = cs10.sigma( ek*cd ) * cs;
                    double s4 = csf.sigma( ek ) * cf;

                    double p1 = s1*ngvdt;
                    double p2 = p1 + s2*ngvdt;
                    double p3 = p2 + s3*ngvdt;
                    double p4 = p3 + s4*ngvdt;
					if ( p4 >= 1 ) throw runtime_error( "too large time step" );
                    
                    if ( r1 < p1 ) {    // 1: H+ + H2 -> fast H (charge exchange)
                        p[j]->set_event( 1 );
						p[j]->set_q( 0 );
                        p[j]->record();
                        event[1]++;
                    }
                    else if ( r1 >= p1 && r1 < p2 ) {   // 2: H+ + H2 -> H_alpha
                        p[j]->set_event( 2 );
                    	p[j]->set_q( 0 );
                        p[j]->record();
                        // record position and velocity with Ha emission
                        pa.push_back( 2 );
                        xa.push_back( p[j]->x() );
                        va.push_back( p[j]->vx() );
                        event[2]++;
                    }
                    else if ( r1 >= p2 && r1 < p3 ) {	// 3: H+ + H2 -> momentum transfer
						p[j]->set_event( 3 );
                        p[j]->record();
                        event[3]++;
                        break;	// stop
                    }
                    else if ( r1 >= p3 && r1 < p4 ) {    // 22: D-D fusion
                        p[j]->set_event( 22 );
                        p[j]->record();
                        event[22]++;
                        // record position and velocity with fusion
                        pf.push_back( 22 );
                        xf.push_back( p[j]->x() );
                        vf.push_back( p[j]->vx() );
                        break;    // stop
                    }
                }
                
                else if ( p[j]->m() == 2*md && p[j]->q() == 1 ) {  // H2+
                    double ek = p[j]->ek();
                    double s1 = cs12.sigma( ek*cd );
                    double s2 = cs14.sigma( ek*cd );
                    double s3 = cs16.sigma( ek*cd );
                    double s4 = csf.sigma( ek*0.5 ) * cf;
                    double s5 = cs11.sigma( ek*cd );
                    
                    double p1 = s1*ngvdt;
                    double p2 = p1 + s2*ngvdt;
                    double p3 = p2 + s3*ngvdt;
                    double p4 = p3 + s4*ngvdt;
                    double p5 = p4 + s5*ngvdt;
					if ( p5 >= 1 ) throw runtime_error( "too large time step" );

                    if ( r1 < p1 ) {    // 4: H2+ + H2 -> slow H2+
                        p[j]->set_event( 4 );
                        p[j]->set_q( 0 );
                        p[j]->record();
                        event[4]++;
                    }
                    else if ( r1 >= p1 && r1 < p2 ) {	// 5: H2+ + H2 -> fast H+
                        p[j]->set_event( 5 );
                        p[j]->set_m( md );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[5]++;
                    }
                    else if ( r1 >= p2 && r1 < p3 ) {   // 6: H2+ + H2 -> H_alpha
                        p[j]->set_event( 6 );
						p[j]->set_m( md );
                        p[j]->set_q( 0 );
                        p[j]->record();
                        event[6]++;
                        // record position and velocity with Ha emission
                        pa.push_back( 6 );
                        xa.push_back( p[j]->x() );
                        va.push_back( p[j]->vx() );
                    }
                    else if ( r1 >= p3 && r1 < p4 ) {    // 23: D-D fusion
						p[j]->set_event( 23 );
                        p[j]->record();
                        event[23]++;
                        // record position and velocity with fusion
                        pf.push_back( 23 );
                        xf.push_back( p[j]->x() );
                        vf.push_back( p[j]->vx() );
                        break;    // stop
                    }
                    else if ( r1 >= p4 && r1 < p5 ) {   // 21: H2+ + H2 -> H3+ + H
                        p[j]->set_event( 21 );
                        p[j]->set_m( 3*md );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[21]++;
                    }
                }

                else if ( p[j]->m() == 3*md && p[j]->q() == 1 ) {  // H3+
                    double ek = p[j]->ek();
                    double s1 = cs18.sigma( ek*cd );
                    double s2 = cs19.sigma( ek*cd );
                    double s3 = cs20.sigma( ek*cd );
                    double s4 = cs21.sigma( ek*cd );
                    double s5 = cs23.sigma( ek*cd );
                    double s6 = cs25.sigma( ek*cd ) * cs;
                    double s7 = csf.sigma( ek*0.333 ) * cf;

                    double p1 = s1*ngvdt;
                    double p2 = p1 + s2*ngvdt;
                    double p3 = p2 + s3*ngvdt;
                    double p4 = p3 + s4*ngvdt;
                    double p5 = p4 + s5*ngvdt;
                	double p6 = p5 + s6*ngvdt;
                    double p7 = p6 + s7*ngvdt;
					if ( p7 >= 1 ) throw runtime_error( "too large time step" );

                    if ( r1 < p1 ) {    // 7: H3+ + H2 -> fast H+
                        p[j]->set_event( 7 );
                        p[j]->set_m( md );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[7]++;
                    }
                    else if ( r1 >= p1 && r1 < p2 ) {    // 8: H3+ + H2 -> fast H2+
                        p[j]->set_event( 8 );
                        p[j]->set_m( 2*md );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[8]++;
                    }
                    else if ( r1 >= p2 && r1 < p3 ) {   // 9: H3+ + H2 -> fast H
                        p[j]->set_event( 9 );
                        p[j]->set_m( md );
                        p[j]->set_q( 0 );
                        p[j]->record();
                        event[9]++;
                    }
                    else if ( r1 >= p3 && r1 < p4 ) {	// 10: H3+ + H2 -> fast H2
                        p[j]->set_event( 10 );
                        p[j]->set_m( 2*md );
						p[j]->set_q( 0 );
                        p[j]->record();
                        event[10]++;
                    }
                    else if ( r1 >= p4 && r1 < p5 ) {	// 11: H3+ + H2 -> H_alpha
                        p[j]->set_event( 11 );
                        p[j]->set_m( md );
                        p[j]->set_q( 0 );
                        p[j]->record();
                        event[11]++;
                        // record position and velocity
                        pa.push_back( 11 );
                        xa.push_back( p[j]->x() );
                        va.push_back( p[j]->vx() );
                    }
                    else if ( r1 >= p5 && r1 < p6 ) {    // 12: H3+ + H2 -> momentum transfer
                        p[j]->set_event( 12 );
                        p[j]->record();
                        event[12]++;
                        break;	// stop
                    }
                    else if ( r1 >= p6 && r1 < p7 ) {    // 24: D-D fusion
                        p[j]->set_event( 24 );
                        p[j]->record();
                        event[24]++;
                        // record position and velocity
                        pf.push_back( 24 );
                        xf.push_back( p[j]->x() );
                        vf.push_back( p[j]->vx() );
                        break;    // stop
                    }
                }
                
                else if ( p[j]->m() == md && p[j]->q() == 0 ) {  // H
                    double ek = p[j]->ek();
                    double s1 = cs30.sigma( ek*cd );
                    double s2 = cs31.sigma( ek*cd );
                    double s3 = cs33.sigma( ek*cd );
                    double s4 = cs36.sigma( ek*cd ) * cs;
                    double s5 = csf.sigma( ek ) * cf;
                    
                    double p1 = s1*ngvdt;
                    double p2 = p1 + s2*ngvdt;
                    double p3 = p2 + s3*ngvdt;
                    double p4 = p3 + s4*ngvdt;
                    double p5 = p4 + s5*ngvdt;
					if ( p5 >= 1 ) throw runtime_error( "too large time step" );

                    if ( r1 < p1 ) {    // 13: H + H2 -> slow H2+
                        p[j]->set_event( 13 );
                        p[j]->record();
                        event[13]++;
                    }
                    else if ( r1 >= p1 && r1 < p2 ) {    // 14: H + H2 -> fast H+
                        p[j]->set_event( 14 );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[14]++;
                    }
                    else if ( r1 >= p2 && r1 < p3 ) {   // 15: H + H2 -> H_alpha
                        p[j]->set_event( 15 );
                        p[j]->record();
                        event[15]++;
                        // record position and velocity
                        pa.push_back( 15 );
                        xa.push_back( p[j]->x() );
                        va.push_back( p[j]->vx() );
                    }
                    else if ( r1 >= p3 && r1 < p4 ) {    // 16: H + H2 -> momentum transfer
                        p[j]->set_event( 16 );
                        p[j]->record();
                        event[16]++;
                        break;	// stop
                    }
                    else if ( r1 >= p4 && r1 < p5 ) {    // 25: D-D fusion
                        p[j]->set_event( 25 );
                        p[j]->record();
                        event[25]++;
                        // record position and velocity
                        pf.push_back( 25 );
                        xf.push_back( p[j]->x() );
                        vf.push_back( p[j]->vx() );
                        break;    // stop
                    }
                }
                else if ( p[j]->m() == 2*md && p[j]->q() == 0 ) {  // H2
                    double ek = p[j]->ek();
                    double s1 = cs38.sigma( ek*cd );
                    double s2 = cs43.sigma( ek*cd );
                    double s3 = cs44.sigma( ek*cd );
                    double s4 = cs46.sigma( ek*cd ) * cs;
                    double s5 = csf.sigma( ek*0.5 ) * cf;
                    
                    double p1 = s1*ngvdt;
                    double p2 = p1 + s2*ngvdt;
                    double p3 = p2 + s3*ngvdt;
                    double p4 = p3 + s4*ngvdt;
                    double p5 = p4 + s5*ngvdt;
					if ( p5 >= 1 ) throw runtime_error( "too large time step" );

                    
                    if ( r1 < p1 ) {    // 17: H2 + H2 -> fast H2+
                        p[j]->set_event( 17 );
                        p[j]->set_q( 1 );
                        p[j]->record();
                        event[17]++;
                    }
                    else if ( r1 >= p1 && r1 < p2 ) {    // 18: H2 + H2 -> fast H+
                        p[j]->set_event( 18 );
                        p[j]->set_m( md );
						p[j]->set_q( 1 );
                        p[j]->record();
                        event[18]++;
                    }
                    else if ( r1 >= p2 && r1 < p3 ) {   // 19: H2 + H2 -> H_alpha
                        p[j]->set_event( 19 );
                        p[j]->record();
                        event[19]++;
                        // record position and velocity
                        pa.push_back( 19 );
                        xa.push_back( p[j]->x() );
                        va.push_back( p[j]->vx() );
                    }
                    else if ( r1 >= p3 && r1 < p4 ) {    // 20: H2 + H2 -> momentum transfer
                        p[j]->set_event( 20 );
                        p[j]->record();
                        event[20]++;
                        break;    // stop
                    }
                    else if ( r1 >= p4 && r1 < p5 ) {    // 26: D-D fusion
                        p[j]->set_event( 26 );
                        p[j]->record();
                        event[26]++;
                        // record position and velocity
                        pf.push_back( 26 );
                        xf.push_back( p[j]->x() );
                        vf.push_back( p[j]->vx() );
                        break;    // stop
                    }
                }
                
                // if no collision, proceed the particle
                p[j]->move( dt );
                
                // if out of range, terminate
                if ( p[j]->x() < xmin || p[j]->x() > xmax ) {
                    p[j]->set_event( 99 );
                    p[j]->record();
                    break;  // stop
                }
                
                // acceleration
                if ( p[j]->q() != 0 ) {
                    p[j]->accelerate( dt );
                    if ( p[j]->ek() > ek_max*1.1 ) {
                        cout << '!';
                        p[j]->set_event( 101 );
                        p[j]->record();
                        pt.push_back( 101 );
                        xt.push_back( p[j]->x() );
                        vt.push_back( p[j]->vx() );
                        break;
                    }
				}
				
				// eliminate particles taking too much time (20190729)
				if ( k > 1e7 ) {
					p[j]->set_event( 100 );
					p[j]->record();
					cout << '*';
					// record position and velocity of terminated particle
					pt.push_back( 100 );
					xt.push_back( p[j]->x() );
					vt.push_back( p[j]->vx() );
					break;  // stop
				}
				
            }
        }
        
        // record history of particles
        for ( long j = 0; j < np; j++ ) {
            long ne = p[j]->ne();
            for ( long k = 0; k < ne; k++ )
                fout1 << p[j]->event(k) << '\t' << p[j]->hm(k) << '\t' << p[j]->hq(k) << '\t' << p[j]->hx(k) << '\t' << p[j]->hvx(k) << '\n';
        }
    }
    cout << endl;

    cout << event[0] << " ions are generated\n";
    for ( int i = 1; i < num_of_events; i++ )
        cout << "Event #" << i << " : " << event[i] << '\n';
    cout << pf.size() << " fusion reactions\n";
	cout << pt.size() << " forced terminations\n";

    // record positions and velocities of particles when Ha emissions occur
    for ( long i = 0; i < pa.size(); i++ )
        fout2 << pa[i] << '\t' << xa[i] << '\t' << va[i] << '\n';

    // record positions and velocities of particles when D-D fusion reactions occur
    for ( long i = 0; i < pf.size(); i++ )
        fout3 << pf[i] << '\t' << xf[i] << '\t' << vf[i] << '\n';
    
    // record positions and velocities of particles when errors occur
    for ( long i = 0; i < pt.size(); i++ )
        fout4 << pt[i] << '\t' << xt[i] << '\t' << vt[i] << '\n';

    return 0;
}
