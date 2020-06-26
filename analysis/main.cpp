//
//  main.cpp
//  analysis
//
//  Created by Jun on 2017/12/24.
//  Copyright © 2017年 Tokyo Tech. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "phys_const.h"

int main()
{
    using namespace std;

    // prepare an input stream for all particle data
    ifstream fin( "out.txt", ios::in );
    
    double phi;    // cathode potential [V]
    cout << "cathode voltage [-kV] = ";
    cin >> phi;
    phi *= -1e3;    // -kV -> V

    double xo;
    cout << "observation position [m] = ";
    cin >> xo;
    
    double dxo;
    cout << "half width of observation region [m] = ";
    cin >> dxo;
    
    int n;
    cout << "number of data points in output spectrum = ";
    cin >> n;
    
    // define ranges of space and velocity
    double xmin = -0.21;	// m
    double xmax = 0.21;	// m
    double dx = xmax/n;
    double vmax = sqrt(2*fabs(phi)*si::e/(1*si::amu))*1.1;
    double dv = vmax/n;
    double de = fabs(phi)/(2*n);
    
    int num_of_events = 27;
    vector<double> fx[num_of_events], fv[num_of_events], fek[num_of_events];    // for total distribution functions
    vector<double> fx1[num_of_events], fv1[num_of_events], fek1[num_of_events]; // for local distribution functions

    // set counters to zero
    for ( int i = 0; i < num_of_events; i++ ) {
        fx[i].resize(2*n,0);
        fv[i].resize(2*n,0);
        fek[i].resize(2*n,0);
        fx1[i].resize(2*n,0);
        fv1[i].resize(2*n,0);
        fek1[i].resize(2*n,0);
    }
    
    // load particle data and count the number of each event
    for ( long i = 0; ; i++ ) {
        int evt;
        double m, q, x, v;
        
        fin >> evt >> m >> q >> x >> v;
        if ( fin.eof() ) break;
		if ( evt == 99 ) continue;	// 190724
		if ( evt == 100 ) continue;	// 190730
		if ( evt == 101 ) continue;

        double ek = 0.5*m*si::amu*v*v/si::e;    // kinetic energy [eV]
        
        // record position, velocity, and kinetic energy of the particle
        long j = floor(x/dx) + n;
        long k = floor(v/dv) + n;
        long l = floor(ek/de);
        fx[evt].at(j)++;
        fv[evt].at(k)++;
        fek[evt].at(l)++;
        
        // if the event occurs in the region of interest, record position, velocity, and kinetic energy
        if ( x >= xo-dxo && x <= xo+dxo ) {
            long j = floor(x/dx) + n;
            long k = floor(v/dv) + n;
            long l = floor(ek/de);
            fx1[evt].at(j)++;
            fv1[evt].at(k)++;
            fek1[evt].at(l)++;
        }
    }
    
    // output distribution functions
    for ( int i = 0; i < num_of_events; i++ ) {
        
        stringstream fn( "" );
        fn << "fx" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout1( fn.str(), ios::out );   // total f(x)
        
        fn.str( "" );    // clear
        fn << "fv" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout2( fn.str(), ios::out );   // total f(v)
        
        fn.str( "" );    // clear
        fn << "fek" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout3( fn.str(), ios::out );   // total f(E)
        
        fn.str( "" );    // clear
        fn << "fxl" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout4( fn.str(), ios::out );   // local f(x)
        
        fn.str( "" );    // clear
        fn << "fvl" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout5( fn.str(), ios::out );   // local f(v)
        
        fn.str( "" );    // clear
        fn << "fekl" << setw( 2 ) << setfill( '0' ) << i << ".txt";
        ofstream fout6( fn.str(), ios::out );   // local f(E)
        
        // output distribution functions
        for ( long j = 0; j < 2*n; j++ ) {
            fout1 << xmin + (j+0.5)*dx << '\t' << fx[i].at(j) << '\n';
            fout2 << -vmax + (j+0.5)*dv << '\t' << fv[i].at(j) << '\n';
            fout3 << (j+0.5)*de << '\t' << fek[i].at(j) << '\n';
            fout4 << xmin + (j+0.5)*dx << '\t' << fx1[i].at(j) << '\n';
            fout5 << -vmax + (j+0.5)*dv << '\t' << fv1[i].at(j) << '\n';
            fout6 << (j+0.5)*de << '\t' << fek1[i].at(j) << '\n';
        }
    }
    
    // prepare an input stream for fusion reactions
    ifstream fin2( "fusion.txt", ios::in );
    
    // reset counters
    fill( fx[0].begin(), fx[0].end(), 0 );
    fill( fv[0].begin(), fv[0].end(), 0 );
    fill( fek[0].begin(), fek[0].end(), 0 );
    
    // load fusion reaction data and count the number of fusion event
    for ( long i = 0;; i++ ) {
        int evt;
        double x, v;
        
        fin2 >> evt >> x >> v;
        if ( fin2.eof() ) break;
        
        double ek = 0.5*si::amu*v*v/si::e;    // eV/u
        
        long j = floor(x/dx) + n;
        long k = floor(v/dv) + n;
        long l = floor(ek/de);
        fx[0].at(j)++;
        fv[0].at(k)++;
        fek[0].at(l)++;
    }
    
    // output distribution functions of particles causing fusion reactions
    ofstream fout7( "fx_f.txt", ios::out );
    ofstream fout8( "fv_f.txt", ios::out );
    ofstream fout9( "fek_f.txt", ios::out );

    for ( long i = 0; i < 2*n; i++ ) {
        fout7 << xmin + (i+0.5)*dx << '\t' << fx[0].at(i) << '\n';
        fout8 << -vmax + (i+0.5)*dv << '\t' << fv[0].at(i) << '\n';
        fout9 << (i+0.5)*de << '\t' << fek[0].at(i) << '\n';
    }
    
    // prepare an input stream for Ha emission data
    ifstream fin3( "ha.txt", ios::in );
    
    // reset counters
    fill( fx[0].begin(), fx[0].end(), 0 );
    fill( fv[0].begin(), fv[0].end(), 0 );
    
    for ( long i = 0;; i++ ) {
        int evt;
        double x, v;
        
        fin3 >> evt >> x >> v;
        if ( fin3.eof() ) break;
        
        long j = floor(x/dx) + n;
        long k = floor(v/dv) + n;
        fx[0].at(j)++;
        fv[0].at(k)++;
    }
    
    // output distribution functions of particles causing Ha emission
    ofstream fout10( "fx_ha.txt", ios::out );
    ofstream fout11( "fv_ha.txt", ios::out );

    for ( long i = 0; i < 2*n; i++ ) {
        fout10 << xmin + (i+0.5)*dx << '\t' << fx[0].at(i) << '\n';
        fout11 << -vmax + (i+0.5)*dv << '\t' << fv[0].at(i) << '\n';
    }

}
