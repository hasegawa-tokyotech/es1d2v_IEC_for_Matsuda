//
//  main.cpp
//  pic2mc
//
//  Created by Jun on 2019/07/19.
//  Copyright © 2019 Jun. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>

int main()
{
    using namespace std;
	
    int i;
    cout << "field file number = ";
    cin >> i;
	
    // load normalization factors
    ifstream fin0( "norm_params.txt", ios::in );
    if ( !fin0 ) throw runtime_error( "file not open" );
    char buf[255];
    double x_norm, phi_norm, ef_norm;
    fin0 >> buf >> x_norm
         >> buf >> buf
		 >> buf >> phi_norm;
	fin0.getline( buf, 255 );
	fin0.getline( buf, 255 );
	fin0.getline( buf, 255 );
	fin0.getline( buf, 255 );
	fin0 >> buf >> ef_norm;
    
    // load field data
    char fn[255];
    vector<double> x, phi, ex;
	sprintf( fn, "fi%01d_%03d.txt", 0, i );
	ifstream fin1( fn, ios::in );
    if ( !fin1 ) throw runtime_error( "file not open" );
    fin1.getline( buf, 255 );   // skip first line
    for (;;) {
        double xi, phii, exi;
        fin1 >> xi >> buf >> phii >> exi;
		fin1.getline( buf, 255 );

		if ( fin1.eof() ) break;

        // store values in SI unit
        x.push_back( xi*x_norm );
        phi.push_back( phii*phi_norm );
        ex.push_back( exi*ef_norm );
        
    }
	long n = x.size();	// number of loaded data points
	cout << n << " data points were loaded\n";
	
	// store field data in the format for MC code
	double lx = x.at(n-1);
	double v0 = phi.at(n-1);	// cathode voltage
	vector<double> xx, pp, ee;
	for ( long i = 0; i <= n-2; i++ ) {
		xx.push_back(x.at(i) - lx);
		pp.push_back(phi.at(i));
		ee.push_back(ex.at(i));
	}
	xx.push_back(x.at(n-1)-lx);
	pp.push_back(phi.at(n-1));
	ee.push_back(ex.at(0.0));
	for ( long i = n-2; i >= 0; i-- ) {
		xx.push_back(lx - x.at(i));
		pp.push_back(phi.at(i));
		ee.push_back(-ex.at(i));
	}
	long nn = xx.size();	// number of data points
    //cout << "xx.size() = " << xx.size() << '\n';
    //cout << "pp.size() = " << pp.size() << '\n';
    //cout << "ee.size() = " << ee.size() << '\n';

	// smoothing
	int m;
	cout << "number of neighboring data points for smoothing = ";
	cin >> m;
	
	x.clear();
	phi.clear();
	ex.clear();
	for ( long i = 0; i < m; i++ ) {
		double psum = 0, esum = 0;
		for ( long j = 0; j <= 2*i; j++ ) {
			psum += pp.at(j);
			esum += ee.at(j);
			//cout << '[' << j << ']';
		}
		x.push_back(xx.at(i));
		phi.push_back(psum/(2*i+1));
		ex.push_back(esum/(2*i+1));
        //cout << " : " << 2*i+1 << '\n';
	}
	for ( long i = m; i < nn-m; i++ ) {
		double psum = 0, esum = 0;
		for ( long j = i-m; j <= i+m; j++ ) {
			psum += pp.at(j);
			esum += ee.at(j);
			//cout << '[' << j << ']';
		}
		x.push_back(xx.at(i));
		phi.push_back(psum/(2*m+1));
		ex.push_back(esum/(2*m+1));
        //cout << " : " << 2*m+1 << '\n';
	}
    for ( long i = m-1; i >= 0; i-- ) {
        double psum = 0, esum = 0;
        for ( long j = 2*i; j >= 0; j-- ) {
            psum += pp.at(nn-1-j);
            esum += ee.at(nn-1-j);
            //cout << '[' << n-1-j << ']';
        }
        x.push_back(xx.at(nn-1-i));
        phi.push_back(psum/(2*i+1));
        ex.push_back(esum/(2*i+1));
        //cout << " : " << 2*i+1 << '\n';
    }
	
	// get collection factor for cathode voltage
	double vc = v0/phi.at(n-1);
	
	// get E field from phi (for debug)
	/*
	vector<double> exx;
	exx.push_back( -(phi.at(1)-phi.at(0))/(x.at(1)-x.at(0)) );
	for ( long i = 1; i < nn-1; i++ ) {
		double dx = x.at(i+1)-x.at(i-1);
		double dphi = phi.at(i+1)-phi.at(i-1);
		exx.push_back(-dphi/dx);
	}
	exx.push_back( -(phi.at(nn-1)-phi.at(nn-2))/(x.at(nn-1)-x.at(nn-2)) );
	*/
	
    // output field data
    ofstream fout( "pic_field.txt", ios::out );
    fout << "#x(m)\tphi(V)\tex(V/m)\n";
    for ( long i = 0; i < nn; i++ ) {
        fout << x.at(i) << '\t' << phi.at(i)*vc << '\t' << ex.at(i)*vc << '\n';
	}
    
	return 0;
}
