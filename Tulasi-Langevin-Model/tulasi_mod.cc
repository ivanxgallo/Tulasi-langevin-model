#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <fstream>
#include <complex>
#include "randonailo.h"
#include "vectores.h"
#include "matrices.h"

using namespace std ;


int main()

{
	/*------ Cambia la semilla del Random -------*/
	srand((int)time(NULL)) ;

	/*------ IDENTIDAD 3x3 -------*/
	complex<double> i(0,1.) ;

	/*------- parametros para el paso de tiempo --------*/
	double dt = 0.001, t = 0, T = 1e3, m ;

	/*------- parametros naturales del sistema --------*/
	double om = 1., l = 0.1, g ;

	/*-------- variables para el RK4 --------*/
	complex<double> u = 0, ut = u, u1, u2, u3, u4 ;

	//----------------------------------------------------------------------------//
	//-------------------- archivos a escribir (ofstream) -----------------------//

	ofstream vel("u_k.dat") ;

	vel << "velocidad" << endl ;

 	/*---------------------------------- DINAMICA ----------------------------------------*/
	g = normailon_complex() ;
	while(t<T)
	{

		g = memonailon_complex(g,m) ;
		u = ut ;
		u1 = i*om*u - 0.5*abs(u)*u/l + pow(l,-1./3.)*g ;

		u = ut + u1*dt/2. ;
		u2 = i*om*u - 0.5*abs(u)*u/l ;

		u = ut + u2*dt/2. ;
		u3 = i*om*u - 0.5*abs(u)*u/l ;

		g = memonailon_complex(g,m) ;
		u = ut + u3*dt ;
		u4 = i*om*u - 0.5*abs(u)*u/l + pow(l,-1./3.)*g ;

		m = (u1 + 2*u2 + 2*u3 + u4)/6. ;

		ut += m*dt ;


		vel << real(u) << endl ;



		t += dt ;
		cout << (t/T)*100. << " % " << endl ;

	}

    vel.close() ;

	return 0;
}
