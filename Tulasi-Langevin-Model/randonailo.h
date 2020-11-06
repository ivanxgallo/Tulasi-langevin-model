#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "vectores.h"

using namespace std;



// genera un numero random entre ini y fin, preci indica cuantos decimales de pricisión requerimos
double randonilo(int ini=0, int fin=1, int preci=4)
{
    int p=pow(10,preci);
    return (rand() % ((fin-ini)*p+1) + ini*p)/(p*1.0);
}

double randoneilo(double xi=0., double xf=1., int preci =4)
{
    return randonilo(0,1,preci)*(xf - xi) + xi;
}


double normanailon(int preci = 4)
{
    double min = pow(10,-preci) ;
    double xx = randoneilo(min,1), yy = randoneilo(min,1) ;

    return cos(2*M_PI*yy)*sqrt(-2*log(xx)) ;
} 


complex<double> normanailon_complex(int preci = 4)
{
    double min = pow(10,-preci) ;
    double xx = randoneilo(min,1), yy = randoneilo(min,1) ;
    complex<double> x(cos(2*M_PI*yy)*sqrt(-2*log(xx)),sin(2*M_PI*yy)*sqrt(-2*log(xx))) ;

    return x ;
} 


int sgnrand(int xi= 0, int xf=1)
{
    int s = randonilo(xi,xf,0)*2 - 1;


    return s;
}

Vector<double> gauss_vector(void)
{
    Vector<double> v(3);

    double min = 0.0001;
    double xx = randoneilo(min,1), yy = randoneilo(min,1), zz = randoneilo(min,1), ww = randoneilo(min,1) ;
    v(0) = cos(2*M_PI*yy)*sqrt(-2*log(xx)); v(1) = sin(2*M_PI*yy)*sqrt(-2*log(xx)); v(2) = cos(2*M_PI*zz)*sqrt(-2*log(ww));

    return  v ;
}

Vector< complex<double> > complex_gauss_vector(void)
{
    Vector< complex<double> > v(3);

    double min = 0.0001;
    double xx = randoneilo(min,1), yy = randoneilo(min,1) ;
    double zz = randoneilo(min,1), ww = randoneilo(min,1) ;
    double uu = randoneilo(min,1), vv = randoneilo(min,1) ;

    complex<double> cx(cos(2*M_PI*yy)*sqrt(-2*log(xx)),sin(2*M_PI*yy)*sqrt(-2*log(xx)));
    complex<double> cy(cos(2*M_PI*zz)*sqrt(-2*log(ww)),sin(2*M_PI*zz)*sqrt(-2*log(ww)));
    complex<double> cz(cos(2*M_PI*uu)*sqrt(-2*log(vv)),sin(2*M_PI*uu)*sqrt(-2*log(vv)));

    v(0) = cx; v(1) = cy; v(2) = cz;

    return  v;
}

Vector<double> random_unitario(void)
{
    Vector<double> v(3);

    double min = 0.0001;
    double theta = randoneilo(min,M_PI), phi = randoneilo(min,2*M_PI) ;
    v(0) = sin(theta)*cos(phi) ; v(1) = sin(theta)*cos(phi) ; v(2) = cos(theta) ;

    return  v ;
}


/*
Vector<double> gauss_vector(int preci = 4)
{
    Vector<double> v(3);
    double min = pow(10,-preci);

    double xx = randoneilo(min,1,preci), yy = randoneilo(min,1,preci), zz = randoneilo(min,1,preci), ww = randoneilo(min,1,preci) ;
    v(0) = cos(2*M_PI*yy)*sqrt(-2*log(xx)); v(1) = sin(2*M_PI*yy)*sqrt(-2*log(xx)); v(2) = cos(2*M_PI*zz)*sqrt(-2*log(ww));

    return  v ;
}
*/

