#ifndef CLASE_VECTOR
#define CLASE_VECTOR

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

using namespace std;

template<class T>
class Vector
{
private:
	int largo;
	vector<T> vec;
public:
	Vector(); //-----------------------------------------------------> constructor por defecto
	Vector(int); //--------------------------------------------------> constructor por argumento de largo
	Vector(const Vector<T> &);//-------------------------------------> constructor de copia blanda
	Vector<T> & operator = (const Vector<T> &); //-------------------> constructor de copia dura
	~Vector(); //----------------------------------------------------> destructor
	T & operator () (int); //----------------------------------------> operador de seleccion (se toca)
	T const & operator () (int) const; //----------------------------------> operador de seleccion (se mira)

	int Largo() const;

	void push_back(T=T(0)) const;
};

// FUNCIONES EXTERNAS //
template<class T>
Vector<T> operator + (const Vector<T> & v1,const Vector<T> & v2);

template<class T>
Vector<T> operator - (const Vector<T> & v1,const Vector<T> & v2);

template<class T, class K>
T operator * (const Vector<K> & v1,const Vector<T> & v2);

template<class T, class K>
Vector<T> operator * (const K & alpha,const Vector<T> & v2);

template<class T, class K>
Vector<K> operator * (const Vector<T> & v2, const K & alpha);

template<class T>
Vector<T> operator / (const Vector<T> & v,const T & alpha);

template<class T, class K>
Vector<T> operator % (const Vector<K> & v1,const Vector<T> & v2);

template<class T>
T mod(const Vector<T> & v);

template<class T>
ostream & operator << (ostream & out, const Vector<T> & v);



//DERIVADAS (nuevo)
template<class T>
T df_centrada(T (*f)(T), T, T, T);




// IMPLEMENTACION //
template<class T>
Vector<T> :: Vector()
{
	largo = 0;
	vec.resize(largo);
}


template<class T>
Vector<T> :: Vector(int large)
{
	largo = large;
	vec.resize(largo);
}


template<class T>
Vector<T> :: Vector(const Vector<T> & v)
{
	largo = v.largo;
	vec = v.vec;
}


template<class T>
Vector<T> & Vector<T> :: operator = (const Vector<T> & v)
{
	largo = v.largo;
	vec = v.vec;

	return *this;
}


template<class T>
Vector<T> :: ~Vector()
{
	//OH SENIOOR
}


template<class T>
T & Vector<T> ::  operator () (int pos)
{
	return vec[pos];
}


template<class T>
T const & Vector<T> ::  operator () (int pos) const
{
	return vec[pos];
}


template<class T>
int Vector<T> :: Largo() const
{
	return largo;
}

template<class T>
void Vector<T>::push_back(T a) const
{
	vec.push_back(a);
}



template<class T>
Vector<T> operator + (const Vector<T> & v1,const Vector<T> & v2)
{
	int largo1 = v1.Largo() ;
	int largo2 = v2.Largo() ;

	Vector<T> SUMA(largo1);

	if(largo1 != largo2)
	{
		cout << "no es posible sumar vectores de distinto rango" << endl;
	}

	else
	{
		for(int i=0;i<largo1;i++)
			SUMA(i) = v1(i) + v2(i);
	}

	return SUMA;
}



template<class T>
Vector<T> operator - (const Vector<T> & v1,const Vector<T> & v2)
{
	int largo1 = v1.Largo() ;
	int largo2 = v2.Largo() ;

	Vector<T> RESTA(largo1);

	if(largo1 != largo2)
	{
		cout << "no es posible restar vectores de distinto rango" << endl;
	}

	else
	{
		for(int i=0;i<largo1;i++)
			RESTA(i) = v1(i) - v2(i);
	}

	return RESTA;
}


template<class T, class K>
T operator * (const Vector<K> & v1,const Vector<T> & v2)
{
	int largo1 = v1.Largo() ;
	int largo2 = v2.Largo() ;

	T PUNTO = T(0);

	for(int i=0;i<largo1;i++)
		PUNTO += v1(i) * v2(i);

	return PUNTO;
}


template<class T, class K>
Vector<T> operator * (const K & alpha,const Vector<T> & v)
{
	int largo = v.Largo() ;

	Vector<T> POND(largo);

	for(int i=0;i<largo;i++)
		POND(i) = alpha * v(i);

	return POND;
}

template<class T, class K>
Vector<K> operator * (const Vector<T> & v, const K & alpha)
{
	int largo = v.Largo() ;

	Vector<T> POND(largo);

	for(int i=0;i<largo;i++)
		POND(i) = alpha * v(i);

	return POND;
}


template<class T>
Vector<T> operator / (const Vector<T> & v,const T & alpha)
{
	int largo = v.Largo() ;

	Vector<T> DIV(largo);

	for(int i=0;i<largo;i++)
		DIV(i) = v(i)/alpha;

	return DIV;
}


template<class T, class K>
Vector<T> operator % (const Vector<K> & v1,const Vector<T> & v2)
{
	int largo1 = v1.Largo() ;
	int largo2 = v2.Largo() ;

	Vector<T> CRUZ(largo1);

	for(int i=0;i<largo1;i++)
		CRUZ(i) = v1((i+1)%3)*v2((i+2)%3) - v1((i+2)%3)*v2((i+1)%3);

	return CRUZ;
}


template<class T>
T mod(const Vector<T> & v)
{
	T mod = pow(v*v,.5);
	return mod;
}




template<class T>
ostream & operator << (ostream & out,const Vector<T> & v)
{
	int large = v.Largo();
	//out << "(";
	for(int i=0;i<large;i++)
	{
		(i != large-1) ? out << v(i) << "\t" : out << v(i) ;
	}

	return out;
}


//DERIVADAS

template<class T>
T df_centrada(T (*f)(T), T x_0, T dx)
{
	T df = (f(x_0+dx)-f(x_0-dx))/(2*dx);
  	return df;
}



Vector< complex<double> > conj(const Vector< complex<double> > & V)
{
	int len = V.Largo();
	Vector< complex<double> > C(len) ;
	for(int i=0;i<len;i++)
		C(i) = conj(V(i));

	return C ;
}

Vector<double> real(const Vector< complex<double> > & V)
{
	int len = V.Largo();
	Vector<double> C(len) ;
	for(int i=0;i<len;i++)
		C(i) = real(V(i));

	return C ;
}

Vector<double> imag(const Vector< complex<double> > & V)
{
	int len = V.Largo();
	Vector<double> C(len) ;
	for(int i=0;i<len;i++)
		C(i) = imag(V(i));

	return C ;
}


#endif
