#ifndef MATRIZ_H  //Para poder implementarlo en otro codigo
#define MATRIZ_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "vectores.h"

using namespace std;

//Creamos la clase template//
template <class T>
class Matriz
{
private:
	//definimos las variables intrinsecas de una matriz//
	int filas,columnas;
	vector< vector<T> > matrix;
public:
	Matriz();  //---------------------------------> constructor por defecto
	Matriz(string); //----------------------------> constructor por archivo
	Matriz(int,int);  //--------------------------> constructor por filas y columnas
	Matriz(const Matriz<T> &); //-----------------> constructor de copia blanda
  	~Matriz();  //--------------------------------> destructor

  	T & operator () (int,int); //-----------------> operador de seleccion
  	T const & operator () (int,int) const; //----->
	Matriz<T> operator = (const Matriz<T> &); //--> constructor de copia dura

// Funciones internas //

  	int Filas() const;
  	int Columnas() const;
  	void cambio_f(int,int);
  	void cambio_c(int,int);
  	void operacion_fila(T,T,int,int);
  	void agregar_fila(int);
  	void agregar_columna(int);
  	void ponderar_fila(int,T);
  	void ponderar_columna(int,T);
  	T traza();
};


//declaracion de operadores//

template <class T>
Matriz<T> operator + (const Matriz<T> &, const Matriz<T> &);

template <class T>
Matriz<T> operator - (const Matriz<T> &, const Matriz<T> &);

template <class T>
Matriz<T> operator * (const Matriz<T> &, const Matriz<T> &);

template <class T, class K>
Matriz<T> operator * (K, const Matriz<T> &);

template <class T, class K>
Vector<T> operator * (const Matriz<T> &, const Vector<K> &);

template <class T>
ostream &operator << (ostream &, const Matriz<T> );

// Funciones externas //



template<class T>
Matriz<T> Adjuntos(const Matriz<T> &);

template <class T>
Matriz<T> Inversa(const Matriz<T> &);

template <class T>
Matriz<T> Traspuesta(const Matriz<T> &);

template<class T>
Matriz<T> Adjunta(const Matriz<T> &, int, int);

template <class T, class K>
Matriz<T> Diad(const Vector<T> &, const Vector<K> &);

template<class T>
Matriz<T> rCruz(const Matriz<T> &);

template<class T>
Matriz<T> Cruz(const Matriz<T> &);

template<class T>
T max(const Matriz<T> &);



//--------------- implementacion------------------//

//Una matriz por defecto//
template <class T>
Matriz<T> :: Matriz()
{
	filas = 0;
	columnas = 0;
	matrix.resize(filas);
	matrix[0].resize(columnas);
}




//Una matriz por archivo//
template <class T>
Matriz<T> :: Matriz(string archivo)
{
	string archivo1=archivo;
	ifstream entrada(archivo1.c_str());
	if(!entrada.is_open())
	{
		cout<< "Su archivo " << archivo1 << " no existe." << endl;
	}

	entrada >> filas >> columnas;

	matrix.resize(filas);
	for(int k=0;k<filas;k++)
	{
		matrix[k].resize(columnas);
	}

	for(int i=0;i<filas;i++)
	{
		for (int j=0;j<columnas;j++)
		{
			entrada >> matrix[i][j];
		}
	}
	entrada.close();
}




//Una matriz dadas sus filas y columnas//
template <class T>
Matriz<T> :: Matriz(int f,int c)
{
	filas = f;
	columnas = c;
	matrix.resize(filas);
	for(int i=0;i<filas;i++)
	{
		matrix[i].resize(columnas);
	}
}




// Destructor //
template <class T>
Matriz<T> :: ~Matriz()
{
	// Oh, señiooor //
}



// Constructor de copia blando
template <class T>
Matriz<T> :: Matriz(const Matriz<T> &A)
{
	filas = A.filas;
	columnas = A.columnas;
	matrix = A.matrix;
}




// Filas //
template <class T>
int Matriz<T> :: Filas() const
{
	return filas;
}



// Columnas //
template <class T>
int Matriz<T> :: Columnas() const
{
	return columnas;
}



// Operador de algo para modificar //
template <class T>
T &Matriz<T> :: operator () (int i,int j)
{
	return matrix[i][j];
}



// Operador para asignar //
template <class T>
T const &Matriz<T> :: operator () (int i,int j) const
{
	return matrix[i][j];
}



// Constructor de copia duro //
template <class T>
Matriz<T> Matriz<T> :: operator = (const Matriz<T> &A)
{
	filas = A.filas;
	columnas = A.columnas;
	matrix = A.matrix;
	return *this;
}



// Funcion para cambiar filas //
template <class T>
void Matriz<T> :: cambio_f(int i1,int i2)
{

	T auxiliar; //cero para templates no es la mejor opcion
	for(int j=0;j<columnas;j++) //esto deberia recorrer las columnas
	{ //los indices estan mal, verifica cada uno de ellos.
		auxiliar = matrix[i1][j];
		matrix[i1][j] = matrix[i2][j];
		matrix[i2][j] = auxiliar;
	}
}



// Funcion para cambiar columnas //
template <class T>
void Matriz<T> :: cambio_c(int j1,int j2)
{
	T auxiliar;
	for(int i=0;i<filas;i++)
	{
		auxiliar = matrix[i][j1];
		matrix[i][j1] = matrix[i][j2];
		matrix[i][j2] = auxiliar;
	}
}



// Funcion para operar filas //
template <class T>
void Matriz<T> :: operacion_fila(T alpha,T beta,int i1,int i2)
{
	for(int j=0;j<columnas;j++)
	{
		matrix[i2][j] = beta*matrix[i2][j] ;
		matrix[i1][j] = alpha*matrix[i1][j] ;
		matrix[i2][j] = matrix[i2][j]-matrix[i1][j] ;
	}
}


template <class T>
void Matriz<T> :: ponderar_fila(int i,T alpha)
{
	for (int j=0; j<columnas; j++) //recorremos las columnas multiplicando por r
	{
		matrix[i][j]=alpha*matrix[i][j];
	}
}


template <class T>
void Matriz<T> :: ponderar_columna(int i,T alpha)
{
	for (int j=0; j<filas; j++)
	{
		matrix[j][i]=alpha*matrix[j][i];
	}
}


// Funcion que agrega fila //
template<class T>
void Matriz<T> :: agregar_fila(int i) // esto podria tener argumento por defecto (i=1), pero mi a mi compilador le molesta
{
	filas = filas + i;
	matrix.resize(filas);
	for (int i=0;i<filas;i++)
	{
		matrix[i].resize(columnas);
	}
}


// Funcion para agregar columna //
template<class T>
void Matriz<T> :: agregar_columna(int i)
{
	columnas = columnas + i;
	matrix.resize(columnas);
	for (int i=0;i<filas;i++)
	{
		matrix[i].resize(columnas);
	}
}


// TRAZA //
template<class T>
T Matriz<T> :: traza()
{
	T traza;
	if(filas != columnas)
	{
		cout<< "Su matriz no es cuadrada "<<endl;
	}
	else
	{
		for(int i=0;i<filas;i++)
		{
			traza += matrix[i][i];
		}
	}
	return traza;
}



// Implementacion de operadores basicos //

// SUMA //
template <class T>
Matriz<T> operator + (const Matriz<T> &A, const Matriz<T> &B)
{
	int Af = A.Filas();
	int Ac = A.Columnas();
	int Bf = B.Filas();
	int Bc = B.Columnas();
	Matriz<T> SUMA(Af,Bc);

	if(!(Af == Bf) || !(Ac == Bc))
		cout << "Sus matrices no son de las mismas dimensiones" << endl;
	else
		for (int i=0;i<Af; i++)
			for (int j=0;j<Bc;j++)
				SUMA(i,j) = A(i,j) + B(i,j);
	return SUMA;
}


// RESTA //
template <class T>
Matriz<T> operator - (const Matriz<T> &A, const Matriz<T> &B)
{
	int Af = A.Filas();
	int Ac = A.Columnas();
	int Bf = B.Filas();
	int Bc = B.Columnas();
	Matriz<T> RESTA(Af,Bc);

	if(!(Af == Bf) || !(Ac == Bc))
		cout << "Sus matrices no son de las mismas dimensiones" << endl;
	else
		for (int i=0;i<Af; i++)
			for(int j=0;j<Bc;j++)
				RESTA(i,j) = A(i,j) - B(i,j);

	return RESTA;
}


// MULTIPLICACION //
template <class T>
Matriz<T> operator * (const Matriz<T> &A, const Matriz<T> &B)
{
	int Af = A.Filas();
	int Ac = A.Columnas();
	int Bf = B.Filas();
	int Bc = B.Columnas();
	Matriz<T> MULT(Af,Bc);

	if(!(Ac == Bf))
		cout << "Sus matrices no tienen dimensiones donde este definida la multiplicacion" << endl;
	else
		for (int i=0;i<Af; i++)
			for(int j=0;j<Bc;j++)
				for(int k=0;k<Bc;k++)
					MULT(i,j) += A(i,k)*B(k,j);




	return MULT;
}


// PONDERACION POR ESCALAR //
template <class T, class K>
Matriz<T> operator * (K alpha, const Matriz<T> &A)
{
	int Af = A.Filas();
	int Ac = A.Columnas();

	Matriz<T> PONDERACION(Af,Ac);
	for (int i=0;i<Af; i++)
		for(int j=0;j<Ac;j++)
			PONDERACION(i,j) = alpha*A(i,j);

	return PONDERACION;
}



// MATRIZ POR VECTOR //

template <class T, class K>
Vector<T> operator * (const Matriz<T> &A, const Vector<K> &v)
{
	int len = v.Largo();

	Vector<T> POND(len);
	for (int i=0;i < A.Filas(); i++)
	{
		for(int j=0;j < A.Columnas();j++)
			POND(i) += A(i,j)*v(j) ;
	}
	return POND;
}




// MODO GRAFICO //
template <class T>
ostream &operator << (ostream &out, const Matriz<T> A)
{
	int Af = A.Filas();
	int Ac = A.Columnas();

	for(int i=0;i<Af;i++)
	{
		for(int j=0;j<Ac;j++)
			out << A(i,j) << "\t";
		out << endl;
	}
	return out;
}


// MATRIZ DE ADJUNTOS //
template<class T>
Matriz<T> Adjuntos(const Matriz<T> &M)
{
	int f = M.Filas();
	int c = M.Columnas();
	Matriz<T> ADJ(f,c) ;
	for (int i=0; i<f ; i++)
		for (int j=0; j<c ; j++)
			ADJ(i,j) = Det_Recursivo(Adjunta(M,i,j))*pow(-1,i+j);
	return ADJ;

}




// MATRIZ TRASPUESTA //
template<class T>
Matriz<T> Traspuesta(const Matriz<T> &A)
{
	int f = A.Filas();
	int c = A.Columnas();
	Matriz<T> TRASP(c,f);
	for(int i=0; i<c; i++)
       	for(int j=0; j<f; j++)
           	TRASP(i,j) = A(j,i);    //cambiamos filas por columnas

    return TRASP;
}



// MATRIZ ADJUNTA //
template<class T>
Matriz<T> Adjunta(const Matriz<T> &M, int fil, int col)
{
  int f = M.Filas();
  int c = M.Columnas();
  Matriz<T> ADJ(f-1,c-1);
  for (int i=0; i<f; i++)
    for (int j=0; j<c; j++)
    	if (i != fil && j != col)
	    {
			if (i<fil && j<col)
		 		ADJ(i,j) = M(i,j);
			else if (i>fil && j<col)
		  		ADJ(i-1,j) = M(i,j);
			else if (i<fil && j>col)
		  		ADJ(i,j-1) = M(i,j);
			else
		  		ADJ(i-1,j-1) = M(i,j);
	    }

  return ADJ;
}


// DETERMINANTE DE METODO RECURSIVO //
template<class T>
T Det_Recursivo(const Matriz<T> &M, bool fil = true , int n = 0)
{
	T det = T(0);

	if (M.Filas() == 2)
	{
		det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
	}
	else
	{
		if (fil == true)
			for (int j = 0; j < M.Columnas(); j++) // if (M(n,j) != 0)
				{
					det += M(n,j)*Det_Recursivo(Adjunta(M,n,j))*pow(-1,j+n);
					//cout << Adjunta(M,n,j) << 9 << endl;
				}
		else
			for (int i = 0; i < M.Filas(); i++)
				if ( M(i,n) != T(0) )
					det += M(i,n)*Det_Recursivo(Adjunta(M,i,n))*pow(-1,i+n);
	}

	return det ;
}



// MATRIZ INVERSA //
template<class T>
Matriz<T> Inversa(const Matriz<T> &A)
{
	int f = A.Filas();
	int c = A.Columnas();
	T detA = Det_Recursivo(A);
	Matriz<T> INV(f,c);


	INV = (1./detA)*Traspuesta(Adjuntos(A));

	return INV ;
}

template <class T, class K>
Matriz<T> Diad(const Vector<T> &A, const Vector<K> &B)
{
	int lenA = A.Largo();
	int lenB = B.Largo();
	Matriz<T> DIAD(lenA,lenB);

	for(int i=0;i<lenA;i++)
		for(int j=0;j<lenB;j++)
			DIAD(i,j) = A(i)*B(j);

	return DIAD ;
}

// GENERA LA MATRIZ EQUIVALENTE A UN PRODUCTO CRUZ
template<class T>
Matriz<T> rCruz(const Vector<T> &V)
{
	Matriz<T> CRUZ(3,3) ;

	CRUZ(0,1) = V(2) ; CRUZ(0,2) = -V(1)  ; CRUZ(1,2) = V(0) ;
	CRUZ(1,0) = -V(2) ; CRUZ(2,0) = V(1) ; CRUZ(2,1) = -V(0) ;

	return CRUZ ;
}

template<class T>
Matriz<T> Cruz(const Vector<T> &V)
{
	return -1.*rCruz(V) ;
}


/*---------- Funciones en estado beta -----------*/


template<class T>
T max(const Matriz<T> &M)
{
	T MAX = T(0);
	int N = M.Filas();
	Matriz<T> A(N,N);
	A = M;

	for(int i=0;i<N;i++)
    	for(int j=0;j<N;j++)
    	{
    		if(j < N-1)
    			if(A(i,j) > A(i,j+1))
    				A(i,j+1) = A(i,j);
    	}

    for(int i=0;i<N;i++)
    {
    	if(i < N-1)
    		if(A(i,N-1) > A(i+1,N-1))
    			A(i+1,N-1) = A(i,N-1);
    }

    MAX = A(N-1,N-1);

	return MAX;
}







#endif
