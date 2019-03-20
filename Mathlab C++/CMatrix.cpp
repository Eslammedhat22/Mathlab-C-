/////////////////////////////////////////////////* Functions Definition *///////////////////////////////////////////
#include "CMatrix.h"
#include <iostream>
#include <exception>
#include <string>
#include <string.h>
#include <stdarg.h>
#include <algorithm>
#include <cstdio>
#include <math.h>
////////////////////////////////////////////////* default Constructor *///////////////////////////////////
CMatrix::CMatrix()
{
	nR = 0;
	nC = 0;
	values = NULL;
}
///////////////////////////////////////////////////* Constructor */////////////////////////////////////////////////
/*CMatrix::CMatrix(int nR, int nC)
{
if (nR <= 0 || nC <= 0){ throw(" error:you try to define a matrix have a negtive dimensions or zero dimensions "); }
else
{
this->nR = nR;
this->nC = nC;
values = new double*[nR];
for (int row = 0; row < nR; row++)
{
values[row] = new double[nC];
for (int col = 0; col < nC; col++)
values[row][col] = 0;
}
}
}*/
//////////////////////////////////////////////////* Copy Constructor *//////////////////////////////////////////
CMatrix::CMatrix(const CMatrix &m)
{
	nR = 0;
	nC = 0;
	values = NULL;
	copy(m);
}

//////////////////////////////////////////////////////* Destructor *///////////////////////////////////////////
CMatrix::~CMatrix()
{
	reset();
}
////////////////////////////////////////////////////////* Reset *////////////////////////////////////////////
void CMatrix::reset()
{
	if (values)
	{
		for (int i = 0; i<nR; i++)
			delete[] values[i];
		delete[] values;
	}
	nR = nC = 0;
	values = NULL;
}
/////////////////////////////////////* to get an element with given coordinates from a matrix *////////////////////////////////
double CMatrix::getElement(int r, int c)
{
	if (r >= nR || r < 0 || c >= nC || c <0){ throw(" error:you try to access element out of matrix by getElement() function "); }
	else return values[r][c];
}
////////////////////////////////////* to set an element with given coordinates & value to a matrix *///////////////////////////////

void CMatrix::setElement(int r, int c, double v)
{
	if (r >= nR || r < 0 || c >= nC || c <0){ throw(" error:you try to set element out of matrix by setElement() function "); }
	else values[r][c] = v;
}
////////////////////////////////////////////////////////* Copy function *////////////////////////////////////////////
void CMatrix::copy(const CMatrix& m)
{
	reset();
	this->nR = m.nR;
	this->nC = m.nC;
	if ((nR*nC) == 0)
	{
		values = NULL; return;
	}
	values = new double*[nR];
	for (int iR = 0; iR<nR; iR++)
	{
		values[iR] = new double[nC];
		for (int iC = 0; iC<nC; iC++)
		{
			values[iR][iC] = m.values[iR][iC];
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nR, int nC, int initialization, double initializationValue)
{
	this->nR = nR;
	this->nC = nC;
	if ((nR*nC) == 0)
	{
		values = NULL;
		return;
	}
	values = new double*[nR];
	for (int iR = 0; iR<nR; iR++)
	{
		values[iR] = new double[nC];
		for (int iC = 0; iC<nC; iC++)
		{
			switch (initialization)
			{
			case MI_ZEROS: values[iR][iC] = 0; break;
			case MI_ONES: values[iR][iC] = 1; break;
			case MI_EYE: values[iR][iC] = (iR == iC) ? 1 : 0; break;
			case MI_RAND: values[iR][iC] = (rand() % 1000000) / 1000000.0; break;
			case MI_VALUE: values[iR][iC] = initializationValue; break;
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(int nR, int nC, double first, ...)
{
	this->nR = nR;
	this->nC = nC;
	if ((nR*nC) == 0){ values = NULL; return; }
	values = new double*[nR];
	va_list va;
	va_start(va, first);
	for (int iR = 0; iR<nR; iR++)
	{
		values[iR] = new double[nC];
		for (int iC = 0; iC<nC; iC++)
		{
			values[iR][iC] = (iC == 0 && iR == 0) ? first : va_arg(va, double);
		}
	}
	va_end(va);
}

/////////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(string s)
{
	nR = nC = 0;
	values = NULL;
	copy(s);
}
////////////////////////////////////////////////////////////////////////
CMatrix::CMatrix(double d)
{
	nR = nC = 1;
	values = new double*[nR];
	for (int row = 0; row < nR; row++)
	{
		values[row] = new double[nC];
		for (int col = 0; col < nC; col++)
			values[row][col] = d;
	}

}
///////////////////////////////////////////* to get the Minor of an element *///////////////////////////////////////
CMatrix CMatrix::getMinor(int r, int c)
{
	if (nR == 1 || nC == 1) throw(" error:Invalid matrix dimension(1x1) in getMinor()function ");
	if (nR == 0 && nC == 0)throw(" error:undefined matrix(0x0) in getMinor()function ");
	if ((nR == 1 && nC > 1) || (nR > 1 && nC == 1))throw("error:Can't get a minor to a vector (must be 2D dimension) ");
	if (r >= nR || r < 0 || c >= nC || c <0){ throw("error:you try to getminor of element out of matrix by getMinor() function "); }
	else
	{
		CMatrix m(nR - 1, nC - 1);
		for (int iR = 0; iR<m.nR; iR++)
		for (int iC = 0; iC<m.nC; iC++)
		{
			int sR = (iR<r) ? iR : iR + 1;
			int sC = (iC<c) ? iC : iC + 1;
			m.values[iR][iC] = values[sR][sC];
		}
		return m;
	}
}
/////////////////////////////////////////* to calculate the DETERMINANT of the matrix *////////////////////////////
int CMatrix::checkzero(int r, int c)
{
	int flag=0;
	if(values[r][c]==0)
	{
		for(int x=r+1;x<nR;x++)
		{
			if(values[x][c]!=0)
			{
				for(int k=0;k<nC;k++)
				{
					values[r][k]=values[r][k]+values[x][k];
					flag=1;

				}
			}
			if(flag==1){return 1 ; break;}
		}

	}
	else
	{
		return 1;
	}
	if (flag==0)return 0;
}
double CMatrix::getDeterminant()
{
	if (nR != nC)throw(" error:Invalid matrix dimension(not square matrix) in getDeterminant() function ");
	if (nR == 0 && nC == 0)throw(" error:undefined matrix(0x0) in getDeterminant()nfunction ");
	else
	{
		CMatrix m=*this;
		for(int i=0;i<(m.nR-1);i++)
		{
			if(m.values[i][i]==0)
			{
			   	int v=m.checkzero(i,i);
				if(v==0) return 0;
			}
			double k;
			for(int j=i+1;j<m.nR;j++)
			{
				k=(-1*m.values[j][i])/m.values[i][i];
				for(int x=0;x<m.nC;x++)
				{
					m.values[j][x]=m.values[j][x]+(k*m.values[i][x]);
				}
			}
		}

		/*for(int i=1;i<m.nR;i++)
		{
			for(int j=0;j<i;j++)
			{
				double k;
				if(m.values[j][j]==0)
				{
					int v=m.checkzero(j,j);
					if(v==0) return 0;
				}
				k=(-1*m.values[i][j])/m.values[j][j];
				for(int x=0;x<m.nC;x++)
				{
					m.values[i][x]=m.values[i][x]+(k*m.values[j][x]);
				}

			}
		}*/
		double det =1.0;
		for (int u=0;u<m.nR;u++)
		{
			det=det*m.values[u][u];
		}
		return det;

	}
}

////////////////////////////////////////////////* to get the COFACTOR of the matrix *////////////////////////////////////
CMatrix CMatrix::getCofactor()
{
	if (nR != nC)throw(" error:Invalid matrix dimension(not square matrix) in getCofactor() function ");
	if (nR == 1 && nC == 1)throw(" error:Invalid matrix dimension(1*1) in getCofactor()function ");
	if (nR == 0 && nC == 0)throw(" error:undefined matrix(0x0) in getCofactor() function ");
	else
	{
		CMatrix m(nR, nC);
		int k = 1;
		int v;
		for (int row = 0; row < nR; row++)
		{
			for (int col = 0; col < nC; col++)
			{
				if (col == 0) v = k;
				m.values[row][col] = v *getMinor(row, col).getDeterminant();
				v *= -1;
			}
			k *= -1;
		}
		return m;
	}
}

/////////////////////////////////////////////* to get the total number of ELEMENTS inside the matrix *///////////////////////////
int CMatrix::getn()
{
	return nR*nC;
}
///////////////////////////////////////////* to get the number of ROWS inside the matrix *///////////////////////////////////////
int CMatrix::getnR()
{
	return nR;
}
////////////////////////////////////////////* to get the number of COLOUMS inside the matrix *////////////////////////////////
int CMatrix::getnC()
{
	return nC;
}
////////////////////////////////////////////* to get the TRANSPOSE of the matrix */////////////////////////////////////
CMatrix CMatrix::getTranspose()
{
	if (nR == 0 && nC == 0)throw(" error:undefined matrix(0x0) in getTranspose()function ");
	else
	{
		CMatrix m(nC, nR);
		for (int i = 0; i<nR; i++)
		{
			for (int j = 0; j<nC; j++)
			{
				m.values[j][i] = values[i][j];
			}
		}
		return m;
	}
}
/////////////////////////////////////////////////* to get the INVERSE of the matrix *//////////////////////////////
CMatrix CMatrix::getInverse()
{
	if (nR != nC) throw(" error:Invalid matrix dimension(not square matrix) in getInverse() function ");
	if (nR == 0 && nC == 0)throw(" error:undefined matrix(0x0) in getInverse() function ");
	double det = getDeterminant();
	CMatrix m(nR, nC);
	if (det == 0){ throw(" error:can't get inverse to a singluar matrix "); }
	else
	{
		m = getCofactor();
		m = m.getTranspose();
		for (int i = 0; i < nR; i++)
		{
			for (int j = 0; j < nC; j++)
			{
				m.values[i][j] /= det;
			}
		}
		return m;
	}

}
////////////////////////////////////////////////////* Equal ("=") operator *////////////////////////////////////
CMatrix CMatrix::operator=(const CMatrix& m)
{
	copy(m);
	return *this;
}
////////////////////////////////////////////* []operator to get an element with a given index *////////////////
double CMatrix:: operator [] (int index)
{
	if (nR == 0 && nC == 0)throw(" error: you try to access element of undefined matrix(0x0) by [] operator ");
	if (index >= nR*nC || index < 0){ throw(" error:you try to access element out of matrix by [] operator"); }
	else
	{

		int r, c;
		r = index / nC;
		c = index % nC;
		return values[r][c];
	}
}
/////////////////////////////////////////* ()operator to get an element with a given Coordinates *////////////////////
double CMatrix:: operator () (int r, int c)
{
	if (nR == 0 && nC == 0)throw(" error: you try to access element of undefined matrix(0x0) by () operator ");
	if (r >= nR || r < 0 || c >= nC || c <0){ throw("error:you try to access element out of matrix by() operator"); }
	else
	{

		return values[r][c];
	}
}
///////////////////////////////////////* ()operator to get an element with a given index */////////////////////////
double CMatrix:: operator () (int index)
{
	if (index >= nR*nC || index < 0){ throw("error:you try to access element out of matrix by () operator"); }
	else
	{

		int r, c;
		r = index / nC;
		c = index % nC;
		return values[r][c];
	}
}
///////////////////////////////////////////////////////////////////////////////////////

istream& operator >> (istream &is, CMatrix& m)
{
	string s;
	getline(is, s, ']');
	s += "]";
	m = CMatrix(s);
	return is;
}
//////////////////////////////////////////////////////////////////////////////////////
ostream& operator << (ostream &os, CMatrix& m)
{
	os << m.getString();
	return os;
}

/////////////////////////////////////////////////////////////////////


void CMatrix::add(const CMatrix& m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to add undefined matrix(0x0) by add() function ");
	if ((nR != m.nR) || (nC != m.nC))throw ("error:you try to add two matrices have different size by add()function");
	else
	{
		for (int iR = 0; iR < nR; iR++)
		{
			for (int iC = 0; iC < nC; iC++)
			{
				values[iR][iC] += m.values[iR][iC];
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
void CMatrix::add(double c)
{
	if (nR == 0 && nC == 0)throw(" error: you try to add undefined matrix(0x0) by add() function ");

	else
	{
		for (int iR = 0; iR < nR; iR++)
		{
			for (int iC = 0; iC < nC; iC++)
			{
				values[iR][iC] += c;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////
void CMatrix::operator+=(const CMatrix& m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to add undefined matrix(0x0) by  += operator ");
	if ((nR != m.nR) || (nC != m.nC)) throw ("error:you try to add two matrices have different size by += operator");
	else add(m);
}
/////////////////////////////////////////////////////////////////////
void CMatrix::operator+=(double c)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to add undefined matrix(0x0) by  += operator ");
	else add(c);

}
/////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator+(const CMatrix &m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to add undefined matrix(0x0) by  + operator ");
	if ((nR != m.nR) || (nC != m.nC))throw ("error:you try to add two matrices have different size by + operator");
	else
	{
		CMatrix s = *this;
		s += m;
		return s;
	}
}
////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator+(double c)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to add undefined matrix(0x0) by + operator ");
	else
	{
		CMatrix s = *this;
		s += c;
		return s;
	}
}
///////////////////////////////////////////////////////////////////
void CMatrix::sub(const CMatrix& m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to sub undefined matrix(0x0) by sub() function ");
	if ((nR != m.nR) || (nC != m.nC))throw ("error:you try to sub two matrices have different size by sub()function");

	for (int iR = 0; iR < nR; iR++)
	{
		for (int iC = 0; iC < nC; iC++)
			values[iR][iC] -= m.values[iR][iC];
	}
}
//////////////////////////////////////////////////////////////////
void CMatrix::sub(double c)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to sub undefined matrix(0x0) by sub() function ");


	for (int iR = 0; iR < nR; iR++)
	{
		for (int iC = 0; iC < nC; iC++)
			values[iR][iC] -= c;
	}
}
//////////////////////////////////////////////////////////////////
void CMatrix::operator-=(const CMatrix& m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to sub undefined matrix(0x0) by  -= operator ");
	if ((nR != m.nR) || (nC != m.nC)) throw ("error:you try to sub two matrices have different size by -= operator");
	sub(m);
}
//////////////////////////////////////////////////////////////////
void CMatrix::operator-=(double c)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to sub undefined matrix(0x0) by  -= operator ");
	else sub(c);
}
/////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator-(CMatrix m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to sub undefined matrix(0x0) by  - operator ");
	if ((nR != m.nR) || (nC != m.nC))throw ("error:you try to sub two matrices have different size by - operator");
	else
	{
		CMatrix s = *this;
		s -= m;
		return s;
	}
}
////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator-(double c)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to sub undefined matrix(0x0) by - operator ");
	else
	{
		CMatrix s = *this;
		s -= c;
		return s;
	}

}
////////////////////////////////////////////////////////////////
void CMatrix::mul(const CMatrix &m)
{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to multiply undefined matrix(0x0) by mul() function ");
	if ((nC != m.nR))throw ("error:you try to multiply two matrices have invalid dimension by mul() function");
	else
	{
		CMatrix Result(nR, m.nC);
		for (int iR = 0; iR<Result.nR; iR++)

		for (int iC = 0; iC<Result.nC; iC++)

		{

			Result.values[iR][iC] = 0;

			for (int k = 0; k<m.nC; k++)

				Result.values[iR][iC] += values[iR][k] * m.values[k][iC];

		}
		copy(Result);
	}
}
//////////////////////////////////////////////////////////////////////////
void CMatrix::mul(double c)
{
	if (nR == 0 && nC == 0)throw(" error: you try to multiply undefined matrix(0x0) by mul() function ");

	else
	{
		for (int iR = 0; iR < nR; iR++)
		{
			for (int iC = 0; iC < nC; iC++)
			{
				values[iR][iC] *= c;
			}
		}
	}

}
//////////////////////////////////////////////////////////////////////////

void CMatrix::operator*=(const CMatrix& m)

{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to multiply undefined matrix(0x0) by  *= operator ");
	if ((nC != m.nR))throw ("error:you try to multiply two matrices have invalid dimension by *= operator");
	else mul(m);
}
/////////////////////////////////////////////////////////////////////////
void CMatrix::operator*=(double d)
{
	if ((nR == 0 && nC == 0))throw(" error: you try to multiply undefined matrix(0x0) by  *= operator ");
	else
	{
		mul(d);
	}


}
/////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator*(const CMatrix& m)

{
	if ((nR == 0 && nC == 0) || (m.nR == 0 && m.nC == 0))throw(" error: you try to multiply undefined matrix(0x0) by * operator ");
	if ((nC != m.nR))throw ("error:you try to multiply two matrices have invalid dimension by * operator");
	else
	{
		CMatrix r = *this;

		r *= m;
		return r;
	}

}
////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator*(double d)

{
	if ((nR == 0 && nC == 0))throw(" error: you try to multiply undefined matrix(0x0) by  * operator ");
	else
	{
		CMatrix r = *this;

		r *= d;

		return r;

	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////

void CMatrix::addColumn(CMatrix& m)
{
	CMatrix n(max(nR, m.nR), nC + m.nC);
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(0, nC, m);
	*this = n;
}
/////////////////////////////////////////////////////////////////////////////
void CMatrix::addRow(CMatrix& m)
{
	CMatrix n(nR + m.nR, max(nC, m.nC));
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(nR, 0, m);
	*this = n;
}
/////////////////////////////////////////////////////////////////////////////
void CMatrix::copy(double d)
{
	for (int i = 0; i < nR; i++)
	for (int j = 0; j < nC; j++)
		values[i][j] = d;
}
////////////////////////////////////////////////////////////////////////////////

void CMatrix::copy(string s)
{
	reset();

	char* buffer = new char[s.length() + 1];
	strcpy(buffer,s.c_str());
	char* lineContext;
	char* lineSeparators = ";\r\n";
	char* line = strtok_r(buffer, lineSeparators, &lineContext);
	char* context;
	char* separators = " []";
	while (line)
	{
		CMatrix row;
		char* token = strtok_r(line, separators, &context);
		while (token)
		{
			CMatrix item = atof(token);
			row.addColumn(item);
			token = strtok_r(NULL, separators, &context);
		}
		if (row.nC>0 && (row.nC == nC || nR == 0)) addRow(row);
		line = strtok_r(NULL, lineSeparators, &lineContext);
		if ((row.nC != nC)&&(row.nC!=0))    ///////// edit/////////////////
			{
				reset();
				throw("error:invalid matrix");
			}	
	}
	delete[] buffer;
}
//////////////////////////////////////////////////////////////////////////////////
string CMatrix::getString()
{
	string s;
	for (int iR = 0; iR<nR; iR++)
	{
		for (int iC = 0; iC<nC; iC++)
		{
			char buffer[50];
			snprintf(buffer,50, "%g\t", values[iR][iC]);
			s += buffer;
		} s += "\n";
	}

	if (nR>0 && nC>0 )return s;				///////// edit/////////////////
	else throw("empty matrix");
}

////////////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator=(double d)
{
	copy(d);
	return *this;
}
////////////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator=(string s)
{
	copy(s);
	return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////// basma & mohamed
CMatrix CMatrix::operator++()        //pre-increment -> ++a
{
	add(1.0);
	return *this;
}

///////////////////////////////////////////////////////////////////////////


CMatrix CMatrix::operator++(int)     //post-increment -> a++
{
	CMatrix C = *this;
	add(1.0);
	return C;
}

///////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator--()      //pre-decrement -> --a
{
	add(-1.0);
	return *this;
}

//////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator--(int)     //post-decrement -> a--
{
	CMatrix r = *this;
	add(-1.0);
	return r;
}

///////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator-()      //negative
{
	for (int iR = 0; iR<nR; iR++)
	{
		for (int iC = 0; iC<nC; iC++)
			values[iR][iC] = -values[iR][iC];
	}
	return *this;
}

///////////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::operator+()        //positive
{
	return *this;
}

////////////////////////////////////////////////////////////////////////// setSubMatrix
void CMatrix::setSubMatrix(int r, int c, CMatrix& m)
{
	if ((r + m.nR)>nR || (c + m.nC)>nC)throw("Invalid matrix dimension");
	for (int iR = 0; iR<m.nR; iR++)
	for (int iC = 0; iC<m.nC; iC++)
		values[r + iR][c + iC] = m.values[iR][iC];
}

////////////////////////////////////////////////////////////////////  getSubMatrix
CMatrix CMatrix::getSubMatrix(int r, int c, int nr, int nc)
{
	if ((r + nr)>nR || (c + nc)>nC)throw("Invalid matrix dimension");
	CMatrix m(nr, nc);
	for (int iR = 0; iR<m.nR; iR++)
	for (int iC = 0; iC<m.nC; iC++)
		m.values[iR][iC] = values[r + iR][c + iC];
	return m;
}

/////////////////////////////////////////////////////////////////// division
void CMatrix::div(const CMatrix& m)
{
	if (m.nR != m.nC)
		throw("Invalid matrix dimension");
	if (nC != m.nR)throw("Invalid matrix dimension");
	else
	{

		CMatrix r = m;
		r = r.getInverse();
		mul(r);
	}

}

/////////////////////////////////////////////////////////////// division over matrix
void CMatrix::operator/=(const CMatrix& m)
{
	div(m);
}

///////////////////////////////////////////////////////////// division over value
void CMatrix::operator/=(double d)
{
	for (int iR = 0; iR<nR; iR++)
	for (int iC = 0; iC<nC; iC++)
		values[iR][iC] /= d;
}

//////////////////////////////////////////////////////////// division over matrix returns matrix
CMatrix CMatrix::operator/(const CMatrix& m)
{
	CMatrix r = *this;
	r /= m;
	return r;
}

//////////////////////////////////////////////////////// division over value returns matrix
CMatrix CMatrix::operator/(double d)
{
	CMatrix r = *this;
	r /= d;
	return r;

}
/////////////////////////////////////////////////////////////////////////////////
CMatrix CMatrix::elementDiv(double d)
{
	CMatrix m(nR, nC);
	for (int i = 0; i<nR; i++)
	for (int j = 0; j < nC; j++)
	{
		double v = getElement(i, j);
		m.setElement(i, j, (d / v));
	}
	return m;

}


////////////////////////////////mohammed w basma w rbna ystor ///////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////// trigonometric fns

CMatrix sin(CMatrix &m)                  ////////////////////// gets sin matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = sin(m.values[i][j]);
	 return n;
}

CMatrix asin(CMatrix &m)                  ////////////////////// gets asin matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = asin(m.values[i][j]);
	 return n;
}
CMatrix acos(CMatrix &m)                  ////////////////////// gets asin matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = acos(m.values[i][j]);
	 return n;
}
CMatrix atan(CMatrix &m)                  ////////////////////// gets asin matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = atan(m.values[i][j]);
	 return n;
}



CMatrix cos(CMatrix &m)                  ////////////////////// gets cos matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = cos(m.values[i][j]);
	 return n;
}

CMatrix tan(CMatrix &m)                  ////////////////////// gets tan matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = tan(m.values[i][j]);
	 return n;
}

CMatrix sinh(CMatrix &m)                  ////////////////////// gets sinh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = sinh(m.values[i][j]);
	 return n;
}

CMatrix cosh(CMatrix &m)                  ////////////////////// gets cosh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = cosh(m.values[i][j]);
	 return n;
}

CMatrix tanh(CMatrix &m)                  ////////////////////// gets tanh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = tanh(m.values[i][j]);
	 return n;
}
CMatrix sech(CMatrix &m)                  ////////////////////// gets sinh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = 1/cosh(m.values[i][j]);
	 return n;
}
CMatrix csch(CMatrix &m)                  ////////////////////// gets sinh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = 1/sinh(m.values[i][j]);
	 return n;
}
CMatrix coth(CMatrix &m)                  ////////////////////// gets sinh matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = 1/tanh(m.values[i][j]);
	 return n;
}

CMatrix sec(CMatrix &m)                  ////////////////////// gets sec matrix
{	
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			 if(cos(m.values[i][j])==0)
				 throw ("Marh Error .. Invalid matrix");
			 n.values[i][j] = 1/(cos(m.values[i][j]));
		 }
	 return n;
}

CMatrix csc(CMatrix &m)                  ////////////////////// gets csc matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			 if (sin(m.values[i][j]) ==0)
				 throw ("Marh Error .. Invalid matrix");
			 n.values[i][j] = 1/(sin(m.values[i][j]));
		 }
	 return n;
}

CMatrix cot(CMatrix &m)                  ////////////////////// gets cot matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			 if (tan(m.values[i][j]) ==0)
				 throw ("Marh Error .. Invalid matrix");
			 n.values[i][j] = 1/(tan(m.values[i][j]));
		 }
	 return n;
}

/////////////////////////////////////////////////////////////////////////////////////////// Roots

CMatrix sqrt(CMatrix &m)                  ////////////////////// gets sqrt matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			 if (m.values[i][j] <=0)
				 throw ("Marh Error .. Invalid matrix");
			 n.values[i][j] = sqrt(m.values[i][j]);
		 }
	 return n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////powers

CMatrix CMatrix :: powElement(double p)                      ///////////////////power each element
{
	CMatrix r = *this;
	 for (int iR=0 ; iR<nR ; iR++) 
		 for (int iC=0 ; iC<nC ; iC++)
		 {
			 r.values[iR][iC]=pow(r.values[iR][iC],p);
		 }
		 return r;
}

CMatrix CMatrix :: pow_1x1_matrix(double p)                  /////////////////// fraction power for 1x1 matrix
{
	CMatrix r = *this;
	r.values[0][0] = pow(values[0][0] , p) ;
	return r;
}

int binary(int number ,int location)
{
  if((1 &(number>>location))==0) return 0;
  else return 1;
}




CMatrix CMatrix :: pow_matrix(int p)                                    /////////////////// power for matrix
{
	CMatrix result (nR , nC ,2,0);
	CMatrix y = *this ;
	int flag=0;
	
	if (p<0)
	{
		p*= -1; flag =1 ;
	}
	double h = log10(2) ;
	int k=((log10(p))/h) +1;
	for ( int i=0 ; i<k ; i++ )          
	{
		if (binary (p,i)==1)
			result = result * y;
		y = y*y;
	}
	if (flag == 1)
		return result.getInverse();
	
	return result ; 
}



/////////////////////////////////////////////////////////////////////////////////////////logarithmic

CMatrix log(CMatrix &m)                  ////////////////////// gets log matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			  if (m.values[i][j] <= 0)
				 throw ("Math Error .. Invalid matrix ");
			 n.values[i][j] = log10 (m.values[i][j]);
		 }
	 return n;
}

CMatrix ln(CMatrix &m)                  ////////////////////// gets ln matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
		 {
			 if (m.values[i][j] <= 0)
				 throw ("Math Error .. Invalid matrix ");
			 n.values[i][j] = log(m.values[i][j]);
		 }
	 return n;
}

CMatrix exp(CMatrix &m)                  ////////////////////// gets exp matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = exp(m.values[i][j]);
	 return n;
}

//////////////////////////////////////////////////////////////////////////////////

CMatrix asec(CMatrix &m)                  ////////////////////// gets asec matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = acos(1/m.values[i][j]);
	 return n;
}

CMatrix acsc(CMatrix &m)                  ////////////////////// gets acsc matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = asin(1/m.values[i][j]);
	 return n;
}


CMatrix acot(CMatrix &m)                  ////////////////////// gets acot matrix
{
	 CMatrix n(m.nR, m.nC);
	 for (int i = 0; i < n.nR;i++)
		 for (int j = 0; j < n.nC; j++)
			 n.values[i][j] = atan(1/m.values[i][j]) ;
	 return n;
}





////////////////////////////////////medhat &amira &aya///////////////////////////////////
string CMatrix::sendString()
{
	string s = "[";
	for (int iR = 0; iR < nR; iR++)
	{
		for (int iC = 0; iC < nC; iC++)
		{
			char buffer[50];
			if (iC < nC - 1)
			{
				snprintf(buffer, 50, "%g ", values[iR][iC]);
				s += buffer;
			}
			else if (iC == nC - 1)
			{
				snprintf(buffer, 50, "%g", values[iR][iC]);
				s += buffer;
			}
		}
		if (iR != nR - 1)
		{
			s += ";";
		}
	} s += "]";
	if (nR>0 && nC>0 )return s;				
	else throw("empty matrix sendstring");             
}
//////////////////////////////////////////////////////////////////////////////////////
void CMatrix::concatinate(CMatrix& m)
{
	if (nR != m.nR)throw("Invalid matrix dimension");
	CMatrix n(nR, nC + m.nC);
	n.setSubMatrix(0, 0, *this);
	int r = 0;
	int c = nC;
	n.setSubMatrix(r, c, m);
	copy(n);
}
