/////////////////////////////////////////////////* Functions Declaration *///////////////////////////////////////////
#ifndef CMATRIX_H
#define CMATRIX_H
#include <iostream>
#include <string>
#define PI 3.142857142857143
using namespace std;
 


class CMatrix
{

	int nR, nC;
	double** values;
	int checkzero(int r,int c);
public:

	CMatrix();                                      /* default Constructor */
	/*CMatrix(int nR, int nC); */                       /* constructor */
	CMatrix(const  CMatrix &m);                            /* Copy constructor */
	~CMatrix();						                /* Destructor */
	void reset();					                /* Reset */
	void setElement(int nR, int nC, double v);      /* to set an element with given coordinates & value to a matrix */
	double getElement(int r, int c);                /* to get an element with given coordinates from a matrix */
	void copy(const CMatrix& m);                          /* Copy Function */
	CMatrix getMinor(int r, int c);                 /* to get minor of a an element */
	double getDeterminant();		                /* to calculate the determinant of the matrix */
	CMatrix  getCofactor();			                /* to get cofactor of the matrix*/
	int getn();                                     /* to get the total number of ELEMENTS inside the matrix */
	int getnR();									/* to get the number of ROWS inside the matrix */
	int getnC();									/* to get the number of COLOUMS inside the matrix */
	CMatrix getTranspose();							/* to get the TRANSPOSE of the matrix */
	CMatrix  getInverse();							/* to get the INVERSE of the matrix */
	CMatrix operator = (const CMatrix& m);                /* equal ("=") opertor */
	double  operator [] (int index);				/* []operator to get an element with a given index */
	double  operator () (int index);				/* ()operator to get an element with a given index */
	double  operator () (int r, int c);				/* ()operator to get an element with a given Coordinates */
	friend istream& operator >> (istream &is, CMatrix& C);    /* Out Stream*/
	friend ostream& operator << (ostream &os, CMatrix& C);    /* In Stream Void*/
	///////////////////////////////////////////////////////////////////////////////////////////////
	void add(const CMatrix& m);
	void add(double c);
	void operator+=(const CMatrix& m);
	void operator+=(double c);
	CMatrix operator+(const CMatrix &m);
	CMatrix operator+(double c);
	void sub(const CMatrix& m);
	void sub(double c);
	void operator-=(const CMatrix& m);
	void operator-=(double c);
	CMatrix operator-(CMatrix m);
	CMatrix operator-(double c);
	void mul(const CMatrix& m);
	void mul(double c);
	void operator*=(const CMatrix& m);
	void operator*=(double d);
	CMatrix operator*(const CMatrix& m);
	CMatrix operator*(double d);
	///////////////////////////////////////////////////////////////////////////////////////////////////
	enum MI{ MI_ZEROS, MI_ONES, MI_EYE, MI_RAND, MI_VALUE };
	CMatrix(int nR, int nC, int initialization = MI_ZEROS, double initializationValue = 0.0);
	CMatrix(int nR, int nC, double first, ...);
	CMatrix(double d);
	CMatrix(string s);
	void copy(double d);
	void copy(string s);
	string getString();
	CMatrix operator=(double d);
	CMatrix operator=(string s);
	void addColumn(CMatrix& m);
	void addRow(CMatrix& m);
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void div(const CMatrix& m);                    //division  
	void operator/=(const CMatrix& m);             // division matrix over matrix internal return
	void operator/=(double d);               //division matrix over value internal return
	CMatrix operator/(const CMatrix& m);           // division matrix over matrix returns matrix           
	CMatrix operator/(double d);             // division matrix over value returns matrix       
	CMatrix operator++();                    //Pre Increment
	CMatrix operator++(int);                 //Post Increment, int is not used 
	CMatrix operator--();                    //Pre Increment
	CMatrix operator--(int);                 //Post Increment, int is not used
	CMatrix operator-();                     //positive
	CMatrix operator+();                     //negative
	void setSubMatrix(int iR, int iC, CMatrix& m);                  /*setSubMatrix*/
	CMatrix getSubMatrix(int r, int c, int nr, int nc);                /*getSubMatrix*/
	CMatrix elementDiv(double d =1.0);                                           /*return matrix of element division*/

	///////////////////////////////mohamed w basma w rbna ystor tany o.O/////////////////////////







/*/////////////////////////////////////////////////////////////////////////////////// trigonometric fns*/

friend CMatrix sin(CMatrix &m);                ///////////////////// sine matrix
friend CMatrix cos(CMatrix &m);                ///////////////////// cosine matrix
friend CMatrix tan(CMatrix &m);                ///////////////////// tan matrix
friend CMatrix sec(CMatrix &m);                ///////////////////// sec
friend CMatrix csc(CMatrix &m);                ///////////////////// cosec
friend CMatrix cot(CMatrix &m);                ///////////////////// cot 


friend CMatrix asin(CMatrix &m);               ///////////////////// asine matrix
friend CMatrix acos(CMatrix &m);               ///////////////////// asine matrix
friend CMatrix atan(CMatrix &m);               ///////////////////// asine matrix
friend CMatrix asec(CMatrix &m);               ///////////////////// asec matrix
friend CMatrix acsc(CMatrix &m);               ///////////////////// acosec matrix
friend CMatrix acot(CMatrix &m);               ///////////////////// acot matrix

friend CMatrix sinh(CMatrix &m);               ///////////////////// sinh matrix
friend CMatrix cosh(CMatrix &m);               ///////////////////// cosh matrix
friend CMatrix tanh(CMatrix &m);               ///////////////////// tanh matrix
friend CMatrix csch(CMatrix &m);               ///////////////////// sinh matrix
friend CMatrix sech(CMatrix &m);               ///////////////////// cosh matrix
friend CMatrix coth(CMatrix &m);               ///////////////////// tanh matrix


 
//////////////////////////////////////////////////////////////////////////////////// Roots

friend CMatrix sqrt(CMatrix &m);               /////// sqrt matrix

//////////////////////////////////////////////////////////////////////////////////// Powers

CMatrix powElement(double p) ;                 /////////////////// power each element
CMatrix pow_1x1_matrix(double p) ;             /////////////////// fraction power for 1x1 matrix
CMatrix pow_matrix(int p) ;                    /////////////////// power for matrix

//////////////////////////////////////////////////////////////////////////////////// Logarithmic

friend CMatrix log(CMatrix &m);               /////// log matrix
friend CMatrix ln(CMatrix &m);                /////// ln matrix
friend CMatrix exp(CMatrix &m);               /////// exp matrix


//////////////////////////////////////////////////////////////////////////////////// zeros,ones,random,eye


///////////////////////////////////////////////////////////////////////////
void concatinate(CMatrix& m);
string sendString();


};


#endif // CMATRIX_H


