/*
gcomplex - a simple encapsulation of gsl gsl_complex

All member functions:
	1. Properties of complex numbers
		creal,cimag,real,imag,arg,abs,abs2,logabs,conj
	2. Complex arithmetic operators
		+ - * / = += -= *= /= with gcomplex and with double,= with gsl_complex
	3. Elementary Complex Functions
		sqrt,pow (with gcomplex and with double),exp,ln,log10,log

All friend functions
	1. Properties of complex numbers
		real,imag,arg,abs,abs2,logabs,conj
	2. Complex arithmetic operators
		+ - * / with double
	3. Elementary Complex Functions
		sqrt,gsqrt (sqrt with double returning gcomplex),pow (with gcomplex and with double),exp,ln,log10,log
	4. Complex Trigonometric Functions, Inverse Complex Trigonometric Functions
		sin,cos,tan,sec,csc,cot,asin,acos,atan,asec,acsc,acot
	5. Complex Hyperbolic Functions, Inverse Complex Hyperbolic Functions
		sinh,cosh,tanh,sech,csch,coth,asinh,acosh,atanh,asech,acsch,acoth

NOTE 
	1. The memory of gcomplex should be the same as gsl_complex so that gcomplex[] can be compatible with gsl_complex[]. 
	   Thus you should never add any virtual functions or non-static variables.

TODO 

*/

#ifndef _GCOMPLEX_H_
#define _GCOMPLEX_H_

#include <iostream>
#ifdef __cplusplus
extern "C" {
#endif

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef __cplusplus
}
#endif

class gcomplex {
private:
	gsl_complex _cp;
public:
	//gsl提供的数学函数中，仅特性函数、初等函数有成员函数和友元函数两个版本，其余如三角函数只有友元函数版本
	gcomplex();
	gcomplex(const gcomplex&);
	gcomplex(const gsl_complex&);
	gcomplex(const gsl_complex*);
	gcomplex(double x = 0, double y = 0,bool type=0);//type=0，rect 否则 polar
	~gcomplex() {};
	//特性成员函数
	double creal()const;
	double cimag()const;
	double &real();
	double &imag();
	double arg()const;
	double abs()const;
	double abs2()const;
	double logabs()const;

	gcomplex conj()const;
	//重载运算符
	gcomplex operator +(const gcomplex&)const;
	gcomplex operator -(const gcomplex&)const;
	gcomplex operator *(const gcomplex&)const;
	gcomplex operator /(const gcomplex&)const;
	gcomplex operator +(double)const;
	gcomplex operator -(double)const;
	gcomplex operator *(double)const;
	gcomplex operator /(double)const;

	gcomplex& operator +=(const gcomplex&);
	gcomplex& operator -=(const gcomplex&);
	gcomplex& operator *=(const gcomplex&);
	gcomplex& operator /=(const gcomplex&);
	gcomplex& operator +=(double);
	gcomplex& operator -=(double);
	gcomplex& operator *=(double);
	gcomplex& operator /=(double);

	gcomplex& operator =(const gcomplex&);
	gcomplex& operator =(const gsl_complex&);
	gcomplex& operator =(double);


	//初等运算成员函数
	//gsl_complex_sqrt_real接收double类型，故不作为成员函数提供
	gcomplex sqrt()const;
	gcomplex pow(const gcomplex&)const;
	gcomplex pow(double)const;
	gcomplex exp()const;
	gcomplex ln()const;
	gcomplex log10()const;
	gcomplex log(gcomplex)const;

	//友元特性函数
	friend double& real(gcomplex&);
	friend double& imag(gcomplex&);
	friend double abs(const gcomplex&);
	friend double arg(const gcomplex&);
	friend double abs2(const gcomplex&);
	friend double logabs(const gcomplex&);
	friend gcomplex conj(const gcomplex &);

	//友元运算符重载函数
	friend gcomplex operator+(double, const gcomplex&);
	friend gcomplex operator-(double, const gcomplex&);
	friend gcomplex operator*(double, const gcomplex&);
	friend gcomplex operator/(double, const gcomplex&);
	

	//友元初等函数
	//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
	friend gcomplex sqrt(const gcomplex&);
	friend gcomplex gsqrt(double);
	friend gcomplex pow(const gcomplex&, const gcomplex&);
	friend gcomplex pow(const gcomplex&, double);
	friend gcomplex exp(const gcomplex&);
	friend gcomplex ln(const gcomplex&);
	friend gcomplex log10(const gcomplex&);
	friend gcomplex log(const gcomplex&, const gcomplex&);
	
	//友元三角、反三角、双曲函数
	friend gcomplex sin(const gcomplex&);
	friend gcomplex cos(const gcomplex&);
	friend gcomplex tan(const gcomplex&);
	friend gcomplex sec(const gcomplex&);
	friend gcomplex csc(const gcomplex&);
	friend gcomplex cot(const gcomplex&);

	friend gcomplex asin(const gcomplex&);
	friend gcomplex acos(const gcomplex&);
	friend gcomplex atan(const gcomplex&);
	friend gcomplex asec(const gcomplex&);
	friend gcomplex acsc(const gcomplex&);
	friend gcomplex acot(const gcomplex&);

	friend gcomplex sinh(const gcomplex&);
	friend gcomplex cosh(const gcomplex&);
	friend gcomplex tanh(const gcomplex&);
	friend gcomplex sech(const gcomplex&);
	friend gcomplex csch(const gcomplex&);
	friend gcomplex coth(const gcomplex&);

	friend gcomplex asinh(const gcomplex&);
	friend gcomplex acosh(const gcomplex&);
	friend gcomplex atanh(const gcomplex&);
	friend gcomplex asech(const gcomplex&);
	friend gcomplex acsch(const gcomplex&);
	friend gcomplex acoth(const gcomplex&);
	//其他友元函数
	friend std::ostream& operator<<(std::ostream& o, const gcomplex&);

};

//全局函数声明
//友元特性函数
double& real(gcomplex&);
double& imag(gcomplex&);
double abs(const gcomplex&);
double arg(const gcomplex&);
double abs2(const gcomplex&);
double logabs(const gcomplex&);
gcomplex conj(const gcomplex&);

//友元运算符重载函数
gcomplex operator+(double, const gcomplex&);
gcomplex operator-(double, const gcomplex&);
gcomplex operator*(double, const gcomplex&);
gcomplex operator/(double, const gcomplex&);

//友元初等函数
//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
gcomplex sqrt(const gcomplex&);
gcomplex gsqrt(double);
gcomplex pow(const gcomplex&, const gcomplex&);
gcomplex pow(const gcomplex&, double);
gcomplex exp(const gcomplex&);
gcomplex ln(const gcomplex&);
gcomplex log10(const gcomplex&);
gcomplex log(const gcomplex&, const gcomplex&);

//友元三角、反三角、双曲函数
gcomplex sin(const gcomplex&);
gcomplex cos(const gcomplex&);
gcomplex tan(const gcomplex&);
gcomplex sec(const gcomplex&);
gcomplex csc(const gcomplex&);
gcomplex cot(const gcomplex&);

gcomplex asin(const gcomplex&);
gcomplex acos(const gcomplex&);
gcomplex atan(const gcomplex&);
gcomplex asec(const gcomplex&);
gcomplex acsc(const gcomplex&);
gcomplex acot(const gcomplex&);

gcomplex sinh(const gcomplex&);
gcomplex cosh(const gcomplex&);
gcomplex tanh(const gcomplex&);
gcomplex sech(const gcomplex&);
gcomplex csch(const gcomplex&);
gcomplex coth(const gcomplex&);

gcomplex asinh(const gcomplex&);
gcomplex acosh(const gcomplex&);
gcomplex atanh(const gcomplex&);
gcomplex asech(const gcomplex&);
gcomplex acsch(const gcomplex&);
gcomplex acoth(const gcomplex&);

//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex&);


extern const gcomplex I;
extern const gcomplex J;



#endif