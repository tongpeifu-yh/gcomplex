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
	gcomplex operator +(gcomplex)const;
	gcomplex operator -(gcomplex)const;
	gcomplex operator *(gcomplex)const;
	gcomplex operator /(gcomplex)const;
	gcomplex operator +(double)const;
	gcomplex operator -(double)const;
	gcomplex operator *(double)const;
	gcomplex operator /(double)const;

	//初等运算成员函数
	//gsl_complex_sqrt_real接收double类型，故不作为成员函数提供
	gcomplex sqrt()const;
	gcomplex pow(gcomplex)const;
	gcomplex pow(double)const;
	gcomplex exp()const;
	gcomplex ln()const;
	gcomplex log10()const;
	gcomplex log(gcomplex)const;

	//友元特性函数
	friend double& real(gcomplex&);
	friend double& imag(gcomplex&);
	friend double abs(gcomplex);
	friend double arg(gcomplex);
	friend double abs2(gcomplex);
	friend double logabs(gcomplex);

	//友元运算符重载函数
	friend gcomplex operator+(double, gcomplex);
	friend gcomplex operator-(double, gcomplex);
	friend gcomplex operator*(double, gcomplex);
	friend gcomplex operator/(double, gcomplex);

	//友元初等函数
	//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
	friend gcomplex sqrt(gcomplex);
	friend gcomplex gsqrt(double);
	friend gcomplex pow(gcomplex, gcomplex);
	friend gcomplex pow(gcomplex, double);
	friend gcomplex exp(gcomplex);
	friend gcomplex ln(gcomplex);
	friend gcomplex log10(gcomplex);
	friend gcomplex log(gcomplex,gcomplex);
	

	//其他友元函数
	friend std::ostream& operator<<(std::ostream& o, const gcomplex&);

};

//全局函数声明
//友元特性函数
double& real(gcomplex&);
double& imag(gcomplex&);
double abs(gcomplex);
double arg(gcomplex);
double abs2(gcomplex);
double logabs(gcomplex);

//友元运算符重载函数
gcomplex operator+(double, gcomplex);
gcomplex operator-(double, gcomplex);
gcomplex operator*(double, gcomplex);
gcomplex operator/(double, gcomplex);

//友元初等函数
//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
gcomplex sqrt(gcomplex);
gcomplex gsqrt(double);
gcomplex pow(gcomplex, gcomplex);
gcomplex pow(gcomplex, double);
gcomplex exp(gcomplex);
gcomplex ln(gcomplex);
gcomplex log10(gcomplex);
gcomplex log(gcomplex, gcomplex);


//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex&);


extern const gcomplex I;
extern const gcomplex J;



#endif