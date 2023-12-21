#ifndef _GCOMPLEX_HPP_
#define _GCOMPLEX_HPP_
#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex_float_math.h>
#include <gsl/gsl_complex_long_double_math.h>

template <class _ty>
struct _common_c_complex {
	_ty dat[2];
};

template<class _ty>
class _basic_complex {
public:
	_ty _cp;
};

template <class _ty>
class gcomplex :private _basic_complex<_common_c_complex<_ty>> {
public:
	gcomplex() {
		this->_cp.dat[0]=0;
		this->_cp.dat[1]=0;
	}

};

template <>
class gcomplex<float> :private _basic_complex<gsl_complex_float> {
public:
	gcomplex()
	{
		_cp = gsl_complex_float_rect(0, 0);
	}
	gcomplex(const gcomplex& z)
	{
		this->_cp = z._cp;
	}
	gcomplex(const gsl_complex_float& z)
	{
		this->_cp = z;
	}
	gcomplex(const gsl_complex_float* z)
	{
		this->_cp = *z;
	}
	gcomplex(float x = 0, float y = 0, bool type = 0)
	{
		if (type)
			_cp = gsl_complex_float_polar(x, y);
		else
			_cp = gsl_complex_float_rect(x, y);
	}
	//成员特性函数
	float creal()const
	{
		return GSL_REAL(_cp);
	}
	float cimag()const
	{
		return GSL_IMAG(_cp);
	}
	float& real()
	{
		return GSL_REAL(_cp);
	}
	float& imag()
	{
		return GSL_IMAG(_cp);
	}
	float arg()const
	{
		return gsl_complex_float_arg(_cp);
	}
	float abs()const
	{
		return gsl_complex_float_abs(_cp);
	}
	float abs2()const
	{
		return gsl_complex_float_abs2(_cp);
	}
	float logabs()const
	{
		return gsl_complex_float_logabs(_cp);
	}

	gcomplex conj()const
	{
		return gcomplex(gsl_complex_float_conjugate(_cp));
	}

	//运算符重载
	gcomplex operator +(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_float_add(_cp, n._cp));
		return t;
	}
	gcomplex operator -(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_float_sub(_cp, n._cp));
		return t;
	}
	gcomplex operator *(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_float_mul(_cp, n._cp));
		return t;
	}
	gcomplex operator /(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_float_div(_cp, n._cp));
		return t;
	}

	gcomplex operator +(float n)const
	{
		gcomplex t(gsl_complex_float_add_real(_cp, n));
		return t;
	}
	gcomplex operator -(float n)const
	{
		gcomplex t(gsl_complex_float_sub_real(_cp, n));
		return t;
	}
	gcomplex operator *(float n)const
	{
		gcomplex t(gsl_complex_float_mul_real(_cp, n));
		return t;
	}
	gcomplex operator /(float n)const
	{
		gcomplex t(gsl_complex_float_div_real(_cp, n));
		return t;
	}

	gcomplex& operator +=(const gcomplex& z)
	{
		this->_cp = gsl_complex_float_add(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator -=(const gcomplex& z)
	{
		this->_cp = gsl_complex_float_sub(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator *=(const gcomplex& z)
	{
		this->_cp = gsl_complex_float_mul(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator /=(const gcomplex& z)
	{
		this->_cp = gsl_complex_float_div(this->_cp, z._cp);
		return *this;
	}

	gcomplex& operator +=(float x)
	{
		this->_cp = gsl_complex_float_add_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator -=(float x)
	{
		this->_cp = gsl_complex_float_sub_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator *=(float x)
	{
		this->_cp = gsl_complex_float_mul_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator /=(float x)
	{
		this->_cp = gsl_complex_float_div_real(this->_cp, x);
		return *this;
	}

	gcomplex& operator =(const gcomplex& z)
	{
		this->_cp = z._cp;
		return *this;
	}
	gcomplex& operator =(const gsl_complex_float& z)
	{
		this->_cp = z;
		return *this;
	}
	gcomplex& operator =(float a)
	{
		GSL_REAL(this->_cp) = a;
		GSL_IMAG(this->_cp) = 0;
		return *this;
	}

	//初等运算成员函数
	gcomplex sqrt()const
	{
		return gcomplex(gsl_complex_float_sqrt(_cp));
	}
	gcomplex pow(const gcomplex& a)const
	{
		return gcomplex(gsl_complex_float_pow(_cp, a._cp));
	}
	gcomplex pow(float a)const
	{
		return gcomplex(gsl_complex_float_pow_real(_cp, a));
	}
	gcomplex exp()const
	{
		return gcomplex(gsl_complex_float_exp(_cp));
	}
	gcomplex ln()const
	{
		return gcomplex(gsl_complex_float_log(_cp));
	}
	gcomplex log10()const
	{
		return gcomplex(gsl_complex_float_log10(_cp));
	}
	gcomplex log(gcomplex a)const
	{
		return gcomplex(gsl_complex_float_log_b(_cp, a._cp));
	}

	template <class _ty1>
	operator gcomplex<_ty1>()
	{
		gcomplex<_ty1> t;
		t._cp.dat[0] = (_ty1)(this->_cp.dat[0]);
		t._cp.dat[1] = (_ty1)(this->_cp.dat[1]);
		return *this;
	}
	//友元特性函数
	friend float& real(gcomplex&);
	friend float& imag(gcomplex&);
	friend float abs(const gcomplex&);
	friend float arg(const gcomplex&);
	friend float abs2(const gcomplex&);
	friend float logabs(const gcomplex&);
	friend gcomplex conj(const gcomplex&);

	//友元运算符重载函数
	friend gcomplex operator+(float, const gcomplex&);
	friend gcomplex operator-(float, const gcomplex&);
	friend gcomplex operator*(float, const gcomplex&);
	friend gcomplex operator/(float, const gcomplex&);


	//友元初等函数
	//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
	friend gcomplex sqrt(const gcomplex&);
	friend gcomplex gsqrt(float);
	friend gcomplex pow(const gcomplex&, const gcomplex&);
	friend gcomplex pow(const gcomplex&, float);
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

template <>
class gcomplex<double> :private _basic_complex<gsl_complex> {
public:
	gcomplex()
	{
		_cp = gsl_complex_rect(0, 0);
	}
	gcomplex(const gcomplex& z)
	{
		this->_cp = z._cp;
	}
	gcomplex(const gsl_complex& z)
	{
		this->_cp = z;
	}
	gcomplex(const gsl_complex* z)
	{
		this->_cp = *z;
	}
	gcomplex(double x = 0, double y = 0, bool type = 0)
	{
		if (type)
			_cp = gsl_complex_polar(x, y);
		else
			_cp = gsl_complex_rect(x, y);
	}
	//成员特性函数
	double creal()const
	{
		return GSL_REAL(_cp);
	}
	double cimag()const
	{
		return GSL_IMAG(_cp);
	}
	double& real()
	{
		return GSL_REAL(_cp);
	}
	double& imag()
	{
		return GSL_IMAG(_cp);
	}
	double arg()const
	{
		return gsl_complex_arg(_cp);
	}
	double abs()const
	{
		return gsl_complex_abs(_cp);
	}
	double abs2()const
	{
		return gsl_complex_abs2(_cp);
	}
	double logabs()const
	{
		return gsl_complex_logabs(_cp);
	}

	gcomplex conj()const
	{
		return gcomplex(gsl_complex_conjugate(_cp));
	}

	//运算符重载
	gcomplex operator +(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_add(_cp, n._cp));
		return t;
	}
	gcomplex operator -(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_sub(_cp, n._cp));
		return t;
	}
	gcomplex operator *(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_mul(_cp, n._cp));
		return t;
	}
	gcomplex operator /(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_div(_cp, n._cp));
		return t;
	}

	gcomplex operator +(double n)const
	{
		gcomplex t(gsl_complex_add_real(_cp, n));
		return t;
	}
	gcomplex operator -(double n)const
	{
		gcomplex t(gsl_complex_sub_real(_cp, n));
		return t;
	}
	gcomplex operator *(double n)const
	{
		gcomplex t(gsl_complex_mul_real(_cp, n));
		return t;
	}
	gcomplex operator /(double n)const
	{
		gcomplex t(gsl_complex_div_real(_cp, n));
		return t;
	}

	gcomplex& operator +=(const gcomplex& z)
	{
		this->_cp = gsl_complex_add(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator -=(const gcomplex& z)
	{
		this->_cp = gsl_complex_sub(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator *=(const gcomplex& z)
	{
		this->_cp = gsl_complex_mul(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator /=(const gcomplex& z)
	{
		this->_cp = gsl_complex_div(this->_cp, z._cp);
		return *this;
	}

	gcomplex& operator +=(double x)
	{
		this->_cp = gsl_complex_add_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator -=(double x)
	{
		this->_cp = gsl_complex_sub_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator *=(double x)
	{
		this->_cp = gsl_complex_mul_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator /=(double x)
	{
		this->_cp = gsl_complex_div_real(this->_cp, x);
		return *this;
	}

	gcomplex& operator =(const gcomplex& z)
	{
		this->_cp = z._cp;
		return *this;
	}
	gcomplex& operator =(const gsl_complex& z)
	{
		this->_cp = z;
		return *this;
	}
	gcomplex& operator =(double a)
	{
		GSL_REAL(this->_cp) = a;
		GSL_IMAG(this->_cp) = 0;
		return *this;
	}

	template <class _ty1>
	operator gcomplex<_ty1>()
	{
		gcomplex<_ty1> t;
		t._cp.dat[0] = (_ty1)(this->_cp.dat[0]);
		t._cp.dat[1] = (_ty1)(this->_cp.dat[1]);
		return *this;
	}
	//初等运算成员函数
	gcomplex sqrt()const
	{
		return gcomplex(gsl_complex_sqrt(_cp));
	}
	gcomplex pow(const gcomplex& a)const
	{
		return gcomplex(gsl_complex_pow(_cp, a._cp));
	}
	gcomplex pow(double a)const
	{
		return gcomplex(gsl_complex_pow_real(_cp, a));
	}
	gcomplex exp()const
	{
		return gcomplex(gsl_complex_exp(_cp));
	}
	gcomplex ln()const
	{
		return gcomplex(gsl_complex_log(_cp));
	}
	gcomplex log10()const
	{
		return gcomplex(gsl_complex_log10(_cp));
	}
	gcomplex log(gcomplex a)const
	{
		return gcomplex(gsl_complex_log_b(_cp, a._cp));
	}


	//友元特性函数
	friend double& real(gcomplex&);
	friend double& imag(gcomplex&);
	friend double abs(const gcomplex&);
	friend double arg(const gcomplex&);
	friend double abs2(const gcomplex&);
	friend double logabs(const gcomplex&);
	friend gcomplex conj(const gcomplex&);

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

template <>
class gcomplex<long double> :private _basic_complex<gsl_complex_long_double> {
public:
	gcomplex()
	{
		_cp = gsl_complex_long_double_rect(0, 0);
	}
	gcomplex(const gcomplex& z)
	{
		this->_cp = z._cp;
	}
	gcomplex(const gsl_complex_long_double& z)
	{
		this->_cp = z;
	}
	gcomplex(const gsl_complex_long_double* z)
	{
		this->_cp = *z;
	}
	gcomplex(long double x = 0, long double y = 0, bool type = 0)
	{
		if (type)
			_cp = gsl_complex_long_double_polar(x, y);
		else
			_cp = gsl_complex_long_double_rect(x, y);
	}
	//成员特性函数
	long double creal()const
	{
		return GSL_REAL(_cp);
	}
	long double cimag()const
	{
		return GSL_IMAG(_cp);
	}
	long double& real()
	{
		return GSL_REAL(_cp);
	}
	long double& imag()
	{
		return GSL_IMAG(_cp);
	}
	long double arg()const
	{
		return gsl_complex_long_double_arg(_cp);
	}
	long double abs()const
	{
		return gsl_complex_long_double_abs(_cp);
	}
	long double abs2()const
	{
		return gsl_complex_long_double_abs2(_cp);
	}
	long double logabs()const
	{
		return gsl_complex_long_double_logabs(_cp);
	}

	gcomplex conj()const
	{
		return gcomplex(gsl_complex_long_double_conjugate(_cp));
	}

	//运算符重载
	gcomplex operator +(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_long_double_add(_cp, n._cp));
		return t;
	}
	gcomplex operator -(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_long_double_sub(_cp, n._cp));
		return t;
	}
	gcomplex operator *(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_long_double_mul(_cp, n._cp));
		return t;
	}
	gcomplex operator /(const gcomplex& n)const
	{
		gcomplex t(gsl_complex_long_double_div(_cp, n._cp));
		return t;
	}

	gcomplex operator +(long double n)const
	{
		gcomplex t(gsl_complex_long_double_add_real(_cp, n));
		return t;
	}
	gcomplex operator -(long double n)const
	{
		gcomplex t(gsl_complex_long_double_sub_real(_cp, n));
		return t;
	}
	gcomplex operator *(long double n)const
	{
		gcomplex t(gsl_complex_long_double_mul_real(_cp, n));
		return t;
	}
	gcomplex operator /(long double n)const
	{
		gcomplex t(gsl_complex_long_double_div_real(_cp, n));
		return t;
	}

	gcomplex& operator +=(const gcomplex& z)
	{
		this->_cp = gsl_complex_long_double_add(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator -=(const gcomplex& z)
	{
		this->_cp = gsl_complex_long_double_sub(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator *=(const gcomplex& z)
	{
		this->_cp = gsl_complex_long_double_mul(this->_cp, z._cp);
		return *this;
	}
	gcomplex& operator /=(const gcomplex& z)
	{
		this->_cp = gsl_complex_long_double_div(this->_cp, z._cp);
		return *this;
	}

	gcomplex& operator +=(long double x)
	{
		this->_cp = gsl_complex_long_double_add_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator -=(long double x)
	{
		this->_cp = gsl_complex_long_double_sub_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator *=(long double x)
	{
		this->_cp = gsl_complex_long_double_mul_real(this->_cp, x);
		return *this;
	}
	gcomplex& operator /=(long double x)
	{
		this->_cp = gsl_complex_long_double_div_real(this->_cp, x);
		return *this;
	}

	gcomplex& operator =(const gcomplex& z)
	{
		this->_cp = z._cp;
		return *this;
	}
	gcomplex& operator =(const gsl_complex_long_double& z)
	{
		this->_cp = z;
		return *this;
	}
	gcomplex& operator =(long double a)
	{
		GSL_REAL(this->_cp) = a;
		GSL_IMAG(this->_cp) = 0;
		return *this;
	}
	template <class _ty1>
	operator gcomplex<_ty1>()
	{
		gcomplex<_ty1> t;
		t._cp.dat[0] = (_ty1)(this->_cp.dat[0]);
		t._cp.dat[1] = (_ty1)(this->_cp.dat[1]);
		return *this;
	}
	//初等运算成员函数
	gcomplex sqrt()const
	{
		return gcomplex(gsl_complex_long_double_sqrt(_cp));
	}
	gcomplex pow(const gcomplex& a)const
	{
		return gcomplex(gsl_complex_long_double_pow(_cp, a._cp));
	}
	gcomplex pow(long double a)const
	{
		return gcomplex(gsl_complex_long_double_pow_real(_cp, a));
	}
	gcomplex exp()const
	{
		return gcomplex(gsl_complex_long_double_exp(_cp));
	}
	gcomplex ln()const
	{
		return gcomplex(gsl_complex_long_double_log(_cp));
	}
	gcomplex log10()const
	{
		return gcomplex(gsl_complex_long_double_log10(_cp));
	}
	gcomplex log(gcomplex a)const
	{
		return gcomplex(gsl_complex_long_double_log_b(_cp, a._cp));
	}

	//友元特性函数
	friend long double& real(gcomplex&);
	friend long double& imag(gcomplex&);
	friend long double abs(const gcomplex&);
	friend long double arg(const gcomplex&);
	friend long double abs2(const gcomplex&);
	friend long double logabs(const gcomplex&);
	friend gcomplex conj(const gcomplex&);

	//友元运算符重载函数
	friend gcomplex operator+(long double, const gcomplex&);
	friend gcomplex operator-(long double, const gcomplex&);
	friend gcomplex operator*(long double, const gcomplex&);
	friend gcomplex operator/(long double, const gcomplex&);

	//友元初等函数
	//gsl_complex_sqrt_real若封装为sqrt则与cmath库冲突，因此换个名字
	friend gcomplex sqrt(const gcomplex&);
	friend gcomplex gsqrt(long double);
	friend gcomplex pow(const gcomplex&, const gcomplex&);
	friend gcomplex pow(const gcomplex&, long double);
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

extern const gcomplex<double> I;
extern const gcomplex<double> J;
#endif // _GCOMPLEX_HPP_
