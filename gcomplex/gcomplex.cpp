#include "gcomplex.h"
gcomplex::gcomplex()
{
	_cp = gsl_complex_rect(0, 0);
}
gcomplex::gcomplex(const gcomplex& z)
{
	this->_cp = z._cp;
}
gcomplex::gcomplex(const gsl_complex& z)
{
	this->_cp = z;
}
gcomplex::gcomplex(const gsl_complex* z)
{
	this->_cp = *z;
}
gcomplex::gcomplex(double x /* = 0 */, double y /* = 0 */,bool type)
{
	if (type)
		_cp = gsl_complex_polar(x, y);
	else
		_cp = gsl_complex_rect(x, y);
}
//成员特性函数
double gcomplex::creal()const
{
	return GSL_REAL(_cp);
}
double gcomplex::cimag()const
{
	return GSL_IMAG(_cp);
}
double& gcomplex::real()
{
	return GSL_REAL(_cp);
}
double& gcomplex::imag()
{
	return GSL_IMAG(_cp);
}
double gcomplex::arg()const
{
	return gsl_complex_arg(_cp);
}
double gcomplex::abs()const
{
	return gsl_complex_abs(_cp);
}
double gcomplex::abs2()const
{
	return gsl_complex_abs2(_cp);
}
double gcomplex::logabs()const
{
	return gsl_complex_logabs(_cp);
}

gcomplex gcomplex::conj()const
{
	return gcomplex(gsl_complex_conjugate(_cp));
}

//运算符重载
gcomplex gcomplex::operator +(const gcomplex& n)const
{
	gcomplex t(gsl_complex_add(_cp, n._cp));
	return t;
}
gcomplex gcomplex::operator -(const gcomplex& n)const
{
	gcomplex t(gsl_complex_sub(_cp, n._cp));
	return t;
}
gcomplex gcomplex::operator *(const gcomplex& n)const
{
	gcomplex t(gsl_complex_mul(_cp, n._cp));
	return t;
}
gcomplex gcomplex::operator /(const gcomplex& n)const
{
	gcomplex t(gsl_complex_div(_cp, n._cp));
	return t;
}

gcomplex gcomplex::operator +(double n)const
{
	gcomplex t(gsl_complex_add_real(_cp, n));
	return t;
}
gcomplex gcomplex::operator -(double n)const
{
	gcomplex t(gsl_complex_sub_real(_cp, n));
	return t;
}
gcomplex gcomplex::operator *(double n)const
{
	gcomplex t(gsl_complex_mul_real(_cp,n));
	return t;
}
gcomplex gcomplex::operator /(double n)const
{
	gcomplex t(gsl_complex_div_real(_cp, n));
	return t;
}

gcomplex& gcomplex::operator =(const gcomplex&z)
{
	this->_cp = z._cp;
	return *this;
}
gcomplex& gcomplex::operator =(double a)
{
	GSL_REAL(this->_cp) = a;
	GSL_IMAG(this->_cp) = 0;
	return *this;
}



//初等运算成员函数
gcomplex gcomplex::sqrt()const
{
	return gcomplex(gsl_complex_sqrt(_cp));
}
gcomplex gcomplex::pow(const gcomplex& a)const
{
	return gcomplex(gsl_complex_pow(_cp, a._cp));
}
gcomplex gcomplex::pow(double a)const
{
	return gcomplex(gsl_complex_pow_real(_cp, a));
}
gcomplex gcomplex::exp()const
{
	return gcomplex(gsl_complex_exp(_cp));
}
gcomplex gcomplex::ln()const
{
	return gcomplex(gsl_complex_log(_cp));
}
gcomplex gcomplex::log10()const
{
	return gcomplex(gsl_complex_log10(_cp));
}
gcomplex gcomplex::log(gcomplex a)const
{
	return gcomplex(gsl_complex_log_b(_cp, a._cp));
}


//友元特性函数
double& real(gcomplex &z)
{
	return GSL_REAL(z._cp);
}
double& imag(gcomplex &z)
{
	return GSL_IMAG(z._cp);
}
double abs(const gcomplex& z)
{
	return gsl_complex_abs(z._cp);
}
double arg(const gcomplex& z)
{
	return gsl_complex_arg(z._cp);
}
double abs2(const gcomplex& z)
{
	return gsl_complex_abs2(z._cp);
}
double logabs(const gcomplex& z)
{
	return gsl_complex_logabs(z._cp);
}
gcomplex conj(const gcomplex&z)
{
	return gcomplex(gsl_complex_conjugate(z._cp));
}


//友元运算符重载函数
//加、乘直接交换顺序即可，减、除不可
gcomplex operator+(double a, const gcomplex& b)
{
	gcomplex t(gsl_complex_add_real(b._cp, a));
	return t;
}
gcomplex operator-(double a, const gcomplex& b)
{
	gcomplex t(gsl_complex_add_real(gsl_complex_negative(b._cp), a));//a-b=-b+a
	return t;
}
gcomplex operator*(double a, const gcomplex& b)
{
	gcomplex t(gsl_complex_mul_real(b._cp, a));
	return t;
}
gcomplex operator/(double a, const gcomplex& b)
{
	gcomplex t(gsl_complex_div_real(gsl_complex_inverse(b._cp), a));// a/b=1/b*a
	return t;
}

//友元初等函数
gcomplex sqrt(const gcomplex& z)
{
	return gcomplex(gsl_complex_sqrt(z._cp));
}
gcomplex gsqrt(double a)
{
	return gcomplex(gsl_complex_sqrt_real(a));
}
gcomplex pow(const gcomplex& z, const gcomplex& a)
{
	return gcomplex(gsl_complex_pow(z._cp, a._cp));
}
gcomplex pow(const gcomplex& z, double a)
{
	return gcomplex(gsl_complex_pow_real(z._cp, a));
}
gcomplex exp(const gcomplex& z)
{
	return gcomplex(gsl_complex_exp(z._cp));
}
gcomplex ln(const gcomplex& z)
{
	return gcomplex(gsl_complex_log(z._cp));
}
gcomplex log10(const gcomplex& z)
{
	return gcomplex(gsl_complex_log10(z._cp));
}
gcomplex log(const gcomplex& a, const gcomplex& b)
{
	return gcomplex(gsl_complex_log_b(b._cp, a._cp));
}


//友元三角、反三角、双曲函数
gcomplex sin(const gcomplex& z)
{
	return gcomplex(gsl_complex_sin(z._cp));
}
gcomplex cos(const gcomplex& z)
{
	return gcomplex(gsl_complex_cos(z._cp));
}
gcomplex tan(const gcomplex& z)
{
	return gcomplex(gsl_complex_tan(z._cp));
}
gcomplex sec(const gcomplex& z)
{
	return gcomplex(gsl_complex_sec(z._cp));
}
gcomplex csc(const gcomplex& z)
{
	return gcomplex(gsl_complex_csc(z._cp));
}
gcomplex cot(const gcomplex& z)
{
	return gcomplex(gsl_complex_cot(z._cp));
}


gcomplex asin(const gcomplex& z)
{
	return gcomplex(gsl_complex_arcsin(z._cp));
}
gcomplex acos(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccos(z._cp));
}
gcomplex atan(const gcomplex& z)
{
	return gcomplex(gsl_complex_arctan(z._cp));
}
gcomplex asec(const gcomplex& z)
{
	return gcomplex(gsl_complex_arcsec(z._cp));
}
gcomplex acsc(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccsc(z._cp));
}
gcomplex acot(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccot(z._cp));
}


gcomplex sinh(const gcomplex& z)
{
	return gcomplex(gsl_complex_sinh(z._cp));
}
gcomplex cosh(const gcomplex& z)
{
	return gcomplex(gsl_complex_cosh(z._cp));
}
gcomplex tanh(const gcomplex& z)
{
	return gcomplex(gsl_complex_tanh(z._cp));
}
gcomplex sech(const gcomplex& z)
{
	return gcomplex(gsl_complex_sech(z._cp));
}
gcomplex csch(const gcomplex& z)
{
	return gcomplex(gsl_complex_csch(z._cp));
}
gcomplex coth(const gcomplex& z)
{
	return gcomplex(gsl_complex_coth(z._cp));
}

gcomplex asinh(const gcomplex& z)
{
	return gcomplex(gsl_complex_arcsinh(z._cp));
}
gcomplex acosh(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccosh(z._cp));
}
gcomplex atanh(const gcomplex& z)
{
	return gcomplex(gsl_complex_arctanh(z._cp));
}
gcomplex asech(const gcomplex& z)
{
	return gcomplex(gsl_complex_arcsech(z._cp));
}
gcomplex acsch(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccsch(z._cp));
}
gcomplex acoth(const gcomplex& z)
{
	return gcomplex(gsl_complex_arccoth(z._cp));
}

//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex& z)
{
	if (z.creal() == 0 && z.cimag() == 0)
		;
	//o << z.creal() << '+' << z.cimag() << 'i';
	else if (z.creal() == 0)
		o << z.cimag() << 'i';
	else if (z.cimag() == 0)
		o << z.creal();
	else if (z.cimag() > 0)
		o << z.creal() << '+' << z.cimag() << 'i';
	else
		o << z.creal() << z.cimag() << 'i';
	return o;
}

//常数i、j，单位纯虚数
const gcomplex I(0, 1);
const gcomplex J(0, 1);
