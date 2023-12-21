#include "gcomplex.hpp"
// friend functions for gcomplex<float>
//友元特性函数
float& real(gcomplex<float>& z)
{
	return GSL_REAL(z._cp);
}
float& imag(gcomplex<float>& z)
{
	return GSL_IMAG(z._cp);
}
float abs(const gcomplex<float>& z)
{
	return gsl_complex_float_abs(z._cp);
}
float arg(const gcomplex<float>& z)
{
	return gsl_complex_float_arg(z._cp);
}
float abs2(const gcomplex<float>& z)
{
	return gsl_complex_float_abs2(z._cp);
}
float logabs(const gcomplex<float>& z)
{
	return gsl_complex_float_logabs(z._cp);
}
gcomplex<float> conj(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_conjugate(z._cp));
}


//友元运算符重载函数
//加、乘直接交换顺序即可，减、除不可
gcomplex<float> operator+(float a, const gcomplex<float>& b)
{
	gcomplex<float> t(gsl_complex_float_add_real(b._cp, a));
	return t;
}
gcomplex<float> operator-(float a, const gcomplex<float>& b)
{
	gcomplex<float> t(gsl_complex_float_add_real(gsl_complex_float_negative(b._cp), a));//a-b=-b+a
	return t;
}
gcomplex<float> operator*(float a, const gcomplex<float>& b)
{
	gcomplex<float> t(gsl_complex_float_mul_real(b._cp, a));
	return t;
}
gcomplex<float> operator/(float a, const gcomplex<float>& b)
{
	gcomplex<float> t(gsl_complex_float_div_real(gsl_complex_float_inverse(b._cp), a));// a/b=1/b*a
	return t;
}

//友元初等函数
gcomplex<float> sqrt(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_sqrt(z._cp));
}
gcomplex<float> gsqrt(float a)
{
	return gcomplex<float>(gsl_complex_float_sqrt_real(a));
}
gcomplex<float> pow(const gcomplex<float>& z, const gcomplex<float>& a)
{
	return gcomplex<float>(gsl_complex_float_pow(z._cp, a._cp));
}
gcomplex<float> pow(const gcomplex<float>& z, float a)
{
	return gcomplex<float>(gsl_complex_float_pow_real(z._cp, a));
}
gcomplex<float> exp(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_exp(z._cp));
}
gcomplex<float> ln(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_log(z._cp));
}
gcomplex<float> log10(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_log10(z._cp));
}
gcomplex<float> log(const gcomplex<float>& a, const gcomplex<float>& b)
{
	return gcomplex<float>(gsl_complex_float_log_b(b._cp, a._cp));
}


//友元三角、反三角、双曲函数
gcomplex<float> sin(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_sin(z._cp));
}
gcomplex<float> cos(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_cos(z._cp));
}
gcomplex<float> tan(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_tan(z._cp));
}
gcomplex<float> sec(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_sec(z._cp));
}
gcomplex<float> csc(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_csc(z._cp));
}
gcomplex<float> cot(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_cot(z._cp));
}


gcomplex<float> asin(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arcsin(z._cp));
}
gcomplex<float> acos(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccos(z._cp));
}
gcomplex<float> atan(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arctan(z._cp));
}
gcomplex<float> asec(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arcsec(z._cp));
}
gcomplex<float> acsc(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccsc(z._cp));
}
gcomplex<float> acot(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccot(z._cp));
}


gcomplex<float> sinh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_sinh(z._cp));
}
gcomplex<float> cosh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_cosh(z._cp));
}
gcomplex<float> tanh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_tanh(z._cp));
}
gcomplex<float> sech(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_sech(z._cp));
}
gcomplex<float> csch(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_csch(z._cp));
}
gcomplex<float> coth(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_coth(z._cp));
}

gcomplex<float> asinh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arcsinh(z._cp));
}
gcomplex<float> acosh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccosh(z._cp));
}
gcomplex<float> atanh(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arctanh(z._cp));
}
gcomplex<float> asech(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arcsech(z._cp));
}
gcomplex<float> acsch(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccsch(z._cp));
}
gcomplex<float> acoth(const gcomplex<float>& z)
{
	return gcomplex<float>(gsl_complex_float_arccoth(z._cp));
}

//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex<float>& z)
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

// friend functions for gcomplex<double>
//友元特性函数
double& real(gcomplex<double>& z)
{
	return GSL_REAL(z._cp);
}
double& imag(gcomplex<double>& z)
{
	return GSL_IMAG(z._cp);
}
double abs(const gcomplex<double>& z)
{
	return gsl_complex_abs(z._cp);
}
double arg(const gcomplex<double>& z)
{
	return gsl_complex_arg(z._cp);
}
double abs2(const gcomplex<double>& z)
{
	return gsl_complex_abs2(z._cp);
}
double logabs(const gcomplex<double>& z)
{
	return gsl_complex_logabs(z._cp);
}
gcomplex<double> conj(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_conjugate(z._cp));
}


//友元运算符重载函数
//加、乘直接交换顺序即可，减、除不可
gcomplex<double> operator+(double a, const gcomplex<double>& b)
{
	gcomplex<double> t(gsl_complex_add_real(b._cp, a));
	return t;
}
gcomplex<double> operator-(double a, const gcomplex<double>& b)
{
	gcomplex<double> t(gsl_complex_add_real(gsl_complex_negative(b._cp), a));//a-b=-b+a
	return t;
}
gcomplex<double> operator*(double a, const gcomplex<double>& b)
{
	gcomplex<double> t(gsl_complex_mul_real(b._cp, a));
	return t;
}
gcomplex<double> operator/(double a, const gcomplex<double>& b)
{
	gcomplex<double> t(gsl_complex_div_real(gsl_complex_inverse(b._cp), a));// a/b=1/b*a
	return t;
}

//友元初等函数
gcomplex<double> sqrt(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_sqrt(z._cp));
}
gcomplex<double> gsqrt(double a)
{
	return gcomplex<double>(gsl_complex_sqrt_real(a));
}
gcomplex<double> pow(const gcomplex<double>& z, const gcomplex<double>& a)
{
	return gcomplex<double>(gsl_complex_pow(z._cp, a._cp));
}
gcomplex<double> pow(const gcomplex<double>& z, double a)
{
	return gcomplex<double>(gsl_complex_pow_real(z._cp, a));
}
gcomplex<double> exp(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_exp(z._cp));
}
gcomplex<double> ln(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_log(z._cp));
}
gcomplex<double> log10(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_log10(z._cp));
}
gcomplex<double> log(const gcomplex<double>& a, const gcomplex<double>& b)
{
	return gcomplex<double>(gsl_complex_log_b(b._cp, a._cp));
}


//友元三角、反三角、双曲函数
gcomplex<double> sin(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_sin(z._cp));
}
gcomplex<double> cos(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_cos(z._cp));
}
gcomplex<double> tan(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_tan(z._cp));
}
gcomplex<double> sec(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_sec(z._cp));
}
gcomplex<double> csc(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_csc(z._cp));
}
gcomplex<double> cot(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_cot(z._cp));
}


gcomplex<double> asin(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arcsin(z._cp));
}
gcomplex<double> acos(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccos(z._cp));
}
gcomplex<double> atan(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arctan(z._cp));
}
gcomplex<double> asec(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arcsec(z._cp));
}
gcomplex<double> acsc(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccsc(z._cp));
}
gcomplex<double> acot(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccot(z._cp));
}


gcomplex<double> sinh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_sinh(z._cp));
}
gcomplex<double> cosh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_cosh(z._cp));
}
gcomplex<double> tanh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_tanh(z._cp));
}
gcomplex<double> sech(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_sech(z._cp));
}
gcomplex<double> csch(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_csch(z._cp));
}
gcomplex<double> coth(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_coth(z._cp));
}

gcomplex<double> asinh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arcsinh(z._cp));
}
gcomplex<double> acosh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccosh(z._cp));
}
gcomplex<double> atanh(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arctanh(z._cp));
}
gcomplex<double> asech(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arcsech(z._cp));
}
gcomplex<double> acsch(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccsch(z._cp));
}
gcomplex<double> acoth(const gcomplex<double>& z)
{
	return gcomplex<double>(gsl_complex_arccoth(z._cp));
}

//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex<double>& z)
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


// friend functions for gcomplex<long double>
//友元特性函数
long double& real(gcomplex<long double>& z)
{
	return GSL_REAL(z._cp);
}
long double& imag(gcomplex<long double>& z)
{
	return GSL_IMAG(z._cp);
}
long double abs(const gcomplex<long double>& z)
{
	return gsl_complex_long_double_abs(z._cp);
}
long double arg(const gcomplex<long double>& z)
{
	return gsl_complex_long_double_arg(z._cp);
}
long double abs2(const gcomplex<long double>& z)
{
	return gsl_complex_long_double_abs2(z._cp);
}
long double logabs(const gcomplex<long double>& z)
{
	return gsl_complex_long_double_logabs(z._cp);
}
gcomplex<long double> conj(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_conjugate(z._cp));
}


//友元运算符重载函数
//加、乘直接交换顺序即可，减、除不可
gcomplex<long double> operator+(long double a, const gcomplex<long double>& b)
{
	gcomplex<long double> t(gsl_complex_long_double_add_real(b._cp, a));
	return t;
}
gcomplex<long double> operator-(long double a, const gcomplex<long double>& b)
{
	gcomplex<long double> t(gsl_complex_long_double_add_real(gsl_complex_long_double_negative(b._cp), a));//a-b=-b+a
	return t;
}
gcomplex<long double> operator*(long double a, const gcomplex<long double>& b)
{
	gcomplex<long double> t(gsl_complex_long_double_mul_real(b._cp, a));
	return t;
}
gcomplex<long double> operator/(long double a, const gcomplex<long double>& b)
{
	gcomplex<long double> t(gsl_complex_long_double_div_real(gsl_complex_long_double_inverse(b._cp), a));// a/b=1/b*a
	return t;
}

//友元初等函数
gcomplex<long double> sqrt(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_sqrt(z._cp));
}
gcomplex<long double> gsqrt(long double a)
{
	return gcomplex<long double>(gsl_complex_long_double_sqrt_real(a));
}
gcomplex<long double> pow(const gcomplex<long double>& z, const gcomplex<long double>& a)
{
	return gcomplex<long double>(gsl_complex_long_double_pow(z._cp, a._cp));
}
gcomplex<long double> pow(const gcomplex<long double>& z, long double a)
{
	return gcomplex<long double>(gsl_complex_long_double_pow_real(z._cp, a));
}
gcomplex<long double> exp(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_exp(z._cp));
}
gcomplex<long double> ln(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_log(z._cp));
}
gcomplex<long double> log10(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_log10(z._cp));
}
gcomplex<long double> log(const gcomplex<long double>& a, const gcomplex<long double>& b)
{
	return gcomplex<long double>(gsl_complex_long_double_log_b(b._cp, a._cp));
}


//友元三角、反三角、双曲函数
gcomplex<long double> sin(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_sin(z._cp));
}
gcomplex<long double> cos(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_cos(z._cp));
}
gcomplex<long double> tan(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_tan(z._cp));
}
gcomplex<long double> sec(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_sec(z._cp));
}
gcomplex<long double> csc(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_csc(z._cp));
}
gcomplex<long double> cot(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_cot(z._cp));
}


gcomplex<long double> asin(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arcsin(z._cp));
}
gcomplex<long double> acos(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccos(z._cp));
}
gcomplex<long double> atan(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arctan(z._cp));
}
gcomplex<long double> asec(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arcsec(z._cp));
}
gcomplex<long double> acsc(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccsc(z._cp));
}
gcomplex<long double> acot(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccot(z._cp));
}


gcomplex<long double> sinh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_sinh(z._cp));
}
gcomplex<long double> cosh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_cosh(z._cp));
}
gcomplex<long double> tanh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_tanh(z._cp));
}
gcomplex<long double> sech(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_sech(z._cp));
}
gcomplex<long double> csch(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_csch(z._cp));
}
gcomplex<long double> coth(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_coth(z._cp));
}

gcomplex<long double> asinh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arcsinh(z._cp));
}
gcomplex<long double> acosh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccosh(z._cp));
}
gcomplex<long double> atanh(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arctanh(z._cp));
}
gcomplex<long double> asech(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arcsech(z._cp));
}
gcomplex<long double> acsch(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccsch(z._cp));
}
gcomplex<long double> acoth(const gcomplex<long double>& z)
{
	return gcomplex<long double>(gsl_complex_long_double_arccoth(z._cp));
}

//其他友元函数
std::ostream& operator<<(std::ostream& o, const gcomplex<long double>& z)
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


const gcomplex<double> I(0, 1);
const gcomplex<double> J(0, 1);

