#include "ccomplex.h"

std::complex<double> operator +(int a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}
std::complex<double> operator +(long a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}
std::complex<double> operator +(long long a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}

std::complex<double> operator +(unsigned int a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}
std::complex<double> operator +(unsigned long a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}
std::complex<double> operator +(unsigned long long a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}

std::complex<double> operator +(float a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}
std::complex<double> operator +(long double a, const std::complex<double>& b)
{
	std::complex<double> r(b.real() + (double)a, b.imag());
	return r;
}