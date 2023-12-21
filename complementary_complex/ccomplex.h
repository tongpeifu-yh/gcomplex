#pragma once
#ifndef _CCOMPLEX_H_
#define _CCOMPLEX_H_

#include <complex>

std::complex<double> operator +(int a, const std::complex<double>& b);
std::complex<double> operator +(long a, const std::complex<double>& b);
std::complex<double> operator +(long long a, const std::complex<double>& b);

std::complex<double> operator +(unsigned int a, const std::complex<double>& b);
std::complex<double> operator +(unsigned long a, const std::complex<double>& b);
std::complex<double> operator +(unsigned long long a, const std::complex<double>& b);

std::complex<double> operator +(float a, const std::complex<double>& b);
std::complex<double> operator +(long double a, const std::complex<double>& b);
#endif // _CCOMPLEX_H_