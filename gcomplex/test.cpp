
#include <iostream>
#include <string>
#include "gcomplex.h"

using namespace std;
int main()
{
	std::cout << "Hello World!\n";
	gcomplex a = I * 3 + 2;
	gcomplex b = 3 + 4 * I;

	std::cout << "abs(b)=" << b.abs() << endl;
	string s1, s2;
	string s3 = s1 + s2;
	gcomplex c = 5;
	gcomplex d{ 2,8 };
	gcomplex e = c / d;
	cout << e.creal() << "\t" << e.cimag() << endl;
	gcomplex f = sqrt(gcomplex(-100));
	cout << f.creal() << "\t" << f.cimag() << endl;
	cout << abs(f)<<endl;
	cout << f;
	real(f) = 10;
	cout << "\t" << f << endl;
	cout << gsqrt(-5)<<endl;
	cout << ln(f) << endl;
	cout << pow(10 + 6 * I, 2)<<endl;
	cout << 2 + 3 * I << endl;
	cout << conj(2 + 3 * I)<<endl;
	s1 = "25";
	s2 = s1;
	//to_integer(s1);
	cout << abs(3 + 1 * I + 3 * I) << endl<<endl;
	cout << sin(3.1415926 / 2 + 0 * I) << endl;
	cout << cos(3.1415926 / 2 + 0 * I) << endl;
	cout << tan(3.1415926 / 2 + 0 * I) << endl;
	cout << sec(3.1415926 / 2 + 0 * I) << endl;
	cout << csc(3.1415926 / 2 + 0 * I) << endl;
	cout << cot(3.1415926 / 2 + 0 * I) << endl;
	cout << "==========\n";
	cout << sin(3.1415926 / 2) << endl;
	cout << cos(3.1415926 / 2) << endl;
	cout << tan(3.1415926 / 2) << endl;
	cout << sec(3.1415926 / 2) << endl;
	cout << csc(3.1415926 / 2) << endl;
	cout << cot(3.1415926 / 2) << endl;

	f += 2 * I;
	f /= 2;
	f *= 3 + 0 * I;
	f /= I;
	f -= 3;
	f -= 3 * I;
	cout << "f=" << f << endl;
	cout << sizeof(gcomplex)<<" "<<sizeof(gsl_complex)<<endl;
}