
#include <iostream>
#include "gcomplex.h"
using namespace std;
int main()
{
	std::cout << "Hello World!\n";
	gcomplex a = i * 3 + 2;
	gcomplex b = 3 + 4 * i;

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
	cout << ln(f);
}