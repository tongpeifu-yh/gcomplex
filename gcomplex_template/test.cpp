// gcomplex_template.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <complex>
#include "gcomplex.hpp"
int main()
{
    using namespace std;

    std::cout << "Hello World!\n";

    gcomplex<double> a=2.5+3*I;
    cout << a;
    gcomplex <float> b = (5 + 6 * I);
    b = a;
    gcomplex <float> c = 5 + 4 * I;
    c = b;
    //complex b = 10;
}

