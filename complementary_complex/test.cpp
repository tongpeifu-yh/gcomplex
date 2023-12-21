// complementary_complex.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "ccomplex.h"
int main()
{
    using namespace std;
    complex<double> a = 3 + 4i;
    a = a + 5i;
    a = 6.5f + 3i;
    a = 6.3 + 4i;
    a = 5u + 3i;
    a = '2' + 3i;
}

