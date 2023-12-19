# gcomplex - a simple encapsulation of gsl gsl_complex
## Use the project

Until 2023.12.19, only gcomplex.h and gcompex.cpp are needed after setting up the gsl library. Include them in your project and you can use `a+b*I` to get a complex number. 

PS: I really hate it that some languages such as matlab use log to refer to natural logarithm, especially in the case that the more intuitive ln exists. So in this library, `ln` is used instead of `log`.

