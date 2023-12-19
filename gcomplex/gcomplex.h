#ifndef _GCOMPLEX_H_
#define _GCOMPLEX_H_

#include <iostream>
#ifdef __cplusplus
extern "C" {
#endif

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef __cplusplus
}
#endif

class gcomplex {
private:
	gsl_complex _cp;
public:
	//gsl�ṩ����ѧ�����У������Ժ��������Ⱥ����г�Ա��������Ԫ���������汾�����������Ǻ���ֻ����Ԫ�����汾
	gcomplex();
	gcomplex(const gcomplex&);
	gcomplex(const gsl_complex&);
	gcomplex(const gsl_complex*);
	gcomplex(double x = 0, double y = 0,bool type=0);//type=0��rect ���� polar
	~gcomplex() {};
	//���Գ�Ա����
	double creal()const;
	double cimag()const;
	double &real();
	double &imag();
	double arg()const;
	double abs()const;
	double abs2()const;
	double logabs()const;

	gcomplex conj()const;
	//���������
	gcomplex operator +(gcomplex)const;
	gcomplex operator -(gcomplex)const;
	gcomplex operator *(gcomplex)const;
	gcomplex operator /(gcomplex)const;
	gcomplex operator +(double)const;
	gcomplex operator -(double)const;
	gcomplex operator *(double)const;
	gcomplex operator /(double)const;

	//���������Ա����
	//gsl_complex_sqrt_real����double���ͣ��ʲ���Ϊ��Ա�����ṩ
	gcomplex sqrt()const;
	gcomplex pow(gcomplex)const;
	gcomplex pow(double)const;
	gcomplex exp()const;
	gcomplex ln()const;
	gcomplex log10()const;
	gcomplex log(gcomplex)const;

	//��Ԫ���Ժ���
	friend double& real(gcomplex&);
	friend double& imag(gcomplex&);
	friend double abs(gcomplex);
	friend double arg(gcomplex);
	friend double abs2(gcomplex);
	friend double logabs(gcomplex);

	//��Ԫ��������غ���
	friend gcomplex operator+(double, gcomplex);
	friend gcomplex operator-(double, gcomplex);
	friend gcomplex operator*(double, gcomplex);
	friend gcomplex operator/(double, gcomplex);

	//��Ԫ���Ⱥ���
	//gsl_complex_sqrt_real����װΪsqrt����cmath���ͻ����˻�������
	friend gcomplex sqrt(gcomplex);
	friend gcomplex gsqrt(double);
	friend gcomplex pow(gcomplex, gcomplex);
	friend gcomplex pow(gcomplex, double);
	friend gcomplex exp(gcomplex);
	friend gcomplex ln(gcomplex);
	friend gcomplex log10(gcomplex);
	friend gcomplex log(gcomplex,gcomplex);
	

	//������Ԫ����
	friend std::ostream& operator<<(std::ostream& o, const gcomplex&);

};

//ȫ�ֺ�������
//��Ԫ���Ժ���
double& real(gcomplex&);
double& imag(gcomplex&);
double abs(gcomplex);
double arg(gcomplex);
double abs2(gcomplex);
double logabs(gcomplex);

//��Ԫ��������غ���
gcomplex operator+(double, gcomplex);
gcomplex operator-(double, gcomplex);
gcomplex operator*(double, gcomplex);
gcomplex operator/(double, gcomplex);

//��Ԫ���Ⱥ���
//gsl_complex_sqrt_real����װΪsqrt����cmath���ͻ����˻�������
gcomplex sqrt(gcomplex);
gcomplex gsqrt(double);
gcomplex pow(gcomplex, gcomplex);
gcomplex pow(gcomplex, double);
gcomplex exp(gcomplex);
gcomplex ln(gcomplex);
gcomplex log10(gcomplex);
gcomplex log(gcomplex, gcomplex);


//������Ԫ����
std::ostream& operator<<(std::ostream& o, const gcomplex&);


extern const gcomplex I;
extern const gcomplex J;



#endif