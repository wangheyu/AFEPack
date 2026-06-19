/*************************************************************************
* 	rectangle.RT.1.bas_fun.cpp : by R.Lie
*/

#include <cmath>
#include <Miscellaneous.h>

/**
 * 本文件定义的是矩形单元上的零阶 Raviart-Thomas 单元的基函数
 * 及其梯度。这种单元定义在矩形上，每个单元上的四个自由度位于
 * 四条边上，基函数的表达式是由一维的线性元扩展到二维得到的，
 * 用 \f$ (\xi, \eta) \f$ 表示点的参数坐标，这四个基函数分别
 * 为：
 *
 * \f[
 *   \begin{array}{rcl}
 *     \lambda_1 &=& (0, \eta),      \mathrm{at\ } ( 0, -1) \\
 *     \lambda_2 &=& (\xi,  0),      \mathrm{at\ } ( 1,  0) \\
 *     \lambda_3 &=& (0, 1 - \eta),  \mathrm{at\ } ( 0,  1) \\
 *     \lambda_4 &=& (1 - \xi,  0),  \mathrm{at\ } (-1,  0)
 *   \end{array}
 * \f]
 * 
 * 它们的梯度可以相应的根据这些表达式计算出来。
 *
 */

#define GAUSS_ELIMINATION 								\
	for (i = 0;i < 3;i ++) {							\
		k = i;											\
		for (j = i+1;j < 4;j ++)						\
			if (fabs(m[j][i]) > fabs(m[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 4;j ++) {					\
				tmp = m[i][j];							\
				m[i][j] = m[k][j];						\
				m[k][j] = tmp;							\
			}											\
			tmp = a[i][0];								\
			a[i][0] = a[k][0];							\
			a[k][0] = tmp;								\
			tmp = a[i][1];								\
			a[i][1] = a[k][1];							\
			a[k][1] = tmp;								\
		}												\
		for (j = i+1;j < 4;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i+1;k < 4;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];						\
		}												\
	}													\
	a[3][0] /= m[3][3];									\
	a[3][1] /= m[3][3];									\
	for (i = 2;i >= 0;i --) {							\
		for (j = i+1;j < 4;j ++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];					\
		}												\
		a[i][0] /= m[i][i];								\
		a[i][1] /= m[i][i];								\
	}

#ifdef __cplusplus
extern "C" {
#endif

#define vector_length 2
#define vt nVector<vector_length,double>

void lambda_1(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt& val = *((vt *)value);
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = 0.0; val[1] = 1.0 - eta;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt& val = *((vt *)value);
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = xi; val[1] = 0.0;
}

void lambda_3(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt& val = *((vt *)value);
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = 0; val[1] = eta;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt& val = *((vt *)value);
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = 1.0 - xi; val[1] = 0;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt * val = (vt *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0][0] = 0; val[0][1] = -(a[1][1] + a[3][1]*p[1]);
	val[1][0] = 0; val[1][1] = -(a[2][1] + a[3][1]*p[0]);
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt * val = (vt *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0][0] = (a[1][0] + a[3][0]*p[1]); val[0][1] = 0;
	val[1][0] = (a[2][0] + a[3][0]*p[0]); val[1][1] = 0;
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt * val = (vt *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0][0] = 0; val[0][1] = (a[1][1] + a[3][1]*p[1]);
	val[1][0] = 0; val[1][1] = (a[2][1] + a[3][1]*p[0]);	
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	vt * val = (vt *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0][0] = -(a[1][0] + a[3][0]*p[1]); val[0][1] = 0;
	val[1][0] = -(a[2][0] + a[3][0]*p[0]); val[1][1] = 0;
}

#ifdef __cplusplus
}
#endif

/*
*  end of file
**************************************************************************/
