/*************************************************************************
* 	rectangle.1.bas_fun.c : by R.Lie
*/

#include <math.h>

#define GAUSS_ELIMINATION 								\
	for (i = 0;i < 3;i ++) {							\
		k = i;									\
		for (j = i+1;j < 4;j ++)						\
			if (fabs(m[j][i]) > fabs(m[k][i])) k = j;			\
		if (k != i) {								\
			for (j = i;j < 4;j ++) {					\
				tmp = m[i][j];						\
				m[i][j] = m[k][j];					\
				m[k][j] = tmp;						\
			}								\
			tmp = a[i][0];							\
			a[i][0] = a[k][0];						\
			a[k][0] = tmp;							\
			tmp = a[i][1];							\
			a[i][1] = a[k][1];						\
			a[k][1] = tmp;							\
		}									\
		for (j = i+1;j < 4;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i+1;k < 4;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];						\
		}									\
	}										\
	a[3][0] /= m[3][3];								\
	a[3][1] /= m[3][3];								\
	for (i = 2;i >= 0;i --) {							\
		for (j = i+1;j < 4;j ++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];					\
		}									\
		a[i][0] /= m[i][i];							\
		a[i][1] /= m[i][i];							\
	}


void lambda_1(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = (1.0 - xi)*(1.0 - eta);
}

void lambda_2(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = xi*(1.0 - eta);
}

void lambda_3(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = xi*eta;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = (1.0 - xi)*eta;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = -(1.0 - eta)*(a[1][0] + a[3][0]*p[1]) - (1.0 - xi)*(a[1][1] + a[3][1]*p[1]);
	val[1] = -(1.0 - eta)*(a[2][0] + a[3][0]*p[0]) - (1.0 - xi)*(a[2][1] + a[3][1]*p[0]);
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = (1.0 - eta)*(a[1][0] + a[3][0]*p[1]) - xi*(a[1][1] + a[3][1]*p[1]);
	val[1] = (1.0 - eta)*(a[2][0] + a[3][0]*p[0]) - xi*(a[2][1] + a[3][1]*p[0]);
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = eta*(a[1][0] + a[3][0]*p[1]) + xi*(a[1][1] + a[3][1]*p[1]);
	val[1] = eta*(a[2][0] + a[3][0]*p[0]) + xi*(a[2][1] + a[3][1]*p[0]);
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	int i, j, k;
	double m[4][4], tmp, xi, eta;
	double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};
	double * val = (double *)value;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = v[i][0];
		m[i][2] = v[i][1];
		m[i][3] = v[i][0]*v[i][1];
	}
	GAUSS_ELIMINATION;
	xi = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];
	eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];
	val[0] = -eta*(a[1][0] + a[3][0]*p[1]) + (1.0 - xi)*(a[1][1] + a[3][1]*p[1]);
	val[1] = -eta*(a[2][0] + a[3][0]*p[0]) + (1.0 - xi)*(a[2][1] + a[3][1]*p[0]);
}

/*
*  end of file
**************************************************************************/
