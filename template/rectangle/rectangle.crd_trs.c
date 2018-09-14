/*************************************************************************
*  rectangle.crd_trs.cpp : by R.Lie
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

void local_to_global(const double * lp, 
		const double ** lv, 
		const double ** gv, 
		double * gp)
{
	int i, j, k;
	double m[4][4], a[4][2], tmp;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = lv[i][0];
		m[i][2] = lv[i][1];
		m[i][3] = lv[i][0]*lv[i][1];
		a[i][0] = gv[i][0];
		a[i][1] = gv[i][1];
	}
	GAUSS_ELIMINATION;
	gp[0] = a[0][0] + a[1][0]*lp[0] + a[2][0]*lp[1] + a[3][0]*lp[0]*lp[1];
	gp[1] = a[0][1] + a[1][1]*lp[0] + a[2][1]*lp[1] + a[3][1]*lp[0]*lp[1];
}

void global_to_local(const double * gp, 
		const double ** lv, 
		const double ** gv, 
		double * lp)
{
	int i, j, k;
	double m[4][4], a[4][2], tmp;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = gv[i][0];
		m[i][2] = gv[i][1];
		m[i][3] = gv[i][0]*gv[i][1];
		a[i][0] = lv[i][0];
		a[i][1] = lv[i][1];
	}
	GAUSS_ELIMINATION;
	lp[0] = a[0][0] + a[1][0]*gp[0] + a[2][0]*gp[1] + a[3][0]*gp[0]*gp[1];
	lp[1] = a[0][1] + a[1][1]*gp[0] + a[2][1]*gp[1] + a[3][1]*gp[0]*gp[1];
}

double local_to_global_jacobian(const double * lp, 
		const double ** lv, 
		const double ** gv)
{
	int i, j, k;
	double m[4][4], a[4][2], tmp;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = lv[i][0];
		m[i][2] = lv[i][1];
		m[i][3] = lv[i][0]*lv[i][1];
		a[i][0] = gv[i][0];
		a[i][1] = gv[i][1];
	}
	GAUSS_ELIMINATION;
	return (a[1][0] + a[3][0]*lp[1])*(a[2][1] + a[3][1]*lp[0])
		- (a[2][0] + a[3][0]*lp[0])*(a[1][1] + a[3][1]*lp[1]);
}

double global_to_local_jacobian(const double * gp, 
		const double ** lv,
		const double ** gv)
{
	int i, j, k;
	double m[4][4], a[4][2], tmp;
	for (i = 0;i < 4;i ++) {
		m[i][0] = 1.0;
		m[i][1] = gv[i][0];
		m[i][2] = gv[i][1];
		m[i][3] = gv[i][0]*gv[i][1];
		a[i][0] = lv[i][0];
		a[i][1] = lv[i][1];
	}
	GAUSS_ELIMINATION;
	return (a[1][0] + a[3][0]*gp[1])*(a[2][1] + a[3][1]*gp[0])
		- (a[2][0] + a[3][0]*gp[0])*(a[1][1] + a[3][1]*gp[1]);
}

/*
*  end of file
**************************************************************************/
