/*************************************************************************
 * 	rectangle.2.bas_fun.cpp : by R.Lie
 */

#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

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

#define COMMON_PART                                                                     \
    int i, j, k;                                                                        \
    double m[4][4], tmp, xi, eta;                                                       \
    double a[4][2] = {{0.0, 0.0},{1.0, 0.0},{1.0, 1.0},{0.0, 1.0}};                     \
    for (i = 0;i < 4;i ++) {                                                            \
      m[i][0] = 1.0;                                                                    \
      m[i][1] = v[i][0];                                                                \
      m[i][2] = v[i][1];                                                                \
      m[i][3] = v[i][0]*v[i][1];                                                        \
    }                                                                                   \
    GAUSS_ELIMINATION;                                                                  \
    xi  = a[0][0] + a[1][0]*p[0] + a[2][0]*p[1] + a[3][0]*p[0]*p[1];                    \
    eta = a[0][1] + a[1][1]*p[0] + a[2][1]*p[1] + a[3][1]*p[0]*p[1];

#define D_XI_DX  ( a[1][0] + a[3][0]*p[1] )
#define D_XI_DY  ( a[2][0] + a[3][0]*p[0] )
#define D_ETA_DX ( a[1][1] + a[3][1]*p[1] )
#define D_ETA_DY ( a[2][1] + a[3][1]*p[0] )

  void lambda_mm(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 - xi)*(1.0 - eta)*(1.0 - eta);
  }

  void gradient_lambda_mm(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi - 2)*(1.0 - eta)*(1.0 - eta)*D_XI_DX
           + (1.0 - xi)*(1.0 - xi)*(2*eta - 2)*D_ETA_DX;
    val[1] = (2*xi - 2)*(1.0 - eta)*(1.0 - eta)*D_XI_DY
           + (1.0 - xi)*(1.0 - xi)*(2*eta - 2)*D_ETA_DY;
  }

  void lambda_0m(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*(1.0 - xi)*(1.0 - eta)*(1.0 - eta);
  }

  void gradient_lambda_0m(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1 - 2*xi)*(1.0 - eta)*(1.0 - eta)*D_XI_DX
           + xi*(1.0 - xi)*(2*eta - 2)*D_ETA_DX;
    val[1] = (1 - 2*xi)*(1.0 - eta)*(1.0 - eta)*D_XI_DY
           + xi*(1.0 - xi)*(2*eta - 2)*D_ETA_DY;
  }

  void lambda_pm(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*xi*(1.0 - eta)*(1.0 - eta);
  }

  void gradient_lambda_pm(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi)*(1.0 - eta)*(1.0 - eta)*D_XI_DX
           + xi*xi*(2*eta - 2)*D_ETA_DX;
    val[1] = (2*xi)*(1.0 - eta)*(1.0 - eta)*D_XI_DY
           + xi*xi*(2*eta - 2)*D_ETA_DY;
  }

  void lambda_m0(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 - xi)*eta*(1.0 - eta);
  }

  void gradient_lambda_m0(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi - 2)*eta*(1.0 - eta)*D_XI_DX
           + (1.0 - xi)*(1.0 - xi)*(1 - 2*eta)*D_ETA_DX;
    val[1] = (2*xi - 2)*eta*(1.0 - eta)*D_XI_DY
           + (1.0 - xi)*(1.0 - xi)*(1 - 2*eta)*D_ETA_DY;
  }

  void lambda_00(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*(1.0 - xi)*eta*(1.0 - eta);
  }

  void gradient_lambda_00(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1 - 2*xi)*eta*(1.0 - eta)*D_XI_DX
           + xi*(1.0 - xi)*(1 - 2*eta)*D_ETA_DX;
    val[1] = (1 - 2*xi)*eta*(1.0 - eta)*D_XI_DY
           + xi*(1.0 - xi)*(1 - 2*eta)*D_ETA_DY;
  }

  void lambda_p0(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*xi*eta*(1.0 - eta);
  }

  void gradient_lambda_p0(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi)*eta*(1.0 - eta)*D_XI_DX
           + xi*xi*(1 - 2*eta)*D_ETA_DX;
    val[1] = (2*xi)*eta*(1.0 - eta)*D_XI_DY
           + xi*xi*(1 - 2*eta)*D_ETA_DY;
  }

  void lambda_mp(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1.0 - xi)*(1.0 - xi)*eta*eta;
  }

  void gradient_lambda_mp(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi - 2)*eta*eta*D_XI_DX
           + (1.0 - xi)*(1.0 - xi)*(2*eta)*D_ETA_DX;
    val[1] = (2*xi - 2)*eta*eta*D_XI_DY
           + (1.0 - xi)*(1.0 - xi)*(2*eta)*D_ETA_DY;
  }

  void lambda_0p(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*(1.0 - xi)*eta*eta;
  }

  void gradient_lambda_0p(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (1 - 2*xi)*eta*eta*D_XI_DX
           + xi*(1.0 - xi)*(2*eta)*D_ETA_DX;
    val[1] = (1 - 2*xi)*eta*eta*D_XI_DY
           + xi*(1.0 - xi)*(2*eta)*D_ETA_DY;
  }

  void lambda_pp(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = xi*xi*eta*eta;
  }

  void gradient_lambda_pp(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    COMMON_PART;
    val[0] = (2*xi)*eta*eta*D_XI_DX
           + xi*xi*(2*eta)*D_ETA_DX;
    val[1] = (2*xi)*eta*eta*D_XI_DY
           + xi*xi*(2*eta)*D_ETA_DY;
  }

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
