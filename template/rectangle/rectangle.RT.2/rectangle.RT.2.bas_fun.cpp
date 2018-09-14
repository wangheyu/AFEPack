/*************************************************************************
 * 本文件为湘潭大学易年余同学提供
 */

#include <cmath>
#include <Miscellaneous.h>

/**
 * 本文件定义的是矩形单元上的一阶 Raviart-Thomas 单元的基函数
 * 及其梯度。
 *
 */

#define GAUSS_ELIMINATION                       \
  for (i = 0;i < 3;i ++) {                      \
    k = i;                                      \
    for (j = i+1;j < 4;j ++)                    \
      if (fabs(m[j][i]) > fabs(m[k][i])) k = j;	\
    if (k != i) {                               \
      for (j = i;j < 4;j ++) {                  \
        tmp = m[i][j];                          \
        m[i][j] = m[k][j];                      \
        m[k][j] = tmp;                          \
      }                                         \
      tmp = a[i][0];                            \
      a[i][0] = a[k][0];                        \
      a[k][0] = tmp;                            \
      tmp = a[i][1];                            \
      a[i][1] = a[k][1];                        \
      a[k][1] = tmp;                            \
    }                                           \
    for (j = i+1;j < 4;j ++) {                  \
      tmp = m[j][i]/m[i][i];                    \
      for (k = i+1;k < 4;k ++)                  \
        m[j][k] -= tmp*m[i][k];                 \
      a[j][0] -= tmp*a[i][0];                   \
      a[j][1] -= tmp*a[i][1];                   \
    }                                           \
  }                                             \
  a[3][0] /= m[3][3];                           \
  a[3][1] /= m[3][3];                           \
  for (i = 2;i >= 0;i --) {                     \
    for (j = i+1;j < 4;j ++) {                  \
      a[i][0] -= m[i][j]*a[j][0];               \
      a[i][1] -= m[i][j]*a[j][1];               \
    }                                           \
    a[i][0] /= m[i][i];                         \
    a[i][1] /= m[i][i];                         \
  }

#ifdef __cplusplus
extern "C" {
#endif

#define vector_length 2
#define vt nVector<vector_length,double>

  void lambda_01(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = -1.0/8 + 125.0/333*xi-1.0/4*eta+250.0/333*xi*eta+3.0/8*eta*eta-125.0/111*xi*eta*eta;
  }

  void lambda_02(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = -1.0/8 - 125.0/333*xi-1.0/4*eta-250.0/333*xi*eta+3.0/8*eta*eta+125.0/111*xi*eta*eta;
  }

  void lambda_11(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = -1.0/8 + 1.0/4*xi+125.0/333*eta-250.0/333*xi*eta+3.0/8*xi*xi-125.0/111*xi*xi*eta; val[1] = 0.0;
  }

  void lambda_12(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = -1.0/8 + 1.0/4*xi-125.0/333*eta+250.0/333*xi*eta+3.0/8*xi*xi+125.0/111*xi*xi*eta; val[1] = 0.0;
  }

  void lambda_21(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = -1.0/8 - 125.0/333*xi+1.0/4*eta+250.0/333*xi*eta+3.0/8*eta*eta+125.0/111*xi*eta*eta;
  }

  void lambda_22(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = -1.0/8 + 125.0/333*xi+1.0/4*eta-250.0/333*xi*eta+3.0/8*eta*eta-125.0/111*xi*eta*eta;
  }

  void lambda_31(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = -1.0/8 - 1.0/4*xi-125.0/333*eta-250.0/333*xi*eta+3.0/8*xi*xi+125.0/111*xi*xi*eta; val[1] = 0.0;
  }

  void lambda_32(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = -1.0/8 - 1.0/4*xi+125.0/333*eta+250.0/333*xi*eta+3.0/8*xi*xi-125.0/111*xi*xi*eta; val[1] = 0.0;
  }

  void lambda_41(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 3.0/8 - 3.0/8*xi*xi; val[1] = 0.0;
  }

  void lambda_42(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 9.0/8*eta - 9.0/8*xi*xi*eta; val[1] = 0.0;
  }

  void lambda_43(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = 3.0/8 - 3.0/8*eta*eta;
  }

  void lambda_44(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    val[0] = 0.0; val[1] = 9.0/8*xi - 9.0/8*xi*eta*eta;
  }

  void gradient_lambda_01(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = -(-125.0/333-250.0/333*eta+125.0/111*eta*eta)*dxi_dx
      -(1.0/4-250.0/333*xi-6.0/8*eta+250.0/111*xi*eta)*deta_dx;
    val[1][1] = -(-125.0/333-250.0/333*eta+125.0/111*eta*eta)*dxi_dy
      -(1.0/4-250.0/333*xi-6.0/8*eta+250.0/111*xi*eta)*deta_dy;
  }

  void gradient_lambda_02(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = -(125.0/333+250.0/333*eta-125.0/111*eta*eta)*dxi_dx
      -(1.0/4+250.0/333*xi-6.0/8*eta-250.0/111*xi*eta)*deta_dx;
    val[1][1] = -(125.0/333+250.0/333*eta-125.0/111*eta*eta)*dxi_dy
      -(1.0/4+250.0/333*xi-6.0/8*eta-250.0/111*xi*eta)*deta_dy;
  }

  void gradient_lambda_11(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = (1.0/4-250.0/333*eta+6.0/8*xi-250.0/111*xi*eta)*dxi_dx
      +(125.0/333-250.0/333*xi-125.0/111*xi*xi)*deta_dx; 
    val[1][0] = (1.0/4-250.0/333*eta+6.0/8*xi-250.0/111*xi*eta)*dxi_dy
      +(125.0/333-250.0/333*xi-125.0/111*xi*xi)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_12(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = (1.0/4+250.0/333*eta+6.0/8*xi+250.0/111*xi*eta)*dxi_dx
      +(-125.0/333+250.0/333*xi+125.0/111*xi*xi)*deta_dx; 
    val[1][0] = (1.0/4+250.0/333*eta+6.0/8*xi+250.0/111*xi*eta)*dxi_dy
      +(-125.0/333+250.0/333*xi+125.0/111*xi*xi)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_21(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = (-125.0/333+250.0/333*eta+125.0/111*eta*eta)*dxi_dx
      +(1.0/4+250.0/333*xi+6.0/8*eta+250.0/111*xi*eta)*deta_dx;
    val[1][1] = (-125.0/333+250.0/333*eta+125.0/111*eta*eta)*dxi_dy
      +(1.0/4+250.0/333*xi+6.0/8*eta+250.0/111*xi*eta)*deta_dy;
  }

  void gradient_lambda_22(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = (125.0/333-250.0/333*eta-125.0/111*eta*eta)*dxi_dx
      +(1.0/4-250.0/333*xi+6.0/8*eta-250.0/111*xi*eta)*deta_dx;
    val[1][1] = (125.0/333-250.0/333*eta-125.0/111*eta*eta)*dxi_dy
      +(1.0/4-250.0/333*xi+6.0/8*eta-250.0/111*xi*eta)*deta_dy;
  }

  void gradient_lambda_31(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = -(1.0/4+250.0/333*eta-6.0/8*xi-250.0/111*xi*eta)*dxi_dx
      -(125.0/333+250.0/333*xi-125.0/111*xi*xi)*deta_dx; 
    val[1][0] = -(1.0/4+250.0/333*eta-6.0/8*xi-250.0/111*xi*eta)*dxi_dy
      -(125.0/333+250.0/333*xi-125.0/111*xi*xi)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_32(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = -(1.0/4-250.0/333*eta-6.0/8*xi+250.0/111*xi*eta)*dxi_dx
      -(-125.0/333-250.0/333*xi+125.0/111*xi*xi)*deta_dx; 
    val[1][0] = -(1.0/4-250.0/333*eta-6.0/8*xi+250.0/111*xi*eta)*dxi_dy
      -(-125.0/333-250.0/333*xi+125.0/111*xi*xi)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_41(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = (-6.0/8*xi)*dxi_dx
      +(0.0)*deta_dx; 
    val[1][0] = (-6.0/8*xi)*dxi_dy
      +(0.0)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_42(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = (-18.0/8*xi*eta)*dxi_dx
      +(9.0/8 - 9.0/8*xi*xi)*deta_dx; 
    val[1][0] = (-18.0/8*xi*eta)*dxi_dy
      +(9.0/8 - 9.0/8*xi*xi)*deta_dy;
    val[0][1] = 0.0;
    val[1][1] = 0.0;
  }

  void gradient_lambda_43(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = (0.0)*dxi_dx
      +(-6.0/8*eta)*deta_dx;
    val[1][1] = (0.0)*dxi_dy
      +(-6.0/8*eta)*deta_dy;
  }

  void gradient_lambda_44(const double * p, const double ** v, void * value)
  {
    int i, j, k;
    double m[4][4], tmp, xi, eta;
    double a[4][2] = {{-1.0, -1.0},{1.0, -1.0},{1.0, 1.0},{-1.0, 1.0}};
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
    double dxi_dx=a[1][0] + a[3][0]*p[1];
    double dxi_dy=a[2][0] + a[3][0]*p[0];
    double deta_dx=a[1][1] + a[3][1]*p[1];
    double deta_dy=a[2][1] + a[3][1]*p[0];
    val[0][0] = 0; 
    val[1][0] = 0;
    val[0][1] = (9.0/8 - 9.0/8*eta*eta)*dxi_dx
      +(-18.0/8*xi*eta)*deta_dx;
    val[1][1] = (9.0/8 - 9.0/8*eta*eta)*dxi_dy
      +(-18.0/8*xi*eta)*deta_dy;
  }

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
