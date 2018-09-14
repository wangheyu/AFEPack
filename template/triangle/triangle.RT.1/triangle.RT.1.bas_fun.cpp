/**
 * tirangle.RT.1.bas_fun.cpp : by R.Lie
 * 
 */

#include <cmath>
#include <Miscellaneous.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * 三角形上 0 阶的 Raviart-Thomas 元的基函数及其导数的定义。其中
   * 插值点在第 \f$ i \f$ 条边上的基函数的表达式为
   *
   * \f[
   *    \vec{\lambda}_i = \left(\vec{x} - \vec{x}_i\right) l_i/2|K|
   * \f]
   *
   * 其中 \f$ \vec{x}_i \f$ 是第 \f$ i \f$ 条边对面的顶点的坐标，
   * \f$ |K| \f$ 是单元的面积，\f$ l_i \f$ 是第 \f$ i \f$ 条边的长
   * 度。这样的基函数的导数为
   *
   * \f[
   *    \nabla\vec{\lambda}_i = \left(\begin{array}{cc}
   *       \frac{l_i}{2|K|} & 0 \\
   *       0 & \frac{l_i}{2|K|} \end{array}\right)
   * \f]
   *
   * 我们这里将这种单元的描述文件命名为 triangle.RT.1.* 而不是
   * triangle.RT.0.* ，是因为这些基函数事实上是线性函数。
   * 
   * 因为这个文件中用到了类 nVector<2,double> ，必须使用 C++编译器
   * 进行编译。
   *
   */

#define vector_length 2
#define vt nVector<vector_length,double>

#define GET_AREA                                                               \
  double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])                        \
              - (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
#define GET_L(i, j)                                                            \
  double l0 = (v[i][0] - v[j][0]);                                             \
  double l1 = (v[i][1] - v[j][1]);                                             \
  double l = sqrt(l0*l0 + l1*l1);                                              \
  if (fabs(l0) > fabs(l1)) {                                                   \
    if (l0 < 0) l = -l;                                                        \
  } else {                                                                     \
    if (l1 < 0) l = -l;                                                        \
  }


void lambda_1(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  GET_AREA;
  GET_L(2, 1);
  val[0] = (p[0] - v[0][0])*l/area;
  val[1] = (p[1] - v[0][1])*l/area;
};

void lambda_2(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  GET_AREA;
  GET_L(0, 2);
  val[0] = (p[0] - v[1][0])*l/area;
  val[1] = (p[1] - v[1][1])*l/area;
};

void lambda_3(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  GET_AREA;
  GET_L(1, 0);
  val[0] = (p[0] - v[2][0])*l/area;
  val[1] = (p[1] - v[2][1])*l/area;
};

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  GET_AREA;
  GET_L(2, 1);
  area = l/area;
  val[0][0] = area; val[0][1] = 0.0;
  val[1][0] = 0.0; val[1][1] = area;
};

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  GET_AREA;
  GET_L(0, 2);
  area = l/area;
  val[0][0] = area; val[0][1] = 0.0;
  val[1][0] = 0.0; val[1][1] = area;
};

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  GET_AREA;
  GET_L(1, 0);
  area = l/area;
  val[0][0] = area; val[0][1] = 0.0;
  val[1][0] = 0.0; val[1][1] = area;
};

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
