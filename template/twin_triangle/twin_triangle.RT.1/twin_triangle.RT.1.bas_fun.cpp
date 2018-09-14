/**
 * @file   twin_triangle.RT.1.bas_fun.cpp
 * @author Robert Lie
 * @date   Fri Oct 15 10:35:10 2004
 * 
 * @brief  basis functions on twin-triangle for R-T 0 element
 * 
 * 
 */

#include <cmath>
#include <Miscellaneous.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * 双生三角形上 0 阶的 Raviart-Thomas 元的基函数及其导数的定义。
   * 
   * 因为这个文件中用到了类 nVector<2,double> ，必须使用 C++编译器
   * 进行编译。
   *
   */

#define vector_length 2
#define vt nVector<vector_length,double>

#define ZERO (1.0e-6)
#define GET_AREA(v0, v1, v2)                                                   \
  ((v1[0] - v0[0])*(v2[1] - v0[1]) - (v1[1] - v0[1])*(v2[0] - v0[0]))
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
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    GET_L(1, 0);
    val[0] = (p[0] - v[2][0])*l/area;
    val[1] = (p[1] - v[2][1])*l/area;
  }
  if (p_area < area) {
    val[0] /= 2;
    val[1] /= 2;
  }
};

void lambda_2(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area ) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    GET_L(2, 1);
    val[0] = (p[0] - v[0][0])*l/area;
    val[1] = (p[1] - v[0][1])*l/area;
  }
  if (p_area < area) {
    val[0] /= 2;
    val[1] /= 2;
  }
};

void lambda_3(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = GET_AREA(v[0], v[2], v[3]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    GET_L(3, 2);
    val[0] = (p[0] - v[0][0])*l/area;
    val[1] = (p[1] - v[0][1])*l/area;
  }
  if (p_area > -area) {
    val[0] /= 2;
    val[1] /= 2;
  }
};

void lambda_4(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = GET_AREA(v[0], v[2], v[3]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    GET_L(0, 3);
    val[0] = (p[0] - v[2][0])*l/area;
    val[1] = (p[1] - v[2][1])*l/area;
  }
  if (p_area > -area) {
    val[0] /= 2;
    val[1] /= 2;
  }
};

void lambda_5(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    GET_L(0, 2);
    val[0] = (p[0] - v[1][0])*l/area;
    val[1] = (p[1] - v[1][1])*l/area;
  }
  else if (p_area < -area) {
    GET_L(2, 0);
    val[0] = (p[0] - v[3][0])*l/area;
    val[1] = (p[1] - v[3][1])*l/area;
  }
  else {
    {
      GET_L(0, 2);
      val[0] = (p[0] - v[1][0])*l/area;
      val[1] = (p[1] - v[1][1])*l/area;
    }
    {
      GET_L(2, 0);
      val[0] += (p[0] - v[3][0])*l/area;
      val[1] += (p[1] - v[3][1])*l/area;
    }
    val[0] /= 2;
    val[1] /= 2;
  }
};

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(1, 0);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  if (p_area < area) {
    val[0][0] /= 2;
    val[1][1] /= 2;
  }
};

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area ) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(2, 1);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  if (p_area < area) {
    val[0][0] /= 2;
    val[1][1] /= 2;
  }
};

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = GET_AREA(v[0], v[2], v[3]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(3, 2);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  if (p_area > -area) {
    val[0][0] /= 2;
    val[1][1] /= 2;
  }
};

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = GET_AREA(v[0], v[2], v[3]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(0, 3);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  if (p_area > -area) {
    val[0][0] /= 2;
    val[1][0] /= 2;
  }
};

void gradient_lambda_5(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = GET_AREA(v[0], v[1], v[2]);
  double p_area = GET_AREA(v[0], p, v[2])/ZERO;
  if (p_area > area) {
    GET_L(0, 2);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  else if (p_area < -area) {
    GET_L(2, 0);
    area = l/area;
    val[0][0] = area; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = area;
  }
  else {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
};

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
