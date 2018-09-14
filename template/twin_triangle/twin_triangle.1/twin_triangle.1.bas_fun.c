/*************************************************************************
 * 	twin_tirangle.1.bas_fun.c : by R.Lie
 */

#include <math.h>

#define ZERO (1.0e-4)
#define AREA(p0, p1, p2) \
  ((p1[0] - p0[0])*(p2[1] - p0[1]) - \
   (p1[1] - p0[1])*(p2[0] - p0[0]))
#define NORMAL(normal, p0, p1, vol) \
  { normal[0] = (p0[1] - p1[1])/vol; \
    normal[1] = (p1[0] - p0[0])/vol; \
  }
#define NORMAL_P(normal, p0, p1, vol) \
  { normal[0] += (p0[1] - p1[1])/vol; \
    normal[1] += (p1[0] - p0[0])/vol; \
  }

#ifdef __cplusplus
extern "C" {
#endif

void lambda_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area, flag = AREA(p, v[2], v[0]);
  if (flag >= 0) {
    val[0] = AREA(p, v[1], v[2]);
    area = AREA(v[0], v[1], v[2]);
    val[0] /= area;
  }
  else {
    val[0] = AREA(p, v[2], v[3]);
    area = AREA(v[0], v[2], v[3]);
    val[0] /= area;
  }
}

void lambda_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  val[0] = AREA(p, v[2], v[0]);
  if (val[0] < 0) val[0] = 0.0;
  val[0] /= area;	
}

void lambda_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area, flag = AREA(p, v[2], v[0]);
  if (flag >= 0) {
    val[0] = AREA(p, v[0], v[1]);
    area = AREA(v[0], v[1], v[2]);
    val[0] /= area;
  }
  else {
    val[0] = AREA(p, v[3], v[0]);
    area = AREA(v[0], v[2], v[3]);
    val[0] /= area;
  }
}

void lambda_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[2], v[3]);
  val[0] = AREA(p, v[0], v[2]);
  if (val[0] < 0) val[0] = 0.0;
  val[0] /= area;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area1, area2 = AREA(v[0], v[2], v[3]);
  double flag = AREA(p, v[2], v[0]);
  if (fabs(flag) < ZERO*area2) {
    area1 = AREA(v[0], v[1], v[2]);
    NORMAL(val, v[1], v[2], area1);
    NORMAL_P(val, v[2], v[3], area2);
    val[0] *= 0.5; val[1] *= 0.5;
  } else if (flag > 0) {
    area1 = AREA(v[0], v[1], v[2]);
    NORMAL(val, v[1], v[2], area1);
  } else {
    NORMAL(val, v[2], v[3], area2);
  }
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  double flag = AREA(p, v[2], v[0]);
  NORMAL(val, v[2], v[0], area);
  if (fabs(flag) <= ZERO*area) {
    val[0] *= 0.5; val[1] *= 0.5;
  } else if (flag < 0) {
    val[0] = 0.0; val[1] = 0.0;
  }
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area1, area2 = AREA(v[0], v[2], v[3]);
  double flag = AREA(p, v[2], v[0]);
  if (fabs(flag) < ZERO*area2) {
    area1 = AREA(v[0], v[1], v[2]);
    NORMAL(val, v[0], v[1], area1);
    NORMAL_P(val, v[3], v[0], area2);
    val[0] *= 0.5; val[1] *= 0.5;
  } else if (flag > 0) {
    area1 = AREA(v[0], v[1], v[2]);
    NORMAL(val, v[0], v[1], area1);
  } else {
    NORMAL(val, v[3], v[0], area2);
  }
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[2], v[3]);
  double flag = AREA(p, v[0], v[2]);
  NORMAL(val, v[0], v[2], area);
  if (fabs(flag) < ZERO*area) {
    val[0] *= 0.5; val[1] *= 0.5;
  } else if (flag < 0) {
    val[0] = 0.0; val[1] = 0.0;
  }
}

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
