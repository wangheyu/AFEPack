/**
 * @file   triangle.1.bas_fun.c
 * @author Robert Lie
 * @date   Tue Oct 23 12:44:21 2007
 * 
 * @brief  
 * 
 * 
 */

#define AREA(p0, p1, p2) (((p1)[0] - (p0)[0])*((p2)[1] - (p0)[1]) - \
                          ((p1)[1] - (p0)[1])*((p2)[0] - (p0)[0]))
#define NORMAL(normal, p0, p1, vol) \
  { \
    (normal)[0] = ((p0)[1] - (p1)[1])/(vol); \
    (normal)[1] = ((p1)[0] - (p0)[0])/(vol); \
  }

void lambda_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  val[0] = AREA(p, v[1], v[2]);
  val[0] /= area;
}

void lambda_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  val[0] = AREA(v[0], p, v[2]);
  val[0] /= area;
}

void lambda_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  val[0] = AREA(v[0], v[1], p);
  val[0] /= area;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  NORMAL(val, v[1], v[2], area);
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  NORMAL(val, v[2], v[0], area);
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0], v[1], v[2]);
  NORMAL(val, v[0], v[1], area);
}

/**
 * end of file
 * 
 */

