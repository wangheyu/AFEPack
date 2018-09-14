/*=============================================================================
#     FileName: hexahedron.2.D.bas_fun.c
#         Desc: 
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2014-1-06 11:41:28
#      History:
=============================================================================*/

#define CENTER_X \
  double x; \
  x = (v[0][0] + v[1][0] + v[2][0] \
	  +v[3][0] + v[4][0] + v[5][0] \
	  +v[6][0] + v[7][0])/8.0 ; 
#define CENTER_Y\
  double y; \
  y = (v[0][1] + v[1][1] + v[2][1] \
	  +v[3][1] + v[4][1] + v[5][1] \
	  + v[6][1] + v[7][1])/8.0 ; 
#define CENTER_Z \
  double z; \
  z = (v[0][2] + v[1][2] + v[2][2] \
	  +v[3][2] + v[4][2] + v[5][2] \
	  +v[6][2] + v[7][2])/8.0 ; 

void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 1.0;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = p[0]-x; 
}

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = p[1]-y;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = p[2]-z;
}

void lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = (p[0]-x)*(p[0]-x);
}

void lambda_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = (p[1]-y)*(p[1]-y);
}

void lambda_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = (p[2]-z)*(p[2]-z);
}

void lambda_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	CENTER_Y
	val[0] = (p[0]-x)*(p[1]-y);
}

void lambda_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	CENTER_Z
	val[0] = (p[1]-y)*(p[2]-z);
}

void lambda_10(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	CENTER_X
	val[0] = (p[2]-z)*(p[0]-x);
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
	val[2] = 0.0;
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 1.0;
	val[1] = 0.0;
	val[2] = 0.0;
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 1.0;
	val[2] = 0.0;
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
	val[2] = 1.0;
}

void gradient_lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = 2.0*(p[0]-x);
	val[1] = 0.0;
	val[2] = 0.0;
}

void gradient_lambda_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = 0.0;
	val[1] = 2.0*(p[1]-y);
	val[2] = 0.0;
}

void gradient_lambda_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = 0.0;
	val[1] = 0.0;
	val[2] = 2.0*(p[2]-z);
}

void gradient_lambda_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	CENTER_Y
	val[0] = p[1]-y;
	val[1] = p[0]-x;
	val[2] = 0.0;
}

void gradient_lambda_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	CENTER_Z
	val[0] = 0.0;
	val[1] = p[2]-z;
	val[2] = p[1]-y;
}

void gradient_lambda_10(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	CENTER_X
	val[0] = p[2]-z;
	val[1] = 0.0;
	val[2] = p[0]-x;
}
/*
*  end of file
**************************************************************************/
