/*************************************************************************
* 	prism.1.D.bas_fun.c : by R.Lie
*/

#define CENTER_X \
  double x; \
  x = (v[0][0] + v[1][0] + v[2][0] \
	  +v[3][0] + v[4][0] + v[5][0])/6.0 ; 
#define CENTER_Y\
  double y; \
  y = (v[0][1] + v[1][1] + v[2][1] \
	  +v[3][1] + v[4][1] + v[5][1])/6.0 ; 
#define CENTER_Z \
  double z; \
  z = (v[0][2] + v[1][2] + v[2][2] \
	  +v[3][2] + v[4][2] + v[5][2])/6.0 ; 

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

/*
*  end of file
**************************************************************************/
