/*************************************************************************
* 	prism.1.D.bas_fun.c : by Tiao Lu 
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

void phi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 1.0;
}

void phi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = p[0]-x; 
}

void phi_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = p[1]-y;
}

void phi_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = p[2]-z;
}

void phi_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = (p[0]-x)*(p[0]-x);
}

void phi_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = (p[1]-y)*(p[1]-y);
}

void phi_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = (p[2]-z)*(p[2]-z);
}

void phi_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	CENTER_Y
	val[0] = (p[0]-x)*(p[1]-y);
}

void phi_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	CENTER_Z
	val[0] = (p[1]-y)*(p[2]-z);
}

void phi_10(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	CENTER_X
	val[0] = (p[2]-z)*(p[0]-x);
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.;
	val[1] = 0.;
	val[2] = 0.;
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 1.;
	val[1] = 0.;
	val[2] = 0.;
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.;
	val[1] = 1.;
	val[2] = 0.;
}

void gradient_phi_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.;
	val[1] = 0.;
	val[2] = 1.;
}

void gradient_phi_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	val[0] = 2.0*(p[0]-x);
	val[1] = 0.0;
	val[2] = 0.0;
}

void gradient_phi_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	val[0] = 0.0;
	val[1] = 2.0*(p[1]-y);
	val[2] = 0.0;
}

void gradient_phi_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Z
	val[0] = 0.0;
	val[1] = 0.0;
	val[2] = 2.0*(p[2]-z);
}

void gradient_phi_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_X
	CENTER_Y
	val[0] = p[1]-y;
	val[1] = p[0]-x;
	val[2] = 0.0;
}

void gradient_phi_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	CENTER_Y
	CENTER_Z
	val[0] = 0.0;
	val[1] = p[2]-z;
	val[2] = p[1]-y;
}

void gradient_phi_10(const double * p, const double ** v, void * value)
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
