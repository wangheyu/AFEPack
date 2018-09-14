/**
 * @file   triangle.2.bas_fun.c
 * @author Robert Lie
 * @date   Tue Oct 23 12:48:47 2007
 * 
 * @brief  
 * 
 * 
 */

#define AREA(p0, p1, p2) (((p1)[0] - (p0)[0])*((p2)[1] - (p0)[1]) - \
                          ((p1)[1] - (p0)[1])*((p2)[0] - (p0)[0]))
#define dlambda0_dx ((v[1][1] - v[2][1])/area)
#define dlambda0_dy ((v[2][0] - v[1][0])/area)
#define dlambda1_dx ((v[2][1] - v[0][1])/area)
#define dlambda1_dy ((v[0][0] - v[2][0])/area)
#define dlambda2_dx ((v[0][1] - v[1][1])/area)
#define dlambda2_dy ((v[1][0] - v[0][0])/area)


void get_lambda(const double * p, const double ** v, double * lambda, double * area)
{
  area[0] = AREA(v[0], v[1], v[2]);
  lambda[0] = AREA(p, v[1], v[2])/area[0];
  lambda[1] = AREA(p, v[2], v[0])/area[0];
  lambda[2] = AREA(p, v[0], v[1])/area[0];
}

void phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[0]*(2*lambda[0] - 1.0);
}

void phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[1]*(2*lambda[1] - 1.0);
}

void phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[2]*(2*lambda[2] - 1.0);
}

void psi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[1]*lambda[2];
}

void psi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[0]*lambda[2];
}

void psi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[0]*lambda[1];
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = (4.0*lambda[0] - 1.0)*dlambda0_dx;
  val[1] = (4.0*lambda[0] - 1.0)*dlambda0_dy;
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = (4.0*lambda[1] - 1.0)*dlambda1_dx;
  val[1] = (4.0*lambda[1] - 1.0)*dlambda1_dy;
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = (4.0*lambda[2] - 1.0)*dlambda2_dx;
  val[1] = (4.0*lambda[2] - 1.0)*dlambda2_dy;
}

void gradient_psi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*(lambda[2]*dlambda1_dx + lambda[1]*dlambda2_dx);
  val[1] = 4.0*(lambda[2]*dlambda1_dy + lambda[1]*dlambda2_dy);
}

void gradient_psi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*(lambda[0]*dlambda2_dx + lambda[2]*dlambda0_dx);
  val[1] = 4.0*(lambda[0]*dlambda2_dy + lambda[2]*dlambda0_dy);
}

void gradient_psi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*(lambda[1]*dlambda0_dx + lambda[0]*dlambda1_dx);
  val[1] = 4.0*(lambda[1]*dlambda0_dy + lambda[0]*dlambda1_dy);
}

/**
 * end of file
 * 
 */

