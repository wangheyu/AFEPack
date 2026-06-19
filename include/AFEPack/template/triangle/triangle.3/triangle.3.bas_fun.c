/**
 * @file   triangle.3.bas_fun.c
 * @author Robert Lie
 * @date   Thu Sep 20 15:41:58 2007
 * 
 * @brief  Wang Duan and Ruo Li
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
  val[0] = lambda[0]*(3*lambda[0] - 1.0)*(3*lambda[0] - 2.0)/2.0;
}

void phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[1]*(3*lambda[1] - 1.0)*(3*lambda[1] - 2.0)/2.0;
}

void phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[2]*(3*lambda[2] - 1.0)*(3*lambda[2] - 2.0)/2.0;
}

void phi_12(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[0] - 1.0)/2.0;
}

void phi_21(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[1] - 1.0)/2.0;
}

void phi_23(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[1] - 1.0)/2.0;
}

void phi_32(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
}

void phi_13(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[0] - 1.0)/2.0;
}

void phi_31(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
}

void phi_123(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 27.0*lambda[0]*lambda[1]*lambda[2];
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = lambda[0]*(3*lambda[0] - 1.0)*(3*lambda[0] - 2.0)/2.0;
  val[0] = (dlambda0_dx*(3*lambda[0] - 1.0)*(3*lambda[0] - 2.0)/2.0 +
            lambda[0]*3*dlambda0_dx*(3*lambda[0] - 2.0)/2.0 +
            lambda[0]*(3*lambda[0] - 1.0)*3*dlambda0_dx/2.0);
  val[1] = (dlambda0_dy*(3*lambda[0] - 1.0)*(3*lambda[0] - 2.0)/2.0 +
            lambda[0]*3*dlambda0_dy*(3*lambda[0] - 2.0)/2.0 +
            lambda[0]*(3*lambda[0] - 1.0)*3*dlambda0_dy/2.0);
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = lambda[1]*(3*lambda[1] - 1.0)*(3*lambda[1] - 2.0)/2.0;
  val[0] = (dlambda1_dx*(3*lambda[1] - 1.0)*(3*lambda[1] - 2.0)/2.0 + 
            lambda[1]*3*dlambda1_dx*(3*lambda[1] - 2.0)/2.0 +
            lambda[1]*(3*lambda[1] - 1.0)*3*dlambda1_dx/2.0);
  val[1] = (dlambda1_dy*(3*lambda[1] - 1.0)*(3*lambda[1] - 2.0)/2.0 +
            lambda[1]*3*dlambda1_dy*(3*lambda[1] - 2.0)/2.0 +
            lambda[1]*(3*lambda[1] - 1.0)*3*dlambda1_dy/2.0);
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = lambda[2]*(3*lambda[2] - 1.0)*(3*lambda[2] - 2.0)/2.0;
  val[0] = (dlambda2_dx*(3*lambda[2] - 1.0)*(3*lambda[2] - 2.0)/2.0 +
            lambda[2]*3*dlambda2_dx*(3*lambda[2] - 2.0)/2.0 +
            lambda[2]*(3*lambda[2] - 1.0)*3*dlambda2_dx/2.0);
  val[1] = (dlambda2_dy*(3*lambda[2] - 1.0)*(3*lambda[2] - 2.0)/2.0 + 
            lambda[2]*3*dlambda2_dy*(3*lambda[2] - 2.0)/2.0 +
            lambda[2]*(3*lambda[2] - 1.0)*3*dlambda2_dy/2.0);
}

void gradient_phi_12(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[0] - 1.0)/2.0;
  val[0]  = (9.0*dlambda0_dx*lambda[1]*(3*lambda[0] - 1.0)/2.0 + 
             9.0*lambda[0]*dlambda1_dx*(3*lambda[0] - 1.0)/2.0 +
             9.0*lambda[0]*lambda[1]*3*dlambda0_dx/2.0);
  val[1]  = (9.0*dlambda0_dy*lambda[1]*(3*lambda[0] - 1.0)/2.0 + 
             9.0*lambda[0]*dlambda1_dy*(3*lambda[0] - 1.0)/2.0 +
             9.0*lambda[0]*lambda[1]*3*dlambda0_dy/2.0);
}

void gradient_phi_21(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[1] - 1.0)/2.0;
  val[0]  = (9.0*dlambda0_dx*lambda[1]*(3*lambda[1] - 1.0)/2.0 +
             9.0*lambda[0]*dlambda1_dx*(3*lambda[1] - 1.0)/2.0 +
             9.0*lambda[0]*lambda[1]*3*dlambda1_dx/2.0);
  val[1]  = (9.0*dlambda0_dy*lambda[1]*(3*lambda[1] - 1.0)/2.0 + 
             9.0*lambda[0]*dlambda1_dy*(3*lambda[1] - 1.0)/2.0 +
             9.0*lambda[0]*lambda[1]*3*dlambda1_dy/2.0);
}

void gradient_phi_23(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[1] - 1.0)/2.0;
  val[0] = (9.0*dlambda1_dx*lambda[2]*(3*lambda[1] - 1.0)/2.0 +
            9.0*lambda[1]*dlambda2_dx*(3*lambda[1] - 1.0)/2.0 +
            9.0*lambda[1]*lambda[2]*3*dlambda1_dx/2.0);
  val[1] = (9.0*dlambda1_dy*lambda[2]*(3*lambda[1] - 1.0)/2.0 + 
            9.0*lambda[1]*dlambda2_dy*(3*lambda[1] - 1.0)/2.0 +
            9.0*lambda[1]*lambda[2]*3*dlambda1_dy/2.0);
}

void gradient_phi_32(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
  val[0] = (9.0*dlambda1_dx*lambda[2]*(3*lambda[2] - 1.0)/2.0 + 
            9.0*lambda[1]*dlambda2_dx*(3*lambda[2] - 1.0)/2.0 +
            9.0*lambda[1]*lambda[2]*3*dlambda2_dx/2.0);
  val[1] = (9.0*dlambda1_dy*lambda[2]*(3*lambda[2] - 1.0)/2.0 +
            9.0*lambda[1]*dlambda2_dy*(3*lambda[2] - 1.0)/2.0 +
            9.0*lambda[1]*lambda[2]*3*dlambda2_dy/2.0);

}

void gradient_phi_13(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  /* val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[0] - 1.0)/2.0; */
  val[0] = (9.0*dlambda0_dx*lambda[2]*(3*lambda[0] - 1.0)/2.0 +
            9.0*dlambda2_dx*lambda[0]*(3*lambda[0] - 1.0)/2.0 +
            9.0*lambda[0]*lambda[2]*3*dlambda0_dx/2.0);
  val[1] = (9.0*dlambda0_dy*lambda[2]*(3*lambda[0] - 1.0)/2.0 +
            9.0*dlambda2_dy*lambda[0]*(3*lambda[0] - 1.0)/2.0 +
            9.0*lambda[0]*lambda[2]*3*dlambda0_dy/2.0);
}

void gradient_phi_31(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //	val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
  val[0] = (9.0*dlambda0_dx*lambda[2]*(3*lambda[2] - 1.0)/2.0 +
            9.0*dlambda2_dx*lambda[0]*(3*lambda[2] - 1.0)/2.0 +
            9.0*lambda[0]*lambda[2]*3*dlambda2_dx/2.0);
  val[1] = (9.0*dlambda0_dy*lambda[2]*(3*lambda[2] - 1.0)/2.0 +
            9.0*dlambda2_dy*lambda[0]*(3*lambda[2] - 1.0)/2.0 +
            9.0*lambda[0]*lambda[2]*3*dlambda2_dy/2.0);
}

void gradient_phi_123(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 27.0*(dlambda0_dx*lambda[1]*lambda[2] +
                 dlambda1_dx*lambda[2]*lambda[0] + 
                 dlambda2_dx*lambda[1]*lambda[0]);
  val[1] = 27.0*(dlambda0_dy*lambda[1]*lambda[2] +
                 dlambda1_dy*lambda[2]*lambda[0] + 
                 dlambda2_dy*lambda[1]*lambda[0]);
}

/**
 * end of file
 * 
 */

