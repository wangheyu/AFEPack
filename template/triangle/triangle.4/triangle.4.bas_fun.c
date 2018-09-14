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

void phi_0(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0;
}

void phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0;
}

void phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0;
}

void phi_010(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0);
}

void phi_011(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0;
}

void phi_012(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[0]*lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0;
}

void phi_120(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0);
}

void phi_121(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0;
}

void phi_122(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[1]*lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0;
}

void phi_200(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 4.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0);
}

void phi_201(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0;
}

void phi_202(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 16.0*lambda[2]*lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0;
}

void phi_012_0(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 32.0*(4.*lambda[0] - 1.)*lambda[0]*lambda[1]*lambda[2];
}

void phi_012_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 32.0*(4.*lambda[1] - 1.)*lambda[0]*lambda[1]*lambda[2];
  //val[0]= 32.0*(4.*lambda[1] - 1.)*lambda[0]*lambda[1]*lambda[2];
}

void phi_012_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  val[0] = 32.0*(4.*lambda[2] - 1.)*lambda[0]*lambda[1]*lambda[2];
  //val[0]= 32.0*(4.*lambda[2] - 1.)*lambda[0]*lambda[1]*lambda[2];
}

void gradient_phi_0(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0;
  val[0] = (dlambda0_dx*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0 +
            lambda[0]*4.*dlambda0_dx*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0 +
            lambda[0]*(4.*lambda[0] - 1.0)*2.*dlambda0_dx*(4.*lambda[0] - 3.0)/3.0 +
	    lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*4.*dlambda0_dx/3.0);
  val[1] = (dlambda0_dy*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0 +
            lambda[0]*4.*dlambda0_dy*(2.*lambda[0] - 1.0)*(4.*lambda[0] - 3.0)/3.0 +
            lambda[0]*(4.*lambda[0] - 1.0)*2.*dlambda0_dy*(4.*lambda[0] - 3.0)/3.0 +
	    lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)*4.*dlambda0_dy/3.0);
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] = lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0;
  val[0] = (dlambda1_dx*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0 + 
            lambda[1]*4.*dlambda1_dx*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0 +
            lambda[1]*(4.*lambda[1] - 1.0)*2.*dlambda1_dx*(4.*lambda[1] - 3.0)/3.0 +
	    lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*4.*dlambda1_dx/3.0);
  val[1] = (dlambda1_dy*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0 + 
            lambda[1]*4.*dlambda1_dy*(2.*lambda[1] - 1.0)*(4.*lambda[1] - 3.0)/3.0 +
            lambda[1]*(4.*lambda[1] - 1.0)*2.*dlambda1_dy*(4.*lambda[1] - 3.0)/3.0 +
	    lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)*4.*dlambda1_dy/3.0);

}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0;
  val[0] = ( dlambda2_dx*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0+
             lambda[2]*4.*dlambda2_dx*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0 +
             lambda[2]*(4.*lambda[2] - 1.0)*2.*dlambda2_dx*(4.*lambda[2] - 3.0)/3.0 +
	     lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*4.*dlambda2_dx/3.0);
  val[1] =  ( dlambda2_dy*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0+
             lambda[2]*4.*dlambda2_dy*(2.*lambda[2] - 1.0)*(4.*lambda[2] - 3.0)/3.0 +
             lambda[2]*(4.*lambda[2] - 1.0)*2.*dlambda2_dy*(4.*lambda[2] - 3.0)/3.0 +
	     lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)*4.*dlambda2_dy/3.0);
}

void gradient_phi_010(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  4.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0);
  val[0]  = ( 4.0*dlambda0_dx*lambda[1]*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0) + 
              4.0*lambda[0]*dlambda1_dx*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0) +
              4.0*lambda[0]*lambda[1]*4.*dlambda0_dx*(4.*lambda[1] - 1.0) +
	      4.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*4.*dlambda1_dx);
  val[1]  = ( 4.0*dlambda0_dy*lambda[1]*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0) + 
              4.0*lambda[0]*dlambda1_dy*(4.*lambda[0] - 1.0)*(4.*lambda[1] - 1.0) +
              4.0*lambda[0]*lambda[1]*4.*dlambda0_dy*(4.*lambda[1] - 1.0) +
	      4.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*4.*dlambda1_dy);
}

void gradient_phi_011(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda0_dx*lambda[1]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[0]*dlambda1_dx*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[0]*lambda[1]*4.*dlambda0_dx *(2.*lambda[0] - 1.0)/3.0 +
	      16.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*2.*dlambda0_dx/3.0);
  val[1]  = ( 16.0*dlambda0_dy*lambda[1]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[0]*dlambda1_dy*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[0]*lambda[1]*4.*dlambda0_dy *(2.*lambda[0] - 1.0)/3.0 +
	      16.0*lambda[0]*lambda[1]*(4.*lambda[0] - 1.0)*2.*dlambda0_dy/3.0);
}


void gradient_phi_012(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[0]*lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda0_dx*lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[0]*dlambda1_dx*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[0]*lambda[1]*4.*dlambda1_dx *(2.*lambda[1] - 1.0)/3.0 +
	      16.0*lambda[0]*lambda[1]*(4.*lambda[1] - 1.0)*2.*dlambda1_dx/3.0);
  val[1]  = ( 16.0*dlambda0_dy*lambda[1]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[0]*dlambda1_dy*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[0]*lambda[1]*4.*dlambda1_dy *(2.*lambda[1] - 1.0)/3.0 +
	      16.0*lambda[0]*lambda[1]*(4.*lambda[1] - 1.0)*2.*dlambda1_dy/3.0);
}


void gradient_phi_120(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  4.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0);
  val[0]  = ( 4.0*dlambda1_dx*lambda[2]*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0) + 
              4.0*lambda[1]*dlambda2_dx*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0) +
              4.0*lambda[1]*lambda[2]*4.*dlambda1_dx*(4.*lambda[2] - 1.0) +
	      4.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*4.*dlambda2_dx);
  val[1]  = ( 4.0*dlambda1_dy*lambda[2]*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0) + 
              4.0*lambda[1]*dlambda2_dy*(4.*lambda[1] - 1.0)*(4.*lambda[2] - 1.0) +
              4.0*lambda[1]*lambda[2]*4.*dlambda1_dy*(4.*lambda[2] - 1.0) +
	      4.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*4.*dlambda2_dy);
}

void gradient_phi_121(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda1_dx*lambda[2]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[1]*dlambda2_dx*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[1]*lambda[2]*4.*dlambda1_dx *(2.*lambda[1] - 1.0)/3.0 +
	      16.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*2.*dlambda1_dx/3.0);
  val[1]  = ( 16.0*dlambda1_dy*lambda[2]*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[1]*dlambda2_dy*(4.*lambda[1] - 1.0)*(2.*lambda[1] - 1.0)/3.0 +
              16.0*lambda[1]*lambda[2]*4.*dlambda1_dy *(2.*lambda[1] - 1.0)/3.0 +
	      16.0*lambda[1]*lambda[2]*(4.*lambda[1] - 1.0)*2.*dlambda1_dy/3.0);
}


void gradient_phi_122(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[1]*lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda1_dx*lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[1]*dlambda2_dx*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[1]*lambda[2]*4.*dlambda2_dx *(2.*lambda[2] - 1.0)/3.0 +
	      16.0*lambda[1]*lambda[2]*(4.*lambda[2] - 1.0)*2.*dlambda2_dx/3.0);
  val[1]  = ( 16.0*dlambda1_dy*lambda[2]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[1]*dlambda2_dy*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[1]*lambda[2]*4.*dlambda2_dy *(2.*lambda[2] - 1.0)/3.0 +
	      16.0*lambda[1]*lambda[2]*(4.*lambda[2] - 1.0)*2.*dlambda2_dy/3.0);
}

void gradient_phi_200(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  4.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0);
  val[0]  = ( 4.0*dlambda2_dx*lambda[0]*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0) + 
              4.0*lambda[2]*dlambda0_dx*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0) +
              4.0*lambda[2]*lambda[0]*4.*dlambda2_dx*(4.*lambda[0] - 1.0) +
	      4.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*4.*dlambda0_dx);
  val[1]  = ( 4.0*dlambda2_dy*lambda[0]*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0) + 
              4.0*lambda[2]*dlambda0_dy*(4.*lambda[2] - 1.0)*(4.*lambda[0] - 1.0) +
              4.0*lambda[2]*lambda[0]*4.*dlambda2_dy*(4.*lambda[0] - 1.0) +
	      4.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*4.*dlambda0_dy);
}

void gradient_phi_201(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda2_dx*lambda[0]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[2]*dlambda0_dx*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[2]*lambda[0]*4.*dlambda2_dx *(2.*lambda[2] - 1.0)/3.0 +
	      16.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*2.*dlambda2_dx/3.0);
  val[1]  = ( 16.0*dlambda2_dy*lambda[0]*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[2]*dlambda0_dy*(4.*lambda[2] - 1.0)*(2.*lambda[2] - 1.0)/3.0 +
              16.0*lambda[2]*lambda[0]*4.*dlambda2_dy *(2.*lambda[2] - 1.0)/3.0 +
	      16.0*lambda[2]*lambda[0]*(4.*lambda[2] - 1.0)*2.*dlambda2_dy/3.0);
}


void gradient_phi_202(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0] =  16.0*lambda[2]*lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0;
  val[0]  = ( 16.0*dlambda2_dx*lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[2]*dlambda0_dx*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[2]*lambda[0]*4.*dlambda0_dx *(2.*lambda[0] - 1.0)/3.0 +
	      16.0*lambda[2]*lambda[0]*(4.*lambda[0] - 1.0)*2.*dlambda0_dx/3.0);
  val[1]  = ( 16.0*dlambda2_dy*lambda[0]*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[2]*dlambda0_dy*(4.*lambda[0] - 1.0)*(2.*lambda[0] - 1.0)/3.0 +
              16.0*lambda[2]*lambda[0]*4.*dlambda0_dy *(2.*lambda[0] - 1.0)/3.0 +
	      16.0*lambda[2]*lambda[0]*(4.*lambda[0] - 1.0)*2.*dlambda0_dy/3.0);
}

void gradient_phi_012_0(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0]= 32.0*(4.*lambda[0] - 1.)*lambda[0]*lambda[1]*lambda[2];
  val[0] = 32.0*(4.*dlambda0_dx*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[0] - 1.)*dlambda0_dx*lambda[1]*lambda[2] + 
                 (4.*lambda[0] - 1.)*lambda[0]*dlambda1_dx*lambda[2] +
		 (4.*lambda[0] - 1.)*lambda[0]*lambda[1]*dlambda2_dx);
  val[1] =  32.0*(4.*dlambda0_dy*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[0] - 1.)*dlambda0_dy*lambda[1]*lambda[2] + 
                 (4.*lambda[0] - 1.)*lambda[0]*dlambda1_dy*lambda[2] +
		 (4.*lambda[0] - 1.)*lambda[0]*lambda[1]*dlambda2_dy);
}

void gradient_phi_012_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0]= 32.0*(4.*lambda[1] - 1.)*lambda[0]*lambda[1]*lambda[2];
  val[0] = 32.0*(4.*dlambda1_dx*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[1] - 1.)*dlambda0_dx*lambda[1]*lambda[2] + 
                 (4.*lambda[1] - 1.)*lambda[0]*dlambda1_dx*lambda[2] +
		 (4.*lambda[1] - 1.)*lambda[0]*lambda[1]*dlambda2_dx);
  val[1] =  32.0*(4.*dlambda1_dy*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[1] - 1.)*dlambda0_dy*lambda[1]*lambda[2] + 
                 (4.*lambda[1] - 1.)*lambda[0]*dlambda1_dy*lambda[2] +
		 (4.*lambda[1] - 1.)*lambda[0]*lambda[1]*dlambda2_dy);
}

void gradient_phi_012_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[3], area;
  get_lambda(p, v, lambda, &area);
  //val[0]= 32.0*(4.*lambda[2] - 1.)*lambda[0]*lambda[1]*lambda[2];
  val[0] = 32.0*(4.*dlambda2_dx*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[2] - 1.)*dlambda0_dx*lambda[1]*lambda[2] + 
                 (4.*lambda[2] - 1.)*lambda[0]*dlambda1_dx*lambda[2] +
		 (4.*lambda[2] - 1.)*lambda[0]*lambda[1]*dlambda2_dx);
  val[1] =  32.0*(4.*dlambda2_dy*lambda[0]*lambda[1]*lambda[2] +
                 (4.*lambda[2] - 1.)*dlambda0_dy*lambda[1]*lambda[2] + 
                 (4.*lambda[2] - 1.)*lambda[0]*dlambda1_dy*lambda[2] +
		 (4.*lambda[2] - 1.)*lambda[0]*lambda[1]*dlambda2_dy);
}

/**
 * end of file
 * 
 */

