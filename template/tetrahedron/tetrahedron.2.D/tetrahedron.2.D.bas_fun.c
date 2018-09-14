/*************************************************************************
 * 	tetrahedron.2.bas_fun.c : by D.Wang
 */

double get_volume(const double * v0,
                  const double * v1,
                  const double * v2,
                  const double * v3)
{
  return ((v1[0] - v0[0])*(v2[1] - v0[1])*(v3[2] - v0[2])
          + (v1[1] - v0[1])*(v2[2] - v0[2])*(v3[0] - v0[0])
          + (v1[2] - v0[2])*(v2[0] - v0[0])*(v3[1] - v0[1])
          - (v1[0] - v0[0])*(v2[2] - v0[2])*(v3[1] - v0[1])
          - (v1[1] - v0[1])*(v2[0] - v0[0])*(v3[2] - v0[2])
          - (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]));
}

void get_lambda(const double * p, const double ** v, double * lambda, double * volume)
{
  volume[0] = get_volume(v[0], v[1], v[2], v[3]);
  lambda[0] = get_volume(p, v[1], v[2], v[3])/volume[0];
  lambda[1] = get_volume(v[0], p, v[2], v[3])/volume[0];
  lambda[2] = get_volume(v[0], v[1], p, v[3])/volume[0];
  lambda[3] = get_volume(v[0], v[1], v[2], p)/volume[0];
}

#define co_det(v, m, n) (\
		((m%2==0)?-1.:1.) * (\
		  (v[(m+2)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
		* (v[(m+3)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		- (v[(m+2)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		* (v[(m+3)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
	) \
	)


void phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double volume = get_volume(v[0], v[1], v[2], v[3]);
  val[0] = get_volume(p, v[1], v[2], v[3]);
  val[0] /= volume;
  val[0] = val[0]*(2*val[0] - 1);
}

void phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double volume = get_volume(v[0], v[1], v[2], v[3]);
  val[0] = get_volume(v[0], p, v[2], v[3]);
  val[0] /= volume;
  val[0] = val[0]*(2*val[0] - 1);
}

void phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double volume = get_volume(v[0], v[1], v[2], v[3]);
  val[0] = get_volume(v[0], v[1], p, v[3]);
  val[0] /= volume;
  val[0] = val[0]*(2*val[0] - 1);
}

void phi_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double volume = get_volume(v[0], v[1], v[2], v[3]);
  val[0] = get_volume(v[0], v[1], v[2], p);
  val[0] /= volume;
  val[0] = val[0]*(2*val[0] - 1);
}

void psi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[0]*lambda[1];
	
}

void psi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[0]*lambda[2];
	
}

void psi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[0]*lambda[3];
	
}

void psi_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[1]*lambda[2];
	
}

void psi_5(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[1]*lambda[3];
	
}

void psi_6(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4.0*lambda[2]*lambda[3];
	
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = co_det(v, 0, 0)/volume * (4.*lambda[0] - 1.);
  val[1] = co_det(v, 0, 1)/volume * (4.*lambda[0] - 1.);
  val[2] = co_det(v, 0, 2)/volume * (4.*lambda[0] - 1.);
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = co_det(v, 1, 0)/volume * (4.*lambda[1] - 1.);
  val[1] = co_det(v, 1, 1)/volume * (4.*lambda[1] - 1.);
  val[2] = co_det(v, 1, 2)/volume * (4.*lambda[1] - 1.);
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = co_det(v, 2, 0)/volume * (4.*lambda[2] - 1.);
  val[1] = co_det(v, 2, 1)/volume * (4.*lambda[2] - 1.);
  val[2] = co_det(v, 2, 2)/volume * (4.*lambda[2] - 1.);
}

void gradient_phi_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = co_det(v, 3, 0)/volume * (4.*lambda[3] - 1.);
  val[1] = co_det(v, 3, 1)/volume * (4.*lambda[3] - 1.);
  val[2] = co_det(v, 3, 2)/volume * (4.*lambda[3] - 1.);
}

void gradient_psi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 1, 0)/volume * lambda[0] +
                co_det(v, 0, 0)/volume * lambda[1]);
  val[1] = 4 * (co_det(v, 1, 1)/volume * lambda[0] +
                co_det(v, 0, 1)/volume * lambda[1]);
  val[2] = 4 * (co_det(v, 1, 2)/volume * lambda[0] +
                co_det(v, 0, 2)/volume * lambda[1]);
}

void gradient_psi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 2, 0)/volume * lambda[0] +
                co_det(v, 0, 0)/volume * lambda[2]);
  val[1] = 4 * (co_det(v, 2, 1)/volume * lambda[0] +
                co_det(v, 0, 1)/volume * lambda[2]);
  val[2] = 4 * (co_det(v, 2, 2)/volume * lambda[0] +
                co_det(v, 0, 2)/volume * lambda[2]);
}

void gradient_psi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 3, 0)/volume * lambda[0] +
                co_det(v, 0, 0)/volume * lambda[3]);
  val[1] = 4 * (co_det(v, 3, 1)/volume * lambda[0] +
                co_det(v, 0, 1)/volume * lambda[3]);
  val[2] = 4 * (co_det(v, 3, 2)/volume * lambda[0] +
                co_det(v, 0, 2)/volume * lambda[3]);
}

void gradient_psi_4(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 1, 0)/volume * lambda[2] +
                co_det(v, 2, 0)/volume * lambda[1]);
  val[1] = 4 * (co_det(v, 1, 1)/volume * lambda[2] +
                co_det(v, 2, 1)/volume * lambda[1]);
  val[2] = 4 * (co_det(v, 1, 2)/volume * lambda[2] +
                co_det(v, 2, 2)/volume * lambda[1]);
}

void gradient_psi_5(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 1, 0)/volume * lambda[3] +
                co_det(v, 3, 0)/volume * lambda[1]);
  val[1] = 4 * (co_det(v, 1, 1)/volume * lambda[3] +
                co_det(v, 3, 1)/volume * lambda[1]);
  val[2] = 4 * (co_det(v, 1, 2)/volume * lambda[3] +
                co_det(v, 3, 2)/volume * lambda[1]);
}

void gradient_psi_6(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = 4 * (co_det(v, 2, 0)/volume * lambda[3] +
                co_det(v, 3, 0)/volume * lambda[2]);
  val[1] = 4 * (co_det(v, 2, 1)/volume * lambda[3] +
                co_det(v, 3, 1)/volume * lambda[2]);
  val[2] = 4 * (co_det(v, 2, 2)/volume * lambda[3] +
                co_det(v, 3, 2)/volume * lambda[2]);
}
/*
 *  end of file
 **************************************************************************/
