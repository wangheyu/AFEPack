/**
 * @file   tetrahedron.3.D.bas_fun.c
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Thu Aug 29 12:53:03 2013
 * 
 * @brief  
 * 
 * 
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


/// 四个线性基函数
void lambda_I(const double * p, const double ** v, double * val, int i)
{
  const double * q[4] = {v[0], v[1], v[2], v[3]};
  double volume = get_volume(q[0], q[1], q[2], q[3]);
  q[i] = p;
  val[0] = get_volume(q[0], q[1], q[2], q[3]);
  val[0] /= volume;
}

void grad_lambda_I(const double * p, const double ** v, double * val, int i)
{
  double volume = get_volume(v[0], v[1], v[2], v[3]);
  val[0] = co_det(v, i, 0)/volume;
  val[1] = co_det(v, i, 1)/volume;
  val[2] = co_det(v, i, 2)/volume;
}

#define LAMBDA(I)                                                       \
  void lambda_##I(const double * p, const double ** v, void * value)    \
  {                                                                     \
    double * val = (double *)value;                                     \
    lambda_I(p, v, val, I);                                             \
  }                                                                     \
                                                                        \
  void grad_lambda_##I(const double * p, const double ** v, void * value) \
  {                                                                     \
    double * val = (double *)value;                                     \
    grad_lambda_I(p, v, val, I);                                        \
  }                                                                     \


LAMBDA(0)
LAMBDA(1)
LAMBDA(2)
LAMBDA(3)

/// 六个二次基函数
void phi_IJ(const double * p, const double ** v, void * value, int i, int j)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = lambda[i]*lambda[j];
}

void grad_phi_IJ(const double * p, const double ** v, void * value, int i, int j)
{
  double * val = (double *)value;
  double lambda[4], volume, gi[4], gj[4];
  get_lambda(p, v, lambda, &volume);
  grad_lambda_I(p, v, gi, i);
  grad_lambda_I(p, v, gj, j);
  val[0] = gi[0]*lambda[j] + lambda[i]*gj[0];
  val[1] = gi[1]*lambda[j] + lambda[i]*gj[1];
  val[2] = gi[2]*lambda[j] + lambda[i]*gj[2];
}

#define PHI(I, J)                                                       \
  void phi_##I##J(const double * p, const double ** v, void * value)    \
  {                                                                     \
    double * val = (double *)value;                                     \
    phi_IJ(p, v, val, I, J);                                            \
  }                                                                     \
                                                                        \
  void grad_phi_##I##J(const double * p, const double ** v, void * value) \
  {                                                                     \
    double * val = (double *)value;                                     \
    grad_phi_IJ(p, v, val, I, J);                                       \
  }                                                                     \

PHI(0, 1)
PHI(0, 2)
PHI(0, 3)
PHI(1, 2)
PHI(1, 3)
PHI(2, 3)

/// 六个边上的三次基函数
void psi_IJ(const double * p, const double ** v, void * value, int i, int j)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = lambda[i]*lambda[j]*(lambda[i] - lambda[j]);
}

void grad_psi_IJ(const double * p, const double ** v, void * value, int i, int j)
{
  double * val = (double *)value;
  double lambda[4], volume, gi[4], gj[4];
  get_lambda(p, v, lambda, &volume);
  grad_lambda_I(p, v, gi, i);
  grad_lambda_I(p, v, gj, j);
  val[0] = (gi[0]*lambda[j]*(lambda[i] - lambda[j]) + 
            lambda[i]*gj[0]*(lambda[i] - lambda[j]) +
            lambda[i]*lambda[j]*(gi[0] - gj[0]));
  val[1] = (gi[1]*lambda[j]*(lambda[i] - lambda[j]) + 
            lambda[i]*gj[1]*(lambda[i] - lambda[j]) +
            lambda[i]*lambda[j]*(gi[1] - gj[1]));
  val[2] = (gi[2]*lambda[j]*(lambda[i] - lambda[j]) + 
            lambda[i]*gj[2]*(lambda[i] - lambda[j]) +
            lambda[i]*lambda[j]*(gi[2] - gj[2]));
}

#define PSI(I, J)                                                       \
  void psi_##I##J(const double * p, const double ** v, void * value)    \
  {                                                                     \
    double * val = (double *)value;                                     \
    psi_IJ(p, v, val, I, J);                                            \
  }                                                                     \
                                                                        \
  void grad_psi_##I##J(const double * p, const double ** v, void * value) \
  {                                                                     \
    double * val = (double *)value;                                     \
    grad_psi_IJ(p, v, val, I, J);                                       \
  }                                                                     \

PSI(0, 1)
PSI(0, 2)
PSI(0, 3)
PSI(1, 2)
PSI(1, 3)
PSI(2, 3)

/// 四个面上的三次数
void chi_IJK(const double * p, const double ** v, void * value, int i, int j, int k)
{
  double * val = (double *)value;
  double lambda[4], volume;
  get_lambda(p, v, lambda, &volume);
  val[0] = lambda[i]*lambda[j]*lambda[k];
}

void grad_chi_IJK(const double * p, const double ** v, void * value, int i, int j, int k)
{
  double * val = (double *)value;
  double lambda[4], volume, gi[4], gj[4], gk[4];
  get_lambda(p, v, lambda, &volume);
  grad_lambda_I(p, v, gi, i);
  grad_lambda_I(p, v, gj, j);
  grad_lambda_I(p, v, gk, k);
  val[0] = gi[0]*lambda[j]*lambda[k] + lambda[i]*gj[0]*lambda[k] + lambda[i]*lambda[j]*gk[0];
  val[1] = gi[1]*lambda[j]*lambda[k] + lambda[i]*gj[1]*lambda[k] + lambda[i]*lambda[j]*gk[1];
  val[2] = gi[2]*lambda[j]*lambda[k] + lambda[i]*gj[2]*lambda[k] + lambda[i]*lambda[j]*gk[2];
}

#define CHI(I, J, K)                                                      \
  void chi_##I##J##K(const double * p, const double ** v, void * value)    \
  {                                                                     \
    double * val = (double *)value;                                     \
    chi_IJK(p, v, val, I, J, K);                                           \
  }                                                                     \
                                                                        \
  void grad_chi_##I##J##K(const double * p, const double ** v, void * value) \
  {                                                                     \
    double * val = (double *)value;                                     \
    grad_chi_IJK(p, v, val, I, J, K);                                      \
  }                                                                     \


CHI(0, 1, 2)
CHI(0, 1, 3)
CHI(0, 2, 3)
CHI(1, 2, 3)

/*
 *  end of file
 **************************************************************************/
