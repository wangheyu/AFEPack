/*=============================================================================
#     FileName: oblong.to3d.crd_trs.c
#         Desc:
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-09 23:32:20
#      History:
=============================================================================*/


#include <math.h>
#include <malloc.h>
#include <stdio.h>

double det( const double * v1, 
		const double * v2)
{
	return v1[0]*v2[1] 
		- v1[1]*v2[0];
}

void to3d_local_to_global_partial(const double * lp,
                                  const double ** lv,
                                  const double ** gv,
                                  double ** partial);
void to3d_global_to_local(const double * gp,
                          const double ** lv,
                          const double ** gv,
                          double * lp);

void to3d_local_to_global(const double * lp,
                          const double ** lv,
                          const double ** gv,
                          double * gp)
{
	double a[3];
	a[0] = .5*(-lp[0]-lp[1]);
	a[1] = .5*(1.+lp[0]);
	a[2] = .5*(1.+lp[1]);

	gp[0] = a[0]*gv[0][0] + a[1]*gv[1][0] + a[2]*gv[3][0] ; 
	gp[1] = a[0]*gv[0][1] + a[1]*gv[1][1] + a[2]*gv[3][1] ; 
	gp[2] = a[0]*gv[0][2] + a[1]*gv[1][2] + a[2]*gv[3][2] ; 
}

void to3d_global_to_local(const double * gp,
                          const double ** lv,
                          const double ** gv,
                          double * lp)
{

}

double to3d_local_to_global_jacobian(const double * lp,
                                     const double ** lv,
                                     const double ** gv)
{
  double a[2][3];
  a[0][0]=gv[1][0]-gv[0][0];
  a[0][1]=gv[1][1]-gv[0][1];
  a[0][2]=gv[1][2]-gv[0][2];
  a[1][0]=gv[3][0]-gv[0][0];
  a[1][1]=gv[3][1]-gv[0][1];
  a[1][2]=gv[3][2]-gv[0][2];

  double n[3];
  n[0] = a[0][1]*a[1][2]-a[0][2]*a[1][1];
  n[1] = -a[0][0]*a[1][2]-a[0][2]*a[1][0];
  n[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];

  return sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])/4.0;


}

double to3d_global_to_local_jacobian(const double * gp,
                                     const double ** lv,
                                     const double ** gv)
{
  double a[2][3];
  a[0][0]=gv[1][0]-gv[0][0];
  a[0][1]=gv[1][1]-gv[0][1];
  a[0][2]=gv[1][2]-gv[0][2];
  a[1][0]=gv[3][0]-gv[0][0];
  a[1][1]=gv[3][1]-gv[0][1];
  a[1][2]=gv[3][2]-gv[0][2];

  double n[3];
  n[0] = a[0][1]*a[1][2]-a[0][2]*a[1][1];
  n[1] = -a[0][0]*a[1][2]-a[0][2]*a[1][0];
  n[2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];

  return 4.0/sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);



}


void to3d_local_to_global_partial(const double * lp,
                                  const double ** lv,
                                  const double ** gv,
                                  double ** partial)
{
}


/**
 * end of file
 *
 */

/**测试
 * g++ -lm -o oblong.to3d.crd_trs.exe oblong.to3d.crd_trs.c
 * 测试的目的是保证 lp 映射到 gp, 再映射回 lp 要保持一致
 * 测试结果证明了正确性，
 */ 
/*
int main()
{
  int i,j,k;
  const double ** lv, ** gv;
  double tmp_lv[4][2];
  tmp_lv[0][0] = -1.0;
  tmp_lv[0][1] = -1.0;
  tmp_lv[1][0] = 1.0;
  tmp_lv[1][1] = -1.0;
  tmp_lv[2][0] = 1.0;
  tmp_lv[2][1] = 1.0;
  tmp_lv[3][0] = -1.0;
  tmp_lv[3][1] = 1.0;

  lv = (const double **) malloc(4*sizeof(double*));
  for (i=0; i<4; i++) lv[i] = tmp_lv[i];

  double tmp_gv[4][3];
  tmp_gv[0][0] = 0.0;
  tmp_gv[0][1] = 0.0;
  tmp_gv[0][2] = 0.0;
  tmp_gv[2][0] = 2.0;
  tmp_gv[2][1] = 0.0;
  tmp_gv[2][2] = 0.0;
  tmp_gv[3][0] = 2.0;
  tmp_gv[3][1] = 0.0;
  tmp_gv[3][2] = 2.0;
  tmp_gv[1][0] = 0.0;
  tmp_gv[1][1] = 0.0;
  tmp_gv[1][2] = 2.0;

  gv = (const double **) malloc(4*sizeof(double*));
  for (i=0; i<4; i++) gv[i] = tmp_gv[i];


  double lp[2]= {-0.577350269189626, 0.57735026918962};
  double gp[3]= {0.0,0.0,0.0};

  to3d_local_to_global(lp,lv,gv,gp);
  printf("lp =(%lf,%lf), gp=(%lf,%lf,%lf)\n",lp[0],lp[1],gp[0],gp[1],gp[2]);

  to3d_global_to_local(gp,lv,gv,lp);
  printf("lp =(%lf,%lf), gp=(%lf,%lf,%lf)\n",lp[0],lp[1],gp[0],gp[1],gp[2]);
  double jac_g2l, jac_l2g;
  jac_g2l = to3d_global_to_local_jacobian(gp,lv,gv);
  jac_l2g = to3d_local_to_global_jacobian(lp,lv,gv);
  printf("global_to_local_jacobian =%lf,local_to_global_jacobian=%lf \n", jac_g2l, jac_l2g);


  return 0;

}
*/


