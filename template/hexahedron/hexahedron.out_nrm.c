/*=============================================================================
#     FileName: hexahedron.out_nrm.c
#         Desc: 
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-13 11:38:48
#      History:
=============================================================================*/

#include <assert.h>
#include <math.h>
#include <assert.h>
#include <math.h>
#include <malloc.h>

void to3d_local_to_global_partial(const double * lp,
                                  const double ** lv,
                                  const double ** gv,
                                  double ** partial);
void to3d_global_to_local(const double * gp,
                          const double ** lv,
                          const double ** gv,
                          double * lp);
//p 是global 坐标, v 是六面体的六个点的全局坐标
void out_normal(const double * p, const double ** v, int i, double * n)
{
  double area;
  int fidx[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
  int j,k;
  double tmp_lv[4][2];
  tmp_lv[0][0] = -1.0;
  tmp_lv[0][1] = -1.0;
  tmp_lv[1][0] = 1.0;
  tmp_lv[1][1] = -1.0;
  tmp_lv[2][0] = 1.0;
  tmp_lv[2][1] = 1.0;
  tmp_lv[3][0] = -1.0;
  tmp_lv[3][1] = 1.0;
  const double ** lv = (const double **) malloc(4*sizeof(double*));
  for (j=0; j<4; j++) lv[j] = tmp_lv[j];
  double * lp = (double *)malloc(2*sizeof(double));
  const double ** gv = (const double **) malloc(4*sizeof(double*));
  for (j=0; j<4; j++) gv[j] = v[fidx[i][j]];
  to3d_global_to_local(p, lv, gv, lp);

  double ** partial;
  partial = (double **) malloc(3*sizeof(double*) );
  for ( i = 0; i < 3 ; i++)
  {
    partial[i] = (double *) malloc(2*sizeof(double));
  }
  to3d_local_to_global_partial(lp,lv,gv,partial);
  n[0] = partial[1][0]*partial[2][1] - partial[2][0]*partial[1][1];
  n[1] = partial[2][0]*partial[0][1] - partial[0][0]*partial[2][1];
  n[2] = partial[0][0]*partial[1][1] - partial[1][0]*partial[0][1];
  area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= area;
  n[1] /= area;
  n[2] /= area;

  free(lv);
  free(gv);
  free(lp);
  for (j=0; j<3; j++) free(partial[j]);
  free(partial);

  return;
}


/**
 * end of file
 *
 */

/**
 * g++ -o a.exe hexahedron.out_nrm.c ../quadrilateral/quadrilateral.to3d.crd_tr
s.c -lm -g
编译运行，验证了例子中的外法向计算正确。
int main()
{
  int i,j,k;
  const double ** gv;

  double tmp_gv[8][3]={
	{-1.,-1.,-1},{1,-1,-1},{1,1,-1},{-1,1,-1}
	,{-1,-1,1},   {1,-1,1}, {1,1,1}, {-1,1,1}};
  

  gv = (const double **) malloc(8*sizeof(double*));
  for (i=0; i<8; i++) gv[i] = tmp_gv[i];


  double gp[6][3]={
	{0.2,0.1,-1},
	{0.2,0.1,1},
	{0.0,-1,0.0},
	{1.0, 0.5, 0.0},
	{0.5, 1.0, 0.0},
	{-1., 0, 0.0}
  }
	;
  double n[3] ;

  for(i=0;i<6;i++){
  out_normal(gp[i], gv, i, n);
  printf("face %d: gp =(%lf,%lf,%lf), n=(%lf,%lf,%lf)\n",i,gp[i][0],gp[i][1],gp[i][2],n[0],n[1],n[2]);
  }


  return 0;

}
*/

/**
 * end of file
 *
 */
