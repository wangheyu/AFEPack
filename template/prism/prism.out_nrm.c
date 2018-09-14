/*=============================================================================
#     FileName: prism.out_nrm.c
#         Desc:  现在是求一个三维的global的点
出的面上的外法向，这个需要先把这个点对应的局部坐标找出来。 感觉这个地方
如果不是global 的点就好了， 因为实际上我们的积分点是局部的。
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-10 17:25:33
#      History:
=============================================================================*/

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
//p 是global 坐标, v 是三棱柱的六个点的全局坐标
void out_normal(const double * p, const double ** v, int i, double * n)
{
  double area;
  int fidx[5][4]={{0,2,1,1},{3,4,5,5},{1,2,5,4},{2,0,3,5},{0,1,4,3}};
  if (i==0 || i==1)
  {
    const double * vtx[3] = {v[fidx[i][0]], v[fidx[i][1]], v[fidx[i][2]]};
    double d[2][3] =
    {
      {vtx[1][0] - vtx[0][0], vtx[1][1] - vtx[0][1], vtx[1][2] - vtx[0][2]},
      {vtx[2][0] - vtx[0][0], vtx[2][1] - vtx[0][1], vtx[2][2] - vtx[0][2]}
    };
    n[0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
    n[1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
    n[2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
    area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    n[0] /= area;
    n[1] /= area;
    n[2] /= area;
    return;
  }
  else
  {
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
	for(j=0;j<3;j++) free(partial[j]);
	free(partial);
	
    return;
  }

}


/**
 * end of file
 *
 */

/**
 * g++ -o a.exe prism.out_nrm.c ../quadrilateral/quadrilateral.to3d.crd_tr
s.c -lm -g
编译运行，验证了例子中的外法向计算正确。 

int main()
{
  int i,j,k;
  const double ** gv;

  double tmp_gv[6][3]={{0,0,-1},{1,0,-1},{0,1,-1},
  {0,0,1},{1,0,1},{0,1,1}};

  gv = (const double **) malloc(6*sizeof(double*));
  for (i=0; i<6; i++) gv[i] = tmp_gv[i];


  double gp[5][3]={ 
	{0.2,0.1,-1},
	{0.2,0.1,1},
	{0.5,0.5,0.0},
	{0.0, 0.5, 0.0},
	{0.5, 0.0, 0.0}
  }
	;
  double n[3] ;

  for(i=0;i<5;i++){
  out_normal(gp[i], gv, i, n);
  printf("face %d: gp =(%lf,%lf,%lf), n=(%lf,%lf,%lf)\n",i,gp[i][0],gp[i][1],gp[i][2],n[0],n[1],n[2]);
  }


  return 0;

}
*/
