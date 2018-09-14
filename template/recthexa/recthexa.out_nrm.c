/*=============================================================================
#     FileName: recthexa.out_nrm.c
#         Desc: 
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2014-08-13 11:38:48
#      History:
=============================================================================*/

#include <assert.h>
#include <math.h>
#include <assert.h>
#include <math.h>
#include <malloc.h>

//p 是global 坐标, v 是六面体的六个点的全局坐标
void out_normal(const double * p, const double ** gv, int i, double * n)
{
  double area;
  int fidx[6][4]={{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
  int j,k;

  double a[3],b[3];  
  a[0] = gv[fidx[i][1]][0] - gv[fidx[i][0]][0];  
  a[1] = gv[fidx[i][1]][1] - gv[fidx[i][0]][1];  
  a[2] = gv[fidx[i][1]][2] - gv[fidx[i][0]][2];  
  b[0] = gv[fidx[i][2]][0] - gv[fidx[i][0]][0];  
  b[1] = gv[fidx[i][2]][1] - gv[fidx[i][0]][1];  
  b[2] = gv[fidx[i][2]][2] - gv[fidx[i][0]][2];  

  n[0] = a[1]*b[2] - a[2]*b[1];
  n[1] = -a[0]*b[2] + a[2]*b[0];
  n[2] = a[0]*b[1] - a[1]*b[0];

  area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

  n[0] /= area;
  n[1] /= area;
  n[2] /= area;

  return;
}


/**
 * end of file
 *
 */

/**
 * g++ -o a.exe recthexa.out_nrm.c 
s.c -lm -g
编译运行，验证了例子中的外法向计算正确。

int main()
{
  int i,j,k;
  double gp[3]={0,0,0};
double dx=2, dy=2 , dz=2;
 double gv1[8][3]
	={ 
	  { 0, 0, 0 }, 
	  { dx, 0, 0}, 
	  { dx,dy, 0}, 
	  { 0,dy,0}, 
	  { 0, 0, dz }, 
	  { dx, 0, dz}, 
	  { dx,dy, dz}, 
	  { 0,dy,dz} 
	};

 //绕z 轴旋转一个角度试试看
 double cos_1 = sqrt(2.)/2., sin_1 = sqrt(1-cos_1*cos_1); 
 for(i=0;i<8;i++){
   double tmp[2];
   tmp[0] = gv1[i][0]*cos_1 - gv1[i][1]*sin_1;
   tmp[1] = gv1[i][0]*sin_1 + gv1[i][1]*cos_1;
   gv1[i][0]=tmp[0];
   gv1[i][1]=tmp[1];
 }

  const double **gv;

  gv = (const double**)malloc(sizeof(double*)*8);
 
  for(i=0;i<8;i++){
	gv[i]=gv1[i];
  }



  double n[3] ;

  for(i=0;i<6;i++){
  out_normal(gp, gv, i, n);
  printf("face %d: n=(%lf,%lf,%lf)\n",i,n[0],n[1],n[2]);
  }


  return 0;

}
*/

/**
 * end of file
 *
 */
