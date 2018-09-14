/*=============================================================================
#     FileName: oblong.crd_trs.c
#         Desc:
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-09 16:06:33
#      History:
=============================================================================*/

#include <math.h>
#include <stdio.h>
#include <malloc.h>


#define _FUZHI_A_ \
  a[0][0] = gv[1][0]-gv[0][0]; \
  a[0][1] = gv[1][1]-gv[0][1]; \
  a[1][0] = gv[3][0]-gv[0][0]; \
  a[1][1] = gv[3][1]-gv[0][1]; 


void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial);


double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv);

void local_to_global(const double * lp,
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
};

void global_to_local(const double * gp,
                     const double ** lv,
                     const double ** gv,
                     double * lp)
{
  double det0, det1, det2, det3;
  double a[3][2];
  _FUZHI_A_

  a[2][0] = 2.*gp[0] - (gv[1][0]+gv[3][0]);
  a[2][1] = 2.*gp[1] - (gv[1][1]+gv[3][1]);

  det0 =  det(a[0],a[1]);
  lp[0] = det(a[2],a[1])/det0;
  lp[1] = det(a[0],a[2])/det0;

};

double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv)
{
  double a[2][2];
  _FUZHI_A_

  return det(a[0],a[1])/4.0;
}

double global_to_local_jacobian(const double * gp,
                                const double ** lv,
                                const double ** gv)
{
  double a[2][2];
  _FUZHI_A_
  return 4.0/det(a[0],a[1]);
}

void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial)
{
  partial[0][0] = (gv[1][0] - gv[0][0])/2.0;
  partial[1][0] = (gv[3][0] - gv[0][0])/2.0;
  partial[0][1] = (gv[1][1] - gv[0][1])/2.0;
  partial[1][1] = (gv[3][1] - gv[0][1])/2.0;
};

void global_to_local_partial(const double * gp,
                             const double ** lv,
                             const double ** gv,
                             double **partial)
{
  double a[3][3];
  double det0;
  _FUZHI_A_
  det0 =  det(a[0],a[1]);
  partial[0][0] = 2.0*a[1][1]/det0; 
  partial[0][1] = -2.0*a[1][0]/det0; 
   
  partial[1][0] = -2.0*a[0][1]/det0; 
  partial[1][1] = 2.0*a[0][0]/det0; 

};
/* 
#include <malloc.h>
int main(){

  double gp[3]={0,0,0};
  int i,j;
double dx=2, dy=2 ;
 double gv1[4][2]
	={ 
	  { 0, 0 }, 
	  { dx, 0}, 
	  { dx,dy}, 
	  { 0,dy}, 
	};
 double lv1[4][2]
	={ 
	  { -1,-1}, 
	  { 1,-1}, 
	  { 1,1}, 
	  { -1,1} 
	};

 //绕z 轴旋转一个角度试试看
 
 double cos_1 = sqrt(2.)/2., sin_1 = sqrt(1-cos_1*cos_1); 
 for(i=0;i<4;i++){
   double tmp[2];
   tmp[0] = gv1[i][0]*cos_1 - gv1[i][1]*sin_1;
   tmp[1] = gv1[i][0]*sin_1 + gv1[i][1]*cos_1;
   gv1[i][0]=tmp[0];
   gv1[i][1]=tmp[1];
 }

  double lp[2];

  const double **gv;
  const double **lv; 
  double ** partial;
  partial = (double**)malloc(sizeof(double*)*2);
  for(i=0;i<2;i++){
	partial[i] =(double*)malloc(sizeof(double)*22);
  }

  gv = (const double**)malloc(sizeof(double*)*4);
  lv = (const double**)malloc(sizeof(double*)*4);
 
  for(i=0;i<4;i++){
	gv[i]=gv1[i];
	lv[i]=lv1[i];
  }


  printf("Input a global point gp\n");
  scanf("%lf%lf",&gp[0],&gp[1]);
  printf("input gp = (%lf,%lf)\n",gp[0],gp[1]);
  global_to_local(gp,lv,gv,lp);
  printf("lp = (%lf,%lf)\n",lp[0],lp[1]);
  local_to_global(lp,lv,gv,gp);
  printf("output gp = (%lf,%lf)\n",gp[0],gp[1]);

  printf("local_to_global_jacobian=%lf\n",local_to_global_jacobian(lp,lv,gv));
  printf("global_to_local_jacobian=%lf\n",global_to_local_jacobian(gp,lv,gv));

  global_to_local_partial(gp,lv,gv,partial);
  printf("global_to_local_parital=\n");
  for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	  printf("%lf ",partial[i][j]);
	}
	printf("\n");
  }
  local_to_global_partial(lp,lv,gv,partial);
  printf("local_to_global_parital=\n");
  for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	  printf("%lf ",partial[i][j]);
	}
	printf("\n");
  }
  return 0;

}
*/
//
// end of file
/////////////////
