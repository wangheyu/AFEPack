/*=============================================================================
#     FileName: recthexa.crd_trs.c
#         Desc:
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2014-08-13 11:01:46
#      History:
=============================================================================*/

#include <math.h>
#include <stdio.h>
#include <malloc.h>


#define _FUZHI_A_ \
  a[0][0] = gv[1][0]-gv[0][0]; \
  a[0][1] = gv[1][1]-gv[0][1]; \
  a[0][2] = gv[1][2]-gv[0][2]; \
  a[1][0] = gv[3][0]-gv[0][0]; \
  a[1][1] = gv[3][1]-gv[0][1]; \
  a[1][2] = gv[3][2]-gv[0][2]; \
  a[2][0] = gv[4][0]-gv[0][0]; \
  a[2][1] = gv[4][1]-gv[0][1]; \
  a[2][2] = gv[4][2]-gv[0][2];

double det( const double * v1, 
		const double * v2, 
		const double * v3)
{
	return v1[0]*v2[1]*v3[2] 
		+ v1[1]*v2[2]*v3[0]
		+ v1[2]*v2[0]*v3[1]
		- v1[0]*v2[2]*v3[1]
		- v1[1]*v2[0]*v3[2]
		- v1[2]*v2[1]*v3[0];
}

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
	double a[4];
	a[0] = .5*(-1.-lp[0]-lp[1]-lp[2]);
	a[1] = .5*(1.+lp[0]);
	a[2] = .5*(1.+lp[1]);
	a[3] = .5*(1.+lp[2]);

	gp[0] = a[0]*gv[0][0] + a[1]*gv[1][0] + a[2]*gv[3][0] + a[3]*gv[4][0]; 
	gp[1] = a[0]*gv[0][1] + a[1]*gv[1][1] + a[2]*gv[3][1] + a[3]*gv[4][1]; 
	gp[2] = a[0]*gv[0][2] + a[1]*gv[1][2] + a[2]*gv[3][2] + a[3]*gv[4][2]; 
};

void global_to_local(const double * gp,
                     const double ** lv,
                     const double ** gv,
                     double * lp)
{
  double det0, det1, det2, det3;
  double a[4][3];
  _FUZHI_A_

  a[3][0] = 2.*gp[0] - (gv[1][0]+gv[3][0]+gv[4][0]-gv[0][0]);
  a[3][1] = 2.*gp[1] - (gv[1][1]+gv[3][1]+gv[4][1]-gv[0][1]);
  a[3][2] = 2.*gp[2] - (gv[1][2]+gv[3][2]+gv[4][2]-gv[0][2]);

  det0 =  det(a[0],a[1],a[2]);
  lp[0] = det(a[3],a[1],a[2])/det0;
  lp[1] = det(a[0],a[3],a[2])/det0;
  lp[2] = det(a[0],a[1],a[3])/det0;

};

double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv)
{
  double a[3][3];
  _FUZHI_A_

  return det(a[0],a[1],a[2])/8.0;
}

double global_to_local_jacobian(const double * gp,
                                const double ** lv,
                                const double ** gv)
{
  double a[3][3];
  _FUZHI_A_
  return 8.0/det(a[0],a[1],a[2]);
}

void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial)
{
  partial[0][0] = (gv[1][0] - gv[0][0])/2.0;
  partial[1][0] = (gv[3][0] - gv[0][0])/2.0;
  partial[2][0] = (gv[4][0] - gv[0][0])/2.0;	
  partial[0][1] = (gv[1][1] - gv[0][1])/2.0;
  partial[1][1] = (gv[3][1] - gv[0][1])/2.0;
  partial[2][1] = (gv[4][1] - gv[0][1])/2.0;	
  partial[0][2] = (gv[1][2] - gv[0][2])/2.0;
  partial[1][2] = (gv[3][2] - gv[0][2])/2.0;
  partial[2][2] = (gv[4][2] - gv[0][2])/2.0;	
};

void global_to_local_partial(const double * gp,
                             const double ** lv,
                             const double ** gv,
                             double **partial)
{
  double a[3][3];
  double det0;
  _FUZHI_A_
  det0 =  det(a[0],a[1],a[2]);
  partial[0][0] = 2.0*(a[1][1]*a[2][2]-a[1][2]*a[2][1])/det0; 
  partial[0][1] = -2.0*(a[1][0]*a[2][2]-a[2][0]*a[1][2])/det0; 
  partial[0][2] = 2.0*(a[1][0]*a[2][1]-a[2][0]*a[1][1])/det0; 
   
  partial[1][0] = -2.0*(a[0][1]*a[2][2]-a[0][2]*a[2][1])/det0; 
  partial[1][1] = 2.0*(a[0][0]*a[2][2]-a[2][0]*a[0][2])/det0; 
  partial[1][2] = -2.0*(a[0][0]*a[2][1]-a[2][0]*a[0][1])/det0; 


  partial[2][0] = 2.0*(a[0][1]*a[1][2]-a[1][1]*a[0][2])/det0; 
  partial[2][1] = -2.0*(a[0][0]*a[1][2]-a[1][0]*a[0][2])/det0; 
  partial[2][2] = 2.0*(a[0][0]*a[1][1]-a[1][0]*a[0][1])/det0; 
};
 
/*
#include <malloc.h>
int main(){

  double gp[3]={0,0,0};
  int i,j;
double dx=2, dy=3 , dz=5;
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
 double lv1[8][3]
	={ 
	  { -1,-1,-1}, 
	  { 1,-1,-1}, 
	  { 1,1,-1}, 
	  { -1,1,-1}, 
	  { -1,-1,1}, 
	  { 1,-1,1}, 
	  { 1,1,1}, 
	  { -1,1,1} 
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

  double lp[3];

  const double **gv;
  const double **lv; 
  double ** partial;
  partial = (double**)malloc(sizeof(double*)*3);
  for(i=0;i<3;i++){
	partial[i] =(double*)malloc(sizeof(double)*3);
  }

  gv = (const double**)malloc(sizeof(double*)*8);
  lv = (const double**)malloc(sizeof(double*)*8);
 
  for(i=0;i<8;i++){
	gv[i]=gv1[i];
	lv[i]=lv1[i];
  }


  printf("Input a global point gp\n");
  scanf("%lf%lf%lf",&gp[0],&gp[1],&gp[2]);
  printf("input gp = (%lf,%lf,%lf)\n",gp[0],gp[1],gp[2]);
  global_to_local(gp,lv,gv,lp);
  printf("lp = (%lf,%lf,%lf)\n",lp[0],lp[1],lp[2]);
  local_to_global(lp,lv,gv,gp);
  printf("output gp = (%lf,%lf,%lf)\n",gp[0],gp[1],gp[2]);

  printf("local_to_global_jacobian=%lf\n",local_to_global_jacobian(lp,lv,gv));
  printf("global_to_local_jacobian=%lf\n",global_to_local_jacobian(gp,lv,gv));

  global_to_local_partial(gp,lv,gv,partial);
  printf("global_to_local_parital=\n");
  for(i=0;i<3;i++){
	for(j=0;j<3;j++){
	  printf("%lf ",partial[i][j]);
	}
	printf("\n");
  }
  local_to_global_partial(lp,lv,gv,partial);
  printf("local_to_global_parital=\n");
  for(i=0;i<3;i++){
	for(j=0;j<3;j++){
	  printf("%lf ",partial[i][j]);
	}
	printf("\n");
  }
  return 0;

}
*/
//
// end of file
///////////////////////////////////////////////////////////////////////////
