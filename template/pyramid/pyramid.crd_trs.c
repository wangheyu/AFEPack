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
	double area;
	//lp[2] = 1.0 会出现零做除数的错误,因此返回顶点
	if(lp[2] == 1.0 ){
	  gp[0] = gv[4][0];
	  gp[1] = gv[4][1];
	  gp[2] = gv[4][2];
	  return;
	}

	area = (1.-lp[2])*(1.-lp[2])*4.;

	a[0] = (1.0-lp[2]-lp[0])*(1.0-lp[2]-lp[1])/area;
	a[1] = (1.0-lp[2]+lp[0])*(1.0-lp[2]-lp[1])/area; 
	a[2] = (1.0-lp[2]+lp[0])*(1.0-lp[2]+lp[1])/area;  
	a[3] = (1.0-lp[2]-lp[0])*(1.0-lp[2]+lp[1])/area; 

	gp[0] =
	  a[0]*(lp[2]*gv[4][0]+(1.-lp[2])*gv[0][0])+
	  a[1]*(lp[2]*gv[4][0]+(1.-lp[2])*gv[1][0])+
	  a[2]*(lp[2]*gv[4][0]+(1.-lp[2])*gv[2][0])+
	  a[3]*(lp[2]*gv[4][0]+(1.-lp[2])*gv[3][0]);
	gp[1] =
	  a[0]*(lp[2]*gv[4][1]+(1.-lp[2])*gv[0][1])+
	  a[1]*(lp[2]*gv[4][1]+(1.-lp[2])*gv[1][1])+
	  a[2]*(lp[2]*gv[4][1]+(1.-lp[2])*gv[2][1])+
	  a[3]*(lp[2]*gv[4][1]+(1.-lp[2])*gv[3][1]);
	gp[2] =
	  a[0]*(lp[2]*gv[4][2]+(1.-lp[2])*gv[0][2])+
	  a[1]*(lp[2]*gv[4][2]+(1.-lp[2])*gv[1][2])+
	  a[2]*(lp[2]*gv[4][2]+(1.-lp[2])*gv[2][2])+
	  a[3]*(lp[2]*gv[4][2]+(1.-lp[2])*gv[3][2]);

};

void global_to_local(const double * gp,
                     const double ** lv,
                     const double ** gv,
                     double * lp)
{

};

double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv)
{

  double partial[3][3]; 
  int i;
  if(lp[2] == 1.) return 1.; // 此处变换是有问题的
  for(i=0;i<3;i++){
  partial[i][0] = (
	  (gv[3][i]-gv[2][i]-gv[1][i]+gv[0][i])*lp[2]
	  +(-gv[3][i]+gv[2][i]-gv[1][i]+gv[0][i])*lp[1]
	  +(-gv[3][i]+gv[2][i]+gv[1][i]-gv[0][i])
	  )/(1.0-lp[2])/4.;
  partial[i][1] = -(
	  (gv[3][i]+gv[2][i]-gv[1][i]-gv[0][i])*lp[2]
	  +(gv[3][i]-gv[2][i]+gv[1][i]-gv[0][i])*lp[0]
	  +(-gv[3][i]-gv[2][i]+gv[1][i]+gv[0][i])
	  )/(1.0-lp[2])/4.;
  partial[i][2] = 
	  (4.*gv[4][i]-gv[3][i]-gv[2][i]-gv[1][i]-gv[0][i])/4.0
	+ (-gv[3][i]+gv[2][i]-gv[1][i]+gv[0][i])*lp[0]*lp[1]
	  /(1.0-lp[2])/(1.0-lp[2])/4.0 ;
  }

  return 
	partial[0][0]*partial[1][1]*partial[2][2]
	+partial[0][1]*partial[1][2]*partial[2][0]
	+partial[1][0]*partial[2][1]*partial[0][2]
	-partial[0][2]*partial[1][1]*partial[2][0]
	-partial[0][1]*partial[1][0]*partial[2][2]
	-partial[0][0]*partial[1][2]*partial[2][1]
	;
}

double global_to_local_jacobian(const double * gp,
                                const double ** lv,
                                const double ** gv)
{
}


void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial)
{
  // partial x / partial xi 
  int i;
  for(i=0;i<3;i++){
  partial[i][0] = (
	  (gv[3][i]-gv[2][i]-gv[1][i]+gv[0][i])*lp[2]
	  +(-gv[3][i]+gv[2][i]-gv[1][i]+gv[0][i])*lp[1]
	  +(-gv[3][i]+gv[2][i]+gv[1][i]-gv[0][i])
	  )/(1.0-lp[2])/4.;
  partial[i][1] = -(
	  (gv[3][i]+gv[2][i]-gv[1][i]-gv[0][i])*lp[2]
	  +(gv[3][i]-gv[2][i]+gv[1][i]-gv[0][i])*lp[0]
	  +(-gv[3][i]-gv[2][i]+gv[1][i]+gv[0][i])
	  )/(1.0-lp[2])/4.;
  partial[i][2] = 
	  (4.*gv[4][i]-gv[3][i]-gv[2][i]-gv[1][i]-gv[0][i])/4.0
	+ (-gv[3][i]+gv[2][i]-gv[1][i]+gv[0][i])*lp[0]*lp[1]
	  /(1.0-lp[2])/(1.0-lp[2])/4.0 ;
  }

};

void global_to_local_partial(const double * gp,
                             const double ** lv,
                             const double ** gv,
                             double **partial)
{
};
 
/*
#include <malloc.h>
int main(){

  double gp[3]={0,0,0};
  int i,j;
double dx=2, dy=2 , dz=1;
 double gv1[5][3]
	={ 
	  { 0, 0, 0 }, 
	  { dx, 0, 0}, 
	  { dx,dy, 0}, 
	  { 0,dy,0}, 
	  { 1, 1, dz }, 
	};
 double lv1[5][3]
	={ 
	  { -1,-1,0}, 
	  { 1,-1,0}, 
	  { 1,1,0}, 
	  { -1,1,0}, 
	  { 0,0,1} 
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

  gv = (const double**)malloc(sizeof(double*)*5);
  lv = (const double**)malloc(sizeof(double*)*5);
 
  for(i=0;i<5;i++){
	gv[i]=gv1[i];
	lv[i]=lv1[i];
  }


  printf("Input a local point lp\n");
  scanf("%lf%lf%lf",&lp[0],&lp[1],&lp[2]);
  printf("input lp = (%lf,%lf,%lf)\n",lp[0],lp[1],lp[2]);
  local_to_global(lp,lv,gv,gp);
  printf("output gp = (%lf,%lf,%lf)\n",gp[0],gp[1],gp[2]);

  printf("local_to_global_jacobian=%lf\n",local_to_global_jacobian(lp,lv,gv));
//  printf("global_to_local_jacobian=%lf\n",global_to_local_jacobian(gp,lv,gv));

  return 0;

}
*/
//
// end of file
///////////////////////////////////////////////////////////////////////////
