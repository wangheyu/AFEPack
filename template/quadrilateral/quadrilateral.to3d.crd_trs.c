/*=============================================================================
#     FileName: quadrilateral.to3d.crd_trs.c
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

void to3d_local_to_global_partial(const double * lp,
                                  const double ** lv,
                                  const double ** gv,
                                  double ** partial);
void to3d_global_to_local(const double * gp,
                          const double ** lv,
                          const double ** gv,
                          double * lp);
#define HEXA_MATJJ_GAUSS_ELIMINATION 								\
	for (i = 0;i < 1;i ++) {							\
		k = i;											\
		for (j = i+1;j < 2;j ++)						\
			if (fabs(matjj[j][i]) > fabs(matjj[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 2;j ++) {					\
				tmp = matjj[i][j];							\
				matjj[i][j] = matjj[k][j];						\
				matjj[k][j] = tmp;							\
			}											\
			tmp = b[i];								\
			b[i] = b[k];							\
			b[k] = tmp;								\
		}												\
		for (j = i+1;j < 2;j ++) {						\
			tmp = matjj[j][i]/matjj[i][i];						\
			for (k = i+1;k < 2;k ++)					\
				matjj[j][k] -= tmp*matjj[i][k];					\
			b[j] -= tmp*b[i];						\
		}												\
	}													\
	b[1] /= matjj[1][1];									\
	for (i = 0;i >= 0;i --) {		for (j = i+1;j < 2;j ++) {						\
			b[i] -= matjj[i][j]*b[j];					\
		}												\
		b[i] /= matjj[i][i];								\
	}

#define HEXA_GAUSS_ELIMINATION 								\
	for (i = 0;i < 3;i ++) {							\
		k = i;											\
		for (j = i+1;j < 4;j ++)						\
			if (fabs(m[j][i]) > fabs(m[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 4;j ++) {					\
				tmp = m[i][j];							\
				m[i][j] = m[k][j];						\
				m[k][j] = tmp;							\
			}											\
			tmp = a[i][0];								\
			a[i][0] = a[k][0];							\
			a[k][0] = tmp;								\
			tmp = a[i][1];								\
			a[i][1] = a[k][1];							\
			a[k][1] = tmp;								\
			tmp = a[i][2];								\
			a[i][2] = a[k][2];							\
			a[k][2] = tmp;								\
		}												\
		for (j = i+1;j < 4;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i+1;k < 4;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];						\
			a[j][2] -= tmp*a[i][2];						\
		}												\
	}													\
	a[3][0] /= m[3][3];									\
	a[3][1] /= m[3][3];									\
	a[3][2] /= m[3][3];									\
	for (i = 2;i >= 0;i --) {	for (j = i+1;j < 4;j ++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];					\
			a[i][2] -= m[i][j]*a[j][2];					\
		}												\
		a[i][0] /= m[i][i];								\
		a[i][1] /= m[i][i];								\
		a[i][2] /= m[i][i];								\
	}

# define HEXA_JACOBIA_SET_VALUE \
	jac[0][0] = a[1][0] + a[3][0]*lp[1] ;	\
	jac[1][0] = a[2][0] + a[3][0]*lp[0] ;	\
	jac[0][1] = a[1][1] + a[3][1]*lp[1] ;	\
	jac[1][1] = a[2][1] + a[3][1]*lp[0] ;	\
	jac[0][2] = a[1][2] + a[3][2]*lp[1] ;	\
	jac[1][2] = a[2][2] + a[3][2]*lp[0] ;


void to3d_local_to_global(const double * lp,
                          const double ** lv,
                          const double ** gv,
                          double * gp)
{
  int i, j, k;
  double m[4][4], tmp;
  double a[4][3];

  for (i = 0; i < 4; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][0]*lv[i][1];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
    a[i][2] = gv[i][2];

  }

  HEXA_GAUSS_ELIMINATION
  gp[0] = a[0][0] + a[1][0]*lp[0] + a[2][0]*lp[1]
          + a[3][0]*lp[0]*lp[1] ;
  gp[1] = a[0][1] + a[1][1]*lp[0] + a[2][1]*lp[1]
          + a[3][1]*lp[0]*lp[1] ;
  gp[2] = a[0][2] + a[1][2]*lp[0] + a[2][2]*lp[1]
          + a[3][2]*lp[0]*lp[1] ;
}
//使用了求解非线性方程组的Newton-Raphson迭代，效率比最速下降法高
void to3d_global_to_local(const double * gp,
                          const double ** lv,
                          const double ** gv,
                          double * lp)
{
  double epsilon = 1e-7;
  const int max_cycle = 1000;
  int i,j,k, cycle;
  double x[3] ;
  double * x0;
  double ** partial;
  double vc[2];
  double xi[2];
  double jac[2][3];
  double tmp;
  double a[2];
  partial = (double **) malloc(3*sizeof(double*) );
  x0 = (double * ) malloc (3*sizeof(double));
  for ( i = 0; i < 3 ; i++)
  {
    partial[i] = (double *) malloc(2*sizeof(double));
  }
  for ( i =0 ; i < 3 ; i ++) 	x[i] = gp[i];
  for ( j = 0; j < 2; j ++)
  {
    vc[j]=0.0;
    for ( i = 0; i< 4 ; i++)
    {
      vc[j] +=lv[i][j];
    }
    vc[j] /=	4.0;
  }

  xi[0] = vc[0];
  xi[1] = vc[1];


  for ( cycle = 0; cycle < max_cycle; cycle++)
  {
    lp[0] = xi[0];
    lp[1] = xi[1];
    to3d_local_to_global( lp, lv, gv, x0) ;
    to3d_local_to_global_partial( lp, lv, gv, partial);

    //printf("lp = (%lf, %lf) , gp=(%lf,%lf,%lf)\n",lp[0],lp[1],x0[0],x0[1],x0[2]);



    tmp = 0.;
    for ( i=0; i<3; i++)
    {
      a[i] = x[i]-x0[i];
      tmp += a[i]*a[i];
    };

    if ( xi[0] < -10.1 || xi[0] >10.1 || xi[1] < -10.1 || xi[1] >10.1 )
    {
      printf("error = %e, %e,%e \n", sqrt(tmp), xi[0],xi[1] );
      printf(" \n not inside of the element quadrilateral\n" )    ;

      free(x0);
      for ( i = 0; i < 2 ; i++)
      {
        free(partial[i]);
      }
      free(partial);

      return ;
    }
    if ( tmp < epsilon*epsilon )
    {
      free(x0);
      for ( i = 0; i < 3 ; i++)
      {
        free(partial[i]);
      }
      free(partial);
      return;
    }
	//printf("tmp = %e at cycle %d\n", tmp, cycle);

    double b[2]= {0.0,0.0};
    for(i=0; i<2; i++)
    {
      for(j=0; j<3; j++)
      {
        b[i]+=partial[j][i]*a[j];
      }
    }
    double matjj[2][2]={{0.0,0.0},{0.0,0.0}};
    for(i=0; i<2; i++)
    {
      for(j=0; j<2; j++)
      {
        for(k=0; k<3; k++)
        {
          matjj[i][j]+=partial[k][i]*partial[k][j];
        }
      }
    }

	HEXA_MATJJ_GAUSS_ELIMINATION

    for ( i= 0; i < 2; i++)
    {
      xi[i]  += b[i];
    }

  }
  free(x0);
  for ( i = 0; i < 2 ; i++)
  {
    free(partial[i]);
  }
  free(partial);
  printf("error^2 = %e, max_cycle = %d, epsilon = %e, can't reach epsilon use mac_cycle steps", tmp,max_cycle, epsilon);
}

double to3d_local_to_global_jacobian(const double * lp,
                                     const double ** lv,
                                     const double ** gv)
{
  int i;
  double ** partial;
  double E, G, F;
  partial = (double **) malloc(3*sizeof(double*) );
  for ( i = 0; i < 3 ; i++)
  {
    partial[i] = (double *) malloc(2*sizeof(double));
  }
  to3d_local_to_global_partial(lp,lv,gv,partial);
  E = (partial[0][0]*partial[0][0]+ partial[1][0]*partial[1][0]
       +partial[2][0]*partial[2][0]);
  G = (partial[0][1]*partial[0][1]+ partial[1][1]*partial[1][1]
       +partial[2][1]*partial[2][1]);
  F = (partial[0][0]*partial[0][1]+ partial[1][0]*partial[1][1]
       +partial[2][0]*partial[2][1]);
  return sqrt(E*G-F*F);
}

double to3d_global_to_local_jacobian(const double * gp,
                                     const double ** lv,
                                     const double ** gv)
{

  double * lp ;
  double result;

  lp = (double *) malloc( 2*sizeof(double) );
  to3d_global_to_local( gp, lv, gv, lp);
  result = 1.0/to3d_local_to_global_jacobian(lp, lv, gv);
  free(lp);
  return result;
}


void to3d_local_to_global_partial(const double * lp,
                                  const double ** lv,
                                  const double ** gv,
                                  double ** partial)
{
  int i, j, k;
  double m[4][4], tmp;
  double a[4][3];
  double jac[2][3];
  for (i = 0; i < 4; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][0]*lv[i][1];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
    a[i][2] = gv[i][2];
  }
  HEXA_GAUSS_ELIMINATION
  HEXA_JACOBIA_SET_VALUE;

  partial[0][0] = jac[0][0];  // x_\xi
  partial[1][0] = jac[0][1];  // y_\xi
  partial[2][0] = jac[0][2];  // z_\xi
  partial[0][1] = jac[1][0];  // x_\eta
  partial[1][1] = jac[1][1];  // y_\eta
  partial[2][1] = jac[1][2];  // z_\eta
}


/**
 * end of file
 *
 */

/**测试
 * g++ -lm -o quadrilateral.to3d.crd_trs.exe quadrilateral.to3d.crd_trs.c
 * 测试的目的是保证 lp 映射到 gp, 再映射回 lp 要保持一致
 * 测试结果证明了正确性，另外非线性迭代的一般只要几次即可收敛。 
tlu@tlu-ThinkPad-X230s:~/Downloads/AFEPack/template/quadrilateral$ ./a.exe 
lp =(-0.577350269189626, 0.57735026918962), gp=()
global_to_local_jacobian =1.,local_to_global_jacobian=1 
 */ 
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


