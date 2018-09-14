/*=============================================================================
#     FileName: quadrilateral.crd_trs.c
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

void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial);
void global_to_local(const double * gp,
                     const double ** lv,
                     const double ** gv,
                     double * lp);

#define HEXA_JAC_GAUSS_ELIMINATION 								\
	for (i = 0;i < 1;i ++) {							\
		k = i;											\
		for (j = i+1;j < 2;j ++)						\
			if (fabs(jac[j][i]) > fabs(jac[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 2;j ++) {					\
				tmp = jac[i][j];							\
				jac[i][j] = jac[k][j];						\
				jac[k][j] = tmp;							\
			}											\
			tmp = a[i];								\
			a[i] = a[k];							\
			a[k] = tmp;								\
		}												\
		for (j = i+1;j < 2;j ++) {						\
			tmp = jac[j][i]/jac[i][i];						\
			for (k = i+1;k < 2;k ++)					\
				jac[j][k] -= tmp*jac[i][k];					\
			a[j] -= tmp*a[i];						\
		}												\
	}													\
	a[1] /= jac[1][1];									\
	for (i = 0;i >= 0;i --) {		for (j = i+1;j < 2;j ++) {						\
			a[i] -= jac[i][j]*a[j];					\
		}												\
		a[i] /= jac[i][i];								\
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
		}												\
		for (j = i+1;j < 4;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i+1;k < 4;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];						\
		}												\
	}													\
	a[3][0] /= m[3][3];									\
	a[3][1] /= m[3][3];									\
	for (i = 2;i >= 0;i --) {	for (j = i+1;j < 4;j ++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];					\
		}												\
		a[i][0] /= m[i][i];								\
		a[i][1] /= m[i][i];								\
	}

# define HEXA_JACOBIA_SET_VALUE \
	jac[0][0] = a[1][0] + a[3][0]*lp[1] ;	\
	jac[1][0] = a[2][0] + a[3][0]*lp[0] ;	\
	jac[0][1] = a[1][1] + a[3][1]*lp[1] ;	\
	jac[1][1] = a[2][1] + a[3][1]*lp[0] ;




void local_to_global(const double * lp,
                     const double ** lv,
                     const double ** gv,
                     double * gp)
{
  int i, j, k;
  double m[4][4], tmp;
  double a[4][2];

  for (i = 0; i < 4; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][0]*lv[i][1];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];

  }

  HEXA_GAUSS_ELIMINATION
  gp[0] = a[0][0] + a[1][0]*lp[0] + a[2][0]*lp[1]
          + a[3][0]*lp[0]*lp[1] ;
  gp[1] = a[0][1] + a[1][1]*lp[0] + a[2][1]*lp[1]
          + a[3][1]*lp[0]*lp[1] ;
}

void global_to_local(const double * gp,
                     const double ** lv,
                     const double ** gv,
                     double * lp)
{
  double epsilon = 1e-7;
  const int max_cycle = 1000;
  int i,j,k, cycle;
  double x[2] ;
  double * x0;
  double ** partial;
  double vc[2];
  double xi[2];
  double jac[2][2];
  double tmp;
  double a[2];
  partial = (double **) malloc(2*sizeof(double*) );
  x0 = (double * ) malloc (2*sizeof(double));
  for( i = 0; i < 2 ; i++)
  {
    partial[i] = (double *) malloc(2*sizeof(double));
  }
  for( i =0 ; i < 2 ; i ++) 	x[i] = gp[i];
  for( j = 0; j < 2; j ++)
  {
    vc[j]=0.0;
    for( i = 0; i< 4 ; i++)
    {
      vc[j] +=lv[i][j];
    }
    vc[j] /=	4.0;
  }
//  Newton Rapson Iteration method
  xi[0] = vc[0];
  xi[1] = vc[1];


  for ( cycle = 0; cycle < max_cycle; cycle++)
  {
    lp[0] = xi[0];
    lp[1] = xi[1];
    local_to_global( lp, lv, gv, x0) ;
    local_to_global_partial( lp, lv, gv, partial);
    for( i = 0 ; i < 2; i++)
    {
      for( j =0 ; j <2; j++)
      {
        jac[i][j] = partial[j][i];
      }
    }
    /*
    		printf("%d  xi = %lf %lf \n ", cycle, xi[0], xi[1]);
    		printf("jac = \n ");
            for( i = 0 ; i < 2; i++) {
            	for( j =0 ; j < 2; j++){
            		printf("%10.10lf  ", jac[i][j]);
            	}
    			printf("\n");
            }
    */

    tmp = 0.;
    for ( i=0; i<2; i++)
    {
      a[i] = x0[i] - x[i];
      tmp += a[i]*a[i];
    };

    if( xi[0] < -10.1 || xi[0] >10.1 || xi[1] < -10.1 || xi[1] >10.1 )
    {
      printf("error = %e, %e,%e \n", sqrt(tmp), xi[0],xi[1] );
      printf(" \n ot inside of the element quadrilateral\n" )    ;

      free(x0);
      for( i = 0; i < 2 ; i++)
      {
        free(partial[i]);
      }
      free(partial);

      return ;
    }
    if( tmp < epsilon*epsilon )
    {
      free(x0);
      for( i = 0; i < 2 ; i++)
      {
        free(partial[i]);
      }
      free(partial);
      return;
    }
    HEXA_JAC_GAUSS_ELIMINATION

    for ( i= 0; i < 2; i++)
    {
      xi[i]  -= a[i];
    }

  }
  free(x0);
  for( i = 0; i < 2 ; i++)
  {
    free(partial[i]);
  }
  free(partial);
  printf("error = %e, max_cycle = %d, epsilon = %lf, can't reach epsilon use mac_cycle steps", sqrt(tmp),max_cycle, epsilon);
}




double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv)
{
  int i, j, k;
  double m[4][4], tmp;
  double a[4][2];
  double jac[2][2];
  for (i = 0; i < 4; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][0]*lv[i][1];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];

  }
  HEXA_GAUSS_ELIMINATION

  HEXA_JACOBIA_SET_VALUE;

  return fabs(jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0] );

}


double global_to_local_jacobian(const double * gp,
                                const double ** lv,
                                const double ** gv)
{
  double * lp ;
  double result;

  lp = (double *) malloc( 2*sizeof(double) );
  global_to_local( gp, lv, gv, lp);
  result = 1.0/local_to_global_jacobian(lp, lv, gv);
  free(lp);
  return result;
}



void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial)
{
  int i, j, k;
  double m[4][4], tmp;
  double a[4][2];
  double jac[2][2];
  for (i = 0; i < 4; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][0]*lv[i][1];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
  }
  HEXA_GAUSS_ELIMINATION
  HEXA_JACOBIA_SET_VALUE;

  partial[0][0] = jac[0][0];
  partial[0][1] = jac[0][1];
  partial[1][0] = jac[1][0];
  partial[1][1] = jac[1][1];
}

// here p is a global point
void global_to_local_partial(const double * gp,
                             const double ** lv,
                             const double ** gv,
                             double **partial)
{
  int i, j;
  double jac[2][2];
  double det;
  double * lp;
  lp = (double *) malloc( 2*sizeof(double) );
  global_to_local(gp,lv,gv,lp);
  local_to_global_partial(lp,lv,gv,partial);
  free(lp);

  for ( i =0 ; i < 2; i++)
  {
    for ( j = 0 ; j < 2 ; j ++)
    {
      jac[i][j] = partial[i][j];
    }
  }

  det = jac[0][0]*jac[1][1]-
        jac[0][1]*jac[1][0];

  partial[0][0] = jac[1][1] /det;
  partial[1][0] =- jac[1][0] /det   ;
  partial[0][1] =- jac[0][1] /det;
  partial[1][1] = jac[0][0] /det   ;
}


/**
 * end of file
 *
 */

/**测试
 * g++ -lm -o quadrilateral.crd_trs.exe quadrilateral.crd_trs.c
 * 测试的目的是保证 lp 映射到 gp, 再映射回 lp 要保持一致
 */
/*
int main(){
  int i,j,k;
  const double ** lv, ** gv;
  double tmp_lv[4][2];
  tmp_lv[0][0] = -1.0; tmp_lv[0][1] = -1.0;
  tmp_lv[1][0] = 1.0; tmp_lv[1][1] = -1.0;
  tmp_lv[2][0] = 1.0; tmp_lv[2][1] = 1.0;
  tmp_lv[3][0] = -1.0; tmp_lv[3][1] = 1.0;

  lv = (const double **) malloc(4*sizeof(double*));
  for(i=0;i<4;i++) lv[i] = tmp_lv[i];

  double tmp_gv[4][2];
  tmp_gv[0][0] = -1.0; tmp_gv[0][1] = -10.0;
  tmp_gv[1][0] = 1.0; tmp_gv[1][1] = -1.0;
  tmp_gv[2][0] = 2.0; tmp_gv[2][1] = 2.0;
  tmp_gv[3][0] = -2.0; tmp_gv[3][1] = 2.0;

  gv = (const double **) malloc(4*sizeof(double*));
  for(i=0;i<4;i++) gv[i] = tmp_gv[i];


  double lp[2]={0.1,0.7};
  double gp[2]={0.0,0.0};

  local_to_global(lp,lv,gv,gp);
  printf("lp =(%lf,%lf), gp=(%lf,%lf)\n",lp[0],lp[1],gp[0],gp[1]);
  global_to_local(gp,lv,gv,lp);
  printf("lp =(%lf,%lf), gp=(%lf,%lf)\n",lp[0],lp[1],gp[0],gp[1]);
  double jac_g2l, jac_l2g;
  jac_g2l = global_to_local_jacobian(gp,lv,gv);
  jac_l2g = local_to_global_jacobian(lp,lv,gv);
  printf("global_to_local_jacobian =%lf,local_to_global_jacobian=%lf \n", jac_g2l, jac_l2g);

  return 0;

}
*/
