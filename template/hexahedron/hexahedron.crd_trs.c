/*=============================================================================
#     FileName: hexahedron.crd_trs.c
#         Desc:
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-13 11:01:46
#      History:
=============================================================================*/

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#define HEXA_JAC_GAUSS_ELIMINATION 								\
	for (i = 0;i < 2;i ++) {							\
		k = i;											\
		for (j = i+1;j < 3;j ++)						\
			if (fabs(jac[j][i]) > fabs(jac[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 3;j ++) {					\
				tmp = jac[i][j];							\
				jac[i][j] = jac[k][j];						\
				jac[k][j] = tmp;							\
			}											\
			tmp = a[i];								\
			a[i] = a[k];							\
			a[k] = tmp;								\
		}												\
		for (j = i+1;j < 3;j ++) {						\
			tmp = jac[j][i]/jac[i][i];						\
			for (k = i+1;k < 3;k ++)					\
				jac[j][k] -= tmp*jac[i][k];					\
			a[j] -= tmp*a[i];						\
		}												\
	}													\
	a[2] /= jac[2][2];									\
	for (i = 1;i >= 0;i --) {							\
		for (j = i+1;j < 3;j ++) {						\
			a[i] -= jac[i][j]*a[j];					\
		}												\
		a[i] /= jac[i][i];								\
	}



#define HEXA_GAUSS_ELIMINATION 								\
	for (i = 0;i < 7;i ++) {							\
		k = i;											\
		for (j = i+1;j < 8;j ++)						\
			if (fabs(m[j][i]) > fabs(m[k][i])) k = j;	\
		if (k != i) {									\
			for (j = i;j < 8;j ++) {					\
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
		for (j = i+1;j < 8;j ++) {						\
			tmp = m[j][i]/m[i][i];						\
			for (k = i+1;k < 8;k ++)					\
				m[j][k] -= tmp*m[i][k];					\
			a[j][0] -= tmp*a[i][0];						\
			a[j][1] -= tmp*a[i][1];						\
			a[j][2] -= tmp*a[i][2];                \
		}												\
	}													\
	a[7][0] /= m[7][7];									\
	a[7][1] /= m[7][7];									\
	a[7][2] /= m[7][7];									\
	for (i = 6;i >= 0;i --) {							\
		for (j = i+1;j < 8;j ++) {						\
			a[i][0] -= m[i][j]*a[j][0];					\
			a[i][1] -= m[i][j]*a[j][1];					\
			a[i][2] -= m[i][j]*a[j][2];					\
		}												\
		a[i][0] /= m[i][i];								\
		a[i][1] /= m[i][i];								\
		a[i][2]	/= m[i][i];								\
	}

# define HEXA_JACOBIA_SET_VALUE \
	jac[0][0] = a[1][0] + a[4][0]*lp[1] + a[6][0]*lp[2] + a[7][0]*lp[1]*lp[2];	\
	jac[1][0] = a[2][0] + a[4][0]*lp[0] + a[5][0]*lp[2] + a[7][0]*lp[0]*lp[2];	\
	jac[2][0] = a[3][0] + a[5][0]*lp[1] + a[6][0]*lp[0] + a[7][0]*lp[0]*lp[1];	\
	jac[0][1] = a[1][1] + a[4][1]*lp[1] + a[6][1]*lp[2] + a[7][1]*lp[1]*lp[2];	\
	jac[1][1] = a[2][1] + a[4][1]*lp[0] + a[5][1]*lp[2] + a[7][1]*lp[0]*lp[2];	\
	jac[2][1] = a[3][1] + a[5][1]*lp[1] + a[6][1]*lp[0] + a[7][1]*lp[0]*lp[1];	\
	jac[0][2] = a[1][2] + a[4][2]*lp[1] + a[6][2]*lp[2] + a[7][2]*lp[1]*lp[2];	\
	jac[1][2] = a[2][2] + a[4][2]*lp[0] + a[5][2]*lp[2] + a[7][2]*lp[0]*lp[2];	\
	jac[2][2] = a[3][2] + a[5][2]*lp[1] + a[6][2]*lp[0] + a[7][2]*lp[0]*lp[1];


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
  int i, j, k;
  double m[8][8], tmp;
  double a[8][3];

  for (i = 0; i < 8; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][2];
    m[i][4] = lv[i][0]*lv[i][1];
    m[i][5] = lv[i][1]*lv[i][2];
    m[i][6] = lv[i][2]*lv[i][0];
    m[i][7] = lv[i][0]*lv[i][1]*lv[i][2];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
    a[i][2]	= gv[i][2];

  }

  HEXA_GAUSS_ELIMINATION
  gp[0] = a[0][0] + a[1][0]*lp[0] + a[2][0]*lp[1] + a[3][0]*lp[2]
          + a[4][0]*lp[0]*lp[1] + a[5][0]*lp[1]*lp[2] + a[6][0]*lp[2]*lp[0] + a[7][0]*lp[0]*lp[1]*lp[2] ;
  gp[1] = a[0][1] + a[1][1]*lp[0] + a[2][1]*lp[1] + a[3][1]*lp[2]
          + a[4][1]*lp[0]*lp[1] + a[5][1]*lp[1]*lp[2] + a[6][1]*lp[2]*lp[0] + a[7][1]*lp[0]*lp[1]*lp[2] ;
  gp[2] = a[0][2] + a[1][2]*lp[0] + a[2][2]*lp[1] + a[3][2]*lp[2]
          + a[4][2]*lp[0]*lp[1] + a[5][2]*lp[1]*lp[2] + a[6][2]*lp[2]*lp[0] + a[7][2]*lp[0]*lp[1]*lp[2] ;

};

void global_to_local(const double * gp,
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
  double vc[3];
  double xi[3];
  double jac[3][3];
  double tmp;
  double a[3];
  partial = (double **) malloc(3*sizeof(double*) );
  x0 = (double * ) malloc (3*sizeof(double));
  for ( i = 0; i < 3 ; i++)
  {
    partial[i] = (double *) malloc(3*sizeof(double));
  }
  for ( i =0 ; i < 3 ; i ++) 	x[i] = gp[i];
  for ( j = 0; j < 3; j ++)
  {
    vc[j]=0.0;
    for ( i = 0; i< 8 ; i++)
    {
      vc[j] +=lv[i][j];
    }
    vc[j] /=	8.0;
  }
//  Newton Rapson Iteration method
  xi[0] = vc[0];
  xi[1] = vc[1];
  xi[2] = vc[2];


  for ( cycle = 0; cycle < max_cycle; cycle++)
  {
    lp[0] = xi[0];
    lp[1] = xi[1];
    lp[2] = xi[2];
    local_to_global( lp, lv, gv, x0) ;
    local_to_global_partial( lp, lv, gv, partial);
    for ( i = 0 ; i < 3; i++)
    {
      for ( j =0 ; j <3; j++)
      {
        jac[i][j] = partial[j][i];
      }
    }
//		printf("jac = \n ");
//        for( i = 0 ; i < 3; i++) {
//        	for( j =0 ; j <3; j++){
//        		printf("%10.10lf  ", jac[i][j]);
//        	}
//			printf("\n");
//        }


    tmp = 0.;
    for ( i=0; i<3; i++)
    {
      a[i] = x0[i] - x[i];
      tmp += a[i]*a[i];
    };

    if ( xi[0] < -10.1 || xi[0] >10.1 || xi[1] < -10.1 || xi[1] >10.1 ||xi[2] < -10.1|| xi[2] >10.1)
    {
      printf("error = %e, %e,%e,%e \n", sqrt(tmp), xi[0],xi[1],xi[2] );
      printf(" \n ot inside of the element hexahedral\n" )    ;

      free(x0);
      for ( i = 0; i < 3 ; i++)
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
    HEXA_JAC_GAUSS_ELIMINATION

    for ( i= 0; i < 3; i++)
    {
      xi[i]  -= a[i];
    }

  }
  free(x0);
  for ( i = 0; i < 3 ; i++)
  {
    free(partial[i]);
  }
  free(partial);
  printf("error = %e, max_cycle = %d, epsilon = %lf, can't reach epsilon use mac_cycle steps", sqrt(tmp),max_cycle, epsilon);
};

double local_to_global_jacobian(const double * lp,
                                const double ** lv,
                                const double ** gv)
{
  int i, j, k;
  double m[8][8], tmp;
  double a[8][3];
  double jac[3][3];
  for (i = 0; i < 8; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][2];
    m[i][4] = lv[i][0]*lv[i][1];
    m[i][5] = lv[i][1]*lv[i][2];
    m[i][6] = lv[i][2]*lv[i][0];
    m[i][7] = lv[i][0]*lv[i][1]*lv[i][2];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
    a[i][2]	= gv[i][2];

  }
  HEXA_GAUSS_ELIMINATION

  HEXA_JACOBIA_SET_VALUE;


  return fabs(jac[0][0]*(jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) -
              jac[0][1]*(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]) +
              jac[0][2]*(jac[1][0]*jac[2][1]	- jac[1][1]*jac[2][0]) ) ;
}

double global_to_local_jacobian(const double * gp,
                                const double ** lv,
                                const double ** gv)
{
  double * lp ;
  double result;

  lp = (double *) malloc( 3*sizeof(double) );
  global_to_local( gp ,lv, gv, lp);
  result =	 1.0/local_to_global_jacobian(lp, lv, gv);
  free(lp);
  return result;
}

void local_to_global_partial(const double * lp,
                             const double ** lv,
                             const double ** gv,
                             double ** partial)
{
  int i, j, k;
  double m[8][8], tmp;
  double a[8][3];
  double jac[3][3];
  for (i = 0; i < 8; i ++)
  {
    m[i][0] = 1.0;
    m[i][1] = lv[i][0];
    m[i][2] = lv[i][1];
    m[i][3] = lv[i][2];
    m[i][4] = lv[i][0]*lv[i][1];
    m[i][5] = lv[i][1]*lv[i][2];
    m[i][6] = lv[i][2]*lv[i][0];
    m[i][7] = lv[i][0]*lv[i][1]*lv[i][2];
    a[i][0] = gv[i][0];
    a[i][1] = gv[i][1];
    a[i][2]	= gv[i][2];
  }
  HEXA_GAUSS_ELIMINATION
  HEXA_JACOBIA_SET_VALUE;

  partial[0][0] = jac[0][0];
  partial[0][1] = jac[0][1];
  partial[0][2] = jac[0][2];
  partial[1][0] = jac[1][0];
  partial[1][1] = jac[1][1];
  partial[1][2] = jac[1][2];
  partial[2][0] = jac[2][0];
  partial[2][1] = jac[2][1];
  partial[2][2] = jac[2][2];
};

void global_to_local_partial(const double * gp,
                             const double ** lv,
                             const double ** gv,
                             double **partial)
{
  int i, j;
  double jac[3][3];
  double det;
  double * lp;
  lp = (double *) malloc( 3*sizeof(double) );
  global_to_local(gp,lv,gv,lp);
  local_to_global_partial(lp,lv,gv,partial);
  free(lp);

  for ( i =0 ; i < 3; i++)
  {
    for ( j = 0 ; j < 3 ; j ++)
    {
      jac[i][j] = partial[i][j];
    }
  }

  det = jac[0][0]*(jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) -
        jac[0][1]*(jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0]) +
        jac[0][2]*(jac[1][0]*jac[2][1]	- jac[1][1]*jac[2][0])  ;

  partial[0][0] =( jac[1][1] * jac[2][2] - jac[1][2]*jac[2][1] )/det;
  partial[1][0] =- ( jac[1][0] * jac[2][2] - jac[1][2]*jac[2][0] )/det   ;
  partial[2][0] =( jac[1][0] * jac[2][1] - jac[1][1]*jac[2][0] )/det ;
  partial[0][1] =- ( jac[0][1] * jac[2][2] - jac[0][2]*jac[2][1] )/det;
  partial[1][1] = ( jac[0][0] * jac[2][2] - jac[0][2]*jac[2][0] )/det   ;
  partial[2][1] =- ( jac[0][0] * jac[2][1] - jac[0][1]*jac[2][0] )/det ;
  partial[0][2] =( jac[0][1] * jac[1][2] - jac[0][2]*jac[1][1] )/det;
  partial[1][2] =- ( jac[0][0] * jac[1][2] - jac[0][2]*jac[1][0] )/det   ;
  partial[2][2] =( jac[0][0] * jac[1][1] - jac[0][1]*jac[1][0] )/det ;
};


//
// end of file
///////////////////////////////////////////////////////////////////////////
