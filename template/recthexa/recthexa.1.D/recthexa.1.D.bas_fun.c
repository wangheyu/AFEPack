/*=============================================================================
#     FileName: recthexa.1.D.bas_fun.c
#         Desc: 
#       Author: Tiao Lu
#        Email: tlu@math.pku.edu.cn
#     HomePage: http://dsec.pku.edu.cn/~tlu
#      Version: 0.0.1
#   LastChange: 2013-10-13 11:41:28
#      History:
=============================================================================*/
#include "../recthexa.crd_trs.c"


/* basis functions are the polynomials of the local coodinates
   such as \xi , \eta , etc */
/* note that the parameters lp are the local points to the reference element. */
#define L1(x) x
#define L2(x) (x*x - 1./3.)
#define L3(x) (x*x*x - 0.6*x)
#define L4(x) (x*x*x*x - 6./7.*x*x + 3./35.)
#define L5(x) (63.*x*x*x*x*x - 70*x*x*x +15*x)
#define L6(x) (231.*x*x*x*x*x*x - 315.*x*x*x*x + 105.*x*x -5.)
#define L7(x) (429.*x*x*x*x*x*x*x -693.*x*x*x*x*x + 315.*x*x*x - 35.*x)
#define L8(x) (6435.*x*x*x*x*x*x*x*x - 12012.*x*x*x*x*x*x + 6930.*x*x*x*x -1260*x*x + 35)

/* the derivetives */
#define DL1(x) 1.
#define DL2(x) (2.*x)
#define DL3(x) (3.*x*x - 0.6)
#define DL4(x) (4.*x*x*x - 12./7.*x)
#define DL5(x) (315.*x*x*x*x - 210.*x*x +15.)
#define DL6(x) (1386.*x*x*x*x*x - 1260.*x*x*x + 210.*x)
#define DL7(x) (3003.*x*x*x*x*x*x - 3465*x*x*x*x + 945.*x*x - 35.)
#define DL8(x) (51480.*x*x*x*x*x*x*x - 72072.*x*x*x*x*x + 27720.*x*x*x - 2520.*x)

#define DEFINE_LV_ \
	int k; \
  double lv1[8][3] \
	={  { -1,-1,-1}, { 1,-1,-1}, { 1,1,-1}, { -1,1,-1},  \
	  { -1,-1,1}, { 1,-1,1}, { 1,1,1}, { -1,1,1}  \
	}; \
  const double **lv;  \
  lv = (const double**)malloc(sizeof(double*)*8); \
  for(k=0;k<8;k++){ \
	lv[k]=lv1[k]; \
  } \
  double lp[3]; \
  global_to_local(p,lv,v, lp); 



#define DEFINE_PARTIAL_ \
    double ** partial; \
  partial = (double**)malloc(sizeof(double*)*3);  \
  for(k=0;k<3;k++){ \
	partial[k] =(double*)malloc(sizeof(double)*3); \
  } \
  global_to_local_partial(p,lv,v,partial); \
  val[0] = grad[0]*partial[0][0] + grad[1]*partial[0][1] + grad[2]*partial[0][2]; \
  val[1] = grad[0]*partial[1][0] + grad[1]*partial[1][1] + grad[2]*partial[1][2]; \
  val[2] = grad[0]*partial[2][0] + grad[1]*partial[2][1] + grad[2]*partial[2][2]; 


void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
   val[0] = 1.0;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	val[0] = L1(lp[0]); 
}

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	val[0] = L1(lp[1]); 
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	val[0] = L1(lp[2]); 
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
	val[2] = 0.0;
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	double grad[3];
	grad[0] = DL1(lp[0]);
	grad[1] = 0.;
	grad[2] = 0.;
  
	DEFINE_PARTIAL_
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	double grad[3];
	grad[0] = 0.;
	grad[1] = DL1(lp[1]);
	grad[2] = 0.;

	DEFINE_PARTIAL_
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	DEFINE_LV_
	double grad[3];
	grad[0] = 0.;
	grad[1] = 0.;
	grad[2] = DL1(lp[2]);

	DEFINE_PARTIAL_
}

/*
*  end of file
**************************************************************************/
