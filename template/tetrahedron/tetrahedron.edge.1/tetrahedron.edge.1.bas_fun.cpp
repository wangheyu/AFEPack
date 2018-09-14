/*************************************************************************
* 	tetrahedron.edge.1.bas_fun.c : by D.Wang
*/

#include <cmath>
#include <Miscellaneous.h>

#define vector_length 3
#define vt nVector<vector_length,double>

#define GET_L(i, j)                                                     \
  double l0 = (v[j][0] - v[i][0]);                                      \
  double l1 = (v[j][1] - v[i][1]);                                      \
  double l2 = (v[j][2] - v[i][2]);                                      \
  double l = 1;                                                         \
  if (fabs(l0) > fabs(l1)) {                                            \
    if(fabs(l0) > fabs(l2)){                                            \
      if (l0 < 0) l = -1;                                               \
    }                                                                   \
    else if(l2 < 0) l = -1;                                             \
  }                                                                     \
  else {                                                                \
    if(fabs(l1) > fabs(l2)){                                            \
      if (l1 < 0) l = -1;                                               \
    }                                                                   \
    else if(l2 < 0) l = -1;                                             \
  }                                                                     \


#ifdef __cplusplus
extern "C" {
#endif

double get_volume(const double * v0,
		const double * v1,
		const double * v2,
		const double * v3)
{
	return ((v1[0] - v0[0])*(v2[1] - v0[1])*(v3[2] - v0[2])
		+ (v1[1] - v0[1])*(v2[2] - v0[2])*(v3[0] - v0[0])
		+ (v1[2] - v0[2])*(v2[0] - v0[0])*(v3[1] - v0[1])
		- (v1[0] - v0[0])*(v2[2] - v0[2])*(v3[1] - v0[1])
		- (v1[1] - v0[1])*(v2[0] - v0[0])*(v3[2] - v0[2])
		- (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]));
}


#define co_det(v, m, n) (\
		((m%2==0)?-1.:1.) * (\
		  (v[(m+2)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
		* (v[(m+3)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		- (v[(m+2)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		* (v[(m+3)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
	) \
	)

void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = get_volume(p, v[1], v[2], v[3]);
	val[0] /= volume;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = get_volume(v[0], p, v[2], v[3]);
	val[0] /= volume;
}

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = get_volume(v[0], v[1], p, v[3]);
	val[0] /= volume;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = get_volume(v[0], v[1], v[2], p);
	val[0] /= volume;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = co_det(v, 0, 0)/volume;
	val[1] = co_det(v, 0, 1)/volume;
	val[2] = co_det(v, 0, 2)/volume;
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = co_det(v, 1, 0)/volume;
	val[1] = co_det(v, 1, 1)/volume;
	val[2] = co_det(v, 1, 2)/volume;
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = co_det(v, 2, 0)/volume;
	val[1] = co_det(v, 2, 1)/volume;
	val[2] = co_det(v, 2, 2)/volume;
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = co_det(v, 3, 0)/volume;
	val[1] = co_det(v, 3, 1)/volume;
	val[2] = co_det(v, 3, 2)/volume;
}


void phi_1(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;

    double tmp1[3],tmp2[3];

    lambda_1(p, v, &p1);
    lambda_2(p, v, &p2);
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_2(p, v, tmp2);

    int i;
    GET_L(0,1);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_1(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;

    double tmp1[3],tmp2[3];

    
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_2(p, v, tmp2);

    int i;
    GET_L(0,1);
    for(i = 0; i < 3 ; i++){
      val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}

void phi_2(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;

    double tmp1[3],tmp2[3];

    lambda_1(p, v, &p1);
    lambda_3(p, v, &p2);
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_3(p, v, tmp2);

    int i;
    GET_L(0,2);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_2(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;

    double tmp1[3],tmp2[3];

    
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_3(p, v, tmp2);

    int i;
    GET_L(0,2);
    for(i = 0; i < 3 ; i++){
      val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}

void phi_3(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    lambda_1(p, v, &p1);
    lambda_4(p, v, &p2);
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(0,3);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_3(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    
    gradient_lambda_1(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(0,3);
    for(i = 0; i < 3 ; i++){
       val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}

void phi_4(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    lambda_3(p, v, &p1);
    lambda_4(p, v, &p2);
    gradient_lambda_3(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(2,3);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_4(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    
    gradient_lambda_3(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(2,3);
    for(i = 0; i < 3 ; i++){
       val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}

void phi_5(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;
    
    double tmp1[3],tmp2[3];

    lambda_2(p, v, &p1);
    lambda_4(p, v, &p2);
    gradient_lambda_2(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(1,3);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_5(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;
    
    double tmp1[3],tmp2[3];

    
    gradient_lambda_2(p, v, tmp1);
    gradient_lambda_4(p, v, tmp2);

    int i;
    GET_L(1,3);
    for(i = 0; i < 3 ; i++){
       val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}

void phi_6(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    lambda_2(p, v, &p1);
    lambda_3(p, v, &p2);
    gradient_lambda_2(p, v, tmp1);
    gradient_lambda_3(p, v, tmp2);

    int i;
    GET_L(1,2);
    for(i = 0; i < 3 ; i++){
      val[i] = l*(p1 * tmp2[i] - p2 * tmp1[i]);
    }
}

void rot_phi_6(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
    double p1,p2;

    double tmp1[3],tmp2[3];
    
    
    gradient_lambda_2(p, v, tmp1);
    gradient_lambda_3(p, v, tmp2);

    int i;
    GET_L(1,2);
    for(i = 0; i < 3 ; i++){
       val[0][i] = 2*l*(tmp1[(i+1)%3]*tmp2[(i+2)%3] -  tmp2[(i+1)%3]*tmp1[(i+2)%3]) ;
    }
}


#ifdef __cplusplus
}
#endif

/*
*  end of file
**************************************************************************/
