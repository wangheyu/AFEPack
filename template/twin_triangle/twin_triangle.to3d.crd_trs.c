/**
 * @file   twin_triangle.to3d.crd_trs.c
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Tue Aug  9 19:26:58 2011
 * 
 * @brief  
 * 
 * 
 */
#include <math.h>

#define AREA(a, b, c) (((b)[0] - (a)[0])*((c)[1] - (a)[1]) - \
                       ((b)[1] - (a)[1])*((c)[0] - (a)[0]))
#define INN_PRD(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

void to3d_local_to_global(const double * lp, 
                          const double ** lv, 
                          const double ** gv, 
                          double * gp)
{
  double lambda[3];
  double area = AREA(lv[0], lv[1], lv[3]);
  lambda[0] = AREA(lp, lv[1], lv[3])/area;
  lambda[1] = AREA(lp, lv[3], lv[0])/area;
  lambda[2] = AREA(lp, lv[0], lv[1])/area;
  gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0] + lambda[2]*gv[3][0];
  gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[1][1] + lambda[2]*gv[3][1];
  gp[2] = lambda[0]*gv[0][2] + lambda[1]*gv[1][2] + lambda[2]*gv[3][2];
}

void to3d_global_to_local(const double * gp, 
                          const double ** lv, 
                          const double ** gv, 
                          double * lp)
{
}

double to3d_local_to_global_jacobian(const double * lp, 
                                     const double ** lv, 
                                     const double ** gv)
{
  double larea = fabs(AREA(lv[0], lv[1], lv[3]));
  double n[3];
  n[0] = ((gv[1][1] - gv[0][1])*(gv[3][2] - gv[0][2]) -
	  (gv[1][2] - gv[0][2])*(gv[3][1] - gv[0][1]));
  n[1] = ((gv[1][2] - gv[0][2])*(gv[3][0] - gv[0][0]) -
	  (gv[1][0] - gv[0][0])*(gv[3][2] - gv[0][2]));
  n[2] = ((gv[1][0] - gv[0][0])*(gv[3][1] - gv[0][1]) -
	  (gv[1][1] - gv[0][1])*(gv[3][0] - gv[0][0]));
  double garea = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return garea/larea;
}

double to3d_global_to_local_jacobian(const double * gp, 
                                     const double ** lv,
                                     const double ** gv)
{
  return 0.0;
}

/**
 * end of file
 * 
 */



