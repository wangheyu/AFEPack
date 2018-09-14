/**
 * @file   triangle.to3d.crd_trs.c
 * @author Yana Di, Robert Lie
 * @date   Tue Oct 23 11:47:40 2007
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
  double area = AREA(lv[0], lv[1], lv[2]);
  lambda[0] = AREA(lp, lv[1], lv[2])/area;
  lambda[1] = AREA(lp, lv[2], lv[0])/area;
  lambda[2] = AREA(lp, lv[0], lv[1])/area;
  gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0] + lambda[2]*gv[2][0];
  gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[1][1] + lambda[2]*gv[2][1];
  gp[2] = lambda[0]*gv[0][2] + lambda[1]*gv[1][2] + lambda[2]*gv[2][2];
}

void to3d_global_to_local(const double * gp, 
                          const double ** lv, 
                          const double ** gv, 
                          double * lp)
{
  double lambda[3];
  double n[3];
  n[0] = ((gv[1][1] - gv[0][1])*(gv[2][2] - gv[0][2]) -
	  (gv[1][2] - gv[0][2])*(gv[2][1] - gv[0][1]));
  n[1] = ((gv[1][2] - gv[0][2])*(gv[2][0] - gv[0][0]) -
	  (gv[1][0] - gv[0][0])*(gv[2][2] - gv[0][2]));
  n[2] = ((gv[1][0] - gv[0][0])*(gv[2][1] - gv[0][1]) -
	  (gv[1][1] - gv[0][1])*(gv[2][0] - gv[0][0]));
  double area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] = ((gv[1][1] - gp[1])*(gv[2][2] - gp[2]) -
	  (gv[1][2] - gp[2])*(gv[2][1] - gp[1]));
  n[1] = ((gv[1][2] - gp[2])*(gv[2][0] - gp[0]) -
	  (gv[1][0] - gp[0])*(gv[2][2] - gp[2]));
  n[2] = ((gv[1][0] - gp[0])*(gv[2][1] - gp[1]) -
	  (gv[1][1] - gp[1])*(gv[2][0] - gp[0]));
  lambda[0] = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])/area;
  n[0] = ((gp[1] - gv[0][1])*(gv[2][2] - gv[0][2]) -
	  (gp[2] - gv[0][2])*(gv[2][1] - gv[0][1]));
  n[1] = ((gp[2] - gv[0][2])*(gv[2][0] - gv[0][0]) -
	  (gp[0] - gv[0][0])*(gv[2][2] - gv[0][2]));
  n[2] = ((gp[0] - gv[0][0])*(gv[2][1] - gv[0][1]) -
	  (gp[1] - gv[0][1])*(gv[2][0] - gv[0][0]));
  lambda[1] = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])/area;
  n[0] = ((gv[1][1] - gv[0][1])*(gp[2] - gv[0][2]) -
	  (gv[1][2] - gv[0][2])*(gp[1] - gv[0][1]));
  n[1] = ((gv[1][2] - gv[0][2])*(gp[0] - gv[0][0]) -
	  (gv[1][0] - gv[0][0])*(gp[2] - gv[0][2]));
  n[2] = ((gv[1][0] - gv[0][0])*(gp[1] - gv[0][1]) -
	  (gv[1][1] - gv[0][1])*(gp[0] - gv[0][0]));
  lambda[2] = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])/area;
  lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0] + lambda[2]*lv[2][0];
  lp[1] = lambda[0]*lv[0][1] + lambda[1]*lv[1][1] + lambda[2]*lv[2][1];
}

double to3d_local_to_global_jacobian(const double * lp, 
                                     const double ** lv, 
                                     const double ** gv)
{
  double larea = fabs(AREA(lv[0], lv[1], lv[2]));
  double n[3];
  n[0] = ((gv[1][1] - gv[0][1])*(gv[2][2] - gv[0][2]) -
	  (gv[1][2] - gv[0][2])*(gv[2][1] - gv[0][1]));
  n[1] = ((gv[1][2] - gv[0][2])*(gv[2][0] - gv[0][0]) -
	  (gv[1][0] - gv[0][0])*(gv[2][2] - gv[0][2]));
  n[2] = ((gv[1][0] - gv[0][0])*(gv[2][1] - gv[0][1]) -
	  (gv[1][1] - gv[0][1])*(gv[2][0] - gv[0][0]));
  double garea = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return garea/larea;
}

double to3d_global_to_local_jacobian(const double * gp, 
                                     const double ** lv,
                                     const double ** gv)
{
  double larea = fabs(AREA(lv[0], lv[1], lv[2]));
  double n[3];
  n[0] = ((gv[1][1] - gv[0][1])*(gv[2][2] - gv[0][2]) -
	  (gv[1][2] - gv[0][2])*(gv[2][1] - gv[0][1]));
  n[1] = ((gv[1][2] - gv[0][2])*(gv[2][0] - gv[0][0]) -
	  (gv[1][0] - gv[0][0])*(gv[2][2] - gv[0][2]));
  n[2] = ((gv[1][0] - gv[0][0])*(gv[2][1] - gv[0][1]) -
	  (gv[1][1] - gv[0][1])*(gv[2][0] - gv[0][0]));
  double garea = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return larea/garea;
}

/**
 * end of file
 * 
 */



