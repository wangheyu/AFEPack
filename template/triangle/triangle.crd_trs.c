/**
 * @file   triangle.crd_trs.c
 * @author Robert Lie
 * @date   Tue Oct 23 12:35:50 2007
 * 
 * @brief  
 * 
 * 
 */

#define AREA(a, b, c) (((b)[0] - (a)[0])*((c)[1] - (a)[1]) - \
                       ((b)[1] - (a)[1])*((c)[0] - (a)[0]))
#define INN_PRD(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

void local_to_global(const double * lp, 
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
}

void global_to_local(const double * gp, 
                     const double ** lv, 
                     const double ** gv, 
                     double * lp)
{
  double lambda[3];
  double area = AREA(gv[0], gv[1], gv[2]);
  lambda[0] = AREA(gp, gv[1], gv[2])/area;
  lambda[1] = AREA(gp, gv[2], gv[0])/area;
  lambda[2] = AREA(gp, gv[0], gv[1])/area;
  lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0] + lambda[2]*lv[2][0];
  lp[1] = lambda[0]*lv[0][1] + lambda[1]*lv[1][1] + lambda[2]*lv[2][1];
}

double local_to_global_jacobian(const double * lp, 
                                const double ** lv, 
                                const double ** gv)
{
  double larea = AREA(lv[0], lv[1], lv[2]);
  double garea = AREA(gv[0], gv[1], gv[2]);
  return garea/larea;
}

double global_to_local_jacobian(const double * gp, 
                                const double ** lv,
                                const double ** gv)
{
  double larea = AREA(lv[0], lv[1], lv[2]);
  double garea = AREA(gv[0], gv[1], gv[2]);
  return larea/garea;
}

/**
 * end of file
 * 
 */

