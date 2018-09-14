/**
 * @file   twin_triangle.iso.crd_trs.c
 * @author Robert Lie
 * @date   Thu May 31 11:17:52 2007
 * 
 * @brief  
 * 
 * 
 */

#define AREA(p0, p1, p2) \
  ((p1[0] - p0[0])*(p2[1] - p0[1]) - (p1[1] - p0[1])*(p2[0] - p0[0]))

void iso_local_to_global(const double * lp, 
                         const double ** lv, 
                         const double ** gv, 
                         double * gp)
{
  double lambda[3], area;
  double ind = AREA(lv[0], lp, lv[2]);
  if (ind > 0) {
    area = AREA(lv[0], lv[1], lv[2]);
    lambda[0] = AREA(lp, lv[1], lv[2])/area;
    lambda[1] = AREA(lp, lv[2], lv[0])/area;
    lambda[2] = AREA(lp, lv[0], lv[1])/area;
    gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0] + lambda[2]*gv[2][0];
    gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[1][1] + lambda[2]*gv[2][1];
  } else {
    area = AREA(lv[0], lv[2], lv[3]);
    lambda[0] = AREA(lp, lv[2], lv[3])/area;
    lambda[1] = AREA(lp, lv[3], lv[0])/area;
    lambda[2] = AREA(lp, lv[0], lv[2])/area;
    gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[2][0] + lambda[2]*gv[3][0];
    gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[2][1] + lambda[2]*gv[3][1];
  }
}

void iso_global_to_local(const double * gp, 
                         const double ** lv, 
                         const double ** gv, 
                         double * lp)
{
  double lambda[3], area;
  double ind = AREA(gv[0], gp, gv[2]);
  if (ind > 0) {
    area = AREA(gv[0], gv[1], gv[2]);
    lambda[0] = AREA(gp, gv[1], gv[2])/area;
    lambda[1] = AREA(gp, gv[2], gv[0])/area;
    lambda[2] = AREA(gp, gv[0], gv[1])/area;
    lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0] + lambda[2]*lv[2][0];
    lp[1] = lambda[0]*lv[0][1] + lambda[1]*lv[1][1] + lambda[2]*lv[2][1];
  } else {
    area = AREA(gv[0], gv[2], gv[3]);
    lambda[0] = AREA(gp, gv[2], gv[3])/area;
    lambda[1] = AREA(gp, gv[3], gv[0])/area;
    lambda[2] = AREA(gp, gv[0], gv[2])/area;
    lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[2][0] + lambda[2]*lv[3][0];
    lp[1] = lambda[0]*lv[0][1] + lambda[1]*lv[2][1] + lambda[2]*lv[3][1];
  }
}

double iso_local_to_global_jacobian(const double * lp, 
                                    const double ** lv, 
                                    const double ** gv)
{
  double larea, garea;
  double ind = AREA(lv[0], lp, lv[2]);
  if (ind > 0) {
    larea = AREA(lv[0], lv[1], lv[2]);
    garea = AREA(gv[0], gv[1], gv[2]);
  } else {
    larea = AREA(lv[0], lv[2], lv[3]);
    garea = AREA(gv[0], gv[2], gv[3]);
  }
  return garea/larea;
}

double iso_global_to_local_jacobian(const double * gp, 
                                    const double ** lv,
                                    const double ** gv)
{
  double larea, garea;
  double ind = AREA(gv[0], gp, gv[2]);
  if (ind > 0) {
    larea = AREA(lv[0], lv[1], lv[2]);
    garea = AREA(gv[0], gv[1], gv[2]);
  } else {
    larea = AREA(lv[0], lv[2], lv[3]);
    garea = AREA(gv[0], gv[2], gv[3]);
  }
  return larea/garea;
};

/**
 * end of file
 * 
 */

