/**
 * @file   twin_triangle.crd_trs.c
 * @author Ruo Li
 * @date   Tue Mar 22 14:54:58 2005
 * 
 * @brief  
 * 
 * 
 */

#define AREA(p0, p1, p2) \
  ((p1[0] - p0[0])*(p2[1] - p0[1]) - \
   (p1[1] - p0[1])*(p2[0] - p0[0]))

#ifdef __cplusplus
extern "C" {
#endif

void local_to_global(const double * lp, 
		     const double ** lv, 
		     const double ** gv, 
		     double * gp)
{
  static int _[6] = {0, 1, 2, 0, 2, 3};
  double lambda[3], chi = AREA(lp, lv[0], lv[2]);
  int I = (chi > 0)?3:0;
  double area = AREA(lv[_[I]], lv[_[I+1]], lv[_[I+2]]);
  lambda[0] = AREA(lp, lv[_[I+1]], lv[_[I+2]])/area;
  lambda[1] = AREA(lp, lv[_[I+2]], lv[_[I]])/area;
  lambda[2] = AREA(lp, lv[_[I]], lv[_[I+1]])/area;
  gp[0] = lambda[0]*gv[_[I]][0] + lambda[1]*gv[_[I+1]][0] + lambda[2]*gv[_[I+2]][0];
  gp[1] = lambda[0]*gv[_[I]][1] + lambda[1]*gv[_[I+1]][1] + lambda[2]*gv[_[I+2]][1];
};

void global_to_local(const double * gp, 
		     const double ** lv, 
		     const double ** gv, 
		     double * lp)
{
  static int _[6] = {0, 1, 2, 0, 2, 3};
  double lambda[3], chi = AREA(gp, gv[0], gv[2]);
  int I = (chi > 0)?3:0;
  double area = AREA(gv[_[I]], gv[_[I+1]], gv[_[I+2]]);
  lambda[0] = AREA(gp, gv[_[I+1]], gv[_[I+2]])/area;
  lambda[1] = AREA(gp, gv[_[I+2]], gv[_[I]])/area;
  lambda[2] = AREA(gp, gv[_[I]], gv[_[I+1]])/area;
  lp[0] = lambda[0]*lv[_[I]][0] + lambda[1]*lv[_[I+1]][0] + lambda[2]*lv[_[I+2]][0];
  lp[1] = lambda[0]*lv[_[I]][1] + lambda[1]*lv[_[I+1]][1] + lambda[2]*lv[_[I+2]][1];
};

double local_to_global_jacobian(const double * lp, 
				const double ** lv, 
				const double ** gv)
{
  static int _[6] = {0, 1, 2, 0, 2, 3};
  double chi = AREA(lp, lv[0], lv[2]);
  int I = (chi > 0)?3:0;
  double larea = AREA(lv[_[I]], lv[_[I+1]], lv[_[I+2]]);
  double garea = AREA(gv[_[I]], gv[_[I+1]], gv[_[I+2]]);
  return garea/larea;
};

double global_to_local_jacobian(const double * gp, 
				const double ** lv,
				const double ** gv)
{
  static int _[6] = {0, 1, 2, 0, 2, 3};
  double chi = AREA(gp, gv[0], gv[2]);
  int I = (chi > 0)?3:0;
  double larea = AREA(lv[_[I]], lv[_[I+1]], lv[_[I+2]]);
  double garea = AREA(gv[_[I]], gv[_[I+1]], gv[_[I+2]]);
  return larea/garea;
}

#ifdef __cplusplus
}
#endif

/**
 * end of file
 * 
 */

