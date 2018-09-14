/**
 * @file   rectangle.crd_trs.to3d.cpp
 * @author Robert Lie
 * @date   Mon Jul 23 15:05:39 2007
 * 
 * @brief  
 * 
 * 
 */

#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

  void to3d_local_to_global(const double * lp, 
                            const double ** lv, 
                            const double ** gv, 
                            double * gp)
  {
    /**
     * 下面的程序仅仅当参考单元的顶点为 (-1,-1), (1,-1), (1,1) 和
     * (-1,1) 的时候才对。
     */
    double alpha = 0.5*(lp[0] + 1.0);
    double beta  = 0.5*(lp[1] + 1.0);
    gp[0] = ((1 - alpha)*(1 - beta)*gv[0][0] +
             alpha      *(1 - beta)*gv[1][0] +
             alpha      *beta      *gv[2][0] +
             (1 - alpha)*beta      *gv[3][0]);
    gp[1] = ((1 - alpha)*(1 - beta)*gv[0][1] +
             alpha      *(1 - beta)*gv[1][1] +
             alpha      *beta      *gv[2][1] +
             (1 - alpha)*beta      *gv[3][1]);
    gp[2] = ((1 - alpha)*(1 - beta)*gv[0][2] +
             alpha      *(1 - beta)*gv[1][2] +
             alpha      *beta      *gv[2][2] +
             (1 - alpha)*beta      *gv[3][2]);
  }

  void to3d_global_to_local(const double * gp, 
                            const double ** lv, 
                            const double ** gv, 
                            double * lp)
  {
    /// not implemented!
  }

  double to3d_local_to_global_jacobian(const double * lp, 
                                       const double ** lv, 
                                       const double ** gv)
  {
    /**
     * 下面的程序仅仅当参考单元的顶点为 (-1,-1), (1,-1), (1,1) 和
     * (-1,1) 的时候才对。
     */
    double alpha = 0.5*(lp[0] + 1.0);
    double beta  = 0.5*(lp[1] + 1.0);
    double lambda[2][3] = {
      {
        (1 - beta)*(gv[1][0] - gv[0][0]) + beta*(gv[2][0] - gv[2][0]),
        (1 - beta)*(gv[1][1] - gv[0][1]) + beta*(gv[2][1] - gv[2][1]),
        (1 - beta)*(gv[1][2] - gv[0][2]) + beta*(gv[2][2] - gv[2][2])
      }, {
        (1 - alpha)*(gv[3][0] - gv[0][0]) + alpha*(gv[2][0] - gv[1][0]),
        (1 - alpha)*(gv[3][1] - gv[0][1]) + alpha*(gv[2][1] - gv[1][1]),
        (1 - alpha)*(gv[3][2] - gv[0][2]) + alpha*(gv[2][2] - gv[1][2])
      }};
    double J[2][2] = {
      {
        lambda[0][0]*lambda[0][0] + lambda[0][1]*lambda[0][1] + lambda[0][2]*lambda[0][2],
        lambda[0][0]*lambda[1][0] + lambda[0][1]*lambda[1][1] + lambda[0][2]*lambda[1][2],
      }, {
        lambda[1][0]*lambda[0][0] + lambda[1][1]*lambda[0][1] + lambda[1][2]*lambda[0][2],
        lambda[1][0]*lambda[1][0] + lambda[1][1]*lambda[1][1] + lambda[1][2]*lambda[1][2]
      }};
    return 0.5*sqrt(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
  }

  double to3d_global_to_local_jacobian(const double * gp, 
                                       const double ** lv,
                                       const double ** gv)
  {
    /// not implemented!
  }

#ifdef __cplusplus
}
#endif

/**
 * end of file
 * 
 */

