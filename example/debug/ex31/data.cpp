/**
 * @file   data.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Sat Oct 24 15:46:56 2009
 * 
 * @brief  
 * 
 * 
 */

#include "SCL.h"
#define PI M_PIl

double initial_value(const double * x, const double& t)
{
  return sin(2*PI*(x[0] - t))*cos(2*PI*(x[1] - t));
}

std::vector<double> flux(const AFEPack::Point<DIM>& x, const double& u)
{
  std::vector<double> val(DIM);
  val[0] = u*u; val[1] = exp(u);
  return val;
}

double velocity(const AFEPack::Point<DIM>& x, const double& u)
{
  std::vector<double> val(DIM);
  val[0] = 2*u; val[1] = exp(u);
  return sqrt(val[0]*val[0] + val[1]*val[1]);
}

void boundary_value(const AFEPack::Point<DIM>& x,
		    double& u,
                    const std::vector<double>& n,
		    const double& t,
		    const int& bm)
{
  u = initial_value(x, t);
}

void init_data(SCL& app,
               int argc,
               char * argv[])
{
  app.meshfile = "../C";
  if (argc >= 2) app.meshfile = argv[1];
  app.end_t = 1.0;
  if (argc >= 3) app.end_t = atof(argv[2]);

  app.t = 0.0;
  app._u_0_ = &initial_value;
  app._f_ = &flux;
  app._v_ = &velocity;
  app._u_b_ = &boundary_value;
}

/**
 * end of file
 * 
 */
