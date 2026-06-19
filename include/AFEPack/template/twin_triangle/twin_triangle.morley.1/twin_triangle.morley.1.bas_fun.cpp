/**
 * @file   twin_triangle.morley.1.bas_fun.c
 * @author Zhao Weibo <weibo611@163.com>
 * @date   Mon Dec 15 10:20:59 2008
 * 
 * @brief  
 * 
 * 
 */

#ifdef __cplusplus
extern "C" {
#endif

#define AREA(p0, p1, p2) (((p1)[0] - (p0)[0])*((p2)[1] - (p0)[1]) -     \
                          ((p1)[1] - (p0)[1])*((p2)[0] - (p0)[0]))

  void lambda_1(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[1],v[2]);
    val[0] = AREA(p,v[0],v[1]);
    val[0] /= area;
    if(AREA(p,v[2],v[0])>0)
      val[0] = 1.0 - 2*val[0];
    else
      val[0] = 0;
  }

  void lambda_2(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[1],v[2]);
    val[0] = AREA(p,v[1],v[2]);
    val[0] /= area;
    if(AREA(p,v[2],v[0])>0)
      val[0] = 1.0 - 2*val[0];
    else
      val[0] = 0;
  }

  void lambda_3(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[2],v[3]);
    val[0] = AREA(p,v[2],v[3]);
    val[0] /= area;
    if(AREA(p,v[0],v[2])>0)
      val[0] = 1.0 - 2*val[0];
    else
      val[0] = 0;
  }

  void lambda_4(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[2],v[3]);
    val[0] = AREA(p,v[3],v[0]);
    val[0] /= area;
    if(AREA(p,v[0],v[2])>0)
      val[0] = 1.0 - 2*val[0];
    else
      val[0] = 0;
  }

  void lambda_5(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double lambda_tmp = AREA(v[0],p,v[2]);
    if(lambda_tmp > 0.0){
      val[0] = lambda_tmp;
      double area = AREA(v[0],v[1],v[2]);
      val[0] /= area;
      val[0] = 1.0 - 2*val[0];
    }
    else{
      val[0] = -lambda_tmp;
      double area = AREA(v[0],v[2],v[3]);
      val[0] /= area;
      val[0] = 1.0 - 2*val[0];
    }
  }

  void gradient_lambda_1(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[1],v[2]);
    if(AREA(v[0],p,v[2])>0){
      val[0] = -2.0*(v[0][1] - v[1][1])/area;
      val[1] = -2.0*(v[1][0] - v[0][0])/area;
    }
    else{
      val[0] = 0;
      val[1] = 0;
    }
  }

  void gradient_lambda_2(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[1],v[2]);
    if(AREA(v[0],p,v[2])>0){
      val[0] = -2.0*(v[1][1] - v[2][1])/area;
      val[1] = -2.0*(v[2][0] - v[1][0])/area;
    }
    else{
      val[0] = 0;
      val[1] = 0;
    }
  }

  void gradient_lambda_3(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[2],v[3]);
    if(AREA(v[2],p,v[0])>0){
      val[0] = -2.0*(v[2][1] - v[3][1])/area;
      val[1] = -2.0*(v[3][0] - v[2][0])/area;
    }
    else{
      val[0] = 0;
      val[1] = 0;
    }
  }

  void gradient_lambda_4(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double area = AREA(v[0],v[1],v[2]);
    if(AREA(v[2],p,v[0])>0){
      val[0] = -2.0*(v[3][1] - v[0][1])/area;
      val[1] = -2.0*(v[0][0] - v[3][0])/area;
    }
    else{
      val[0] = 0;
      val[1] = 0;
    }
  }

  void gradient_lambda_5(const double * p, const double ** v, void * value)
  {
    double * val = (double *)value;
    double lambda_val = AREA(v[0],p,v[2]);
    if(lambda_val > 0.0){
      double area = AREA(v[0],v[1],v[2]);
      val[0] = -2.0*(v[2][1] - v[0][1])/area;
      val[1] = -2.0*(v[0][0] - v[2][0])/area;
    }
    else{
      double area = AREA(v[0],v[2],v[3]);
      val[0] = -2.0*(v[0][1] - v[2][1])/area;
      val[1] = -2.0*(v[2][0] - v[0][0])/area;
    }
  }

#ifdef __cplusplus
}
#endif

/**
 * end of file
 * 
 */


