/**
 * @file   triangle.RT.2.bas_fun.cpp
 * @author Huayi Wei <weihy1984@163.com>
 * @date   Mon Dec 17 08:36:15 2007
 * 
 * @brief  
 * 
 * 
 */

#include <cmath>
#include <Miscellaneous.h>

#ifdef __cplusplus
extern "C" {
#endif


  /**
   * 三角形上 1 阶的 Raviart-Thomas 元的基函数及其导数的定义。
   *
   *
   * 我们这里将这种单元的描述文件命名为 triangle.RT.2.* 而不是
   * triangle.RT.1.* ，是因为这些基函数事实上是二次函数。
   * 
   * 因为这个文件中用到了类 nVector<2,double> ，必须使用 C++编译器
   * 进行编译。
   *
   */

#define vector_length 2
#define vt nVector<vector_length,double>
#define r2 sqrt(2)

#define detJ ((v[1][0]-v[0][0])*(v[2][1]-v[0][1])-(v[2][0]-v[0][0])*(v[1][1]-v[0][1]))

#define GET_J                                                       \
       double B_k[2][2];                                            \
       B_k[0][0]=(v[1][0]-v[0][0]);                                 \
       B_k[0][1]=(v[2][0]-v[0][0]);                                 \
       B_k[1][0]=(v[1][1]-v[0][1]);                                 \
       B_k[1][1]=(v[2][1]-v[0][1]);                                 \

#define  GET_BASIS                                                  \
        val[0]=(B_k[0][0]*temp[0]+B_k[0][1]*temp[1])*l/detJ;        \
        val[1]=(B_k[1][0]*temp[0]+B_k[1][1]*temp[1])*l/detJ;        \


  /*计算Jacobi矩阵的逆矩阵*/

#define GET_JMIT                                            \
       double Jm[2][2],temp1[2][2];                         \
       Jm[0][0]=( v[2][1]-v[0][1])/detJ;                    \
       Jm[0][1]=(-v[2][0]+v[0][0])/detJ;                    \
       Jm[1][0]=(-v[1][1]+v[0][1])/detJ;                    \
       Jm[1][1]=( v[1][0]-v[0][0])/detJ;                    \
       temp1[0][0]=temp[0][0]*Jm[0][0]+temp[0][1]*Jm[1][0]; \
       temp1[0][1]=temp[0][0]*Jm[0][1]+temp[0][1]*Jm[1][1]; \
       temp1[1][0]=temp[1][0]*Jm[0][0]+temp[1][1]*Jm[1][0]; \
       temp1[1][1]=temp[1][0]*Jm[0][1]+temp[1][1]*Jm[1][1]; \

#define GET_GRAD                                                  \
 val[0][0]=(B_k[0][0]*temp1[0][0]+B_k[0][1]*temp1[1][0])*l/detJ;  \
 val[0][1]=(B_k[0][0]*temp1[0][1]+B_k[0][1]*temp1[1][1])*l/detJ;  \
 val[1][0]=(B_k[1][0]*temp1[0][0]+B_k[1][1]*temp1[1][0])*l/detJ;  \
 val[1][1]=(B_k[1][0]*temp1[0][1]+B_k[1][1]*temp1[1][1])*l/detJ;  \


#define GET_L(i, j)                                                            \
  double l0 = (v[i][0] - v[j][0]);                                             \
  double l1 = (v[i][1] - v[j][1]);                                             \
  double l =1.0;                                                               \
  if (fabs(l0) > fabs(l1)) {                                                   \
    if (l0 < 0) l = -l;                                                        \
  } else {                                                                     \
    if (l1 < 0) l = -l;                                                        \
  }




void get_lambda(const double * p, const double ** v, double * lambda)
{
      
	double area =(v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
        
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;

};



void lambda_01(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];

  get_lambda(p,v,lambda);

  GET_L(1, 0);

  GET_J;

 
  temp[0]= 3.0*lambda[1]-4.0*lambda[1]*(lambda[1]+lambda[2]);
  temp[1]=-2.0+3.0*lambda[1]+6.0*lambda[2]-4.0*lambda[2]*(lambda[1]+lambda[2]);

  GET_BASIS;
 
};

void lambda_02(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];
  
  get_lambda(p,v,lambda);

  GET_L(1, 0);

  GET_J;

  temp[0]= -6.0*lambda[1]+8.0*lambda[1]*(lambda[1]+lambda[2]);
  temp[1]=1.0-3.0*lambda[1]+3.0*lambda[2] - 4.0*lambda[2]*(lambda[1]+lambda[2]);


  GET_BASIS;
 
};

void lambda_11(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];

  get_lambda(p,v,lambda);

  GET_L(2, 1);

  GET_J;
  temp[0]= (-6.0*lambda[1]+8.0*lambda[1]*(lambda[1]+lambda[2]));
  temp[1]= (3.0*lambda[2]-4.0*lambda[2]*(lambda[1]+lambda[2]));



  GET_BASIS;
 
};

void lambda_12(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];

  get_lambda(p,v,lambda);

  GET_L(2, 1);

  GET_J;
  temp[0]=( 3.0*lambda[1]-4.0*lambda[1]*(lambda[1]+lambda[2]));
  temp[1]=(-6.0*lambda[2]+8.0*lambda[2]*(lambda[1]+lambda[2]));
 



  GET_BASIS;
 
};
void lambda_21(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];

  get_lambda(p,v,lambda);

  GET_L(0, 2);

  GET_J;


  temp[0]=1.0+3.0*lambda[1]-3.0*lambda[2]-4.0*lambda[1]*(lambda[1]+lambda[2]);
  temp[1]=-6.0*lambda[2]+8.0*lambda[2]*(lambda[1]+lambda[2]);
  GET_BASIS;
 
};

void lambda_22(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];

  get_lambda(p,v,lambda);

  GET_L(0, 2);

  GET_J;

  temp[0]=-2.0+6.0*lambda[1]+3.0*lambda[2]-4.0*lambda[1]*(lambda[1]+lambda[2]);
  temp[1]=3.0*lambda[2]-4.0*lambda[2]*(lambda[1]+lambda[2]);
 
  GET_BASIS;
  
};

void lambda_31(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3],temp[2];
  double l=1.0;
  get_lambda(p,v,lambda);

  GET_J;

  temp[0]=24.0*lambda[1]-24.0*lambda[1]*(lambda[1]+lambda[2]);
  temp[1]=0.0;

  GET_BASIS;


 };

void lambda_32(const double * p, const double ** v, void * value)
{
   vt& val = *((vt *)value);
  double lambda[3],temp[2];
  double l=1.0; 
  get_lambda(p,v,lambda);

  GET_J;
  temp[0]=0.0;
  temp[1]=24.0*lambda[2]-24.0*lambda[2]*(lambda[1]+lambda[2]);
 
  GET_BASIS;

  };


void gradient_lambda_01(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(1, 0);

  GET_J;

  temp[0][0]=3.0-8.0*lambda[1]-4.0*lambda[2];     
  temp[0][1]=-4.0*lambda[1];
  temp[1][0]= 3.0-4.0*lambda[2];
  temp[1][1]=6.0-4.0*lambda[1]-8.0*lambda[2];
  GET_JMIT; 
  GET_GRAD;
};


void gradient_lambda_02(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(1, 0);

  GET_J;

  temp[0][0]=-6.0+16.0*lambda[1]+8.0*lambda[2];
  temp[0][1]=8.0*lambda[1];
  temp[1][0]=-3.0-4.0*lambda[2];
  temp[1][1]=3.0-4.0*lambda[1]-8.0*lambda[2];
  GET_JMIT; 
  GET_GRAD;
};




void gradient_lambda_11(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(2, 1);

  GET_J;
  temp[0][0]=(-6.0+16.0*lambda[1]+8.0*lambda[2]);
  temp[0][1]=(8.0*lambda[1]);
  temp[1][0]=(-4.0*lambda[2]);
  temp[1][1]=(3.0-4.0*lambda[1]-8.0*lambda[2]);
  GET_JMIT; 
  GET_GRAD;
};



void gradient_lambda_12(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(2, 1);

  GET_J;
  temp[0][0]=(3.0-8.0*lambda[1]-4.0*lambda[2]);
  temp[0][1]=(-4.0*lambda[1]);
  temp[1][0]=(8.0*lambda[2]);
  temp[1][1]=(-6.0+8.0*lambda[1]+16.0*lambda[2]);


  GET_JMIT; 
  GET_GRAD;
};


void gradient_lambda_21(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(0, 2);

  GET_J;

	temp[0][0]= 3.0-8.0*lambda[1]-4.0*lambda[2];
	temp[0][1]=-3.0-4.0*lambda[1];
	temp[1][0]=  8.0*lambda[2];
	temp[1][1]=-6.0+8.0*lambda[1]+16.0*lambda[2];


  GET_JMIT; 
  GET_GRAD;
};


void gradient_lambda_22(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  get_lambda(p,v,lambda);

  GET_L(0, 2);

  GET_J;

	temp[0][0]= 6.0-8.0*lambda[1]-4.0*lambda[2];
	temp[0][1]= 3.0-4.0*lambda[1];
	temp[1][0]= -4.0*lambda[2];
	temp[1][1]=3.0-4.0*lambda[1]-8.0*lambda[2];
  GET_JMIT; 
  GET_GRAD;
};

void gradient_lambda_31(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];
  double l=1.0;
  get_lambda(p,v,lambda);

  GET_J;

  temp[0][0]= 24.0-48.0*lambda[1]-24.0*lambda[2];
  temp[0][1]= -24.0*lambda[1];
  temp[1][0]= 0.0;
  temp[1][1]=0.0;
  GET_JMIT; 
  GET_GRAD;
};
void gradient_lambda_32(const double * p, const double ** v,void * value)
{
  vt * val = (vt *)value;
  double lambda[3],temp[2][2];

  double l=1.0;

  get_lambda(p,v,lambda);
  GET_J;

  temp[0][0]= 0.0;
  temp[0][1]= 0.0;
  temp[1][0]= -24.0*lambda[2];
  temp[1][1]=24.0-24.0*lambda[1]-48.0*lambda[2];

  GET_JMIT; 
  GET_GRAD;
};





























#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
