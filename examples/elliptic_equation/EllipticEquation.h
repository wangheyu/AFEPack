/**
 * @file   EllipticEquation.h
 * @author Robert Lie
 * @date   Wed Feb 21 11:45:03 2007
 *
 * @brief  类 EllipticEquation 的声明，附带的我们还声明了刚度矩阵的类
 *         型 Matrix。
 *
 */

#ifndef __EllipticEquation_h__
#define __EllipticEquation_h__

#include <AFEPack/AMGSolver.h>        /// 代数多重网格求解器
#include <AFEPack/TemplateElement.h>  /// 参考单元
#include <AFEPack/FEMSpace.h>         /// 有限元空间
#include <AFEPack/Operator.h>         /// 算子
#include <AFEPack/BilinearOperator.h> /// 双线性算子
#include <AFEPack/Functional.h>       /// 泛函
#include <AFEPack/EasyMesh.h>         /// EasyMesh 接口

#define DIM 2
#define PI 3.14159265358979323846

/**
 * 类 EllipticEquation 的声明，可以看到，在这个类的声明中，我们基本上
 * 只是把前一节的程序中的很多变量，搬到了这里作为这个类的成员变量。注
 * 意我们使得本类是从 EasyMesh 上派生出来的，从而本类的对象本身可以作
 * 为一个网格来使用，同时这是为了下面扩展到移动网格方法的方便。
 */
class EllipticEquation : public EasyMesh
{
 private:
  /**
   * 下面的几个变量用来构造参考单元，已经在上一节进行过详细解释。
   */
  TemplateGeometry<DIM> triangle_template_geometry;
  CoordTransform<DIM,DIM> triangle_coord_transform;
  TemplateDOF<DIM> triangle_template_dof;
  BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function;

  std::vector<TemplateElement<double,DIM,DIM> > template_element;

  FEMSpace<double,DIM> fem_space; /// 有限元空间
  FEMFunction<double,DIM> u_h;    /// 有限元函数

 public:
  /**
   * 这个函数进行初始化的工作，包括读入网格数据文件，以及构造参考单元，
   * 为建立有限元空间做准备。
   */
  void initialize(std::string&);
  /**
   * 这个函数将构造有限元空间。
   */
  void buildSpace();
  /**
   * 在构造的有限元空间上计算刚度矩阵并求解获得的线性系统，最终完成求
   * 解近似解。
   */
  void solve();
};

/**
 * 这个类是用来计算现在的变系数二次算子的刚度矩阵用的。它从库中的刚度
 * 矩阵类 StiffMatrix 派生出来，唯一需要进行重载的函数就是计算单元刚度
 * 矩阵的函数 getElementMatrix。
 */
class Matrix : public StiffMatrix<DIM,double>
{
 public:
 Matrix(FEMSpace<double,DIM>& sp) :
  StiffMatrix<DIM,double>(sp) {};
  virtual ~Matrix() {};
 public:
  virtual void
    getElementMatrix(const Element<double,2>& ele0,
                     const Element<double,2>& ele1,
                     const ActiveElementPairIterator<DIM>::State state);
};

#endif // __EllipticEquation_h__

/**
 * end of file
 *
 */
