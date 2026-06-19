/**
 * @file   EllipticEquation.cpp
 * @author Robert Lie
 * @date   Wed Feb 21 15:26:48 2007
 *
 * @brief  类 EllipticEquation 和 Matrix 的实现文件。
 *
 */

#include "EllipticEquation.h"

/**
 * 计算二次项系数的函数，表达式为
 *                4, x^2 > y^2;
 *    a(x, y) = {
 *                1, x^2 < y^2.
 *
 * @param p 点的坐标： x = p[0], y = p[1]
 *
 * @return 返回函数 a(x, y) 的值
 *
 */
double _a_(const double * p)
{
  if (p[0]*p[0] > p[1]*p[1])
    return 400.;
  else
    return 1.;
}

double _u_(const double * p)
{
  return cos(2.0 * PI * p[0]) * sin(3.0 * PI * p[1]);
}

double _f_(const double * p)
{
  return _a_(p) * 13 * PI * PI * _u_(p);
}

/**
 * 计算矩阵 Matrix 的单元刚度矩阵。这个函数本来是为了计算双线性型 b(u,
 * v) 而设计的，其中 u 和 v 可以在不同的有限元空间上。因此，这个函数接
 * 受两个单元的参数，分别是 u 和 v 所在的有限元空间的单元，在我们现在
 * 的情况下，这两个参数是同一个单元。第三个参数是为了实现库中的所谓多
 * 套网格方法而预留的，在这里我们可以先忽略这个参数的作用。
 *
 * @param ele0 有限元空间的单元
 * @param ele1 这里是和 ele0 相同的单元
 *
 * 下面的这段程序是比较标准的进行单元上的数值积分的程序的套路，弄清楚
 * 其作用对于掌握 AFEPack 的程序设计是至关重要的。计算一个函数 f(x) 在
 * 一个单元上的数值积分可以使用通用的公式
 *
 *    sum |J| w_l f(x_l) |e|
 *
 * 其中 |e| 是参考单元的体积，w_l 是积分公式中在第 l 个积分点上的权重，
 * 而 |J| 是从参考单元到实际单元进行变换的雅可比行列式，x_l 是第 l 个
 * 积分点。在 AFEPack 中，积分点是首先给在参考单元上的，然后通过坐标变
 * 换映射到实际单元上，成为实际单元上的积分点。下面的程序基本上就是在
 * 实现这样的一个数值积分公式。
 *
 */
void
Matrix::getElementMatrix(const Element<double,DIM>& ele0,
                         const Element<double,DIM>& ele1,
                         const ActiveElementPairIterator<DIM>::State state)
{
  /// 计算参考单元的体积
  double volume = ele0.templateElement().volume();

  /// 取出参考单元上的数值积分公式
  const QuadratureInfo<DIM>& quad_info =
    ele0.findQuadratureInfo(algebricAccuracy());

  /// 计算在积分点上的坐标变换的雅可比矩阵
  std::vector<double> jacobian =
    ele0.local_to_global_jacobian(quad_info.quadraturePoint());

  /// 这是积分点的个数
  int n_quadrature_point = quad_info.n_quadraturePoint();

  /// 将参考单元上的积分点变换到网格中的实际单元上去
  std::vector<AFEPack::Point<DIM> > q_point =
    ele0.local_to_global(quad_info.quadraturePoint());

  /// 计算单元上的基函数的梯度值
  std::vector<std::vector<std::vector<double> > > basis_grad =
    ele0.basis_function_gradient(q_point);

  int n_ele_dof = ele0.dof().size(); /// 单元上的自由度的个数

  /// 对积分点做求和
  for (int l = 0;l < n_quadrature_point;l ++) {
    /// Jxw 就是积分公式中的项 |e| |J| w_l
    double Jxw = quad_info.weight(l)*jacobian[l]*volume;

    /// 计算函数 a(x, y) 在积分点的值
    double a_val = _a_(q_point[l]);

    for (int j = 0;j < n_ele_dof;j ++) { /// 对自由度做循环
      for (int k = 0;k < n_ele_dof;k ++) { /// 对自由度做循环

        /// 将积分点上的值往单元刚度矩阵上进行累加
        elementMatrix(j,k) += Jxw*a_val*innerProduct(basis_grad[j][l],
                                                     basis_grad[k][l]);
      }
    }
  }
}

/**
 * 初始化函数中准备了参考单元以及读入了网格数据文件。
 */
void EllipticEquation::initialize(std::string& filename)
{
  /// 读入参考单元的几何信息及坐标变换
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_coord_transform.readData("triangle.crd_trs");

  /// 在参考单元的几何信息更新以后，模板单元上的自由度分布需要重新初始
  /// 化一下，对于下面的基函数，也是相似的情况。
  triangle_template_dof.reinit(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  triangle_basis_function.reinit(triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");

  /// 用上面读入的信息合成参考单元
  template_element.resize(1);
  template_element[0].reinit(triangle_template_geometry,
                             triangle_template_dof,
                             triangle_coord_transform, 
                             triangle_basis_function);

  /// 读入网格数据文件
  this->readData(filename);
}

/**
 * 构造有限元空间，这里的代码是完全从前一节中搬过来的。
 */
void EllipticEquation::buildSpace()
{
  /// 对有限元空间的对象重新初始化一下
  fem_space.reinit(*this, template_element);
        
  int n_element = this->n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++) {
    fem_space.element(i).reinit(fem_space,i,0);
  }

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();
}

/**
 * 构造线性系统，并求解，得到需要的有限元函数。这一段的程序也基本上是
 * 由前一节搬过来的。
 */
void EllipticEquation::solve()
{
  /// 注意这里我们使用 Matrix 来当我们的刚度矩阵
  Matrix stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 3;
  stiff_matrix.build();

  /// 有限元函数 u_h 初始化
  u_h.reinit(fem_space);

  /**
   * 向量 f_h 用来计算右端项 f 在离散过后得到的向量。
   *
   * 注意：我们需要使用前一节中的函数表达式 _f_ 和 _u_，为了节省篇幅，
   * 我们这个文件中没有重复写这两个函数的表达式。
   */
  Vector<double> f_h;
  Operator::L2Discretize(&_f_, fem_space, f_h, 3);

  /// 使用狄氏边值条件。
  BoundaryFunction<double,DIM> boundary(BoundaryConditionInfo::DIRICHLET,
                                        1, &_u_);
  BoundaryConditionAdmin<double,DIM> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, u_h, f_h);

  /// 调用 AFEPack 中的代数多重网格求解器，求解线性系统。
  AMGSolver solver(stiff_matrix);
  solver.solve(u_h, f_h);      

  /// 将解输出。
  u_h.writeOpenDXData("u.dx");
}

/**
 * end of file
 *
 */
