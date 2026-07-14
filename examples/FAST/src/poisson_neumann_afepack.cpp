/**
 * @file poisson_neumann_afepack.cpp
 * @brief AFEPack Poisson 方程算例：纯 Neumann 问题。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 纯 Neumann 问题：处理常数核空间和兼容性条件，验证压力/势函数零均值约束思想。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct NeumannResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    double       rhs_mean_before = 0.0;
    double       rhs_mean_after = 0.0;
    SolveInfo    cg;
    SolveInfo    jacobi_pcg;
    double       cg_l2_error = 0.0;
    double       cg_h1_error = 0.0;
    double       pcg_l2_error = 0.0;
    double       pcg_h1_error = 0.0;
  };

  double exact_solution(const double *p)
  {
    return std::cos(pi * p[0]) + std::cos(pi * p[1]);
  }

  std::vector<double> exact_gradient(const double *p)
  {
    return {-pi * std::sin(pi * p[0]), -pi * std::sin(pi * p[1])};
  }

  double right_hand_side(const double *p)
  {
    return pi * pi * exact_solution(p);
  }

  std::string default_template_path()
  {
    if (const char *env = std::getenv("AFEPACK_PATH"))
      return std::string(env) + "/template/triangle";
    return "/home/steve/local/AFEPack/include/AFEPack/template/triangle";
  }

  void ensure_template_path()
  {
    if (std::getenv("AFEPACK_TEMPLATE_PATH") != nullptr)
      return;

    const std::string path = default_template_path();
    if (::setenv("AFEPACK_TEMPLATE_PATH", path.c_str(), 0) != 0)
      throw std::runtime_error("failed to set AFEPACK_TEMPLATE_PATH");
  }

  std::vector<TemplateElement<double, dim, dim>>
  make_triangle_template(TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle.1.tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle.1.bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
  }

  FEMSpace<double, dim> build_space(
    EasyMesh &mesh,
    std::vector<TemplateElement<double, dim, dim>> &template_element)
  {
    FEMSpace<double, dim> fem_space;
    fem_space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      fem_space.element(i).reinit(fem_space, i, 0);

    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    return fem_space;
  }

  double arithmetic_mean(const Vector<double> &v)
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i)
      sum += v(i);
    return sum / static_cast<double>(v.size());
  }

  void project_arithmetic_mean_zero(Vector<double> &v)
  {
    const double mean = arithmetic_mean(v);
    for (std::size_t i = 0; i < v.size(); ++i)
      v(i) -= mean;
  }

  void apply_preconditioner(const SparseMatrix<double> &matrix,
                            const Vector<double> &src,
                            Vector<double> &dst,
                            const bool use_jacobi)
  {
    dst.reinit(src.size(), false);
    if (!use_jacobi)
      {
        dst = src;
        return;
      }

    for (std::size_t i = 0; i < src.size(); ++i)
      {
        const double diagonal = matrix.diag_element(i);
        dst(i) = std::abs(diagonal) > 1.0e-14 ? src(i) / diagonal : src(i);
      }
    project_arithmetic_mean_zero(dst);
  }

  SolveInfo solve_projected_cg(const SparseMatrix<double> &matrix,
                               Vector<double> &x,
                               const Vector<double> &compatible_rhs,
                               const bool use_jacobi,
                               const double tolerance,
                               const int max_iterations)
  {
    const double rhs_norm = std::max(compatible_rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(compatible_rhs.size());
    matrix.vmult(ax, x);

    Vector<double> r = compatible_rhs;
    r.add(-1.0, ax);
    project_arithmetic_mean_zero(r);

    double relative_residual = r.l2_norm() / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    Vector<double> z(compatible_rhs.size());
    apply_preconditioner(matrix, r, z, use_jacobi);
    Vector<double> p = z;
    project_arithmetic_mean_zero(p);
    Vector<double> ap(compatible_rhs.size());

    double rz_old = r.dot(z);
    if (std::abs(rz_old) < 1.0e-30)
      return {false, 0, relative_residual};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        matrix.vmult(ap, p);
        project_arithmetic_mean_zero(ap);

        const double denominator = p.dot(ap);
        if (std::abs(denominator) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        const double alpha = rz_old / denominator;
        x.add(alpha, p);
        project_arithmetic_mean_zero(x);

        r.add(-alpha, ap);
        project_arithmetic_mean_zero(r);

        relative_residual = r.l2_norm() / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, iteration, relative_residual};

        apply_preconditioner(matrix, r, z, use_jacobi);
        const double rz_new = r.dot(z);
        if (std::abs(rz_old) < 1.0e-30)
          return {false, iteration, relative_residual};

        const double beta = rz_new / rz_old;
        p.scale(beta);
        p.add(1.0, z);
        project_arithmetic_mean_zero(p);
        rz_old = rz_new;
      }

    return {false, max_iterations, relative_residual};
  }

  double finite_element_mean(const FEMFunction<double, dim> &solution)
  {
    const auto &space = solution.femSpace();
    double integral = 0.0;
    double area = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info = element.findQuadratureInfo(4);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = solution.value(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            integral += jxw * values[q];
            area += jxw;
          }
      }

    return integral / std::max(area, 1.0e-30);
  }

  void subtract_continuous_mean(FEMFunction<double, dim> &solution)
  {
    const double mean = finite_element_mean(solution);
    for (std::size_t i = 0; i < solution.size(); ++i)
      solution(i) -= mean;
  }

  void copy_to_function(const Vector<double> &values,
                        FEMFunction<double, dim> &solution)
  {
    for (std::size_t i = 0; i < values.size(); ++i)
      solution(i) = values(i);
  }

  NeumannResult solve_neumann_case(const std::string &mesh_file,
                                   const double h)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto fem_space = build_space(mesh, template_element);

    StiffMatrix<dim, double> stiffness(fem_space);
    stiffness.algebricAccuracy() = 4;
    stiffness.build();

    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    const double rhs_mean_before = arithmetic_mean(rhs);
    project_arithmetic_mean_zero(rhs);
    const double rhs_mean_after = arithmetic_mean(rhs);

    Vector<double> cg_solution(rhs.size());
    cg_solution = 0.0;
    const SolveInfo cg =
      solve_projected_cg(stiffness, cg_solution, rhs, false, 1.0e-10, 2000);
    if (!cg.converged)
      throw std::runtime_error("projected CG failed");

    Vector<double> pcg_solution(rhs.size());
    pcg_solution = 0.0;
    const SolveInfo jacobi_pcg =
      solve_projected_cg(stiffness, pcg_solution, rhs, true, 1.0e-10, 2000);
    if (!jacobi_pcg.converged)
      throw std::runtime_error("projected Jacobi-PCG failed");

    FEMFunction<double, dim> cg_function(fem_space);
    copy_to_function(cg_solution, cg_function);
    subtract_continuous_mean(cg_function);

    FEMFunction<double, dim> pcg_function(fem_space);
    copy_to_function(pcg_solution, pcg_function);
    subtract_continuous_mean(pcg_function);

    FunctionFunction<double> exact(&exact_solution, &exact_gradient);

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            rhs_mean_before,
            rhs_mean_after,
            cg,
            jacobi_pcg,
            Functional::L2Error(cg_function, exact, 4),
            Functional::H1SemiError(cg_function, exact, 4),
            Functional::L2Error(pcg_function, exact, 4),
            Functional::H1SemiError(pcg_function, exact, 4)};
  }

  double convergence_rate(const double coarse_error,
                          const double fine_error,
                          const double coarse_h,
                          const double fine_h)
  {
    if (coarse_error <= 0.0 || fine_error <= 0.0)
      return 0.0;
    return std::log(coarse_error / fine_error) / std::log(coarse_h / fine_h);
  }
}

int main()
{
  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack pure Neumann Poisson projected CG test\n";
      std::cout << "u=cos(pi x)+cos(pi y), homogeneous Neumann boundary, "
                   "zero-mean projected algebraic solve.\n";
      std::cout << "h cells dofs rhs_mean_before rhs_mean_after CG CG_relres "
                   "JPCG JPCG_relres L2 H1 L2_rate H1_rate\n";
      std::cout << std::flush;

      const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
      const std::vector<std::pair<double, std::string>> meshes = {
        {0.20, mesh_dir + "/unit_square_h0p20"},
        {0.10, mesh_dir + "/unit_square_h0p10"},
        {0.05, mesh_dir + "/unit_square_h0p05"},
        {0.025, mesh_dir + "/unit_square_h0p025"}};

      std::vector<NeumannResult> results;
      results.reserve(meshes.size());
      for (const auto &[h, mesh_file] : meshes)
        results.push_back(solve_neumann_case(mesh_file, h));

      for (std::size_t i = 0; i < results.size(); ++i)
        {
          const auto &result = results[i];
          std::cout << result.h << ' '
                    << result.cells << ' '
                    << result.dofs << ' '
                    << result.rhs_mean_before << ' '
                    << result.rhs_mean_after << ' '
                    << result.cg.iterations << ' '
                    << result.cg.relative_residual << ' '
                    << result.jacobi_pcg.iterations << ' '
                    << result.jacobi_pcg.relative_residual << ' '
                    << result.pcg_l2_error << ' '
                    << result.pcg_h1_error << ' ';
          if (i == 0)
            std::cout << "- -\n";
          else
            std::cout << convergence_rate(results[i - 1].pcg_l2_error,
                                          result.pcg_l2_error,
                                          results[i - 1].h,
                                          result.h)
                      << ' '
                      << convergence_rate(results[i - 1].pcg_h1_error,
                                          result.pcg_h1_error,
                                          results[i - 1].h,
                                          result.h)
                      << '\n';
        }

      if (!(results.back().pcg_l2_error < results.front().pcg_l2_error &&
            results.back().pcg_h1_error < results.front().pcg_h1_error))
        throw std::runtime_error("projected Neumann errors should decrease");
      if (!(std::abs(results.back().rhs_mean_after) < 1.0e-13))
        throw std::runtime_error("projected Neumann RHS is not mean-free");
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
