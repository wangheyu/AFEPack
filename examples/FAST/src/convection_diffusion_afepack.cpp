/**
 * @file convection_diffusion_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：基础离散流程。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 基础离散流程：给出模型问题从网格读取、有限元空间构造、线性系统装配到结果输出的完整流程。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/BilinearOperator.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

#ifndef AFEPACK_DEFAULT_OUTPUT_DIR
#  define AFEPACK_DEFAULT_OUTPUT_DIR "build"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;
  constexpr double beta_x = 1.0;
  constexpr double beta_y = 1.0;

  double current_epsilon = 1.0e-2;

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct CaseResult
  {
    double       epsilon = 0.0;
    bool         supg = false;
    int          cells = 0;
    unsigned int dofs = 0;
    int          iterations = 0;
    double       relative_residual = 0.0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
  };

  double exact_solution(const double *p)
  {
    return std::sin(pi * p[0]) * std::sin(pi * p[1]);
  }

  std::vector<double> exact_gradient(const double *p)
  {
    return {pi * std::cos(pi * p[0]) * std::sin(pi * p[1]),
            pi * std::sin(pi * p[0]) * std::cos(pi * p[1])};
  }

  double right_hand_side(const double *p)
  {
    const double u = exact_solution(p);
    const auto   grad_u = exact_gradient(p);
    return 2.0 * current_epsilon * pi * pi * u +
           beta_x * grad_u[0] + beta_y * grad_u[1];
  }

  double beta_dot(const std::vector<double> &gradient)
  {
    return beta_x * gradient[0] + beta_y * gradient[1];
  }

  double beta_norm()
  {
    return std::sqrt(beta_x * beta_x + beta_y * beta_y);
  }

  double coth(const double x)
  {
    return 1.0 / std::tanh(x);
  }

  double supg_delta(const double epsilon, const double cell_diameter)
  {
    const double norm = beta_norm();
    if (norm == 0.0)
      return 0.0;

    const double peclet = norm * cell_diameter / (2.0 * epsilon);
    if (peclet < 1.0e-8)
      return cell_diameter * cell_diameter / (12.0 * epsilon);

    return cell_diameter / (2.0 * norm) * (coth(peclet) - 1.0 / peclet);
  }

  double cell_diameter(const Element<double, dim> &element)
  {
    std::vector<Point<dim>> vertices;
    element.buildVertexArray(vertices);

    double diameter = 0.0;
    for (std::size_t i = 0; i < vertices.size(); ++i)
      for (std::size_t j = i + 1; j < vertices.size(); ++j)
        diameter = std::max(diameter, distance(vertices[i], vertices[j]));

    return diameter;
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

  class ConvectionDiffusionMatrix : public BilinearOperator<dim, double>
  {
  public:
    ConvectionDiffusionMatrix(FEMSpace<double, dim> &space,
                              const double epsilon,
                              const bool use_supg)
      : BilinearOperator<dim, double>(space, space),
        epsilon(epsilon),
        use_supg(use_supg)
    {}

    void getElementMatrix(
      const Element<double, dim> &ele0,
      const Element<double, dim> &,
      const ActiveElementPairIterator<dim>::State =
        ActiveElementPairIterator<dim>::EQUAL) override
    {
      const double volume = ele0.templateElement().volume();
      const QuadratureInfo<dim> &quad_info =
        ele0.findQuadratureInfo(algebricAccuracy());
      const auto jacobian =
        ele0.local_to_global_jacobian(quad_info.quadraturePoint());
      const auto q_point = ele0.local_to_global(quad_info.quadraturePoint());
      const auto basis_value = ele0.basis_function_value(q_point);
      const auto basis_grad = ele0.basis_function_gradient(q_point);

      const int n_q = quad_info.n_quadraturePoint();
      const int n_dof = ele0.dof().size();
      const double delta =
        use_supg ? supg_delta(epsilon, cell_diameter(ele0)) : 0.0;

      for (int q = 0; q < n_q; ++q)
        {
          const double jxw = quad_info.weight(q) * jacobian[q] * volume;
          for (int i = 0; i < n_dof; ++i)
            {
              const double beta_grad_test = beta_dot(basis_grad[i][q]);
              for (int j = 0; j < n_dof; ++j)
                {
                  const double diffusion =
                    epsilon * innerProduct(basis_grad[i][q],
                                           basis_grad[j][q]);
                  const double convection =
                    beta_dot(basis_grad[j][q]) * basis_value[i][q];
                  const double streamline =
                    delta * beta_dot(basis_grad[j][q]) * beta_grad_test;

                  elementMatrix(i, j) +=
                    jxw * (diffusion + convection + streamline);
                }
            }
        }
    }

  private:
    const double epsilon;
    const bool   use_supg;
  };

  void add_supg_rhs(const FEMSpace<double, dim> &space,
                    Vector<double> &rhs,
                    const double epsilon,
                    const bool use_supg)
  {
    if (!use_supg)
      return;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(4);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto basis_grad = element.basis_function_gradient(q_point);

        const double delta = supg_delta(epsilon, cell_diameter(element));
        const int n_q = quad_info.n_quadraturePoint();
        const int n_dof = element.dof().size();
        const auto &dof = element.dof();

        for (int q = 0; q < n_q; ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double f_value = right_hand_side(q_point[q]);
            for (int i = 0; i < n_dof; ++i)
              rhs(dof[i]) += jxw * delta * f_value *
                             beta_dot(basis_grad[i][q]);
          }
      }
  }

  Vector<double> residual(const SparseMatrix<double> &matrix,
                          const Vector<double> &x,
                          const Vector<double> &rhs)
  {
    Vector<double> result;
    result.reinit(rhs.size(), false);
    matrix.vmult(result, x);
    result.sadd(-1.0, rhs);
    return result;
  }

  void jacobi_apply(const SparseMatrix<double> &matrix,
                    const Vector<double> &src,
                    Vector<double> &dst)
  {
    dst.reinit(src.size(), false);
    for (std::size_t i = 0; i < src.size(); ++i)
      {
        const double diagonal = matrix.el(i, i);
        dst(i) = std::abs(diagonal) > 1.0e-14 ? src(i) / diagonal : src(i);
      }
  }

  SolveInfo solve_bicgstab(const SparseMatrix<double> &matrix,
                           Vector<double> &x,
                           const Vector<double> &rhs,
                           const double tolerance,
                           const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);
    Vector<double> r = residual(matrix, x, rhs);
    const Vector<double> r_hat = r;

    Vector<double> p(rhs.size());
    Vector<double> v(rhs.size());
    Vector<double> s(rhs.size());
    Vector<double> t(rhs.size());
    Vector<double> p_hat(rhs.size());
    Vector<double> s_hat(rhs.size());
    p = 0.0;
    v = 0.0;

    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    double relres = r.l2_norm() / rhs_norm;
    if (relres <= tolerance)
      return {true, 0, relres};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        const double rho_new = r_hat.dot(r);
        if (std::abs(rho_new) < 1.0e-30)
          return {false, iteration - 1, relres};

        if (iteration == 1)
          p = r;
        else
          {
            const double beta =
              (rho_new / rho_old) * (alpha / omega);
            for (std::size_t i = 0; i < p.size(); ++i)
              p(i) = r(i) + beta * (p(i) - omega * v(i));
          }

        jacobi_apply(matrix, p, p_hat);
        matrix.vmult(v, p_hat);

        const double alpha_denominator = r_hat.dot(v);
        if (std::abs(alpha_denominator) < 1.0e-30)
          return {false, iteration - 1, relres};
        alpha = rho_new / alpha_denominator;

        s = r;
        s.add(-alpha, v);

        relres = s.l2_norm() / rhs_norm;
        if (relres <= tolerance)
          {
            x.add(alpha, p_hat);
            return {true, iteration, relres};
          }

        jacobi_apply(matrix, s, s_hat);
        matrix.vmult(t, s_hat);

        const double omega_denominator = t.dot(t);
        if (std::abs(omega_denominator) < 1.0e-30)
          return {false, iteration - 1, relres};
        omega = t.dot(s) / omega_denominator;

        x.add(alpha, p_hat);
        x.add(omega, s_hat);

        r = s;
        r.add(-omega, t);

        relres = r.l2_norm() / rhs_norm;
        if (relres <= tolerance)
          return {true, iteration, relres};
        if (std::abs(omega) < 1.0e-30)
          return {false, iteration, relres};

        rho_old = rho_new;
      }

    return {false, max_iterations, relres};
  }

  void normalize(Vector<double> &vector, const double norm)
  {
    for (std::size_t i = 0; i < vector.size(); ++i)
      vector(i) /= norm;
  }

  void add_scaled(Vector<double> &destination,
                  const double coefficient,
                  const Vector<double> &source)
  {
    for (std::size_t i = 0; i < destination.size(); ++i)
      destination(i) += coefficient * source(i);
  }

  void apply_plane_rotation(double &dx,
                            double &dy,
                            const double cs,
                            const double sn)
  {
    const double temporary = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temporary;
  }

  void generate_plane_rotation(const double dx,
                               const double dy,
                               double &cs,
                               double &sn)
  {
    if (dy == 0.0)
      {
        cs = 1.0;
        sn = 0.0;
      }
    else if (std::abs(dy) > std::abs(dx))
      {
        const double temporary = dx / dy;
        sn = 1.0 / std::sqrt(1.0 + temporary * temporary);
        cs = temporary * sn;
      }
    else
      {
        const double temporary = dy / dx;
        cs = 1.0 / std::sqrt(1.0 + temporary * temporary);
        sn = temporary * cs;
      }
  }

  std::vector<double> solve_upper_hessenberg(
    const std::vector<std::vector<double>> &h,
    const std::vector<double> &g,
    const int size)
  {
    std::vector<double> y(size, 0.0);
    for (int i = size - 1; i >= 0; --i)
      {
        double sum = g[i];
        for (int j = i + 1; j < size; ++j)
          sum -= h[i][j] * y[j];
        y[i] = sum / h[i][i];
      }
    return y;
  }

  SolveInfo solve_gmres(const SparseMatrix<double> &matrix,
                        Vector<double> &x,
                        const Vector<double> &rhs,
                        const double tolerance,
                        const int restart,
                        const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);
    Vector<double> r = residual(matrix, x, rhs);
    double beta = r.l2_norm();
    double relative_residual = beta / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    int total_iterations = 0;
    while (total_iterations < max_iterations)
      {
        std::vector<Vector<double>> v(restart + 1, Vector<double>(rhs.size()));
        std::vector<std::vector<double>> h(
          restart + 1, std::vector<double>(restart, 0.0));
        std::vector<double> cs(restart, 0.0);
        std::vector<double> sn(restart, 0.0);
        std::vector<double> g(restart + 1, 0.0);

        v[0] = r;
        normalize(v[0], beta);
        g[0] = beta;

        int inner_iterations = 0;
        for (; inner_iterations < restart &&
               total_iterations < max_iterations;
             ++inner_iterations, ++total_iterations)
          {
            matrix.vmult(v[inner_iterations + 1], v[inner_iterations]);

            for (int i = 0; i <= inner_iterations; ++i)
              {
                h[i][inner_iterations] =
                  v[inner_iterations + 1].dot(v[i]);
                add_scaled(v[inner_iterations + 1],
                           -h[i][inner_iterations],
                           v[i]);
              }

            h[inner_iterations + 1][inner_iterations] =
              v[inner_iterations + 1].l2_norm();
            if (h[inner_iterations + 1][inner_iterations] != 0.0)
              normalize(v[inner_iterations + 1],
                        h[inner_iterations + 1][inner_iterations]);

            for (int i = 0; i < inner_iterations; ++i)
              apply_plane_rotation(h[i][inner_iterations],
                                   h[i + 1][inner_iterations],
                                   cs[i],
                                   sn[i]);

            generate_plane_rotation(h[inner_iterations][inner_iterations],
                                    h[inner_iterations + 1][inner_iterations],
                                    cs[inner_iterations],
                                    sn[inner_iterations]);
            apply_plane_rotation(h[inner_iterations][inner_iterations],
                                 h[inner_iterations + 1][inner_iterations],
                                 cs[inner_iterations],
                                 sn[inner_iterations]);
            apply_plane_rotation(g[inner_iterations],
                                 g[inner_iterations + 1],
                                 cs[inner_iterations],
                                 sn[inner_iterations]);

            relative_residual =
              std::abs(g[inner_iterations + 1]) / rhs_norm;
            if (relative_residual <= tolerance)
              {
                const int system_size = inner_iterations + 1;
                const auto y =
                  solve_upper_hessenberg(h, g, system_size);
                for (int i = 0; i < system_size; ++i)
                  add_scaled(x, y[i], v[i]);
                return {true, total_iterations + 1, relative_residual};
              }
          }

        const auto y =
          solve_upper_hessenberg(h, g, inner_iterations);
        for (int i = 0; i < inner_iterations; ++i)
          add_scaled(x, y[i], v[i]);

        r = residual(matrix, x, rhs);
        beta = r.l2_norm();
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  CaseResult run_case(const std::string &mesh_file,
                      const std::string &output_file,
                      const double epsilon,
                      const bool use_supg,
                      const bool write_solution)
  {
    current_epsilon = epsilon;

    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);

    auto fem_space = build_space(mesh, template_element);

    ConvectionDiffusionMatrix matrix(fem_space, epsilon, use_supg);
    matrix.algebricAccuracy() = 4;
    matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    add_supg_rhs(fem_space, rhs, epsilon, use_supg);

    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(fem_space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(matrix, solution, rhs);

    for (std::size_t i = 0; i < solution.size(); ++i)
      solution(i) = 0.0;
    const SolveInfo solve_info =
      solve_gmres(matrix, solution, rhs, 1.0e-10, 50, 2000);
    if (!solve_info.converged)
      {
        std::ostringstream message;
        message << "GMRES failed to converge for epsilon=" << epsilon
                << ", method=" << (use_supg ? "SUPG" : "Galerkin")
                << ", iterations=" << solve_info.iterations
                << ", relative residual=" << solve_info.relative_residual;
        throw std::runtime_error(message.str());
      }

    if (write_solution)
      solution.writeOpenDXData(output_file);

    FunctionFunction<double> exact(&exact_solution, &exact_gradient);
    return {epsilon,
            use_supg,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            solve_info.iterations,
            solve_info.relative_residual,
            Functional::L2Error(solution, exact, 4),
            Functional::H1SemiError(solution, exact, 4)};
  }

  void print_results(const std::vector<CaseResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion P1 manufactured solution\n";
    std::cout << std::setw(12) << "epsilon" << std::setw(10) << "method"
              << std::setw(10) << "cells" << std::setw(10) << "dofs"
              << std::setw(8) << "iter" << std::setw(14) << "relres"
              << std::setw(14) << "L2" << std::setw(14) << "H1-semi"
              << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(12) << result.epsilon
                  << std::setw(10) << (result.supg ? "SUPG" : "Galerkin")
                  << std::setw(10) << result.cells
                  << std::setw(10) << result.dofs
                  << std::setw(8) << result.iterations
                  << std::setw(14) << result.relative_residual
                  << std::setw(14) << result.l2_error
                  << std::setw(14) << result.h1_error << '\n';
      }
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p05";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/convection_diffusion_afepack";

  try
    {
      std::vector<CaseResult> results;
      const std::vector<double> epsilons = {1.0e-1, 1.0e-2, 1.0e-3};
      for (const double epsilon : epsilons)
        {
          results.push_back(run_case(mesh_file,
                                     output_prefix + "_galerkin.dx",
                                     epsilon,
                                     false,
                                     false));
          results.push_back(run_case(mesh_file,
                                     output_prefix + "_supg.dx",
                                     epsilon,
                                     true,
                                     epsilon == epsilons.back()));
        }

      print_results(results);
      std::cout << "Wrote " << output_prefix << "_supg.dx\n";
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
