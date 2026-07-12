/**
 * @file poisson_solver_comparison_afepack.cpp
 * @brief AFEPack Poisson 方程算例：求解器对比。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 求解器对比：在同一离散系统上比较不同线性求解器或预处理器的迭代行为。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/AMGSolver.h>
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
#include <stdexcept>
#include <string>
#include <utility>
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

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct SolverComparisonResult
  {
    double       h = 0.0;
    unsigned int cells = 0;
    unsigned int dofs = 0;
    unsigned int nonzeros = 0;
    SolveInfo    cg;
    SolveInfo    jacobi_pcg;
    SolveInfo    amg_bicgstab;
    double       amg_relative_residual = 0.0;
    double       cg_l2_error = 0.0;
    double       jacobi_l2_error = 0.0;
    double       amg_bicgstab_l2_error = 0.0;
    double       amg_l2_error = 0.0;
  };

  double exact_solution(const double *p)
  {
    return std::sin(pi * p[0]) * std::sin(pi * p[1]);
  }

  double right_hand_side(const double *p)
  {
    return 2.0 * pi * pi * exact_solution(p);
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

  void apply_dirichlet_boundary(StiffMatrix<dim, double> &stiff_matrix,
                                FEMFunction<double, dim> &solution,
                                Vector<double> &rhs,
                                const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(stiff_matrix, solution, rhs);
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
  }

  double relative_residual(const SparseMatrix<double> &matrix,
                           const Vector<double> &x,
                           const Vector<double> &rhs)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> residual = rhs;
    residual.add(-1.0, ax);
    return residual.l2_norm() / rhs_norm;
  }

  SolveInfo solve_cg(const SparseMatrix<double> &matrix,
                     Vector<double> &x,
                     const Vector<double> &rhs,
                     const bool use_jacobi,
                     const double tolerance,
                     const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> r = rhs;
    r.add(-1.0, ax);

    double rel_residual = r.l2_norm() / rhs_norm;
    if (rel_residual <= tolerance)
      return {true, 0, rel_residual};

    Vector<double> z(rhs.size());
    apply_preconditioner(matrix, r, z, use_jacobi);
    Vector<double> p = z;
    Vector<double> ap(rhs.size());

    double rz_old = r.dot(z);
    if (std::abs(rz_old) < 1.0e-30)
      return {false, 0, rel_residual};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        matrix.vmult(ap, p);
        const double denominator = p.dot(ap);
        if (std::abs(denominator) < 1.0e-30)
          return {false, iteration - 1, rel_residual};

        const double alpha = rz_old / denominator;
        x.add(alpha, p);
        r.add(-alpha, ap);

        rel_residual = r.l2_norm() / rhs_norm;
        if (rel_residual <= tolerance)
          return {true, iteration, rel_residual};

        apply_preconditioner(matrix, r, z, use_jacobi);
        const double rz_new = r.dot(z);
        if (std::abs(rz_old) < 1.0e-30)
          return {false, iteration, rel_residual};

        const double beta = rz_new / rz_old;
        p.scale(beta);
        p.add(1.0, z);
        rz_old = rz_new;
      }

    return {false, max_iterations, rel_residual};
  }

  SolveInfo solve_amg_bicgstab(const SparseMatrix<double> &matrix,
                               AMGPreconditioner &preconditioner,
                               Vector<double> &x,
                               const Vector<double> &rhs,
                               const double tolerance,
                               const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> r = rhs;
    r.add(-1.0, ax);
    Vector<double> r_hat = r;

    double rel_residual = r.l2_norm() / rhs_norm;
    if (rel_residual <= tolerance)
      return {true, 0, rel_residual};

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

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        const double rho_new = r_hat.dot(r);
        if (std::abs(rho_new) < 1.0e-30)
          return {false, iteration - 1, rel_residual};

        const double beta = (rho_new / rho_old) * (alpha / omega);
        p.scale(beta);
        p.add(1.0, r);
        p.add(-beta * omega, v);

        p_hat = 0.0;
        preconditioner.vmult(p_hat, p);
        matrix.vmult(v, p_hat);

        const double rhat_v = r_hat.dot(v);
        if (std::abs(rhat_v) < 1.0e-30)
          return {false, iteration - 1, rel_residual};

        alpha = rho_new / rhat_v;
        s = r;
        s.add(-alpha, v);

        rel_residual = s.l2_norm() / rhs_norm;
        if (rel_residual <= tolerance)
          {
            x.add(alpha, p_hat);
            return {true, iteration, rel_residual};
          }

        s_hat = 0.0;
        preconditioner.vmult(s_hat, s);
        matrix.vmult(t, s_hat);

        const double tt = t.dot(t);
        if (std::abs(tt) < 1.0e-30)
          return {false, iteration - 1, rel_residual};

        omega = t.dot(s) / tt;
        if (std::abs(omega) < 1.0e-30)
          return {false, iteration - 1, rel_residual};

        x.add(alpha, p_hat);
        x.add(omega, s_hat);

        r = s;
        r.add(-omega, t);
        rel_residual = r.l2_norm() / rhs_norm;
        if (rel_residual <= tolerance)
          return {true, iteration, rel_residual};

        rho_old = rho_new;
      }

    return {false, max_iterations, rel_residual};
  }

  std::vector<std::pair<double, std::string>> default_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  SolverComparisonResult run_mesh_case(const double h,
                                       const std::string &mesh_file,
                                       const bool write_solution)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof(geometry);
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto fem_space = build_space(mesh, template_element);

    StiffMatrix<dim, double> stiff_matrix(fem_space);
    stiff_matrix.algebricAccuracy() = 4;
    stiff_matrix.build();

    FEMFunction<double, dim> boundary_solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    apply_dirichlet_boundary(stiff_matrix, boundary_solution, rhs, fem_space);

    FEMFunction<double, dim> cg_solution(fem_space);
    FEMFunction<double, dim> jacobi_solution(fem_space);
    FEMFunction<double, dim> amg_bicgstab_solution(fem_space);
    FEMFunction<double, dim> amg_solution(fem_space);

    const SolveInfo cg_info =
      solve_cg(stiff_matrix, cg_solution, rhs, false, 1.0e-8, 20000);
    if (!cg_info.converged)
      throw std::runtime_error("plain CG did not converge");

    const SolveInfo jacobi_info =
      solve_cg(stiff_matrix, jacobi_solution, rhs, true, 1.0e-8, 20000);
    if (!jacobi_info.converged)
      throw std::runtime_error("Jacobi-PCG did not converge");

    AMGPreconditioner amg_preconditioner(stiff_matrix, 2, 2);
    const SolveInfo amg_bicgstab_info =
      solve_amg_bicgstab(stiff_matrix,
                         amg_preconditioner,
                         amg_bicgstab_solution,
                         rhs,
                         1.0e-8,
                         2000);
    if (!amg_bicgstab_info.converged)
      throw std::runtime_error("AMG-preconditioned BiCGSTAB did not converge");

    AMGSolver amg_solver(stiff_matrix);
    amg_solver.solve(amg_solution, rhs, 1.0e-10, 400);
    std::cout << '\n';

    if (write_solution)
      {
        cg_solution.writeOpenDXData(
          std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
          "/poisson_solver_comparison_cg_h0p025.dx");
        amg_solution.writeOpenDXData(
          std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
          "/poisson_solver_comparison_amg_h0p025.dx");
      }

    return {h,
            mesh.n_geometry(dim),
            fem_space.n_dof(),
            static_cast<unsigned int>(stiff_matrix.n_nonzero_elements()),
            cg_info,
            jacobi_info,
            amg_bicgstab_info,
            relative_residual(stiff_matrix, amg_solution, rhs),
            Functional::L2Error(cg_solution,
                                FunctionFunction<double>(&exact_solution),
                                4),
            Functional::L2Error(jacobi_solution,
                                FunctionFunction<double>(&exact_solution),
                                4),
            Functional::L2Error(amg_bicgstab_solution,
                                FunctionFunction<double>(&exact_solution),
                                4),
            Functional::L2Error(amg_solution,
                                FunctionFunction<double>(&exact_solution),
                                4)};
  }

  void run_default_study()
  {
    const auto meshes = default_meshes();
    std::vector<SolverComparisonResult> results;
    results.reserve(meshes.size());

    for (std::size_t i = 0; i < meshes.size(); ++i)
      results.push_back(run_mesh_case(meshes[i].first,
                                      meshes[i].second,
                                      i + 1 == meshes.size()));

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson iterative solver comparison\n";
    std::cout << "P1 triangles, manufactured solution, homogeneous Dirichlet boundary.\n";
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(10) << "nnz"
              << std::setw(10) << "CG" << std::setw(14) << "CG res"
              << std::setw(12) << "J-PCG" << std::setw(14)
              << "J-PCG res" << std::setw(12) << "AMG-BiCG"
              << std::setw(14) << "AMG-BiCG res"
              << std::setw(14) << "AMG res"
              << std::setw(14) << "AMG L2" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(10) << result.h << std::setw(10)
                  << result.cells << std::setw(10) << result.dofs
                  << std::setw(10) << result.nonzeros << std::setw(10)
                  << result.cg.iterations << std::setw(14)
                  << result.cg.relative_residual << std::setw(12)
                  << result.jacobi_pcg.iterations << std::setw(14)
                  << result.jacobi_pcg.relative_residual << std::setw(14)
                  << result.amg_bicgstab.iterations << std::setw(14)
                  << result.amg_bicgstab.relative_residual << std::setw(14)
                  << result.amg_relative_residual << std::setw(14)
                  << result.amg_l2_error << '\n';
      }

    std::cout << "Wrote "
              << std::string(AFEPACK_DEFAULT_OUTPUT_DIR)
              << "/poisson_solver_comparison_{cg,amg}_h0p025.dx\n";
  }
}

int main()
{
  try
    {
      run_default_study();
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
