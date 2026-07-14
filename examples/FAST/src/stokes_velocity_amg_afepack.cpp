/**
 * @file stokes_velocity_amg_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：速度块 AMG。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 速度块 AMG：展示速度子问题的代数多重网格预处理在鞍点系统中的使用方式。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
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
#include <sstream>
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
  constexpr double viscosity = 1.0;

  using ScalarFunction = double (*)(const double *);
  using GradientFunction = std::vector<double> (*)(const double *);

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct ComponentResult
  {
    unsigned int nonzeros = 0;
    SolveInfo    cg;
    SolveInfo    jacobi_pcg;
    SolveInfo    amg_bicgstab;
    double       amg_l2_error = 0.0;
    double       amg_h1_error = 0.0;
  };

  struct VelocityBlockResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int scalar_dofs = 0;
    unsigned int velocity_dofs = 0;
    unsigned int scalar_nonzeros = 0;
    int          cg_max_iterations = 0;
    int          cg_total_iterations = 0;
    int          jacobi_max_iterations = 0;
    int          jacobi_total_iterations = 0;
    int          amg_max_iterations = 0;
    int          amg_total_iterations = 0;
    double       cg_relative_residual = 0.0;
    double       jacobi_relative_residual = 0.0;
    double       amg_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    bool         has_rates = false;
    double       velocity_l2_rate = 0.0;
    double       velocity_h1_rate = 0.0;
  };

  double sx(const double *p) { return std::sin(pi * p[0]); }
  double cx(const double *p) { return std::cos(pi * p[0]); }
  double sy(const double *p) { return std::sin(pi * p[1]); }
  double cy(const double *p) { return std::cos(pi * p[1]); }
  double c2x(const double *p) { return std::cos(2.0 * pi * p[0]); }
  double c2y(const double *p) { return std::cos(2.0 * pi * p[1]); }

  double velocity_x_exact(const double *p)
  {
    return 2.0 * pi * sx(p) * sx(p) * sy(p) * cy(p);
  }

  double velocity_y_exact(const double *p)
  {
    return -2.0 * pi * sx(p) * cx(p) * sy(p) * sy(p);
  }

  std::vector<double> velocity_x_gradient(const double *p)
  {
    return {4.0 * pi * pi * sx(p) * cx(p) * sy(p) * cy(p),
            2.0 * pi * pi * sx(p) * sx(p) * c2y(p)};
  }

  std::vector<double> velocity_y_gradient(const double *p)
  {
    return {-2.0 * pi * pi * c2x(p) * sy(p) * sy(p),
            -4.0 * pi * pi * sx(p) * cx(p) * sy(p) * cy(p)};
  }

  double laplace_velocity_x(const double *p)
  {
    const double b = sy(p) * cy(p);
    return 4.0 * pi * pi * pi * b *
           (c2x(p) - 2.0 * sx(p) * sx(p));
  }

  double laplace_velocity_y(const double *p)
  {
    const double a = sx(p) * cx(p);
    return 4.0 * pi * pi * pi * a *
           (2.0 * sy(p) * sy(p) - c2y(p));
  }

  double velocity_x_rhs(const double *p)
  {
    return -viscosity * laplace_velocity_x(p);
  }

  double velocity_y_rhs(const double *p)
  {
    return -viscosity * laplace_velocity_y(p);
  }

  double zero_boundary_value(const double *)
  {
    return 0.0;
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
    dof.readData("triangle.2.tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle.2.bas_fun");

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

  void apply_homogeneous_dirichlet_boundary(
    StiffMatrix<dim, double> &stiff_matrix,
    FEMFunction<double, dim> &solution,
    Vector<double> &rhs,
    const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &zero_boundary_value);

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

    double relative_residual = r.l2_norm() / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    Vector<double> z(rhs.size());
    apply_preconditioner(matrix, r, z, use_jacobi);
    Vector<double> p = z;
    Vector<double> ap(rhs.size());

    double rz_old = r.dot(z);
    if (std::abs(rz_old) < 1.0e-30)
      return {false, 0, relative_residual};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        matrix.vmult(ap, p);
        const double denominator = p.dot(ap);
        if (std::abs(denominator) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        const double alpha = rz_old / denominator;
        x.add(alpha, p);
        r.add(-alpha, ap);

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
        rz_old = rz_new;
      }

    return {false, max_iterations, relative_residual};
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

    double relative_residual = r.l2_norm() / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

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
          return {false, iteration - 1, relative_residual};

        const double beta = (rho_new / rho_old) * (alpha / omega);
        p.scale(beta);
        p.add(1.0, r);
        p.add(-beta * omega, v);

        p_hat = 0.0;
        preconditioner.vmult(p_hat, p);
        matrix.vmult(v, p_hat);

        const double rhat_v = r_hat.dot(v);
        if (std::abs(rhat_v) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        alpha = rho_new / rhat_v;
        s = r;
        s.add(-alpha, v);

        relative_residual = s.l2_norm() / rhs_norm;
        if (relative_residual <= tolerance)
          {
            x.add(alpha, p_hat);
            return {true, iteration, relative_residual};
          }

        s_hat = 0.0;
        preconditioner.vmult(s_hat, s);
        matrix.vmult(t, s_hat);

        const double tt = t.dot(t);
        if (std::abs(tt) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        omega = t.dot(s) / tt;
        if (std::abs(omega) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        x.add(alpha, p_hat);
        x.add(omega, s_hat);

        r = s;
        r.add(-omega, t);
        relative_residual = r.l2_norm() / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, iteration, relative_residual};

        rho_old = rho_new;
      }

    return {false, max_iterations, relative_residual};
  }

  ComponentResult solve_component(FEMSpace<double, dim> &space,
                                  ScalarFunction rhs_function,
                                  ScalarFunction exact_function,
                                  GradientFunction exact_gradient,
                                  const bool write_solution,
                                  const std::string &output_name)
  {
    StiffMatrix<dim, double> stiff_matrix(space);
    stiff_matrix.algebricAccuracy() = 5;
    stiff_matrix.build();

    Vector<double> rhs;
    Operator::L2Discretize(rhs_function, space, rhs, 5);

    FEMFunction<double, dim> boundary_solution(space);
    apply_homogeneous_dirichlet_boundary(stiff_matrix,
                                         boundary_solution,
                                         rhs,
                                         space);

    FEMFunction<double, dim> cg_solution(space);
    FEMFunction<double, dim> jacobi_solution(space);
    FEMFunction<double, dim> amg_solution(space);

    const SolveInfo cg_info =
      solve_cg(stiff_matrix, cg_solution, rhs, false, 1.0e-8, 40000);
    if (!cg_info.converged)
      throw std::runtime_error("velocity CG did not converge");

    const SolveInfo jacobi_info =
      solve_cg(stiff_matrix, jacobi_solution, rhs, true, 1.0e-8, 40000);
    if (!jacobi_info.converged)
      throw std::runtime_error("velocity Jacobi-PCG did not converge");

    AMGPreconditioner amg_preconditioner(stiff_matrix, 2, 2);
    amg_preconditioner.smoothStep() = 2;
    const SolveInfo amg_info =
      solve_amg_bicgstab(stiff_matrix,
                         amg_preconditioner,
                         amg_solution,
                         rhs,
                         1.0e-8,
                         4000);
    if (!amg_info.converged)
      throw std::runtime_error(
        "velocity AMG-preconditioned BiCGSTAB did not converge");

    if (write_solution)
      amg_solution.writeOpenDXData(std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                                   "/" + output_name);

    FunctionFunction<double> exact(exact_function, exact_gradient);
    return {static_cast<unsigned int>(stiff_matrix.n_nonzero_elements()),
            cg_info,
            jacobi_info,
            amg_info,
            Functional::L2Error(amg_solution, exact, 5),
            Functional::H1SemiError(amg_solution, exact, 5)};
  }

  double convergence_rate(const double coarse_h,
                          const double fine_h,
                          const double coarse_error,
                          const double fine_error)
  {
    if (coarse_h <= 0.0 || fine_h <= 0.0 ||
        coarse_error <= 0.0 || fine_error <= 0.0 ||
        coarse_h == fine_h)
      return 0.0;

    return std::log(coarse_error / fine_error) /
           std::log(coarse_h / fine_h);
  }

  std::string format_rate(const bool has_rate, const double rate)
  {
    if (!has_rate)
      return "-";

    std::ostringstream output;
    output << std::fixed << std::setprecision(3) << rate;
    return output.str();
  }

  VelocityBlockResult run_mesh_case(const double h,
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
    auto velocity_space = build_space(mesh, template_element);

    const ComponentResult ux =
      solve_component(velocity_space,
                      &velocity_x_rhs,
                      &velocity_x_exact,
                      &velocity_x_gradient,
                      write_solution,
                      "stokes_velocity_amg_ux_h0p05.dx");
    const ComponentResult uy =
      solve_component(velocity_space,
                      &velocity_y_rhs,
                      &velocity_y_exact,
                      &velocity_y_gradient,
                      write_solution,
                      "stokes_velocity_amg_uy_h0p05.dx");

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            2u * velocity_space.n_dof(),
            std::max(ux.nonzeros, uy.nonzeros),
            std::max(ux.cg.iterations, uy.cg.iterations),
            ux.cg.iterations + uy.cg.iterations,
            std::max(ux.jacobi_pcg.iterations,
                     uy.jacobi_pcg.iterations),
            ux.jacobi_pcg.iterations + uy.jacobi_pcg.iterations,
            std::max(ux.amg_bicgstab.iterations,
                     uy.amg_bicgstab.iterations),
            ux.amg_bicgstab.iterations + uy.amg_bicgstab.iterations,
            std::max(ux.cg.relative_residual,
                     uy.cg.relative_residual),
            std::max(ux.jacobi_pcg.relative_residual,
                     uy.jacobi_pcg.relative_residual),
            std::max(ux.amg_bicgstab.relative_residual,
                     uy.amg_bicgstab.relative_residual),
            std::sqrt(ux.amg_l2_error * ux.amg_l2_error +
                      uy.amg_l2_error * uy.amg_l2_error),
            std::sqrt(ux.amg_h1_error * ux.amg_h1_error +
                      uy.amg_h1_error * uy.amg_h1_error)};
  }

  void add_rates(std::vector<VelocityBlockResult> &results)
  {
    for (std::size_t i = 1; i < results.size(); ++i)
      {
        results[i].has_rates = true;
        results[i].velocity_l2_rate =
          convergence_rate(results[i - 1].h,
                           results[i].h,
                           results[i - 1].velocity_l2_error,
                           results[i].velocity_l2_error);
        results[i].velocity_h1_rate =
          convergence_rate(results[i - 1].h,
                           results[i].h,
                           results[i - 1].velocity_h1_error,
                           results[i].velocity_h1_error);
      }
  }

  void run_default_study(const std::string &mesh_dir)
  {
    const std::vector<std::pair<double, std::string>> mesh_cases = {
      {0.20, mesh_dir + "/unit_square_h0p20"},
      {0.10, mesh_dir + "/unit_square_h0p10"},
      {0.05, mesh_dir + "/unit_square_h0p05"}};

    std::vector<VelocityBlockResult> results;
    for (std::size_t i = 0; i < mesh_cases.size(); ++i)
      results.push_back(run_mesh_case(mesh_cases[i].first,
                                      mesh_cases[i].second,
                                      i + 1 == mesh_cases.size()));

    add_rates(results);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Stokes velocity block AMG diagnostic\n";
    std::cout << "P2 scalar velocity blocks, homogeneous Dirichlet, "
                 "manufactured divergence-free velocity.\n";
    std::cout << "h cells scalar_dofs velocity_dofs scalar_nnz "
                 "CGmax CGsum JPCGmax JPCGsum AMG-BiCGmax AMG-BiCGsum "
                 "CG_res JPCG_res AMG_res velocity_L2 velocity_H1 "
                 "velocity_L2_rate velocity_H1_rate\n";

    for (const auto &result : results)
      std::cout << result.h << ' '
                << result.cells << ' '
                << result.scalar_dofs << ' '
                << result.velocity_dofs << ' '
                << result.scalar_nonzeros << ' '
                << result.cg_max_iterations << ' '
                << result.cg_total_iterations << ' '
                << result.jacobi_max_iterations << ' '
                << result.jacobi_total_iterations << ' '
                << result.amg_max_iterations << ' '
                << result.amg_total_iterations << ' '
                << result.cg_relative_residual << ' '
                << result.jacobi_relative_residual << ' '
                << result.amg_relative_residual << ' '
                << result.velocity_l2_error << ' '
                << result.velocity_h1_error << ' '
                << std::setw(9)
                << format_rate(result.has_rates,
                               result.velocity_l2_rate) << ' '
                << std::setw(9)
                << format_rate(result.has_rates,
                               result.velocity_h1_rate) << '\n';

    if (!(results.back().velocity_l2_error < results.front().velocity_l2_error &&
          results.back().velocity_h1_error < results.front().velocity_h1_error))
      throw std::runtime_error("velocity block errors did not decrease");

    for (const auto &result : results)
      if (result.amg_max_iterations > result.jacobi_max_iterations)
        throw std::runtime_error(
          "AMG-preconditioned velocity solve should not be worse than Jacobi");

    std::cout << "Wrote "
              << std::string(AFEPACK_DEFAULT_OUTPUT_DIR)
              << "/stokes_velocity_amg_{ux,uy}_h0p05.dx\n";
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-directory]\n";
      return 2;
    }

  const std::string mesh_dir =
    argc > 1 ? argv[1] : std::string(AFEPACK_DEFAULT_MESH_DIR);

  try
    {
      run_default_study(mesh_dir);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
