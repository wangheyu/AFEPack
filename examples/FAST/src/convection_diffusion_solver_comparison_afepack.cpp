/**
 * @file convection_diffusion_solver_comparison_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：求解器对比。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
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

#define main convection_diffusion_boundary_layer_default_main
#include "convection_diffusion_boundary_layer_afepack.cpp"
#undef main

namespace
{
  enum class SolverKind
  {
    gmres,
    jacobi_gmres,
    jacobi_bicgstab
  };

  struct SolverStudyResult
  {
    std::string  discretization;
    std::string  solver;
    double       epsilon = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    int          iterations = 0;
    double       relative_residual = 0.0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
  };

  const char *solver_name(const SolverKind kind)
  {
    switch (kind)
      {
        case SolverKind::gmres:
          return "GMRES";
        case SolverKind::jacobi_gmres:
          return "J-GMRES";
        case SolverKind::jacobi_bicgstab:
          return "J-BiCG";
      }
    return "unknown";
  }

  SolveInfo solve_jacobi_gmres(const SparseMatrix<double> &matrix,
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
        std::vector<Vector<double>> z(restart, Vector<double>(rhs.size()));
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
            jacobi_apply(matrix, v[inner_iterations], z[inner_iterations]);
            matrix.vmult(v[inner_iterations + 1], z[inner_iterations]);

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
                  add_scaled(x, y[i], z[i]);
                return {true, total_iterations + 1, relative_residual};
              }
          }

        const auto y =
          solve_upper_hessenberg(h, g, inner_iterations);
        for (int i = 0; i < inner_iterations; ++i)
          add_scaled(x, y[i], z[i]);

        r = residual(matrix, x, rhs);
        beta = r.l2_norm();
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  SolverStudyResult run_solver_case(const std::string &mesh_file,
                                    const double epsilon,
                                    const bool use_supg,
                                    const SolverKind solver_kind)
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
    SolveInfo solve_info;
    switch (solver_kind)
      {
        case SolverKind::gmres:
          solve_info = solve_gmres(matrix, solution, rhs, 1.0e-10, 50, 2000);
          break;
        case SolverKind::jacobi_gmres:
          solve_info =
            solve_jacobi_gmres(matrix, solution, rhs, 1.0e-10, 50, 2000);
          break;
        case SolverKind::jacobi_bicgstab:
          solve_info =
            solve_bicgstab(matrix, solution, rhs, 1.0e-10, 2000);
          break;
      }

    if (!solve_info.converged)
      {
        std::ostringstream message;
        message << solver_name(solver_kind)
                << " failed for method=" << (use_supg ? "SUPG" : "Galerkin")
                << ", iterations=" << solve_info.iterations
                << ", relative residual=" << solve_info.relative_residual;
        throw std::runtime_error(message.str());
      }

    FunctionFunction<double> exact(&exact_solution, &exact_gradient);
    return {use_supg ? "SUPG" : "Galerkin",
            solver_name(solver_kind),
            epsilon,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            solve_info.iterations,
            solve_info.relative_residual,
            Functional::L2Error(solution, exact, 4),
            Functional::H1SemiError(solution, exact, 4)};
  }

  void print_solver_results(const std::vector<SolverStudyResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion solver comparison\n";
    std::cout << "epsilon=1e-3, h=0.05 exponential outflow layer, P1 triangles.\n";
    std::cout << std::setw(12) << "method" << ' ' << std::setw(10)
              << "solver" << ' ' << std::setw(10) << "cells" << ' '
              << std::setw(10) << "dofs" << ' ' << std::setw(8)
              << "iter" << ' ' << std::setw(14) << "relres" << ' '
              << std::setw(14) << "L2" << ' ' << std::setw(14)
              << "H1-semi" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(12) << result.discretization << ' '
                  << std::setw(10) << result.solver << ' ' << std::setw(10)
                  << result.cells << ' ' << std::setw(10) << result.dofs
                  << ' ' << std::setw(8) << result.iterations << ' '
                  << std::setw(14) << result.relative_residual << ' '
                  << std::setw(14) << result.l2_error << ' ' << std::setw(14)
                  << result.h1_error << '\n';
      }
  }
}

int main()
{
  try
    {
      const std::string mesh_file =
        std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p05";
      const double epsilon = 1.0e-3;
      const std::vector<SolverKind> solvers = {SolverKind::gmres,
                                               SolverKind::jacobi_gmres};

      std::vector<SolverStudyResult> results;
      for (const bool use_supg : {false, true})
        for (const SolverKind solver : solvers)
          results.push_back(run_solver_case(mesh_file,
                                            epsilon,
                                            use_supg,
                                            solver));

      print_solver_results(results);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
