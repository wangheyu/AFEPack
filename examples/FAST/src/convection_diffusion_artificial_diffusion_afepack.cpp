/**
 * @file convection_diffusion_artificial_diffusion_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：人工扩散稳定化。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 人工扩散稳定化：比较在对流占优区域加入人工扩散后对振荡抑制和解剖面扩散的影响。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define main convection_diffusion_interior_layer_hidden_main
#include "convection_diffusion_interior_layer_afepack.cpp"
#undef main

namespace
{
  enum class Stabilization
  {
    galerkin,
    artificial_diffusion,
    supg
  };

  struct StabilizationResult
  {
    Stabilization method = Stabilization::galerkin;
    double        physical_epsilon = 0.0;
    double        diffusion_used = 0.0;
    double        physical_peclet = 0.0;
    int           cells = 0;
    unsigned int  dofs = 0;
    int           iterations = 0;
    double        relative_residual = 0.0;
    double        min_value = 0.0;
    double        max_value = 0.0;
    double        undershoot = 0.0;
    double        overshoot = 0.0;
    double        l2_norm = 0.0;
    double        h1_seminorm = 0.0;
    double        transition_fraction = 0.0;
  };

  const char *method_name(const Stabilization method)
  {
    switch (method)
      {
      case Stabilization::galerkin:
        return "Galerkin";
      case Stabilization::artificial_diffusion:
        return "artificial";
      case Stabilization::supg:
        return "SUPG";
      }
    return "unknown";
  }

  double transition_fraction(const FEMFunction<double, dim> &solution)
  {
    int interior_count = 0;
    int transition_count = 0;
    const auto &space = solution.femSpace();
    for (std::size_t i = 0; i < solution.size(); ++i)
      {
        if (space.dofBoundaryMark(i) != 0)
          continue;
        ++interior_count;
        if (solution(i) > 0.1 && solution(i) < 0.9)
          ++transition_count;
      }
    if (interior_count == 0)
      return 0.0;
    return static_cast<double>(transition_count) /
           static_cast<double>(interior_count);
  }

  StabilizationResult run_stabilization_case(const std::string &mesh_file,
                                             const double physical_epsilon,
                                             const double mesh_h,
                                             const Stabilization method)
  {
    current_epsilon = physical_epsilon;

    const bool use_supg = method == Stabilization::supg;
    const double wind_norm = std::sqrt(1.0 + sqrt3 * sqrt3);
    const double diffusion_used =
      method == Stabilization::artificial_diffusion ?
        std::max(physical_epsilon, 0.5 * wind_norm * mesh_h) :
        physical_epsilon;

    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto fem_space = build_space(mesh, template_element);

    ConvectionDiffusionMatrix matrix(fem_space, diffusion_used, use_supg);
    matrix.algebricAccuracy() = 4;
    matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    add_supg_rhs(fem_space, rhs, physical_epsilon, use_supg);

    BoundaryFunction<double, dim> bottom_boundary(
      BoundaryConditionInfo::DIRICHLET, 1, &boundary_bottom);
    BoundaryFunction<double, dim> right_boundary(
      BoundaryConditionInfo::DIRICHLET, 2, &boundary_one);
    BoundaryFunction<double, dim> top_boundary(
      BoundaryConditionInfo::DIRICHLET, 3, &boundary_zero);
    BoundaryFunction<double, dim> left_boundary(
      BoundaryConditionInfo::DIRICHLET, 4, &boundary_zero);

    BoundaryConditionAdmin<double, dim> boundary_admin(fem_space);
    boundary_admin.add(bottom_boundary);
    boundary_admin.add(right_boundary);
    boundary_admin.add(top_boundary);
    boundary_admin.add(left_boundary);
    boundary_admin.apply(matrix, solution, rhs);

    for (std::size_t i = 0; i < solution.size(); ++i)
      solution(i) = 0.0;
    const SolveInfo solve_info =
      solve_gmres(matrix, solution, rhs, 1.0e-9, 50, 2500);
    if (!solve_info.converged)
      {
        std::ostringstream message;
        message << "GMRES failed for " << method_name(method)
                << ", iterations=" << solve_info.iterations
                << ", relative residual=" << solve_info.relative_residual;
        throw std::runtime_error(message.str());
      }

    const SolutionDiagnostics diagnostics = diagnose_solution(solution);

    return {method,
            physical_epsilon,
            diffusion_used,
            wind_norm * mesh_h / (2.0 * physical_epsilon),
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            solve_info.iterations,
            solve_info.relative_residual,
            diagnostics.min_value,
            diagnostics.max_value,
            diagnostics.undershoot,
            diagnostics.overshoot,
            diagnostics.l2_norm,
            diagnostics.h1_seminorm,
            transition_fraction(solution)};
  }

  void print_stabilization_results(
    const std::vector<StabilizationResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion artificial diffusion test\n";
    std::cout << "Book Example 3.1.3 interior/downstream layer, "
                 "epsilon=1/200, h=0.05.\n";
    std::cout << "The artificial row solves a modified isotropic-diffusion "
                 "operator with epsilon_eff=max(epsilon,|beta|h/2).\n";
    std::cout << std::setw(14) << "method" << ' ' << std::setw(12)
              << "eps_used" << ' ' << std::setw(10) << "Pe_h" << ' '
              << std::setw(10) << "cells" << ' ' << std::setw(10)
              << "dofs" << ' ' << std::setw(8) << "iter" << ' '
              << std::setw(14) << "relres" << ' ' << std::setw(14)
              << "min" << ' ' << std::setw(14) << "max" << ' '
              << std::setw(14) << "under" << ' ' << std::setw(14)
              << "over" << ' ' << std::setw(14) << "L2-norm" << ' '
              << std::setw(14) << "H1-semi" << ' ' << std::setw(14)
              << "transition" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(14) << method_name(result.method) << ' '
                  << std::setw(12) << result.diffusion_used << ' '
                  << std::setw(10) << result.physical_peclet << ' '
                  << std::setw(10) << result.cells << ' '
                  << std::setw(10) << result.dofs << ' '
                  << std::setw(8) << result.iterations << ' '
                  << std::setw(14) << result.relative_residual << ' '
                  << std::setw(14) << result.min_value << ' '
                  << std::setw(14) << result.max_value << ' '
                  << std::setw(14) << result.undershoot << ' '
                  << std::setw(14) << result.overshoot << ' '
                  << std::setw(14) << result.l2_norm << ' '
                  << std::setw(14) << result.h1_seminorm << ' '
                  << std::setw(14) << result.transition_fraction << '\n';
      }

    std::cout << "transition is the fraction of interior dofs with "
                 "0.1<u_h<0.9.\n";
  }
}

int main()
{
  try
    {
      const double physical_epsilon = 1.0 / 200.0;
      const double mesh_h = 0.05;
      const std::string mesh_file =
        std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p05";

      std::vector<StabilizationResult> results;
      results.push_back(run_stabilization_case(mesh_file,
                                               physical_epsilon,
                                               mesh_h,
                                               Stabilization::galerkin));
      results.push_back(run_stabilization_case(
        mesh_file,
        physical_epsilon,
        mesh_h,
        Stabilization::artificial_diffusion));
      results.push_back(run_stabilization_case(mesh_file,
                                               physical_epsilon,
                                               mesh_h,
                                               Stabilization::supg));

      print_stabilization_results(results);

      if (!(results[2].overshoot < results[0].overshoot))
        throw std::runtime_error("SUPG should reduce overshoot here");
      if (!(results[1].transition_fraction > results[2].transition_fraction))
        throw std::runtime_error(
          "artificial diffusion should produce a wider transition region here");
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
