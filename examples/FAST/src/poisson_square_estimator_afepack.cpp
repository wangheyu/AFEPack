/**
 * @file poisson_square_estimator_afepack.cpp
 * @brief AFEPack Poisson 方程算例：方形区域误差估计。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 方形区域误差估计：在规则方形区域上给出后验误差估计的基准结果。
 * - 后验误差估计：计算误差指示子并输出可视化数据，为自适应网格加密提供依据。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define POISSON_LSHAPE_ESTIMATOR_NO_MAIN
#include "poisson_lshape_estimator_afepack.cpp"
#undef POISSON_LSHAPE_ESTIMATOR_NO_MAIN

namespace
{
  struct SquareEstimatorResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    double       h1_error = 0.0;
    double       estimator = 0.0;
    double       effectivity = 0.0;
    double       top_ten_percent_share = 0.0;
    double       boundary_band_share = 0.0;
    bool         has_rate = false;
    double       h1_rate = 0.0;
    double       estimator_rate = 0.0;
  };

  double boundary_band_share(const FEMSpace<double, dim> &space,
                             const std::vector<double> &local_eta_square,
                             const double total_square)
  {
    if (total_square <= 0.0)
      return 0.0;

    constexpr double band_width = 0.15;
    double selected = 0.0;
    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto centroid = element_centroid(space, space.element(e));
        const double boundary_distance =
          std::min({centroid[0],
                    centroid[1],
                    1.0 - centroid[0],
                    1.0 - centroid[1]});
        if (boundary_distance <= band_width)
          selected += local_eta_square[e];
      }
    return selected / total_square;
  }

  SquareEstimatorResult solve_square_and_estimate(
    const std::string &mesh_file,
    const double h)
  {
    active_problem = BookProblem::square_harmonic;

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

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    apply_dirichlet_boundary(stiff_matrix, solution, rhs, fem_space);

    AMGSolver solver(stiff_matrix);
    solver.solve(solution, rhs, 1.0e-10, 500);

    const ErrorNorms error = compute_error_norms(solution);
    const auto [estimator, local_eta_square] = jump_estimator(solution);
    const double total_square = estimator * estimator;

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            error.h1_seminorm,
            estimator,
            estimator / std::max(error.h1_seminorm, 1.0e-30),
            top_share(local_eta_square, total_square, 0.10),
            boundary_band_share(fem_space, local_eta_square, total_square)};
  }

  double rate(const double coarse_h,
              const double fine_h,
              const double coarse_value,
              const double fine_value)
  {
    if (coarse_h <= 0.0 || fine_h <= 0.0 ||
        coarse_value <= 0.0 || fine_value <= 0.0 ||
        coarse_h == fine_h)
      return 0.0;
    return std::log(coarse_value / fine_value) /
           std::log(coarse_h / fine_h);
  }

  void add_rates(std::vector<SquareEstimatorResult> &results)
  {
    for (std::size_t i = 1; i < results.size(); ++i)
      {
        results[i].has_rate = true;
        results[i].h1_rate = rate(results[i - 1].h,
                                  results[i].h,
                                  results[i - 1].h1_error,
                                  results[i].h1_error);
        results[i].estimator_rate = rate(results[i - 1].h,
                                         results[i].h,
                                         results[i - 1].estimator,
                                         results[i].estimator);
      }
  }
}

int main()
{
  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Poisson square residual jump estimator test\n";
      std::cout << "Book Example 1.1.3; P1 triangles, harmonic exact "
                   "solution, interior-edge normal-gradient jumps only.\n";
      std::cout << "h cells dofs H1_error eta effectivity "
                   "top10_eta_share boundary_band_eta_share H1_rate "
                   "eta_rate\n";
      std::cout << std::flush;

      const std::vector<std::pair<std::string, double>> meshes = {
        {"unit_square_h0p20", 0.20},
        {"unit_square_h0p10", 0.10},
        {"unit_square_h0p05", 0.05},
        {"unit_square_h0p025", 0.025}};

      std::vector<SquareEstimatorResult> results;
      for (const auto &[mesh_name, h] : meshes)
        results.push_back(solve_square_and_estimate(
          std::string(AFEPACK_DEFAULT_MESH_DIR) + "/" + mesh_name,
          h));
      add_rates(results);

      for (const auto &result : results)
        {
          std::cout << result.h << ' '
                    << result.cells << ' '
                    << result.dofs << ' '
                    << result.h1_error << ' '
                    << result.estimator << ' '
                    << result.effectivity << ' '
                    << result.top_ten_percent_share << ' '
                    << result.boundary_band_share << ' ';
          if (result.has_rate)
            std::cout << result.h1_rate << ' '
                      << result.estimator_rate << '\n';
          else
            std::cout << "- -\n";
        }

      if (!(results.back().h1_error < results.front().h1_error))
        throw std::runtime_error("H1 error should decrease on this smooth test");
      if (!(results.back().effectivity > 0.5 &&
            results.back().effectivity < 10.0))
        throw std::runtime_error("unexpected residual estimator effectivity");
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
