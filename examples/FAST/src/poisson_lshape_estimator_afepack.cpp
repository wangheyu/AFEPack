/**
 * @file poisson_lshape_estimator_afepack.cpp
 * @brief AFEPack Poisson 方程算例：后验误差估计。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 后验误差估计：计算误差指示子并输出可视化数据，为自适应网格加密提供依据。
 * - L 形区域奇异性：考察非凸区域角点奇异性对误差分布和收敛阶的影响。
 *
 * 网格与数据：主要依赖 L 形区域网格、单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define main poisson_book_examples_hidden_main
#include "poisson_book_examples_afepack.cpp"
#undef main

namespace
{
  struct ElementEdge
  {
    int    element = -1;
    double grad_x = 0.0;
    double grad_y = 0.0;
  };

  struct EstimatorResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    double       h1_error = 0.0;
    double       estimator = 0.0;
    double       effectivity = 0.0;
    double       top_ten_percent_share = 0.0;
    double       near_corner_share = 0.0;
  };

  double distance_between_dofs(const FEMSpace<double, dim> &space,
                               const int left,
                               const int right)
  {
    const auto &p = space.dofInfo(left).interp_point;
    const auto &q = space.dofInfo(right).interp_point;
    return std::hypot(p[0] - q[0], p[1] - q[1]);
  }

  Point<dim> element_centroid(const FEMSpace<double, dim> &space,
                              const Element<double, dim> &element)
  {
    Point<dim> centroid;
    centroid[0] = 0.0;
    centroid[1] = 0.0;

    const auto &dofs = element.dof();
    for (const int dof : dofs)
      {
        const auto &point = space.dofInfo(dof).interp_point;
        centroid[0] += point[0];
        centroid[1] += point[1];
      }
    centroid[0] /= static_cast<double>(dofs.size());
    centroid[1] /= static_cast<double>(dofs.size());
    return centroid;
  }

  std::array<double, 2> element_gradient(
    const FEMFunction<double, dim> &solution,
    const Element<double, dim> &element)
  {
    const auto centroid = element_centroid(solution.femSpace(), element);
    const std::vector<Point<dim>> points = {centroid};
    const auto gradients = solution.gradient(points, element);
    return {gradients[0][0], gradients[0][1]};
  }

  std::pair<double, std::vector<double>>
  jump_estimator(const FEMFunction<double, dim> &solution)
  {
    const auto &space = solution.femSpace();
    std::map<std::pair<int, int>, std::vector<ElementEdge>> edge_map;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto &element = space.element(e);
        const auto &dofs = element.dof();
        if (dofs.size() != 3)
          throw std::runtime_error("P1 triangle estimator expects 3 dofs");

        const auto gradient = element_gradient(solution, element);
        const std::array<std::pair<int, int>, 3> edges = {
          std::make_pair(dofs[0], dofs[1]),
          std::make_pair(dofs[1], dofs[2]),
          std::make_pair(dofs[2], dofs[0])};

        for (auto edge : edges)
          {
            if (edge.first > edge.second)
              std::swap(edge.first, edge.second);
            edge_map[edge].push_back({e, gradient[0], gradient[1]});
          }
      }

    double estimator_square = 0.0;
    std::vector<double> local_eta_square(space.n_element(), 0.0);

    for (const auto &[edge, adjacent] : edge_map)
      {
        if (adjacent.size() != 2)
          continue;

        const auto &p = space.dofInfo(edge.first).interp_point;
        const auto &q = space.dofInfo(edge.second).interp_point;
        const double dx = q[0] - p[0];
        const double dy = q[1] - p[1];
        const double length = distance_between_dofs(space,
                                                    edge.first,
                                                    edge.second);
        if (length <= 1.0e-14)
          continue;

        const double normal_x = -dy / length;
        const double normal_y = dx / length;
        const double jump =
          (adjacent[0].grad_x - adjacent[1].grad_x) * normal_x +
          (adjacent[0].grad_y - adjacent[1].grad_y) * normal_y;
        const double contribution = length * length * jump * jump;

        estimator_square += contribution;
        local_eta_square[adjacent[0].element] += 0.5 * contribution;
        local_eta_square[adjacent[1].element] += 0.5 * contribution;
      }

    return {std::sqrt(estimator_square), local_eta_square};
  }

  double top_share(std::vector<double> local_eta_square,
                   const double total_square,
                   const double fraction)
  {
    if (total_square <= 0.0)
      return 0.0;

    std::sort(local_eta_square.begin(),
              local_eta_square.end(),
              std::greater<double>());
    const std::size_t count =
      std::max<std::size_t>(1,
                            static_cast<std::size_t>(
                              std::ceil(fraction *
                                        local_eta_square.size())));
    const double selected =
      std::accumulate(local_eta_square.begin(),
                      local_eta_square.begin() +
                        std::min(count, local_eta_square.size()),
                      0.0);
    return selected / total_square;
  }

  double near_corner_share(const FEMSpace<double, dim> &space,
                           const std::vector<double> &local_eta_square,
                           const double total_square)
  {
    if (total_square <= 0.0)
      return 0.0;

    constexpr double radius = 0.25;
    double selected = 0.0;
    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto centroid = element_centroid(space, space.element(e));
        const double dx = centroid[0] - 0.5;
        const double dy = centroid[1] - 0.5;
        if (std::hypot(dx, dy) <= radius)
          selected += local_eta_square[e];
      }
    return selected / total_square;
  }

  EstimatorResult solve_and_estimate(const std::string &mesh_file,
                                     const double h)
  {
    active_problem = BookProblem::lshape_singular;

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
            near_corner_share(fem_space, local_eta_square, total_square)};
  }
}

#ifndef POISSON_LSHAPE_ESTIMATOR_NO_MAIN
int main()
{
  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Poisson L-shape residual jump estimator test\n";
      std::cout << "Book Example 1.1.4; P1 triangles, f=0, interior-edge "
                   "normal-gradient jumps only.\n";
      std::cout << "h cells dofs H1_error eta effectivity "
                   "top10_eta_share near_corner_eta_share\n";

      const std::vector<std::pair<std::string, double>> meshes = {
        {"lshape_h0p20", 0.20},
        {"lshape_h0p10", 0.10},
        {"lshape_h0p05", 0.05}};

      for (const auto &[mesh_name, h] : meshes)
        {
          const EstimatorResult result =
            solve_and_estimate(std::string(AFEPACK_DEFAULT_MESH_DIR) +
                                 "/" + mesh_name,
                               h);
          std::cout << result.h << ' '
                    << result.cells << ' '
                    << result.dofs << ' '
                    << result.h1_error << ' '
                    << result.estimator << ' '
                    << result.effectivity << ' '
                    << result.top_ten_percent_share << ' '
                    << result.near_corner_share << '\n';
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
#endif
