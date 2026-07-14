/**
 * @file poisson_lshape_stretched_afepack.cpp
 * @brief AFEPack Poisson 方程算例：拉伸网格实验。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 拉伸网格实验：比较各向异性或拉伸网格下误差估计与求解器表现。
 * - L 形区域奇异性：考察非凸区域角点奇异性对误差分布和收敛阶的影响。
 *
 * 网格与数据：主要依赖 L 形区域网格、单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
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
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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
  constexpr int    half_intervals = 8;
  constexpr double pi = 3.141592653589793238462643383279502884;

  struct Node
  {
    double x = 0.0;
    double y = 0.0;
    int    boundary_mark = 0;
  };

  struct ElementData
  {
    std::array<int, 3> node = {};
    std::array<int, 3> neighbor = {{-1, -1, -1}};
    std::array<int, 3> side = {{-1, -1, -1}};
  };

  struct SideData
  {
    int a = -1;
    int b = -1;
    int ea = -1;
    int eb = -1;
    int boundary_mark = 0;
  };

  struct GeneratedMesh
  {
    std::vector<Node>        nodes;
    std::vector<ElementData> elements;
    std::vector<SideData>    sides;
    double                   min_edge = 0.0;
    double                   max_edge = 0.0;
  };

  struct StretchResult
  {
    double       alpha = 1.0;
    std::string  mesh_base;
    std::string  output_file;
    unsigned int cells = 0;
    unsigned int dofs = 0;
    double       min_edge = 0.0;
    double       max_edge = 0.0;
    double       l2_error = 0.0;
    double       h1_seminorm_error = 0.0;
  };

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

  void book_coordinates(const double *p, double &x, double &y)
  {
    x = 2.0 * p[0] - 1.0;
    y = 2.0 * p[1] - 1.0;
  }

  double exact_solution(const double *p)
  {
    double x = 0.0;
    double y = 0.0;
    book_coordinates(p, x, y);

    const double radius = std::hypot(x, y);
    if (radius < 1.0e-14)
      return 0.0;

    const double theta = std::atan2(y, x);
    return std::pow(radius, 2.0 / 3.0) *
           std::sin((2.0 * theta + pi) / 3.0);
  }

  std::vector<double> exact_gradient(const double *p)
  {
    double x = 0.0;
    double y = 0.0;
    book_coordinates(p, x, y);

    const double radius = std::hypot(x, y);
    if (radius < 1.0e-14)
      return {0.0, 0.0};

    const double alpha = 2.0 / 3.0;
    const double theta = std::atan2(y, x);
    const double psi = (2.0 * theta + pi) / 3.0;
    const double scale = alpha * std::pow(radius, alpha - 1.0);

    const double du_dx_book =
      scale * (std::sin(psi) * std::cos(theta) -
               std::cos(psi) * std::sin(theta));
    const double du_dy_book =
      scale * (std::sin(psi) * std::sin(theta) +
               std::cos(psi) * std::cos(theta));

    return {2.0 * du_dx_book, 2.0 * du_dy_book};
  }

  double right_hand_side(const double *p)
  {
    (void)p;
    return 0.0;
  }

  std::vector<double> stretched_coordinates(const double stretch_alpha)
  {
    const int n = half_intervals;
    std::vector<double> coordinates(2 * n + 1, 0.5);
    std::vector<double> lengths(n, 0.0);

    if (std::abs(stretch_alpha - 1.0) < 1.0e-14)
      {
        std::fill(lengths.begin(), lengths.end(), 0.5 / n);
      }
    else
      {
        const double first =
          0.5 * (stretch_alpha - 1.0) /
          (std::pow(stretch_alpha, n) - 1.0);
        for (int i = 0; i < n; ++i)
          lengths[i] = first * std::pow(stretch_alpha, i);
      }

    coordinates[n] = 0.5;
    for (int i = 0; i < n; ++i)
      {
        coordinates[n + i + 1] = coordinates[n + i] + lengths[i];
        coordinates[n - i - 1] = coordinates[n - i] - lengths[i];
      }

    coordinates.front() = 0.0;
    coordinates.back() = 1.0;
    return coordinates;
  }

  bool nearly_equal(const double a, const double b)
  {
    return std::abs(a - b) < 1.0e-12;
  }

  bool in_removed_quadrant(const double x, const double y)
  {
    return x < 0.5 - 1.0e-12 && y < 0.5 - 1.0e-12;
  }

  bool on_lshape_boundary(const double x, const double y)
  {
    return nearly_equal(x, 1.0) || nearly_equal(y, 1.0) ||
           (nearly_equal(x, 0.0) && y >= 0.5 - 1.0e-12) ||
           (nearly_equal(y, 0.0) && x >= 0.5 - 1.0e-12) ||
           (nearly_equal(x, 0.5) && y <= 0.5 + 1.0e-12) ||
           (nearly_equal(y, 0.5) && x <= 0.5 + 1.0e-12);
  }

  int add_node(std::vector<Node> &nodes,
               std::vector<std::vector<int>> &node_id,
               const int i,
               const int j,
               const std::vector<double> &x,
               const std::vector<double> &y)
  {
    if (node_id[j][i] >= 0)
      return node_id[j][i];

    const double px = x[i];
    const double py = y[j];
    const int id = static_cast<int>(nodes.size());
    nodes.push_back({px, py, on_lshape_boundary(px, py) ? 1 : 0});
    node_id[j][i] = id;
    return id;
  }

  double edge_length(const Node &a, const Node &b)
  {
    return std::hypot(a.x - b.x, a.y - b.y);
  }

  GeneratedMesh generate_lshape_mesh(const double stretch_alpha)
  {
    const auto x = stretched_coordinates(stretch_alpha);
    const auto y = stretched_coordinates(stretch_alpha);
    const int n = static_cast<int>(x.size()) - 1;

    GeneratedMesh mesh;
    std::vector<std::vector<int>> node_id(
      x.size(), std::vector<int>(x.size(), -1));

    for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= n; ++i)
        if (!in_removed_quadrant(x[i], y[j]))
          add_node(mesh.nodes, node_id, i, j, x, y);

    auto add_triangle = [&mesh](const int a, const int b, const int c) {
      ElementData element;
      element.node = {{a, b, c}};
      mesh.elements.push_back(element);
    };

    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        {
          const double x_mid = 0.5 * (x[i] + x[i + 1]);
          const double y_mid = 0.5 * (y[j] + y[j + 1]);
          if (in_removed_quadrant(x_mid, y_mid))
            continue;

          const int n00 = node_id[j][i];
          const int n10 = node_id[j][i + 1];
          const int n11 = node_id[j + 1][i + 1];
          const int n01 = node_id[j + 1][i];
          if (n00 < 0 || n10 < 0 || n11 < 0 || n01 < 0)
            throw std::runtime_error("invalid L-shape cell connectivity");

          add_triangle(n00, n10, n11);
          add_triangle(n00, n11, n01);
        }

    std::map<std::pair<int, int>, int> edge_to_side;
    auto register_edge = [&](const int element_index,
                             const int local_side,
                             const int a,
                             const int b) {
      const std::pair<int, int> key = std::minmax(a, b);
      auto [it, inserted] =
        edge_to_side.emplace(key, static_cast<int>(mesh.sides.size()));
      if (inserted)
        {
          SideData side;
          side.a = key.first;
          side.b = key.second;
          side.ea = element_index;
          mesh.sides.push_back(side);
        }
      else
        {
          SideData &side = mesh.sides[it->second];
          if (side.eb >= 0)
            throw std::runtime_error("nonmanifold generated side");
          side.eb = element_index;
        }
      mesh.elements[element_index].side[local_side] = it->second;
    };

    for (int e = 0; e < static_cast<int>(mesh.elements.size()); ++e)
      {
        const auto node = mesh.elements[e].node;
        register_edge(e, 0, node[1], node[2]);
        register_edge(e, 1, node[2], node[0]);
        register_edge(e, 2, node[0], node[1]);
      }

    for (int e = 0; e < static_cast<int>(mesh.elements.size()); ++e)
      for (int s = 0; s < 3; ++s)
        {
          const SideData &side = mesh.sides[mesh.elements[e].side[s]];
          mesh.elements[e].neighbor[s] =
            side.ea == e ? side.eb : side.ea;
        }

    mesh.min_edge = std::numeric_limits<double>::max();
    mesh.max_edge = 0.0;
    for (SideData &side : mesh.sides)
      {
        side.boundary_mark = side.eb < 0 ? 1 : 0;
        const double length =
          edge_length(mesh.nodes[side.a], mesh.nodes[side.b]);
        mesh.min_edge = std::min(mesh.min_edge, length);
        mesh.max_edge = std::max(mesh.max_edge, length);
      }

    return mesh;
  }

  std::string alpha_tag(const double alpha)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << alpha;
    std::string tag = stream.str();
    for (char &c : tag)
      if (c == '.')
        c = 'p';
    while (tag.size() > 1 && tag.back() == '0')
      tag.pop_back();
    if (!tag.empty() && tag.back() == 'p')
      tag.pop_back();
    return tag;
  }

  void write_generated_mesh(const GeneratedMesh &mesh,
                            const std::string &mesh_base)
  {
    std::filesystem::create_directories(
      std::filesystem::path(mesh_base).parent_path());

    {
      std::ofstream out(mesh_base + ".n");
      out << std::setw(6) << mesh.nodes.size() << std::setw(8)
          << mesh.elements.size() << std::setw(8) << mesh.sides.size()
          << " **(Nnd, Nee, Nsd)**\n";
      out << std::scientific << std::setprecision(16);
      for (std::size_t i = 0; i < mesh.nodes.size(); ++i)
        out << std::setw(6) << i << " " << std::setw(24) << mesh.nodes[i].x
            << " " << std::setw(24) << mesh.nodes[i].y << " "
            << std::setw(6) << mesh.nodes[i].boundary_mark << "\n";
    }

    {
      std::ofstream out(mesh_base + ".e");
      out << std::setw(6) << mesh.elements.size() << std::setw(8)
          << mesh.nodes.size() << std::setw(8) << mesh.sides.size()
          << " **(Nee, Nnd, Nsd)**\n";
      for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        {
          const auto &element = mesh.elements[i];
          out << std::setw(6) << i;
          for (const int node : element.node)
            out << std::setw(8) << node;
          for (const int neighbor : element.neighbor)
            out << std::setw(8) << neighbor;
          for (const int side : element.side)
            out << std::setw(8) << side;
          out << "\n";
        }
    }

    {
      std::ofstream out(mesh_base + ".s");
      out << mesh.sides.size() << "\n";
      for (std::size_t i = 0; i < mesh.sides.size(); ++i)
        {
          const auto &side = mesh.sides[i];
          out << std::setw(6) << i << std::setw(8) << side.a
              << std::setw(8) << side.b << std::setw(8) << side.ea
              << std::setw(8) << side.eb << std::setw(8)
              << side.boundary_mark << " \n";
        }
    }
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
    BoundaryFunction<double, dim> boundary(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary);
    boundary_admin.apply(stiff_matrix, solution, rhs);
  }

  std::pair<double, double>
  compute_error_norms(const FEMFunction<double, dim> &solution)
  {
    double l2_square = 0.0;
    double h1_square = 0.0;

    const auto &space = solution.femSpace();
    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info = element.findQuadratureInfo(5);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = solution.value(q_point, element);
        const auto gradients = solution.gradient(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double point[dim] = {q_point[q][0], q_point[q][1]};
            const double value_error = exact_solution(point) - values[q];
            const auto exact_grad = exact_gradient(point);
            const double grad_error_x = exact_grad[0] - gradients[q][0];
            const double grad_error_y = exact_grad[1] - gradients[q][1];
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;

            l2_square += jxw * value_error * value_error;
            h1_square +=
              jxw * (grad_error_x * grad_error_x +
                     grad_error_y * grad_error_y);
          }
      }

    return {2.0 * std::sqrt(l2_square), std::sqrt(h1_square)};
  }

  StretchResult solve_for_alpha(const double stretch_alpha)
  {
    const std::string tag = alpha_tag(stretch_alpha);
    const std::string mesh_base =
      std::string(AFEPACK_DEFAULT_MESH_DIR) +
      "/generated_lshape_stretched_alpha" + tag;
    const std::string output_file =
      std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
      "/poisson_lshape_stretched_alpha" + tag + ".dx";

    const GeneratedMesh generated_mesh =
      generate_lshape_mesh(stretch_alpha);
    write_generated_mesh(generated_mesh, mesh_base);

    EasyMesh mesh;
    mesh.readData(mesh_base);

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
    std::cout << '\n';
    solution.writeOpenDXData(output_file);

    const auto [l2_error, h1_error] = compute_error_norms(solution);

    return {stretch_alpha,
            mesh_base,
            output_file,
            mesh.n_geometry(dim),
            fem_space.n_dof(),
            generated_mesh.min_edge,
            generated_mesh.max_edge,
            l2_error,
            h1_error};
  }
}

int main()
{
  try
    {
      const std::vector<double> alphas = {1.0, 1.25, 1.5, 2.0};
      std::vector<StretchResult> results;
      results.reserve(alphas.size());

      std::cout << std::scientific << std::setprecision(6);
      std::cout
        << "AFEPack Poisson L-shaped singular stretched-grid study\n";
      std::cout
        << "Book Example 1.1.4, P1 triangles, fixed 16x16 tensor grid with lower-left quadrant removed.\n";

      for (const double alpha : alphas)
        results.push_back(solve_for_alpha(alpha));

      std::cout << std::setw(10) << "alpha" << std::setw(10) << "cells"
                << std::setw(10) << "dofs" << std::setw(16)
                << "min edge" << std::setw(16) << "max edge"
                << std::setw(16) << "L2" << std::setw(16)
                << "H1 semi" << '\n';

      for (const auto &result : results)
        {
          std::cout << std::setw(10) << result.alpha << std::setw(10)
                    << result.cells << std::setw(10) << result.dofs
                    << std::setw(16) << result.min_edge
                    << std::setw(16) << result.max_edge << std::setw(16)
                    << result.l2_error << std::setw(16)
                    << result.h1_seminorm_error << '\n';
        }

      const auto best =
        std::min_element(results.begin(),
                         results.end(),
                         [](const StretchResult &a,
                            const StretchResult &b) {
                           return a.h1_seminorm_error <
                                  b.h1_seminorm_error;
                         });
      if (best != results.end())
        std::cout << "Best H1-semi error alpha: " << best->alpha
                  << " error: " << best->h1_seminorm_error << '\n';

      std::cout << "Wrote "
                << std::string(AFEPACK_DEFAULT_OUTPUT_DIR)
                << "/poisson_lshape_stretched_alpha*.dx\n";
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
