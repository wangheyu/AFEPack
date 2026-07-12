/**
 * @file convection_diffusion_shishkin_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：Shishkin 网格实验。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - Shishkin 网格实验：使用分层网格思想处理奇异摄动型对流扩散问题中的边界层。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <array>
#include <filesystem>
#include <fstream>
#include <map>

#define main convection_diffusion_book_3_1_1_default_main
#include "convection_diffusion_book_3_1_1_afepack.cpp"
#undef main

namespace
{
  struct ShishkinNode
  {
    double x = 0.0;
    double y = 0.0;
    int    boundary_mark = 0;
  };

  struct ShishkinSide
  {
    int a = -1;
    int b = -1;
    int ea = -1;
    int eb = -1;
    int boundary_mark = 0;
  };

  struct ShishkinElement
  {
    std::array<int, 3> node = {{-1, -1, -1}};
    std::array<int, 3> neighbor = {{-1, -1, -1}};
    std::array<int, 3> side = {{-1, -1, -1}};
  };

  struct ShishkinMesh
  {
    std::vector<ShishkinNode>    nodes;
    std::vector<ShishkinSide>    sides;
    std::vector<ShishkinElement> elements;
    double                      min_edge = 0.0;
    double                      max_edge = 0.0;
    double                      transition_height = 0.0;
  };

  struct ShishkinResult
  {
    int          n = 0;
    std::string mesh_kind;
    double       transition_height = 0.0;
    double       min_edge = 0.0;
    double       max_edge = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    int          iterations = 0;
    double       relative_residual = 0.0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
    double       lower_strip_h1_error = 0.0;
    double       minimum = 0.0;
    double       maximum = 0.0;
    double       h1_rate = 0.0;
  };

  bool generated_boundary_node(const double x, const double y)
  {
    constexpr double tolerance = 1.0e-12;
    return std::abs(x) < tolerance || std::abs(x - 1.0) < tolerance ||
           std::abs(y) < tolerance || std::abs(y - 1.0) < tolerance;
  }

  std::pair<long long, long long>
  coordinate_key(const double x, const double y)
  {
    constexpr double scale = 1.0e12;
    return {static_cast<long long>(std::llround(scale * x)),
            static_cast<long long>(std::llround(scale * y))};
  }

  double edge_length(const ShishkinNode &a, const ShishkinNode &b)
  {
    return std::hypot(a.x - b.x, a.y - b.y);
  }

  int add_node(ShishkinMesh &mesh,
               std::map<std::pair<long long, long long>, int> &node_id,
               const double x,
               const double y)
  {
    const auto key = coordinate_key(x, y);
    const auto it = node_id.find(key);
    if (it != node_id.end())
      return it->second;

    const int id = static_cast<int>(mesh.nodes.size());
    mesh.nodes.push_back({x, y, generated_boundary_node(x, y) ? 1 : 0});
    node_id.emplace(key, id);
    return id;
  }

  std::vector<double> uniform_coordinates(const int n)
  {
    std::vector<double> coordinates(n + 1);
    for (int i = 0; i <= n; ++i)
      coordinates[i] = static_cast<double>(i) / static_cast<double>(n);
    return coordinates;
  }

  std::vector<double>
  shishkin_y_coordinates(const int n, const double epsilon)
  {
    if (n % 2 != 0)
      throw std::runtime_error("Shishkin grid requires an even interval count");

    const double transition_height =
      std::min(0.5, epsilon * std::log(static_cast<double>(n)));
    const double y_star = 1.0 - transition_height;

    std::vector<double> coordinates(n + 1);
    const int half = n / 2;
    for (int j = 0; j <= half; ++j)
      coordinates[j] = y_star * static_cast<double>(j) /
                       static_cast<double>(half);
    for (int j = half + 1; j <= n; ++j)
      coordinates[j] = y_star +
                       transition_height *
                         static_cast<double>(j - half) /
                         static_cast<double>(half);
    return coordinates;
  }

  void add_rectangular_patch(
    ShishkinMesh &mesh,
    std::map<std::pair<long long, long long>, int> &node_id,
    const std::vector<double> &x_coordinates,
    const std::vector<double> &y_coordinates)
  {
    for (std::size_t j = 0; j + 1 < y_coordinates.size(); ++j)
      for (std::size_t i = 0; i + 1 < x_coordinates.size(); ++i)
        {
          const int n00 = add_node(mesh,
                                   node_id,
                                   x_coordinates[i],
                                   y_coordinates[j]);
          const int n10 = add_node(mesh,
                                   node_id,
                                   x_coordinates[i + 1],
                                   y_coordinates[j]);
          const int n11 = add_node(mesh,
                                   node_id,
                                   x_coordinates[i + 1],
                                   y_coordinates[j + 1]);
          const int n01 = add_node(mesh,
                                   node_id,
                                   x_coordinates[i],
                                   y_coordinates[j + 1]);

          ShishkinElement lower;
          lower.node = {{n00, n10, n11}};
          mesh.elements.push_back(lower);

          ShishkinElement upper;
          upper.node = {{n00, n11, n01}};
          mesh.elements.push_back(upper);
        }
  }

  void finalize_sides(ShishkinMesh &mesh)
  {
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
          ShishkinSide side;
          side.a = key.first;
          side.b = key.second;
          side.ea = element_index;
          mesh.sides.push_back(side);
        }
      else
        {
          ShishkinSide &side = mesh.sides[it->second];
          if (side.eb >= 0)
            throw std::runtime_error("nonmanifold Shishkin side");
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
          const ShishkinSide &side =
            mesh.sides[mesh.elements[e].side[s]];
          mesh.elements[e].neighbor[s] =
            side.ea == e ? side.eb : side.ea;
        }

    mesh.min_edge = std::numeric_limits<double>::max();
    mesh.max_edge = 0.0;
    for (ShishkinSide &side : mesh.sides)
      {
        side.boundary_mark = side.eb < 0 ? 1 : 0;
        const double length =
          edge_length(mesh.nodes[side.a], mesh.nodes[side.b]);
        mesh.min_edge = std::min(mesh.min_edge, length);
        mesh.max_edge = std::max(mesh.max_edge, length);
      }
  }

  ShishkinMesh generate_tensor_mesh(const int n,
                                    const double epsilon,
                                    const bool use_shishkin)
  {
    const std::vector<double> x_coordinates = uniform_coordinates(n);
    const std::vector<double> y_coordinates =
      use_shishkin ? shishkin_y_coordinates(n, epsilon) :
                     uniform_coordinates(n);

    ShishkinMesh mesh;
    mesh.transition_height =
      use_shishkin ? 1.0 - y_coordinates[y_coordinates.size() / 2] : 0.0;

    std::map<std::pair<long long, long long>, int> node_id;
    add_rectangular_patch(mesh, node_id, x_coordinates, y_coordinates);
    finalize_sides(mesh);
    return mesh;
  }

  void write_generated_mesh(const ShishkinMesh &mesh,
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

  std::string mesh_label(const int n, const bool use_shishkin)
  {
    std::ostringstream stream;
    stream << (use_shishkin ? "shishkin" : "uniform") << "_n" << n;
    return stream.str();
  }

  ShishkinResult run_shishkin_case(const std::string &mesh_dir,
                                   const int n,
                                   const double epsilon,
                                   const bool use_shishkin)
  {
    const ShishkinMesh mesh =
      generate_tensor_mesh(n, epsilon, use_shishkin);
    const std::string label = mesh_label(n, use_shishkin);
    const std::string mesh_base = mesh_dir + "/" + label;
    write_generated_mesh(mesh, mesh_base);

    const CaseResult result = run_case(mesh_base,
                                       "shishkin",
                                       label,
                                       1.0 / static_cast<double>(n),
                                       "",
                                       epsilon,
                                       false,
                                       false);

    return {n,
            use_shishkin ? "Shishkin" : "uniform",
            mesh.transition_height,
            mesh.min_edge,
            mesh.max_edge,
            result.cells,
            result.dofs,
            result.iterations,
            result.relative_residual,
            result.l2_error,
            result.h1_error,
            result.lower_strip_h1_error,
            result.minimum,
            result.maximum,
            0.0};
  }

  void compute_rates(std::vector<ShishkinResult> &results)
  {
    std::map<std::string, std::pair<int, double>> previous;
    for (ShishkinResult &result : results)
      {
        const std::string key = result.mesh_kind;
        const auto it = previous.find(key);
        if (it != previous.end())
          {
            const int previous_n = it->second.first;
            const double previous_error = it->second.second;
            result.h1_rate =
              std::log(previous_error / result.h1_error) /
              std::log(static_cast<double>(result.n) /
                       static_cast<double>(previous_n));
          }
        previous[key] = {result.n, result.h1_error};
      }
  }

  void print_shishkin_results(const std::vector<ShishkinResult> &results,
                              const double epsilon)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion Example 3.1.1 Shishkin grids\n";
    std::cout << "epsilon=" << epsilon
              << ", book transition y*=1-2*epsilon*log(n); "
              << "unit-square top layer height=epsilon*log(n).\n";
    std::cout << std::setw(6) << "n" << ' ' << std::setw(10) << "mesh"
              << ' ' << std::setw(12) << "tau_y" << ' '
              << std::setw(12) << "min_h" << ' ' << std::setw(12)
              << "max_h" << ' ' << std::setw(10) << "cells" << ' '
              << std::setw(10)
              << "dofs" << ' ' << std::setw(8) << "GMRES" << ' '
              << std::setw(14) << "relres" << ' ' << std::setw(14)
              << "L2" << ' ' << std::setw(14) << "H1-semi" << ' '
              << std::setw(14) << "H1-low" << ' ' << std::setw(10)
              << "H1-rate" << ' ' << std::setw(14) << "min" << ' '
              << std::setw(14) << "max" << '\n';

    for (const ShishkinResult &result : results)
      {
        std::cout << std::setw(6) << result.n << ' ' << std::setw(10)
                  << result.mesh_kind << ' ' << std::setw(12)
                  << result.transition_height << ' ' << std::setw(12)
                  << result.min_edge << ' ' << std::setw(12)
                  << result.max_edge << ' '
                  << std::setw(10) << result.cells << ' ' << std::setw(10)
                  << result.dofs << ' ' << std::setw(8)
                  << result.iterations << ' ' << std::setw(14)
                  << result.relative_residual << ' ' << std::setw(14)
                  << result.l2_error << ' ' << std::setw(14)
                  << result.h1_error << ' ' << std::setw(14)
                  << result.lower_strip_h1_error << ' ';
        if (result.h1_rate == 0.0)
          std::cout << std::setw(10) << "-";
        else
          std::cout << std::setw(10) << result.h1_rate;
        std::cout << ' ' << std::setw(14) << result.minimum << ' '
                  << std::setw(14) << result.maximum << '\n';
      }
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [generated-mesh-dir]\n";
      return 2;
    }

  const std::string mesh_dir =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/generated_shishkin_meshes";

  try
    {
      const double epsilon = 1.0 / 64.0;
      const std::vector<int> levels = {8, 16, 32, 64};

      std::vector<ShishkinResult> results;
      for (const int n : levels)
        for (const bool use_shishkin : {false, true})
          results.push_back(
            run_shishkin_case(mesh_dir, n, epsilon, use_shishkin));

      compute_rates(results);
      print_shishkin_results(results, epsilon);

      const auto too_large = std::find_if(
        results.begin(),
        results.end(),
        [](const ShishkinResult &result) {
          return result.mesh_kind == "Shishkin" && result.n == 64 &&
                 result.h1_error > 5.0e-1;
        });
      if (too_large != results.end())
        throw std::runtime_error("Shishkin finest-grid H1 error is too large");
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
