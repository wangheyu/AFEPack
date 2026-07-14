/**
 * @file stokes_contracting_channel_geometry_compare_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：收缩通道流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 收缩通道流：在收缩几何中考察入口、出口和壁面边界对流动结构的影响。
 * - 几何设置对比：在多个几何变体之间比较边界处理、网格质量和物理量变化。
 *
 * 网格与数据：主要依赖 收缩通道网格、收缩通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#define main stokes_contracting_channel_profiles_hidden_main
#include "stokes_contracting_channel_profiles_afepack.cpp"
#undef main

namespace
{
  struct CurvedNode
  {
    double x = 0.0;
    double y = 0.0;
    int    boundary_mark = 0;
  };

  struct CurvedElement
  {
    std::array<int, 3> node = {};
    std::array<int, 3> neighbor = {{-1, -1, -1}};
    std::array<int, 3> side = {{-1, -1, -1}};
  };

  struct CurvedSide
  {
    int a = -1;
    int b = -1;
    int ea = -1;
    int eb = -1;
    int boundary_mark = 0;
  };

  struct CurvedMesh
  {
    std::vector<CurvedNode>    nodes;
    std::vector<CurvedElement> elements;
    std::vector<CurvedSide>    sides;
    double                     min_edge = 0.0;
    double                     max_edge = 0.0;
  };

  struct GeometryCase
  {
    std::string                name;
    std::string                mesh_file;
    std::string                output_prefix;
    ContractingChannelGeometry geometry;
    double                     area = 0.0;
    double                     outlet_half_width = 0.0;
  };

  double quadratic_bow_tie_half_width(const double x)
  {
    const double t = 1.0 - x / outlet_x;
    return 0.5 + 0.5 * t * t;
  }

  double curved_edge_length(const CurvedNode &a, const CurvedNode &b)
  {
    return std::hypot(a.x - b.x, a.y - b.y);
  }

  void build_sides(CurvedMesh &mesh)
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
          CurvedSide side;
          side.a = key.first;
          side.b = key.second;
          side.ea = element_index;
          mesh.sides.push_back(side);
        }
      else
        {
          CurvedSide &side = mesh.sides[it->second];
          if (side.eb >= 0)
            throw std::runtime_error("nonmanifold contracting-channel side");
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
          const CurvedSide &side = mesh.sides[mesh.elements[e].side[s]];
          mesh.elements[e].neighbor[s] =
            side.ea == e ? side.eb : side.ea;
        }

    mesh.min_edge = std::numeric_limits<double>::max();
    mesh.max_edge = 0.0;
    for (CurvedSide &side : mesh.sides)
      {
        side.boundary_mark = side.eb < 0 ? 1 : 0;
        const double length =
          curved_edge_length(mesh.nodes[side.a], mesh.nodes[side.b]);
        mesh.min_edge = std::min(mesh.min_edge, length);
        mesh.max_edge = std::max(mesh.max_edge, length);
      }
  }

  CurvedMesh generate_quadratic_bow_tie_proxy_mesh()
  {
    constexpr int nx = 16;
    constexpr int ny = 8;

    CurvedMesh mesh;
    mesh.nodes.reserve((nx + 1) * (ny + 1));

    auto node_id = [](const int i, const int j) {
      return j * (nx + 1) + i;
    };

    for (int j = 0; j <= ny; ++j)
      {
        const double eta = -1.0 + 2.0 * j / ny;
        for (int i = 0; i <= nx; ++i)
          {
            const double x = outlet_x * i / nx;
            const double width = quadratic_bow_tie_half_width(x);
            const bool boundary =
              i == 0 || i == nx || j == 0 || j == ny;
            mesh.nodes.push_back({x, eta * width, boundary ? 1 : 0});
          }
      }

    mesh.elements.reserve(2 * nx * ny);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        {
          const int n00 = node_id(i, j);
          const int n10 = node_id(i + 1, j);
          const int n11 = node_id(i + 1, j + 1);
          const int n01 = node_id(i, j + 1);

          CurvedElement lower;
          lower.node = {{n00, n10, n11}};
          mesh.elements.push_back(lower);

          CurvedElement upper;
          upper.node = {{n00, n11, n01}};
          mesh.elements.push_back(upper);
        }

    build_sides(mesh);
    return mesh;
  }

  void write_curved_mesh(const CurvedMesh &mesh, const std::string &mesh_base)
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

  std::vector<GeometryCase> default_geometry_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;

    const std::string curved_mesh_file =
      mesh_dir + "/generated_contracting_channel_quadratic_bow_tie_h0p25";
    const CurvedMesh curved_mesh = generate_quadratic_bow_tie_proxy_mesh();
    write_curved_mesh(curved_mesh, curved_mesh_file);

    return {{"linear_trapezoid_proxy",
             mesh_dir + "/contracting_channel_h0p25",
             output_dir + "/stokes_contracting_geometry_linear",
             ContractingChannelGeometry::linear_proxy,
             6.0,
             0.5},
            {"quadratic_bow_tie_like",
             curved_mesh_file,
             output_dir + "/stokes_contracting_geometry_quadratic",
             ContractingChannelGeometry::quadratic_bow_tie_proxy,
             outlet_x + outlet_x / 3.0,
             0.5}};
  }

  const VerticalProfile &profile_at(const ContractingChannelResult &result,
                                    const std::size_t index)
  {
    if (index >= result.profiles.size())
      throw std::runtime_error("contracting-channel profile missing");
    return result.profiles[index];
  }

  void print_result_row(const GeometryCase &geometry_case,
                        const std::string &method,
                        const ContractingChannelResult &result)
  {
    const auto &x0 = profile_at(result, 0);
    const auto &x1 = profile_at(result, 1);
    const auto &x4 = profile_at(result, 2);

    std::cout << geometry_case.name << ' '
              << method << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << geometry_case.area << ' '
              << geometry_case.outlet_half_width << ' '
              << x0.flow_rate << ' '
              << x1.flow_rate << ' '
              << x4.flow_rate << ' '
              << result.flux_imbalance << ' '
              << x4.ux_center << ' '
              << x4.ux_max << ' '
              << result.max_speed << ' '
              << result.divergence_l2_norm << ' '
              << result.pressure_min << ' '
              << result.pressure_max << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [output-prefix]\n";
      return 2;
    }

  const std::string output_prefix =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_contracting_channel_geometry_compare";

  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes contracting-channel geometry comparison "
                   "(book Exercise 5.7 / Exercise 1.5 proxy)\n";
      std::cout << "geometry method cells velocity_dofs pressure_dofs area "
                   "outlet_half_width x0_flow x1_flow x4_flow "
                   "flux_imbalance x4_center_ux x4_max_ux "
                   "max_speed divergence_L2 pressure_min pressure_max\n";

      for (const auto &geometry_case : default_geometry_cases())
        {
          set_contracting_channel_geometry(geometry_case.geometry);
          const auto taylor_hood =
            solve_contracting_channel(geometry_case.mesh_file,
                                      output_prefix + "_" +
                                        geometry_case.name + "_p2p1",
                                      2,
                                      1,
                                      false);
          const auto equal_order =
            solve_contracting_channel(geometry_case.mesh_file,
                                      output_prefix + "_" +
                                        geometry_case.name + "_p1p1_lps",
                                      1,
                                      1,
                                      true);

          print_result_row(geometry_case,
                           "P2P1_TaylorHood",
                           taylor_hood);
          print_result_row(geometry_case, "P1P1_LPS", equal_order);
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
