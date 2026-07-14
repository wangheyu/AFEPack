/**
 * @file navier_stokes_blasius_slit_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：Blasius 开缝模型。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - Blasius 开缝模型：在带薄板/开缝几何中近似平板边界层流动。
 * - Blasius 边界层代理：使用不可压 Navier-Stokes 方程构造平板边界层的数值代理模型。
 *
 * 网格与数据：主要依赖 平板边界层系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#define main navier_stokes_blasius_hidden_main
#include "navier_stokes_blasius_afepack.cpp"
#undef main

namespace
{
  enum class SlitPatch
  {
    upstream,
    lower,
    upper
  };

  bool is_slit_boundary_node(const double x, const double y)
  {
    if (near(x, -1.0) || near(x, 5.0) || near(y, -1.0) || near(y, 1.0))
      return true;
    return x > -boundary_tolerance &&
           x < 5.0 + boundary_tolerance &&
           near(y, 0.0);
  }

  int slit_node_tag(const double x, const double y, const SlitPatch patch)
  {
    if (near(x, 0.0))
      return 0;
    if (!near(y, 0.0))
      return 0;
    return patch == SlitPatch::upper ? 1 : 2;
  }

  int add_slit_node(
    GeneratedMesh &mesh,
    std::map<std::tuple<long long, long long, int>, int> &node_id,
    const double x,
    const double y,
    const SlitPatch patch)
  {
    const auto rounded = coordinate_key(x, y);
    const auto key =
      std::make_tuple(rounded.first,
                      rounded.second,
                      slit_node_tag(x, y, patch));
    const auto it = node_id.find(key);
    if (it != node_id.end())
      return it->second;

    const int id = static_cast<int>(mesh.nodes.size());
    mesh.nodes.push_back(
      {x, y, is_slit_boundary_node(x, y) ? blasius_boundary_mark(x, y) : 0});
    node_id.emplace(key, id);
    return id;
  }

  void add_slit_patch(
    GeneratedMesh &mesh,
    std::map<std::tuple<long long, long long, int>, int> &node_id,
    const std::vector<double> &x_coordinates,
    const std::vector<double> &y_coordinates,
    const SlitPatch patch)
  {
    std::vector<std::vector<int>> nodes(
      y_coordinates.size(),
      std::vector<int>(x_coordinates.size(), -1));

    for (std::size_t j = 0; j < y_coordinates.size(); ++j)
      for (std::size_t i = 0; i < x_coordinates.size(); ++i)
        nodes[j][i] = add_slit_node(mesh,
                                    node_id,
                                    x_coordinates[i],
                                    y_coordinates[j],
                                    patch);

    for (std::size_t j = 0; j + 1 < y_coordinates.size(); ++j)
      for (std::size_t i = 0; i + 1 < x_coordinates.size(); ++i)
        {
          const int n00 = nodes[j][i];
          const int n10 = nodes[j][i + 1];
          const int n11 = nodes[j + 1][i + 1];
          const int n01 = nodes[j + 1][i];

          ElementData lower;
          lower.node = {{n00, n10, n11}};
          mesh.elements.push_back(lower);

          ElementData upper;
          upper.node = {{n00, n11, n01}};
          mesh.elements.push_back(upper);
        }
  }

  void build_slit_sides(GeneratedMesh &mesh)
  {
    mesh.sides.clear();
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
            throw std::runtime_error("nonmanifold zero-thickness Blasius side");
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
        if (side.eb < 0)
          {
            const double midpoint_x =
              0.5 * (mesh.nodes[side.a].x + mesh.nodes[side.b].x);
            const double midpoint_y =
              0.5 * (mesh.nodes[side.a].y + mesh.nodes[side.b].y);
            side.boundary_mark = blasius_boundary_mark(midpoint_x,
                                                       midpoint_y);
            if (side.boundary_mark == 0)
              throw std::runtime_error(
                "unclassified Blasius slit boundary side");
          }
        else
          side.boundary_mark = 0;
        const double length =
          edge_length(mesh.nodes[side.a], mesh.nodes[side.b]);
        mesh.min_edge = std::min(mesh.min_edge, length);
        mesh.max_edge = std::max(mesh.max_edge, length);
      }
  }

  GeneratedMesh generate_zero_thickness_slit_mesh()
  {
    const std::vector<double> upstream_x = {-1.0, -0.5, 0.0};
    const std::vector<double> plate_x =
      {0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
    const std::vector<double> lower_y =
      {-1.0, -0.7, -0.45, -0.25, -0.12, -0.04, 0.0};
    const std::vector<double> upper_y =
      {0.0, 0.04, 0.12, 0.25, 0.45, 0.7, 1.0};
    const std::vector<double> full_y =
      {-1.0, -0.7, -0.45, -0.25, -0.12, -0.04, 0.0,
       0.04, 0.12, 0.25, 0.45, 0.7, 1.0};

    GeneratedMesh mesh;
    std::map<std::tuple<long long, long long, int>, int> node_id;
    add_slit_patch(mesh, node_id, upstream_x, full_y, SlitPatch::upstream);
    add_slit_patch(mesh, node_id, plate_x, lower_y, SlitPatch::lower);
    add_slit_patch(mesh, node_id, plate_x, upper_y, SlitPatch::upper);
    build_slit_sides(mesh);
    return mesh;
  }

  void print_slit_row(const BlasiusResult &result)
  {
    const double delta_x1 =
      average_valid_thickness(result.upper_layer_x1,
                              result.lower_layer_x1);
    const double delta_x4 =
      average_valid_thickness(result.upper_layer_x4,
                              result.lower_layer_x4);

    std::cout << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.system_nonzeros << ' '
              << result.nonlinear_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.final_nonlinear_relative_residual << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.max_speed << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << delta_x1 << ' '
              << delta_x4 << ' '
              << normalized_layer_thickness(delta_x1,
                                            result.viscosity,
                                            1.0)
              << ' '
              << normalized_layer_thickness(delta_x4,
                                            result.viscosity,
                                            4.0)
              << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [viscosities]\n";
      return 2;
    }

  const auto viscosities =
    parse_viscosities(argc > 1 ? argv[1] : "0.1,0.02");

  try
    {
      set_blasius_plate_geometry(0.0, 12.0);

      const std::string mesh_file =
        std::string(AFEPACK_DEFAULT_MESH_DIR) +
        "/generated_blasius_zero_thickness_slit";
      const GeneratedMesh slit_mesh = generate_zero_thickness_slit_mesh();
      write_generated_mesh(slit_mesh, mesh_file);

      const std::string output_prefix =
        std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
        "/navier_stokes_blasius_slit";
      const auto results =
        solve_blasius_continuation(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes Blasius zero-thickness slit "
                   "diagnostic (book Example 7.1.4 / Exercise 7.2)\n";
      std::cout << "linear_solver "
                << blasius_linear_solver_name(blasius_linear_solver())
                << " nonlinear_solver "
                << blasius_nonlinear_solver_name(blasius_nonlinear_solver())
                << " relaxation " << blasius_picard_relaxation() << '\n';
      std::cout << "mesh cells nodes min_edge max_edge\n";
      std::cout << "zero_thickness_slit "
                << slit_mesh.elements.size() << ' '
                << slit_mesh.nodes.size() << ' '
                << slit_mesh.min_edge << ' '
                << slit_mesh.max_edge << '\n';
      std::cout << "nu Re cells velocity_dofs pressure_dofs nonzeros "
                   "Nonlinear GMRES final_update final_gmres final_nonlinear "
                   "velocity_L2 velocity_H1 "
                   "max_speed "
                   "flux_imbalance divergence_L2 delta_x1_avg delta_x4_avg "
                   "delta_x1_over_sqrt_nu delta_x4_over_sqrt_4nu\n";
      for (const auto &result : results)
        print_slit_row(result);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
