/**
 * @file stokes_contracting_channel_discontinuous_pressure_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：间断压力空间。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 间断压力空间：比较连续压力和间断压力离散在收缩通道问题中的差异。
 * - 收缩通道流：在收缩几何中考察入口、出口和壁面边界对流动结构的影响。
 *
 * 网格与数据：主要依赖 收缩通道网格、收缩通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define main stokes_contracting_channel_profiles_hidden_main
#include "stokes_contracting_channel_profiles_afepack.cpp"
#undef main

namespace
{
  struct DiscreteSpaceCase
  {
    std::string name;
    int         velocity_order = 0;
    int         pressure_order = 0;
    bool        stabilized = false;
    bool        discontinuous_pressure = false;
    std::string output_suffix;
  };

  struct DiscreteSpaceRow
  {
    DiscreteSpaceCase       space_case;
    ContractingChannelResult result;
  };

  std::vector<DiscreteSpaceCase> default_space_cases()
  {
    return {{"P2P1_TaylorHood", 2, 1, false, false, "p2p1"},
            {"P2P1DG_LPS_discontinuous_pressure", 2, 1, true, true, "p2p1dg_lps"},
            {"P2P0_discontinuous_pressure", 2, 0, false, false, "p2p0"},
            {"P1P1_LPS", 1, 1, true, false, "p1p1_lps"}};
  }

  const VerticalProfile &profile_at(const ContractingChannelResult &result,
                                    const std::size_t index)
  {
    if (index >= result.profiles.size())
      throw std::runtime_error("missing contracting-channel profile");
    return result.profiles[index];
  }

  void print_result_row(const DiscreteSpaceRow &row)
  {
    const auto &space = row.space_case;
    const auto &result = row.result;
    const auto &x0 = profile_at(result, 0);
    const auto &x1 = profile_at(result, 1);
    const auto &x4 = profile_at(result, 2);
    std::cout << space.name << ' '
              << space.velocity_order << ' '
              << space.pressure_order << ' '
              << (space.discontinuous_pressure ? 1 : 0) << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << x0.flow_rate << ' '
              << x1.flow_rate << ' '
              << x4.flow_rate << ' '
              << result.flux_imbalance << ' '
              << x4.ux_center << ' '
              << x4.ux_max << ' '
              << result.max_speed << ' '
              << result.divergence_l2_norm << ' '
              << result.max_cell_divergence_integral << ' '
              << result.cell_mean_divergence_l2_norm << ' '
              << result.pressure_min << ' '
              << result.pressure_max << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) +
                 "/contracting_channel_ifiss_h0p125";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_contracting_channel_discontinuous_pressure";

  try
    {
      set_contracting_channel_geometry(
        ContractingChannelGeometry::ifiss_bow_tie_half);

      std::vector<DiscreteSpaceRow> rows;
      for (const auto &space_case : default_space_cases())
        rows.push_back({space_case,
                        solve_contracting_channel(
                          mesh_file,
                          output_prefix + "_" + space_case.output_suffix,
                          space_case.velocity_order,
                          space_case.pressure_order,
                          space_case.stabilized,
                          space_case.discontinuous_pressure,
                          false)});

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes contracting-channel discontinuous "
                   "pressure diagnostic\n";
      std::cout << "Book Exercise 5.7 on the original IFISS quad_domain, "
                   "quad_flow, and quad_bc geometry/data: triangular "
                   "P2/P1DG_LPS is the "
                   "closest stabilized discontinuous-linear pressure "
                   "analogue of Q2-P_-1 on the current AFEPack triangle "
                   "mesh family; P2/P0 is retained as a lower-order "
                   "discontinuous pressure control.\n";
      std::cout << "reference_inlet_flow_rate "
                << reference_flow_rate() << '\n';
      std::cout << "method velocity_order pressure_order pressure_DG cells "
                   "velocity_dofs pressure_dofs flow_x0 flow_x1 flow_x4 "
                   "flux_imbalance ux_center_x4 ux_max_x4 max_speed "
                   "divergence_L2 max_cell_divergence_integral "
                   "cell_mean_divergence_L2 pressure_min pressure_max\n";
      for (const auto &row : rows)
        print_result_row(row);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
