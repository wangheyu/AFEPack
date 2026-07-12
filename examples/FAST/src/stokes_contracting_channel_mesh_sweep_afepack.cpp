/**
 * @file stokes_contracting_channel_mesh_sweep_afepack.cpp
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
 * - 网格尺度扫描：在多组网格输入上重复求解，比较自由度增长对误差和迭代次数的影响。
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
  struct ContractingMeshCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      h = 0.0;
  };

  struct ContractingSweepRow
  {
    ContractingMeshCase    mesh_case;
    std::string            method;
    ContractingChannelResult result;
  };

  std::vector<ContractingMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {{"contracting_channel_h0p50",
             mesh_dir + "/contracting_channel_h0p50",
             output_dir + "/stokes_contracting_channel_mesh_h0p50",
             0.50},
            {"contracting_channel_h0p25",
             mesh_dir + "/contracting_channel_h0p25",
             output_dir + "/stokes_contracting_channel_mesh_h0p25",
             0.25}};
  }

  const VerticalProfile &profile_at(const ContractingChannelResult &result,
                                    const std::size_t index)
  {
    if (index >= result.profiles.size())
      throw std::runtime_error("missing contracting-channel profile");
    return result.profiles[index];
  }

  void print_result_row(const ContractingSweepRow &row)
  {
    const auto &result = row.result;
    const auto &x0 = profile_at(result, 0);
    const auto &x1 = profile_at(result, 1);
    const auto &x4 = profile_at(result, 2);
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << row.method << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << x0.flow_rate << ' '
              << x1.flow_rate << ' '
              << x4.flow_rate << ' '
              << x0.ux_center << ' '
              << x1.ux_center << ' '
              << x4.ux_center << ' '
              << result.max_speed << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc != 1)
    {
      std::cerr << "usage: " << argv[0] << '\n';
      return 2;
    }

  try
    {
      std::vector<ContractingSweepRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          rows.push_back({mesh_case,
                          "P2P1_TaylorHood",
                          solve_contracting_channel(mesh_case.mesh_file,
                                                    mesh_case.output_prefix +
                                                      "_p2p1",
                                                    2,
                                                    1,
                                                    false)});
          rows.push_back({mesh_case,
                          "P1P1_LPS",
                          solve_contracting_channel(mesh_case.mesh_file,
                                                    mesh_case.output_prefix +
                                                      "_p1p1_lps",
                                                    1,
                                                    1,
                                                    true)});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes contracting-channel mesh sweep "
                   "(book Exercise 5.7 proxy)\n";
      std::cout << "mesh h method cells velocity_dofs pressure_dofs "
                   "flow_x0 flow_x1 flow_x4 ux_center_x0 ux_center_x1 "
                   "ux_center_x4 max_speed flux_imbalance divergence_L2\n";
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
