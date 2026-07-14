/**
 * @file navier_stokes_contracting_channel_mesh_sweep_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：收缩通道流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
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

#define main navier_stokes_contracting_channel_profiles_hidden_main
#include "navier_stokes_contracting_channel_profiles_afepack.cpp"
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

  struct ContractingNSSweepRow
  {
    ContractingMeshCase      mesh_case;
    ContractingChannelResult result;
  };

  std::vector<ContractingMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {{"contracting_channel_h0p50",
             mesh_dir + "/contracting_channel_h0p50",
             output_dir + "/navier_stokes_contracting_channel_mesh_h0p50",
             0.50},
            {"contracting_channel_h0p25",
             mesh_dir + "/contracting_channel_h0p25",
             output_dir + "/navier_stokes_contracting_channel_mesh_h0p25",
             0.25}};
  }

  const VerticalProfile &profile_at(const ContractingChannelResult &result,
                                    const std::size_t index)
  {
    if (index >= result.profiles.size())
      throw std::runtime_error("missing contracting-channel profile");
    return result.profiles[index];
  }

  void print_result_row(const ContractingNSSweepRow &row)
  {
    const auto &result = row.result;
    const auto &x4 = profile_at(result, 2);
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.viscosity << ' '
              << result.inverse_viscosity << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.picard_iterations << ' '
              << result.final_relative_update << ' '
              << x4.flow_rate << ' '
              << x4.ux_center << ' '
              << x4.ux_max << ' '
              << result.max_speed << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << '\n';
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
    parse_viscosities(argc > 1 ? argv[1] : "1,0.01");

  try
    {
      std::vector<ContractingNSSweepRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          const auto results =
            solve_contracting_channel_continuation(mesh_case.mesh_file,
                                                   mesh_case.output_prefix,
                                                   viscosities);
          for (const auto &result : results)
            rows.push_back({mesh_case, result});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes contracting-channel "
                   "mesh/viscosity sweep (book Exercise 7.6 proxy)\n";
      std::cout << "mesh h nu inverse_nu cells velocity_dofs pressure_dofs "
                   "Picard final_update flow_x4 ux_center_x4 ux_max_x4 "
                   "max_speed flux_imbalance divergence_L2\n";
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
