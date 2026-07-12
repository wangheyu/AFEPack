/**
 * @file navier_stokes_step_flow_equal_order_reattachment_grid_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：再附长度网格实验。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 再附长度网格实验：比较网格尺度对后台阶回流区再附长度预测的影响。
 * - 后台阶流：在带突扩几何的通道中模拟回流区和再附长度等典型流动量。
 *
 * 网格与数据：主要依赖 后台阶通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define AFEPACK_NAVIER_STOKES_STEP_FLOW_EQUAL_ORDER_NO_MAIN
#include "navier_stokes_step_flow_equal_order_stabilized_afepack.cpp"
#undef AFEPACK_NAVIER_STOKES_STEP_FLOW_EQUAL_ORDER_NO_MAIN

namespace
{
  struct StepFlowMeshCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      h = 0.0;
  };

  struct EqualOrderReattachmentRow
  {
    StepFlowMeshCase                mesh_case;
    StabilizedNavierStokesStepResult result;
  };

  std::vector<StepFlowMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {{"step_channel_h0p50",
             mesh_dir + "/step_channel_h0p50",
             output_dir + "/navier_stokes_step_flow_equal_order_reattach_h0p50",
             0.50},
            {"step_channel_h0p25",
             mesh_dir + "/step_channel_h0p25",
             output_dir + "/navier_stokes_step_flow_equal_order_reattach_h0p25",
             0.25}};
  }

  void print_result_row(const EqualOrderReattachmentRow &row)
  {
    const auto &result = row.result;
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.picard_iterations << ' '
              << result.final_relative_update << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << result.wall_shear_samples << ' '
              << result.wall_shear_min << ' '
              << result.wall_shear_max << ' '
              << result.reattachment_found << ' '
              << result.reattachment_x << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.max_speed << '\n';
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
    parse_viscosities(argc > 1 ? argv[1] : "0.05,0.02,0.01");

  try
    {
      std::vector<EqualOrderReattachmentRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          const auto results =
            solve_stabilized_step_flow_continuation(mesh_case.mesh_file,
                                                    mesh_case.output_prefix,
                                                    viscosities);
          for (const auto &result : results)
            rows.push_back({mesh_case, result});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes P1/P1 local-projection "
                   "stabilized expansion step-flow reattachment grid/"
                   "Reynolds sweep (book Exercise 7.5)\n";
      std::cout << "pressure local-projection stabilization parameter beta="
                << equal_order_stabilization_parameter << '\n';
      std::cout << "mesh h nu Re cells velocity_dofs pressure_dofs "
                   "Picard final_update flux_imbalance divergence_L2 "
                   "shear_samples shear_min shear_max reattachment_found "
                   "reattachment_x velocity_L2 velocity_H1 max_speed\n";
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
