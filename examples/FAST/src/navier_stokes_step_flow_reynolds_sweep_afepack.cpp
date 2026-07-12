/**
 * @file navier_stokes_step_flow_reynolds_sweep_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：雷诺数扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 雷诺数扫描：批量改变雷诺数并跟踪非线性求解和流动特征量变化。
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

#define main navier_stokes_step_flow_schur_hidden_main
#include "navier_stokes_step_flow_schur_afepack.cpp"
#undef main

namespace
{
  constexpr const char *default_reynolds_viscosities =
    "1,0.5,0.25,0.125,0.05,0.02,0.01";

  void print_result_row(const StepFlowResult &result)
  {
    std::cout << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.picard_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << result.wall_shear_samples << ' '
              << result.wall_shear_min << ' '
              << result.wall_shear_max << ' '
              << result.reattachment_found << ' '
              << result.reattachment_x << '\n';
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
    parse_viscosities(argc > 1 ? argv[1] : default_reynolds_viscosities);

  const std::string mesh_file =
    std::string(AFEPACK_DEFAULT_MESH_DIR) + "/step_channel_h0p25";
  const std::string output_prefix =
    std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
    "/navier_stokes_step_flow_reynolds_sweep_afepack";

  try
    {
      const auto results =
        solve_step_flow_continuation(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes expansion step-flow "
                   "Reynolds sweep (book Example 7.1.2 and "
                   "Exercise 7.5)\n";
      std::cout << "mesh step_channel_h0p25\n";
      std::cout << "preconditioner tri(Aamg2,Mp-pcg)\n";
      std::cout << "nu Re cells velocity_dofs pressure_dofs Picard GMRES "
                   "final_update final_gmres velocity_L2 velocity_H1 "
                   "flux_imbalance divergence_L2 shear_samples shear_min "
                   "shear_max reattachment_found reattachment_x\n";
      for (const auto &result : results)
        print_result_row(result);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
