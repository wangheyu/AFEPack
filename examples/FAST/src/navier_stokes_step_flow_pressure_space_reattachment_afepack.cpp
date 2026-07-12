/**
 * @file navier_stokes_step_flow_pressure_space_reattachment_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：压力空间再附长度对比。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 压力空间再附长度对比：评估不同压力空间选择对后台阶再附长度计算的影响。
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
#include <stdexcept>
#include <string>
#include <vector>

#define main navier_stokes_step_flow_schur_hidden_main
#include "navier_stokes_step_flow_schur_afepack.cpp"
#undef main

namespace
{
  constexpr double discontinuous_pressure_stabilization = 1.0;

  struct PressureMethod
  {
    std::string name;
    std::string template_name;
    double stabilization = 0.0;
    bool pin_pressure = true;
  };

  std::vector<PressureMethod> pressure_methods(const std::string &name)
  {
    const PressureMethod continuous =
      {"P2P1_continuous", "1", 0.0, true};
    const PressureMethod discontinuous =
      {"P2P1DG_LPS", "1.D", discontinuous_pressure_stabilization, false};

    if (name == "both")
      return {continuous, discontinuous};
    if (name == "continuous")
      return {continuous};
    if (name == "discontinuous")
      return {discontinuous};
    throw std::runtime_error(
      "method must be both, continuous, or discontinuous");
  }

  void print_result(const PressureMethod &method,
                    const StepFlowResult &result)
  {
    if (result.wall_shear_sampling.empty())
      throw std::runtime_error("missing wall-shear sampling diagnostic");
    const auto &wall_shear = result.wall_shear_sampling.back();

    std::cout << method.name << ' '
              << (method.template_name == "1.D" ? 1 : 0) << ' '
              << method.stabilization << ' '
              << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.system_nonzeros << ' '
              << result.picard_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.inflow_rate << ' '
              << result.outflow_rate << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << result.max_cell_divergence_integral << ' '
              << result.cell_mean_divergence_l2_norm << ' '
              << wall_shear.first << ' '
              << wall_shear.second.reattachment_found << ' '
              << wall_shear.second.reattachment_x << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 5)
    {
      std::cerr << "usage: " << argv[0]
                << " [mesh-file [output-prefix [method [viscosities]]]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) +
                 "/step_channel_h0p50";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/navier_stokes_step_flow_pressure_space_reattachment";
  const std::string method_name = argc > 3 ? argv[3] : "both";
  const std::string viscosity_text = argc > 4 ? argv[4] : "0.05,0.02";

  try
    {
      const auto methods = pressure_methods(method_name);
      const auto viscosities = parse_viscosities(viscosity_text);
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes step-flow pressure-space "
                   "reattachment comparison (book Exercise 7.5)\n";
      std::cout << "Triangular analogues on bisected square grids: "
                   "P2/P1 continuous pressure is non-locally-conservative; "
                   "P2/P1DG_LPS retains element-constant divergence tests.\n";
      std::cout << "method pressure_DG stabilization nu Re cells "
                   "velocity_dofs pressure_dofs nonzeros Picard GMRES "
                   "final_update final_gmres inflow outflow flux_imbalance "
                   "divergence_L2 max_cell_divergence_integral "
                   "cell_mean_divergence_L2 shear_subsamples reattachment_found "
                   "reattachment_x\n";

      for (const auto &method : methods)
        {
          const auto results = solve_step_flow_continuation(
            mesh_file,
            output_prefix + "_" + method.name,
            viscosities,
            method.template_name,
            method.stabilization,
            method.pin_pressure,
            true);
          for (const auto &result : results)
            print_result(method, result);
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
