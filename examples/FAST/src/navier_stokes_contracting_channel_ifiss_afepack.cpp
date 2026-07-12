/**
 * @file navier_stokes_contracting_channel_ifiss_afepack.cpp
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
 * - IFISS 风格基准：使用与 IFISS 基准问题相近的收缩通道设置，便于跨代码结果对照。
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
#include <string>

#define main navier_stokes_contracting_channel_profiles_hidden_main
#include "navier_stokes_contracting_channel_profiles_afepack.cpp"
#undef main

int main(int argc, char **argv)
{
  if (argc > 4)
    {
      std::cerr << "usage: " << argv[0]
                << " [mesh-file [output-prefix [viscosities]]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) +
                 "/contracting_channel_ifiss_h0p125";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/navier_stokes_contracting_channel_ifiss";
  const auto requested_viscosities =
    parse_viscosities(argc > 3 ? argv[3] : std::string());

  try
    {
      set_contracting_channel_geometry(
        ContractingChannelGeometry::ifiss_bow_tie_half);
      const auto results = solve_contracting_channel_continuation(
        mesh_file,
        output_prefix,
        requested_viscosities,
        "1.D",
        1.0,
        false);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes original IFISS contracting-channel "
                   "profile diagnostic\n";
      std::cout << "Book Exercise 7.6: original quad_domain/quad_flow/quad_bc "
                   "geometry and boundary data, with triangular P2/P1DG_LPS "
                   "as the Q2-P_-1 analogue.\n";
      std::cout << "reference_inlet_flow_rate "
                << reference_flow_rate() << '\n';
      std::cout << "summary nu inverse_nu x y_min y_max samples flow_rate "
                   "flow_rate_minus_inlet ux_min ux_center ux_max\n";
      for (const auto &result : results)
        print_profile_summary(result);

      std::cout << "diagnostics nu inverse_nu cells velocity_dofs "
                   "pressure_dofs Picard rel_update inlet_rate outlet_rate "
                   "flux_imbalance divergence_L2 max_speed\n";
      for (const auto &result : results)
        std::cout << result.viscosity << ' '
                  << result.inverse_viscosity << ' '
                  << result.cells << ' '
                  << result.velocity_dofs << ' '
                  << result.pressure_dofs << ' '
                  << result.picard_iterations << ' '
                  << result.final_relative_update << ' '
                  << result.inlet_rate << ' '
                  << result.outlet_rate << ' '
                  << result.flux_imbalance << ' '
                  << result.divergence_l2_norm << ' '
                  << result.max_speed << '\n';
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
