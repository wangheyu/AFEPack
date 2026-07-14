/**
 * @file stokes_step_flow_profiles_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：后台阶流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 后台阶流：在带突扩几何的通道中模拟回流区和再附长度等典型流动量。
 * - 剖面数据输出：沿指定截线输出速度或压力剖面，便于与基准数据和图形对比。
 *
 * 网格与数据：主要依赖 后台阶通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define AFEPACK_STOKES_STEP_EQUAL_ORDER_NO_MAIN
#include "stokes_step_flow_equal_order_stabilized_afepack.cpp"
#undef AFEPACK_STOKES_STEP_EQUAL_ORDER_NO_MAIN

namespace
{
  constexpr double exact_step_flow_rate = 2.0 / 3.0;

  void print_profile_summary(const std::string &method,
                             const std::vector<VerticalProfile> &profiles)
  {
    for (const auto &profile : profiles)
      std::cout << method << ' '
                << profile.x << ' '
                << profile.y_min << ' '
                << profile.y_max << ' '
                << profile.samples.size() << ' '
                << profile.flow_rate << ' '
                << profile.flow_rate - exact_step_flow_rate << ' '
                << profile.ux_min << ' '
                << profile.ux_center << ' '
                << profile.ux_max << '\n';
  }

  void print_profile_samples(const std::string &method,
                             const std::vector<VerticalProfile> &profiles)
  {
    for (const auto &profile : profiles)
      for (const auto &sample : profile.samples)
        std::cout << method << ' '
                  << profile.x << ' '
                  << sample.y << ' '
                  << sample.ux << '\n';
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
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/step_channel_h0p25";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_step_flow_profiles_afepack";

  try
    {
      const StokesResult taylor_hood =
        solve_stokes(mesh_file, output_prefix + "_taylor_hood");
      const StabilizedStepResult equal_order =
        solve_stabilized_step(mesh_file, output_prefix + "_p1p1_lps", 0.25);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes step-flow vertical profile diagnostic\n";
      std::cout << "Book Exercise 5.6 analogue on Example 5.1.2: "
                   "u_x profiles at x=-1, x=0, and x=5.\n";
      std::cout << "reference_flow_rate " << exact_step_flow_rate << '\n';
      std::cout << "summary method x y_min y_max samples flow_rate "
                   "flow_rate_minus_reference ux_min ux_center ux_max\n";
      print_profile_summary("P2P1_TaylorHood", taylor_hood.profiles);
      print_profile_summary("P1P1_LPS", equal_order.profiles);

      std::cout << "profile_samples method x y ux\n";
      print_profile_samples("P2P1_TaylorHood", taylor_hood.profiles);
      print_profile_samples("P1P1_LPS", equal_order.profiles);

      std::cout << "diagnostics method cells velocity_dofs pressure_dofs "
                   "velocity_L2 velocity_H1 divergence_L2\n";
      std::cout << "P2P1_TaylorHood "
                << taylor_hood.cells << ' '
                << taylor_hood.velocity_dofs << ' '
                << taylor_hood.pressure_dofs << ' '
                << taylor_hood.velocity_l2_norm << ' '
                << taylor_hood.velocity_h1_seminorm << ' '
                << taylor_hood.divergence_l2_norm << '\n';
      std::cout << "P1P1_LPS "
                << equal_order.cells << ' '
                << equal_order.velocity_dofs << ' '
                << equal_order.pressure_dofs << ' '
                << equal_order.velocity_l2_norm << ' '
                << equal_order.velocity_h1_seminorm << ' '
                << equal_order.divergence_l2_norm << '\n';
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
