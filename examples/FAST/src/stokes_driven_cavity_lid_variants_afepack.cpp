/**
 * @file stokes_driven_cavity_lid_variants_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：驱动方腔流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 驱动方腔流：经典顶盖驱动方腔基准问题，用于验证速度-压力耦合和涡结构。
 * - 顶盖边界变体：比较不同顶盖速度边界设置对方腔流场和压力场的影响。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define AFEPACK_STOKES_CAVITY_NO_MAIN
#include "stokes_driven_cavity_afepack.cpp"

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p10";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_driven_cavity_lid_variants_afepack";

  const std::vector<CavityLidVariant> variants = {
    CavityLidVariant::Leaky,
    CavityLidVariant::Watertight,
    CavityLidVariant::Regularized};

  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes driven cavity lid-variant comparison\n";
      std::cout << "variant cells velocity_dofs pressure_dofs velocity_L2 "
                   "velocity_H1 max_speed pressure_min pressure_max "
                   "pressure_mean pressure_centerline_max "
                   "pressure_antisymmetry_RMS divergence_L2\n";

      for (const CavityLidVariant variant : variants)
        {
          const std::string name = lid_variant_name(variant);
          const StokesResult result =
            solve_stokes(mesh_file, output_prefix + "_" + name, variant);

          std::cout << name << ' ' << result.cells << ' '
                    << result.velocity_dofs << ' ' << result.pressure_dofs
                    << ' ' << result.velocity_l2_norm << ' '
                    << result.velocity_h1_seminorm << ' '
                    << result.max_speed << ' ' << result.pressure_min << ' '
                    << result.pressure_max << ' ' << result.pressure_mean
                    << ' ' << result.pressure_centerline_max << ' '
                    << result.pressure_antisymmetry_rms << ' '
                    << result.divergence_l2_norm << '\n';
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
