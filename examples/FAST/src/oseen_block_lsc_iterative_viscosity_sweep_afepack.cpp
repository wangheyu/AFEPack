/**
 * @file oseen_block_lsc_iterative_viscosity_sweep_afepack.cpp
 * @brief AFEPack Oseen 线性化方程算例：黏性系数扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 Oseen 线性化方程 迁移算例，关注Navier-Stokes 线性化后得到的鞍点系统及其 Schur 补预条件器。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 黏性系数扫描：批量改变黏性系数或雷诺数，观察预条件器和非线性迭代的稳健性。
 * - 迭代线性求解：使用 Krylov 迭代和块预条件器替代直接求解，以观察可扩展性。
 * - LSC 预条件器：构造最小二乘交换子 Schur 补近似，比较其稳健性和代价。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define main oseen_block_lsc_iterative_afepack_main
#include "oseen_block_lsc_iterative_afepack.cpp"
#undef main

namespace
{
  std::vector<double> parse_viscosities(const std::string &csv)
  {
    const std::string values =
      csv.empty() ? "1,0.5,0.25,0.125,0.0625" : csv;
    std::stringstream stream(values);
    std::string token;
    std::vector<double> result;

    while (std::getline(stream, token, ','))
      {
        if (token.empty())
          continue;
        const double value = std::stod(token);
        if (value <= 0.0)
          throw std::runtime_error("viscosities must be positive");
        result.push_back(value);
      }

    if (result.empty())
      throw std::runtime_error("viscosity list is empty");

    return result;
  }

  std::string viscosity_suffix(const double value)
  {
    std::ostringstream out;
    out << std::scientific << std::setprecision(0) << value;
    std::string suffix = out.str();
    for (char &c : suffix)
      if (c == '-' || c == '+' || c == '.')
        c = 'p';
    return suffix;
  }

  struct SweepRow
  {
    double viscosity = 0.0;
    OseenLSCIterativeResult result;
  };

  void print_sweep_row(const SweepRow &row)
  {
    std::cout << row.result.preconditioner << ' '
              << row.viscosity << ' '
              << 1.0 / row.viscosity << ' '
              << row.result.cells << ' '
              << row.result.velocity_dofs << ' '
              << row.result.pressure_dofs << ' '
              << row.result.gmres_iterations << ' '
              << (row.result.gmres_converged ? "ok" : "limit") << ' '
              << row.result.gmres_relative_residual << ' '
              << row.result.velocity_l2_error << ' '
              << row.result.velocity_h1_error << ' '
              << row.result.pressure_l2_error << ' '
              << row.result.divergence_l2_norm << '\n';
  }
}

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
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p20";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/oseen_block_lsc_iterative_viscosity_sweep_afepack";
  const auto viscosities =
    parse_viscosities(argc > 3 ? argv[3] : std::string());

  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Oseen block LSC iterative inner solvers "
                   "viscosity sweep Taylor-Hood manufactured solution\n";
      std::cout << "preconditioner nu formal_Re cells velocity_dofs "
                   "pressure_dofs GMRES status relative_residual velocity_L2 "
                   "velocity_H1 pressure_L2 divergence_L2\n";

      for (const double nu : viscosities)
        {
          viscosity = nu;
          const std::string case_prefix =
            output_prefix + "_nu" + viscosity_suffix(nu);

          const auto fbicg_result =
            solve_oseen_block_lsc_iterative(
              mesh_file,
              case_prefix,
              VelocityInverseKind::jacobi_bicgstab);
          print_sweep_row({nu, fbicg_result});

          const auto famg_result =
            solve_oseen_block_lsc_iterative(
              mesh_file,
              case_prefix,
              VelocityInverseKind::amg_diffusion);
          print_sweep_row({nu, famg_result});

          const auto famg_bicg_result =
            solve_oseen_block_lsc_iterative(
              mesh_file,
              case_prefix,
              VelocityInverseKind::amg_preconditioned_oseen);
          print_sweep_row({nu, famg_bicg_result});
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
