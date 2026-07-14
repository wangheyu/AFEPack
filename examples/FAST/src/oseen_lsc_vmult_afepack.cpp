/**
 * @file oseen_lsc_vmult_afepack.cpp
 * @brief AFEPack Oseen 线性化方程算例：矩阵向量乘验证。
 *
 * @details
 * 本文件是 FAST 目录中的 Oseen 线性化方程 迁移算例，关注Navier-Stokes 线性化后得到的鞍点系统及其 Schur 补预条件器。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 矩阵向量乘验证：单独验证预条件器或算子乘法例程，便于排查块方法实现。
 * - LSC 预条件器：构造最小二乘交换子 Schur 补近似，比较其稳健性和代价。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>

#include <lapacke.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;
  constexpr double viscosity = 1.0;
  constexpr double hp_regularization = 1.0e-8;

  class DenseLU
  {
  public:
    DenseLU() = default;

    DenseLU(const std::vector<double> &matrix, const int size)
    {
      factor(matrix, size);
    }

    void factor(const std::vector<double> &matrix, const int size)
    {
      n = size;
      lu = matrix;
      pivots.assign(n, 0);

      const int info =
        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, lu.data(), n, pivots.data());
      if (info != 0)
        throw std::runtime_error("LAPACKE_dgetrf failed with info=" +
                                 std::to_string(info));
    }

    std::vector<double> solve(const std::vector<double> &rhs) const
    {
      if (static_cast<int>(rhs.size()) != n)
        throw std::runtime_error("DenseLU solve size mismatch");

      std::vector<double> solution = rhs;
      const int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,
                                      'N',
                                      n,
                                      1,
                                      lu.data(),
                                      n,
                                      pivots.data(),
                                      solution.data(),
                                      1);
      if (info != 0)
        throw std::runtime_error("LAPACKE_dgetrs failed with info=" +
                                 std::to_string(info));

      return solution;
    }

  private:
    int n = 0;
    std::vector<double> lu;
    std::vector<int> pivots;
  };

  struct LSCResult
  {
    int          cells = 0;
    unsigned int scalar_velocity_dofs = 0;
    unsigned int velocity_block_dofs = 0;
    unsigned int pressure_dofs = 0;
    double       hp_frobenius_norm = 0.0;
    double       gp_frobenius_norm = 0.0;
    double       input_norm = 0.0;
    double       first_hp_inverse_norm = 0.0;
    double       gp_action_norm = 0.0;
    double       lsc_output_norm = 0.0;
  };

  double sx(const double *p) { return std::sin(pi * p[0]); }
  double cx(const double *p) { return std::cos(pi * p[0]); }
  double sy(const double *p) { return std::sin(pi * p[1]); }
  double cy(const double *p) { return std::cos(pi * p[1]); }

  double velocity_x_exact(const double *p)
  {
    return 2.0 * pi * sx(p) * sx(p) * sy(p) * cy(p);
  }

  double velocity_y_exact(const double *p)
  {
    return -2.0 * pi * sx(p) * cx(p) * sy(p) * sy(p);
  }

  std::string default_template_path()
  {
    if (const char *env = std::getenv("AFEPACK_PATH"))
      return std::string(env) + "/template/triangle";
    return "/home/steve/local/AFEPack/include/AFEPack/template/triangle";
  }

  void ensure_template_path()
  {
    if (std::getenv("AFEPACK_TEMPLATE_PATH") != nullptr)
      return;

    const std::string path = default_template_path();
    if (::setenv("AFEPACK_TEMPLATE_PATH", path.c_str(), 0) != 0)
      throw std::runtime_error("failed to set AFEPACK_TEMPLATE_PATH");
  }

  std::vector<TemplateElement<double, dim, dim>>
  make_triangle_template(const int order,
                         TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle." + std::to_string(order) + ".tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle." + std::to_string(order) + ".bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
  }

  FEMSpace<double, dim> build_space(
    EasyMesh &mesh,
    std::vector<TemplateElement<double, dim, dim>> &template_element)
  {
    FEMSpace<double, dim> fem_space;
    fem_space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      fem_space.element(i).reinit(fem_space, i, 0);

    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    return fem_space;
  }

  int velocity_index(const int component,
                     const int dof,
                     const int n_scalar_velocity_dofs)
  {
    return component * n_scalar_velocity_dofs + dof;
  }

  double vector_dot(const std::vector<double> &a,
                    const std::vector<double> &b)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
      result += a[i] * b[i];
    return result;
  }

  double vector_norm(const std::vector<double> &v)
  {
    return std::sqrt(vector_dot(v, v));
  }

  double dense_frobenius_norm(const std::vector<double> &matrix)
  {
    double sum = 0.0;
    for (const double entry : matrix)
      sum += entry * entry;
    return std::sqrt(sum);
  }

  void subtract_arithmetic_mean(std::vector<double> &v)
  {
    double sum = 0.0;
    for (const double value : v)
      sum += value;

    const double mean = sum / static_cast<double>(v.size());
    for (double &value : v)
      value -= mean;
  }

  void add_entry(std::vector<double> &matrix,
                 const int n_cols,
                 const int row,
                 const int col,
                 const double value)
  {
    matrix[static_cast<std::size_t>(row) * n_cols + col] += value;
  }

  std::vector<double> dense_vmult(const std::vector<double> &matrix,
                                  const std::vector<double> &x,
                                  const int n_rows,
                                  const int n_cols)
  {
    std::vector<double> y(n_rows, 0.0);
    for (int i = 0; i < n_rows; ++i)
      for (int j = 0; j < n_cols; ++j)
        y[i] += matrix[static_cast<std::size_t>(i) * n_cols + j] * x[j];
    return y;
  }

  std::vector<double> dense_tvmult(const std::vector<double> &matrix,
                                   const std::vector<double> &x,
                                   const int n_rows,
                                   const int n_cols)
  {
    std::vector<double> y(n_cols, 0.0);
    for (int i = 0; i < n_rows; ++i)
      for (int j = 0; j < n_cols; ++j)
        y[j] += matrix[static_cast<std::size_t>(i) * n_cols + j] * x[i];
    return y;
  }

  void apply_inverse_velocity_mass(std::vector<double> &vector,
                                   const std::vector<double> &mass_diagonal)
  {
    for (std::size_t i = 0; i < vector.size(); ++i)
      vector[i] /= mass_diagonal[i];
  }

  struct AssembledLSCOperators
  {
    int scalar_velocity_dofs = 0;
    int velocity_block_dofs = 0;
    int pressure_dofs = 0;
    std::vector<double> velocity_oseen;
    std::vector<double> divergence;
    std::vector<double> velocity_mass_diagonal;
  };

  AssembledLSCOperators assemble_lsc_operators(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space)
  {
    const int n_scalar_velocity_dofs = velocity_space.n_dof();
    const int n_velocity_block_dofs = dim * n_scalar_velocity_dofs;
    const int n_pressure_dofs = pressure_space.n_dof();

    AssembledLSCOperators operators;
    operators.scalar_velocity_dofs = n_scalar_velocity_dofs;
    operators.velocity_block_dofs = n_velocity_block_dofs;
    operators.pressure_dofs = n_pressure_dofs;
    operators.velocity_oseen.assign(
      static_cast<std::size_t>(n_velocity_block_dofs) *
        n_velocity_block_dofs,
      0.0);
    operators.divergence.assign(
      static_cast<std::size_t>(n_pressure_dofs) *
        n_velocity_block_dofs,
      0.0);
    operators.velocity_mass_diagonal.assign(n_velocity_block_dofs, 0.0);

    for (int e = 0; e < velocity_space.n_element(); ++e)
      {
        const auto &velocity_element = velocity_space.element(e);
        const auto &pressure_element = pressure_space.element(e);
        const double volume = velocity_element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          velocity_element.findQuadratureInfo(5);
        const auto jacobian =
          velocity_element.local_to_global_jacobian(
            quad_info.quadraturePoint());
        const auto q_point =
          velocity_element.local_to_global(quad_info.quadraturePoint());
        const auto velocity_values =
          velocity_element.basis_function_value(q_point);
        const auto velocity_grads =
          velocity_element.basis_function_gradient(q_point);
        const auto pressure_values =
          pressure_element.basis_function_value(q_point);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double wind[dim] = {velocity_x_exact(p),
                                      velocity_y_exact(p)};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                double lumped_mass = 0.0;

                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    const double diffusion =
                      viscosity * innerProduct(velocity_grads[i][q],
                                               velocity_grads[j][q]);
                    const double convection =
                      (wind[0] * velocity_grads[j][q][0] +
                       wind[1] * velocity_grads[j][q][1]) *
                      velocity_values[i][q];
                    const double velocity_entry =
                      jxw * (diffusion + convection);

                    for (int component = 0; component < dim; ++component)
                      add_entry(operators.velocity_oseen,
                                n_velocity_block_dofs,
                                velocity_index(component,
                                               vi,
                                               n_scalar_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_scalar_velocity_dofs),
                                velocity_entry);

                    lumped_mass +=
                      jxw * velocity_values[i][q] * velocity_values[j][q];
                  }

                for (int component = 0; component < dim; ++component)
                  operators.velocity_mass_diagonal[velocity_index(
                    component, vi, n_scalar_velocity_dofs)] += lumped_mass;
              }

            for (int i = 0; i < n_local_pressure_dofs; ++i)
              {
                const int pi_dof = pressure_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(operators.divergence,
                                n_velocity_block_dofs,
                                pi_dof,
                                velocity_index(component,
                                               vj,
                                               n_scalar_velocity_dofs),
                                -jxw * pressure_values[i][q] *
                                  velocity_grads[j][q][component]);
                  }
              }
          }
      }

    for (double &entry : operators.velocity_mass_diagonal)
      if (std::abs(entry) < 1.0e-14)
        entry = 1.0;

    return operators;
  }

  std::vector<double> apply_hp(const AssembledLSCOperators &operators,
                               const std::vector<double> &src)
  {
    auto velocity_tmp =
      dense_tvmult(operators.divergence,
                   src,
                   operators.pressure_dofs,
                   operators.velocity_block_dofs);
    apply_inverse_velocity_mass(velocity_tmp,
                                operators.velocity_mass_diagonal);
    return dense_vmult(operators.divergence,
                       velocity_tmp,
                       operators.pressure_dofs,
                       operators.velocity_block_dofs);
  }

  std::vector<double> apply_gp(const AssembledLSCOperators &operators,
                               const std::vector<double> &src)
  {
    auto velocity_tmp_1 =
      dense_tvmult(operators.divergence,
                   src,
                   operators.pressure_dofs,
                   operators.velocity_block_dofs);
    apply_inverse_velocity_mass(velocity_tmp_1,
                                operators.velocity_mass_diagonal);

    auto velocity_tmp_2 =
      dense_vmult(operators.velocity_oseen,
                  velocity_tmp_1,
                  operators.velocity_block_dofs,
                  operators.velocity_block_dofs);
    apply_inverse_velocity_mass(velocity_tmp_2,
                                operators.velocity_mass_diagonal);

    return dense_vmult(operators.divergence,
                       velocity_tmp_2,
                       operators.pressure_dofs,
                       operators.velocity_block_dofs);
  }

  std::vector<double> build_dense_pressure_operator(
    const AssembledLSCOperators &operators,
    const bool use_gp)
  {
    std::vector<double> matrix(
      static_cast<std::size_t>(operators.pressure_dofs) *
        operators.pressure_dofs,
      0.0);

    for (int j = 0; j < operators.pressure_dofs; ++j)
      {
        std::vector<double> basis_vector(operators.pressure_dofs, 0.0);
        basis_vector[j] = 1.0;

        auto column = use_gp ? apply_gp(operators, basis_vector) :
                               apply_hp(operators, basis_vector);
        for (int i = 0; i < operators.pressure_dofs; ++i)
          matrix[static_cast<std::size_t>(i) * operators.pressure_dofs + j] =
            column[i];
      }

    return matrix;
  }

  LSCResult run_lsc_vmult(const std::string &mesh_file)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> velocity_geometry;
    CoordTransform<dim, dim> velocity_transform;
    TemplateDOF<dim> velocity_dof;
    BasisFunctionAdmin<double, dim, dim> velocity_basis(velocity_dof);
    auto velocity_template = make_triangle_template(2,
                                                    velocity_geometry,
                                                    velocity_transform,
                                                    velocity_dof,
                                                    velocity_basis);

    TemplateGeometry<dim> pressure_geometry;
    CoordTransform<dim, dim> pressure_transform;
    TemplateDOF<dim> pressure_dof;
    BasisFunctionAdmin<double, dim, dim> pressure_basis(pressure_dof);
    auto pressure_template = make_triangle_template(1,
                                                    pressure_geometry,
                                                    pressure_transform,
                                                    pressure_dof,
                                                    pressure_basis);

    auto velocity_space = build_space(mesh, velocity_template);
    auto pressure_space = build_space(mesh, pressure_template);

    auto operators = assemble_lsc_operators(velocity_space, pressure_space);
    auto hp = build_dense_pressure_operator(operators, false);
    auto gp = build_dense_pressure_operator(operators, true);

    for (int i = 0; i < operators.pressure_dofs; ++i)
      hp[static_cast<std::size_t>(i) * operators.pressure_dofs + i] +=
        hp_regularization;

    DenseLU hp_solver(hp, operators.pressure_dofs);

    std::vector<double> src(operators.pressure_dofs, 0.0);
    for (int i = 0; i < operators.pressure_dofs; ++i)
      src[i] = std::sin(0.11 * static_cast<double>(i + 1)) +
               0.2 * std::cos(0.03 * static_cast<double>(i + 1));
    subtract_arithmetic_mean(src);

    auto tmp_1 = hp_solver.solve(src);
    auto tmp_2 = dense_vmult(gp,
                             tmp_1,
                             operators.pressure_dofs,
                             operators.pressure_dofs);
    subtract_arithmetic_mean(tmp_2);
    auto dst = hp_solver.solve(tmp_2);

    return {static_cast<int>(mesh.n_geometry(dim)),
            static_cast<unsigned int>(operators.scalar_velocity_dofs),
            static_cast<unsigned int>(operators.velocity_block_dofs),
            static_cast<unsigned int>(operators.pressure_dofs),
            dense_frobenius_norm(hp),
            dense_frobenius_norm(gp),
            vector_norm(src),
            vector_norm(tmp_1),
            vector_norm(tmp_2),
            vector_norm(dst)};
  }

  void print_result(const LSCResult &result)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Oseen LSC pressure-space vmult demo\n";
    std::cout << "cells: " << result.cells << '\n';
    std::cout << "scalar velocity dofs: "
              << result.scalar_velocity_dofs << '\n';
    std::cout << "velocity block dofs: "
              << result.velocity_block_dofs << '\n';
    std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
    std::cout << "||H_p + gamma I||_F: "
              << result.hp_frobenius_norm << '\n';
    std::cout << "||G_p||_F: " << result.gp_frobenius_norm << '\n';
    std::cout << "H_p regularization gamma: "
              << hp_regularization << '\n';
    std::cout << "input norm: " << result.input_norm << '\n';
    std::cout << "first H_p inverse action norm: "
              << result.first_hp_inverse_norm << '\n';
    std::cout << "G_p action norm: " << result.gp_action_norm << '\n';
    std::cout << "LSC output norm: " << result.lsc_output_norm << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p10";

  try
    {
      print_result(run_lsc_vmult(mesh_file));
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
