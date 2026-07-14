/**
 * @file oseen_pcd_vmult_afepack.cpp
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
 * - PCD 预条件器：构造压力对流扩散型 Schur 补近似，评估 PCD 预条件的迭代效率。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/BilinearOperator.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>

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
  constexpr double viscosity = 1.0e-2;
  constexpr double pressure_laplace_regularization = 1.0e-6;

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct PCDResult
  {
    unsigned int pressure_dofs = 0;
    std::size_t  mass_nonzeros = 0;
    std::size_t  laplace_nonzeros = 0;
    std::size_t  convdiff_nonzeros = 0;
    double       mass_frobenius_norm = 0.0;
    double       laplace_frobenius_norm = 0.0;
    double       convdiff_frobenius_norm = 0.0;
    double       input_norm = 0.0;
    double       mass_inverse_norm = 0.0;
    double       convdiff_action_norm = 0.0;
    double       pcd_output_norm = 0.0;
    SolveInfo    mass_solve;
    SolveInfo    laplace_solve;
  };

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
  make_triangle_template(TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle.1.tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle.1.bas_fun");

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

  std::vector<double> frozen_velocity(const Point<dim> &p)
  {
    return {2.0 * p[1] * (1.0 - p[0] * p[0]),
            -2.0 * p[0] * (1.0 - p[1] * p[1])};
  }

  class PressureMassMatrix : public BilinearOperator<dim, double>
  {
  public:
    explicit PressureMassMatrix(FEMSpace<double, dim> &space)
      : BilinearOperator<dim, double>(space, space)
    {}

    void getElementMatrix(
      const Element<double, dim> &element,
      const Element<double, dim> &,
      const ActiveElementPairIterator<dim>::State =
        ActiveElementPairIterator<dim>::EQUAL) override
    {
      const double volume = element.templateElement().volume();
      const QuadratureInfo<dim> &quad_info =
        element.findQuadratureInfo(algebricAccuracy());
      const auto jacobian =
        element.local_to_global_jacobian(quad_info.quadraturePoint());
      const auto q_point =
        element.local_to_global(quad_info.quadraturePoint());
      const auto basis_value = element.basis_function_value(q_point);

      const int n_q = quad_info.n_quadraturePoint();
      const int n_dof = element.dof().size();

      for (int q = 0; q < n_q; ++q)
        {
          const double jxw = quad_info.weight(q) * jacobian[q] * volume;
          for (int i = 0; i < n_dof; ++i)
            for (int j = 0; j < n_dof; ++j)
              elementMatrix(i, j) +=
                jxw * basis_value[j][q] * basis_value[i][q];
        }
    }
  };

  class PressureLaplaceMatrix : public BilinearOperator<dim, double>
  {
  public:
    explicit PressureLaplaceMatrix(FEMSpace<double, dim> &space)
      : BilinearOperator<dim, double>(space, space)
    {}

    void getElementMatrix(
      const Element<double, dim> &element,
      const Element<double, dim> &,
      const ActiveElementPairIterator<dim>::State =
        ActiveElementPairIterator<dim>::EQUAL) override
    {
      const double volume = element.templateElement().volume();
      const QuadratureInfo<dim> &quad_info =
        element.findQuadratureInfo(algebricAccuracy());
      const auto jacobian =
        element.local_to_global_jacobian(quad_info.quadraturePoint());
      const auto q_point =
        element.local_to_global(quad_info.quadraturePoint());
      const auto basis_value = element.basis_function_value(q_point);
      const auto basis_grad = element.basis_function_gradient(q_point);

      const int n_q = quad_info.n_quadraturePoint();
      const int n_dof = element.dof().size();

      for (int q = 0; q < n_q; ++q)
        {
          const double jxw = quad_info.weight(q) * jacobian[q] * volume;
          for (int i = 0; i < n_dof; ++i)
            for (int j = 0; j < n_dof; ++j)
              {
                const double laplace =
                  innerProduct(basis_grad[i][q], basis_grad[j][q]);
                const double mass =
                  basis_value[j][q] * basis_value[i][q];
                elementMatrix(i, j) +=
                  jxw * (laplace +
                         pressure_laplace_regularization * mass);
              }
        }
    }
  };

  class PressureConvectionDiffusionMatrix
    : public BilinearOperator<dim, double>
  {
  public:
    explicit PressureConvectionDiffusionMatrix(FEMSpace<double, dim> &space)
      : BilinearOperator<dim, double>(space, space)
    {}

    void getElementMatrix(
      const Element<double, dim> &element,
      const Element<double, dim> &,
      const ActiveElementPairIterator<dim>::State =
        ActiveElementPairIterator<dim>::EQUAL) override
    {
      const double volume = element.templateElement().volume();
      const QuadratureInfo<dim> &quad_info =
        element.findQuadratureInfo(algebricAccuracy());
      const auto jacobian =
        element.local_to_global_jacobian(quad_info.quadraturePoint());
      const auto q_point =
        element.local_to_global(quad_info.quadraturePoint());
      const auto basis_value = element.basis_function_value(q_point);
      const auto basis_grad = element.basis_function_gradient(q_point);

      const int n_q = quad_info.n_quadraturePoint();
      const int n_dof = element.dof().size();

      for (int q = 0; q < n_q; ++q)
        {
          const double jxw = quad_info.weight(q) * jacobian[q] * volume;
          const auto wind = frozen_velocity(q_point[q]);
          for (int i = 0; i < n_dof; ++i)
            for (int j = 0; j < n_dof; ++j)
              {
                const double laplace =
                  innerProduct(basis_grad[i][q], basis_grad[j][q]);
                const double convection =
                  (wind[0] * basis_grad[j][q][0] +
                   wind[1] * basis_grad[j][q][1]) *
                  basis_value[i][q];
                elementMatrix(i, j) +=
                  jxw * (viscosity * laplace + convection);
              }
        }
    }
  };

  double sparse_frobenius_norm(const SparseMatrix<double> &matrix)
  {
    double sum = 0.0;
    for (std::size_t k = 0; k < matrix.n_nonzero_elements(); ++k)
      {
        const double entry = matrix.global_entry(k);
        sum += entry * entry;
      }
    return std::sqrt(sum);
  }

  void subtract_arithmetic_mean(Vector<double> &vector)
  {
    double sum = 0.0;
    for (std::size_t i = 0; i < vector.size(); ++i)
      sum += vector(i);

    const double mean = sum / static_cast<double>(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i)
      vector(i) -= mean;
  }

  void jacobi_apply(const SparseMatrix<double> &matrix,
                    const Vector<double> &src,
                    Vector<double> &dst)
  {
    dst.reinit(src.size(), false);
    for (std::size_t i = 0; i < src.size(); ++i)
      {
        const double diagonal = matrix.diag_element(i);
        dst(i) = std::abs(diagonal) > 1.0e-14 ? src(i) / diagonal : src(i);
      }
  }

  SolveInfo solve_pcg(const SparseMatrix<double> &matrix,
                      Vector<double> &x,
                      const Vector<double> &rhs,
                      const double tolerance,
                      const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> r = rhs;
    r.add(-1.0, ax);

    double relative_residual = r.l2_norm() / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    Vector<double> z(rhs.size());
    jacobi_apply(matrix, r, z);
    Vector<double> p = z;
    Vector<double> ap(rhs.size());

    double rz_old = r.dot(z);
    if (std::abs(rz_old) < 1.0e-30)
      return {false, 0, relative_residual};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        matrix.vmult(ap, p);
        const double denominator = p.dot(ap);
        if (std::abs(denominator) < 1.0e-30)
          return {false, iteration - 1, relative_residual};

        const double alpha = rz_old / denominator;
        x.add(alpha, p);
        r.add(-alpha, ap);

        relative_residual = r.l2_norm() / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, iteration, relative_residual};

        jacobi_apply(matrix, r, z);
        const double rz_new = r.dot(z);
        if (std::abs(rz_old) < 1.0e-30)
          return {false, iteration, relative_residual};

        const double beta = rz_new / rz_old;
        p.scale(beta);
        p.add(1.0, z);
        rz_old = rz_new;
      }

    return {false, max_iterations, relative_residual};
  }

  PCDResult run_pcd_vmult(const std::string &mesh_file)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);

    auto pressure_space = build_space(mesh, template_element);

    PressureMassMatrix mass_matrix(pressure_space);
    mass_matrix.algebricAccuracy() = 4;
    mass_matrix.build();

    PressureLaplaceMatrix laplace_matrix(pressure_space);
    laplace_matrix.algebricAccuracy() = 4;
    laplace_matrix.build();

    PressureConvectionDiffusionMatrix convdiff_matrix(pressure_space);
    convdiff_matrix.algebricAccuracy() = 4;
    convdiff_matrix.build();

    Vector<double> src(pressure_space.n_dof());
    for (std::size_t i = 0; i < src.size(); ++i)
      src(i) = std::sin(0.13 * static_cast<double>(i + 1)) +
               0.25 * std::cos(0.07 * static_cast<double>(i + 1));
    subtract_arithmetic_mean(src);

    Vector<double> mass_inverse(src.size());
    Vector<double> convdiff_action(src.size());
    Vector<double> pcd_output(src.size());

    const SolveInfo mass_solve =
      solve_pcg(mass_matrix, mass_inverse, src, 1.0e-10, 2000);
    if (!mass_solve.converged)
      throw std::runtime_error("PCG failed for pressure mass solve");

    convdiff_matrix.vmult(convdiff_action, mass_inverse);
    subtract_arithmetic_mean(convdiff_action);

    const SolveInfo laplace_solve =
      solve_pcg(laplace_matrix,
                pcd_output,
                convdiff_action,
                1.0e-10,
                4000);
    if (!laplace_solve.converged)
      throw std::runtime_error("PCG failed for pressure Laplace solve");

    return {pressure_space.n_dof(),
            mass_matrix.n_nonzero_elements(),
            laplace_matrix.n_nonzero_elements(),
            convdiff_matrix.n_nonzero_elements(),
            sparse_frobenius_norm(mass_matrix),
            sparse_frobenius_norm(laplace_matrix),
            sparse_frobenius_norm(convdiff_matrix),
            src.l2_norm(),
            mass_inverse.l2_norm(),
            convdiff_action.l2_norm(),
            pcd_output.l2_norm(),
            mass_solve,
            laplace_solve};
  }

  void print_result(const PCDResult &result)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Oseen PCD pressure-space vmult demo\n";
    std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
    std::cout << "M_p nonzeros: " << result.mass_nonzeros << '\n';
    std::cout << "K_p nonzeros: " << result.laplace_nonzeros << '\n';
    std::cout << "F_p nonzeros: " << result.convdiff_nonzeros << '\n';
    std::cout << "||M_p||_F: " << result.mass_frobenius_norm << '\n';
    std::cout << "||K_p + gamma M_p||_F: "
              << result.laplace_frobenius_norm << '\n';
    std::cout << "||F_p||_F: " << result.convdiff_frobenius_norm << '\n';
    std::cout << "viscosity: " << viscosity << '\n';
    std::cout << "pressure Laplace regularization gamma: "
              << pressure_laplace_regularization << '\n';
    std::cout << "input norm: " << result.input_norm << '\n';
    std::cout << "M_p inverse action norm: "
              << result.mass_inverse_norm << '\n';
    std::cout << "F_p action norm: "
              << result.convdiff_action_norm << '\n';
    std::cout << "PCD output norm: " << result.pcd_output_norm << '\n';
    std::cout << "M_p PCG iterations: "
              << result.mass_solve.iterations
              << ", relres: " << result.mass_solve.relative_residual
              << '\n';
    std::cout << "K_p PCG iterations: "
              << result.laplace_solve.iterations
              << ", relres: " << result.laplace_solve.relative_residual
              << '\n';
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
      print_result(run_pcd_vmult(mesh_file));
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
