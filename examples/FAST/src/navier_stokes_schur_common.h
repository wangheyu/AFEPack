/**
 * @file navier_stokes_schur_common.h
 * @brief Navier-Stokes/Oseen Schur 补算例共享的稀疏矩阵、迭代求解器和块预条件器工具。
 *
 * @details
 * 该头文件为 FAST 下多个不可压流算例提供轻量级公共实现，避免在 PCD、LSC、
 * Schur 补和高雷诺数扫描程序之间重复维护相同的线性代数辅助代码。这里的矩阵
 * 容器只服务于示例程序中的小规模块操作和预条件器原型，真正的有限元装配仍由
 * AFEPack 完成。
 *
 * 主要内容包括：
 * - 基于行映射的稀疏矩阵封装，用于保存显式块矩阵和执行矩阵向量乘。
 * - Dirichlet 约束、齐次约束和对角预条件向量等装配后处理工具。
 * - 对称正定压力质量矩阵的 Jacobi-PCG 求解器。
 * - 以速度 AMG 和压力质量矩阵近似组成的 Schur 质量预条件器。
 * - 右预条件 GMRES 及其 Givens 旋转辅助函数。
 *
 * @note 这些工具面向教学和算例验证，接口保持简单直接；若用于大型生产问题，
 *       应优先考虑成熟稀疏矩阵库和更完整的预条件器实现。
 */

#pragma once

#include <AFEPack/AMGSolver.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace afepack_examples::navier_stokes_schur
{
  /**
   * @brief 记录迭代求解器的收敛状态。
   *
   * @details
   * 所有共享迭代器都返回该结构，调用方可以统一读取是否收敛、实际迭代次数
   * 和最终相对残差。相对残差始终按右端项二范数归一化。
   */
  struct SolveInfo
  {
    /// 求解器是否在给定容差内收敛。
    bool   converged = false;
    /// 实际完成的迭代步数。
    int    iterations = 0;
    /// 退出时的相对残差范数。
    double relative_residual = 0.0;
  };

  /**
   * @brief 面向算例的简单行压缩稀疏矩阵。
   *
   * @details
   * 每一行使用 `std::map<int, double>` 保存非零元，便于在装配后快速访问任意
   * 行、清除约束行列和构造块预条件器。该类没有做复杂的内存优化，目标是让
   * Schur 补算例中的线性代数步骤清晰可读。
   */
  class RowSparseMatrix
  {
  public:
    /**
     * @brief 构造具有指定行数的空矩阵。
     * @param n 矩阵行数和列数；此工具类只用于方阵。
     */
    explicit RowSparseMatrix(const int n = 0)
      : rows(n)
    {}

    /**
     * @brief 重置矩阵规模并清空所有非零元。
     * @param n 新的方阵规模。
     */
    void reinit(const int n)
    {
      rows.assign(n, {});
    }

    /**
     * @brief 返回矩阵规模。
     * @return 方阵的行数。
     */
    int size() const
    {
      return static_cast<int>(rows.size());
    }

    /**
     * @brief 向指定位置累加一个矩阵项。
     * @param row 行号。
     * @param column 列号。
     * @param value 要累加的数值；零值会被忽略以减少显式存储。
     */
    void add(const int row, const int column, const double value)
    {
      if (value != 0.0)
        rows[row][column] += value;
    }

    /**
     * @brief 清除指定自由度对应的行和列，并把对角线设为 1。
     * @param index 需要施加强约束的自由度编号。
     */
    void clear_row_and_column(const int index)
    {
      for (auto &row : rows)
        row.erase(index);
      rows[index].clear();
      rows[index][index] = 1.0;
    }

    /**
     * @brief 读取矩阵项。
     * @param row 行号。
     * @param column 列号。
     * @return 对应矩阵项；未显式存储时返回 0。
     */
    double entry(const int row, const int column) const
    {
      const auto position = rows[row].find(column);
      return position == rows[row].end() ? 0.0 : position->second;
    }

    /**
     * @brief 返回指定行的非零项集合。
     * @param row 行号。
     * @return 行内列号到数值的映射。
     */
    const std::map<int, double> &row_entries(const int row) const
    {
      return rows[row];
    }

    /**
     * @brief 执行矩阵向量乘。
     * @param source 输入向量。
     * @return `A * source` 的结果向量。
     */
    std::vector<double> vmult(const std::vector<double> &source) const
    {
      std::vector<double> destination(rows.size(), 0.0);
      for (int row = 0; row < size(); ++row)
        for (const auto &[column, value] : rows[row])
          destination[row] += value * source[column];
      return destination;
    }

    /**
     * @brief 读取指定行的对角项。
     * @param index 自由度编号。
     * @return 对角线数值；缺失时返回 0。
     */
    double diagonal(const int index) const
    {
      return entry(index, index);
    }

    /**
     * @brief 统计显式保存的非零元数量。
     * @return 所有行映射中的条目总数。
     */
    unsigned int nonzeros() const
    {
      unsigned int count = 0;
      for (const auto &row : rows)
        count += static_cast<unsigned int>(row.size());
      return count;
    }

  private:
    std::vector<std::map<int, double>> rows;
  };

  /**
   * @brief 计算两个实向量的欧氏内积。
   * @param left 左向量。
   * @param right 右向量。
   * @return `left` 与 `right` 的点积。
   */
  inline double dot(const std::vector<double> &left,
                    const std::vector<double> &right)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
      result += left[i] * right[i];
    return result;
  }

  /**
   * @brief 计算实向量的二范数。
   * @param vector 输入向量。
   * @return `sqrt(dot(vector, vector))`。
   */
  inline double norm(const std::vector<double> &vector)
  {
    return std::sqrt(dot(vector, vector));
  }

  /**
   * @brief 执行向量加法 `destination += coefficient * source`。
   * @param destination 原地更新的目标向量。
   * @param coefficient 缩放系数。
   * @param source 被缩放并累加的源向量。
   */
  inline void add_scaled(std::vector<double> &destination,
                         const double coefficient,
                         const std::vector<double> &source)
  {
    for (std::size_t i = 0; i < destination.size(); ++i)
      destination[i] += coefficient * source[i];
  }

  /**
   * @brief 对线性系统施加强 Dirichlet 约束。
   * @param matrix 系统矩阵；函数会清除对应行列并设置单位对角。
   * @param rhs 右端项；函数会扣除被消去列对其它方程的贡献。
   * @param index 被约束自由度编号。
   * @param value 指定的边界值。
   */
  inline void apply_dirichlet_constraint(RowSparseMatrix &matrix,
                                         std::vector<double> &rhs,
                                         const int index,
                                         const double value)
  {
    for (int row = 0; row < matrix.size(); ++row)
      if (row != index)
        rhs[row] -= matrix.entry(row, index) * value;

    matrix.clear_row_and_column(index);
    rhs[index] = value;
  }

  /**
   * @brief 对矩阵施加齐次约束。
   * @param matrix 系统矩阵；函数会清除对应行列并设置单位对角。
   * @param index 被约束自由度编号。
   */
  inline void apply_homogeneous_constraint(RowSparseMatrix &matrix,
                                           const int index)
  {
    matrix.clear_row_and_column(index);
  }

  /**
   * @brief 构造 Jacobi 预条件器所需的逆对角向量。
   * @param matrix 输入矩阵。
   * @return 每个自由度对应的逆对角项；近零对角项使用 1 作为保护值。
   */
  inline std::vector<double>
  inverse_diagonal(const RowSparseMatrix &matrix)
  {
    std::vector<double> inverse(matrix.size(), 1.0);
    for (int i = 0; i < matrix.size(); ++i)
      {
        const double diagonal = matrix.diagonal(i);
        inverse[i] =
          std::abs(diagonal) > 1.0e-14 ? 1.0 / diagonal : 1.0;
      }
    return inverse;
  }

  /**
   * @brief 使用 Jacobi 预条件共轭梯度法求解对称正定系统。
   * @param matrix 系统矩阵。
   * @param solution 输入初值并在返回时保存近似解。
   * @param rhs 右端项。
   * @param inverse_diagonal Jacobi 预条件器的逆对角向量。
   * @param tolerance 相对残差收敛阈值。
   * @param max_iterations 最大迭代步数。
   * @return 求解器收敛信息。
   */
  inline SolveInfo solve_pcg(const RowSparseMatrix &matrix,
                             std::vector<double> &solution,
                             const std::vector<double> &rhs,
                             const std::vector<double> &inverse_diagonal,
                             const double tolerance,
                             const int max_iterations)
  {
    const double rhs_norm = std::max(norm(rhs), 1.0e-30);
    auto matrix_times_solution = matrix.vmult(solution);
    std::vector<double> residual = rhs;
    add_scaled(residual, -1.0, matrix_times_solution);

    double relative_residual = norm(residual) / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    std::vector<double> preconditioned(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i)
      preconditioned[i] = inverse_diagonal[i] * residual[i];
    std::vector<double> direction = preconditioned;
    double residual_product = dot(residual, preconditioned);

    for (int iteration = 0; iteration < max_iterations; ++iteration)
      {
        const auto matrix_times_direction = matrix.vmult(direction);
        const double denominator = dot(direction, matrix_times_direction);
        if (std::abs(denominator) <= 1.0e-30)
          break;

        const double alpha = residual_product / denominator;
        add_scaled(solution, alpha, direction);
        add_scaled(residual, -alpha, matrix_times_direction);

        relative_residual = norm(residual) / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, iteration + 1, relative_residual};

        for (std::size_t i = 0; i < rhs.size(); ++i)
          preconditioned[i] = inverse_diagonal[i] * residual[i];

        const double new_residual_product = dot(residual, preconditioned);
        if (std::abs(residual_product) <= 1.0e-30)
          break;

        const double beta = new_residual_product / residual_product;
        for (std::size_t i = 0; i < direction.size(); ++i)
          direction[i] = preconditioned[i] + beta * direction[i];
        residual_product = new_residual_product;
      }

    return {false, max_iterations, relative_residual};
  }

  /**
   * @brief Schur 补算例中使用的块质量预条件器。
   *
   * @details
   * 预条件器把速度块近似为两个分量上的 AMG 求解，并用压力质量矩阵的 PCG
   * 求解近似压力 Schur 补。该类实现 `vmult()`，因此可作为右预条件 GMRES 的
   * 模板参数使用。
   */
  class SchurMassPreconditioner
  {
  public:
    /**
     * @brief 构造 Schur 质量预条件器。
     * @param system_matrix 完整速度-压力块系统矩阵。
     * @param pressure_mass_matrix 压力质量矩阵，用作压力 Schur 补近似。
     * @param velocity_amg 速度块 AMG 预条件器。
     * @param viscosity 黏性系数，用于缩放速度块近似。
     * @param n_velocity_dofs 单个速度分量的自由度数。
     * @param n_pressure_dofs 压力自由度数。
     * @param pressure_pin 被固定的压力自由度，用于处理压力常数核空间。
     * @param context 错误信息上下文，便于定位调用该预条件器的算例。
     */
    SchurMassPreconditioner(const RowSparseMatrix &system_matrix,
                            const RowSparseMatrix &pressure_mass_matrix,
                            AMGPreconditioner &velocity_amg,
                            const double viscosity,
                            const int n_velocity_dofs,
                            const int n_pressure_dofs,
                            const int pressure_pin,
                            std::string context)
      : system_matrix(system_matrix)
      , pressure_mass_matrix(pressure_mass_matrix)
      , velocity_amg(velocity_amg)
      , velocity_inverse_scale(1.0 / viscosity)
      , pressure_inverse_diagonal(inverse_diagonal(pressure_mass_matrix))
      , n_velocity_dofs(n_velocity_dofs)
      , n_pressure_dofs(n_pressure_dofs)
      , pressure_pin(pressure_pin)
      , context(std::move(context))
    {}

    /**
     * @brief 应用块预条件器。
     * @param source 输入残差向量，排列方式为两个速度分量后接压力自由度。
     * @return 预条件后的向量。
     *
     * @throws std::runtime_error 当内部压力质量矩阵 PCG 未收敛时抛出。
     */
    std::vector<double> vmult(const std::vector<double> &source) const
    {
      std::vector<double> destination(source.size(), 0.0);
      std::vector<double> velocity_solution(2 * n_velocity_dofs, 0.0);

      for (int component = 0; component < 2; ++component)
        {
          Vector<double> amg_source(n_velocity_dofs);
          Vector<double> amg_destination(n_velocity_dofs);
          for (int i = 0; i < n_velocity_dofs; ++i)
            amg_source(i) = source[velocity_index(component, i)];
          amg_destination = 0.0;
          velocity_amg.vmult(amg_destination, amg_source);
          for (int i = 0; i < n_velocity_dofs; ++i)
            {
              const double value =
                velocity_inverse_scale * amg_destination(i);
              velocity_solution[velocity_index(component, i)] = value;
              destination[velocity_index(component, i)] = value;
            }
        }

      std::vector<double> pressure_rhs(n_pressure_dofs, 0.0);
      for (int i = 0; i < n_pressure_dofs; ++i)
        {
          double divergence_action = 0.0;
          const int row = pressure_index(i);
          for (const auto &[column, value] : system_matrix.row_entries(row))
            if (column < 2 * n_velocity_dofs)
              divergence_action += value * velocity_solution[column];
          pressure_rhs[i] = source[row] - divergence_action;
        }

      std::vector<double> pressure_solution(n_pressure_dofs, 0.0);
      const SolveInfo pressure_info =
        solve_pcg(pressure_mass_matrix,
                  pressure_solution,
                  pressure_rhs,
                  pressure_inverse_diagonal,
                  1.0e-10,
                  160);
      if (!pressure_info.converged)
        throw std::runtime_error("pressure mass PCG failed inside " +
                                 context + " preconditioner");

      for (int i = 0; i < n_pressure_dofs; ++i)
        destination[pressure_index(i)] =
          i == pressure_pin ? pressure_solution[i] : -pressure_solution[i];

      return destination;
    }

  private:
    /**
     * @brief 将速度分量和局部自由度映射到全局块向量下标。
     * @param component 速度分量编号。
     * @param dof 单分量速度自由度编号。
     * @return 块向量中的全局下标。
     */
    int velocity_index(const int component, const int dof) const
    {
      return component * n_velocity_dofs + dof;
    }

    /**
     * @brief 将压力自由度映射到全局块向量下标。
     * @param dof 压力自由度编号。
     * @return 块向量中的全局下标。
     */
    int pressure_index(const int dof) const
    {
      return 2 * n_velocity_dofs + dof;
    }

    const RowSparseMatrix &system_matrix;
    const RowSparseMatrix &pressure_mass_matrix;
    AMGPreconditioner &velocity_amg;
    double velocity_inverse_scale = 1.0;
    std::vector<double> pressure_inverse_diagonal;
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
    int pressure_pin = 0;
    std::string context;
  };

  /**
   * @brief 对两个标量应用 Givens 平面旋转。
   * @param first 第一个标量，原地更新。
   * @param second 第二个标量，原地更新。
   * @param cosine 旋转余弦。
   * @param sine 旋转正弦。
   */
  inline void apply_plane_rotation(double &first,
                                   double &second,
                                   const double cosine,
                                   const double sine)
  {
    const double temporary = cosine * first + sine * second;
    second = -sine * first + cosine * second;
    first = temporary;
  }

  /**
   * @brief 生成用于消去 Hessenberg 矩阵次对角元的 Givens 旋转。
   * @param first 当前对角元。
   * @param second 当前次对角元。
   * @param cosine 输出旋转余弦。
   * @param sine 输出旋转正弦。
   */
  inline void generate_plane_rotation(const double first,
                                      const double second,
                                      double &cosine,
                                      double &sine)
  {
    if (second == 0.0)
      {
        cosine = 1.0;
        sine = 0.0;
      }
    else if (std::abs(second) > std::abs(first))
      {
        const double ratio = first / second;
        sine = 1.0 / std::sqrt(1.0 + ratio * ratio);
        cosine = ratio * sine;
      }
    else
      {
        const double ratio = second / first;
        cosine = 1.0 / std::sqrt(1.0 + ratio * ratio);
        sine = ratio * cosine;
      }
  }

  /**
   * @brief 回代求解上 Hessenberg 最小二乘子问题。
   * @param hessenberg GMRES Arnoldi 过程中形成并经 Givens 旋转后的矩阵。
   * @param rotated_rhs 同步旋转后的右端项。
   * @param size 当前 Krylov 子空间维数。
   * @return 子空间基向量的线性组合系数。
   */
  inline std::vector<double> solve_upper_hessenberg(
    const std::vector<std::vector<double>> &hessenberg,
    const std::vector<double> &rotated_rhs,
    const int size)
  {
    std::vector<double> coefficients(size, 0.0);
    for (int i = size - 1; i >= 0; --i)
      {
        double sum = rotated_rhs[i];
        for (int j = i + 1; j < size; ++j)
          sum -= hessenberg[i][j] * coefficients[j];
        coefficients[i] = sum / hessenberg[i][i];
      }
    return coefficients;
  }

  /**
   * @brief 使用右预条件 restarted GMRES 求解非对称线性系统。
   *
   * @tparam Preconditioner 提供 `vmult(const std::vector<double>&)` 的预条件器类型。
   * @param matrix 系统矩阵。
   * @param solution 输入初值并在返回时保存近似解。
   * @param rhs 右端项。
   * @param preconditioner 右预条件器。
   * @param tolerance 相对残差收敛阈值。
   * @param restart GMRES 重启长度。
   * @param max_iterations 最大总迭代步数。
   * @return 求解器收敛信息。
   *
   * @details
   * 该实现显式保存 Arnoldi 基、预条件基和上 Hessenberg 矩阵，便于在算例中
   * 打印和调试 Schur 补预条件器的收敛行为。返回的 `iterations` 为累计内迭代
   * 步数，不包含外层非线性迭代。
   */
  template <typename Preconditioner>
  SolveInfo solve_right_preconditioned_gmres(
    const RowSparseMatrix &matrix,
    std::vector<double> &solution,
    const std::vector<double> &rhs,
    const Preconditioner &preconditioner,
    const double tolerance,
    const int restart,
    const int max_iterations)
  {
    const int n = matrix.size();
    const double rhs_norm = std::max(norm(rhs), 1.0e-30);
    auto matrix_times_solution = matrix.vmult(solution);
    std::vector<double> residual = rhs;
    add_scaled(residual, -1.0, matrix_times_solution);

    double residual_norm = norm(residual);
    double relative_residual = residual_norm / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    int total_iterations = 0;
    while (total_iterations < max_iterations)
      {
        std::vector<std::vector<double>> basis(
          restart + 1, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> preconditioned_basis(
          restart, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> hessenberg(
          restart + 1, std::vector<double>(restart, 0.0));
        std::vector<double> cosines(restart, 0.0);
        std::vector<double> sines(restart, 0.0);
        std::vector<double> rotated_rhs(restart + 1, 0.0);

        basis[0] = residual;
        for (double &entry : basis[0])
          entry /= residual_norm;
        rotated_rhs[0] = residual_norm;

        int inner_iterations = 0;
        for (; inner_iterations < restart &&
               total_iterations < max_iterations;
             ++inner_iterations, ++total_iterations)
          {
            preconditioned_basis[inner_iterations] =
              preconditioner.vmult(basis[inner_iterations]);
            basis[inner_iterations + 1] =
              matrix.vmult(preconditioned_basis[inner_iterations]);

            for (int i = 0; i <= inner_iterations; ++i)
              {
                hessenberg[i][inner_iterations] =
                  dot(basis[inner_iterations + 1], basis[i]);
                add_scaled(basis[inner_iterations + 1],
                           -hessenberg[i][inner_iterations],
                           basis[i]);
              }

            hessenberg[inner_iterations + 1][inner_iterations] =
              norm(basis[inner_iterations + 1]);
            if (hessenberg[inner_iterations + 1][inner_iterations] != 0.0)
              for (double &entry : basis[inner_iterations + 1])
                entry /= hessenberg[inner_iterations + 1][inner_iterations];

            for (int i = 0; i < inner_iterations; ++i)
              apply_plane_rotation(hessenberg[i][inner_iterations],
                                   hessenberg[i + 1][inner_iterations],
                                   cosines[i],
                                   sines[i]);

            generate_plane_rotation(
              hessenberg[inner_iterations][inner_iterations],
              hessenberg[inner_iterations + 1][inner_iterations],
              cosines[inner_iterations],
              sines[inner_iterations]);
            apply_plane_rotation(
              hessenberg[inner_iterations][inner_iterations],
              hessenberg[inner_iterations + 1][inner_iterations],
              cosines[inner_iterations],
              sines[inner_iterations]);
            apply_plane_rotation(rotated_rhs[inner_iterations],
                                 rotated_rhs[inner_iterations + 1],
                                 cosines[inner_iterations],
                                 sines[inner_iterations]);

            relative_residual =
              std::abs(rotated_rhs[inner_iterations + 1]) / rhs_norm;
            if (relative_residual <= tolerance)
              {
                const int system_size = inner_iterations + 1;
                const auto coefficients = solve_upper_hessenberg(
                  hessenberg, rotated_rhs, system_size);
                for (int i = 0; i < system_size; ++i)
                  add_scaled(solution,
                             coefficients[i],
                             preconditioned_basis[i]);
                return {true, total_iterations + 1, relative_residual};
              }
          }

        const auto coefficients = solve_upper_hessenberg(
          hessenberg, rotated_rhs, inner_iterations);
        for (int i = 0; i < inner_iterations; ++i)
          add_scaled(solution,
                     coefficients[i],
                     preconditioned_basis[i]);

        matrix_times_solution = matrix.vmult(solution);
        residual = rhs;
        add_scaled(residual, -1.0, matrix_times_solution);
        residual_norm = norm(residual);
        relative_residual = residual_norm / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }
}
