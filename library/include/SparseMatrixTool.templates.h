/**
 * SparseMatrixTool.templates.h
 * 
 */

#include "SparseMatrixTool.h"

namespace SparseMatrixTool {

  /**
   * 将两个稀疏矩阵按照水平方式连接成为一个新的稀疏矩阵。
   * 
   * @param m0 左边的矩阵
   * @param m1 右边的矩阵
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */
  template <typename number>
    void hCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok)
    {
      const SparsityPattern& sp0 = m0.get_sparsity_pattern();
      const SparsityPattern& sp1 = m1.get_sparsity_pattern();
      if (! is_pattern_ok)
	hCatSparsityPattern(sp0, sp1, sp);
      m.reinit(sp);
      int i, j, n = sp0.n_cols();
      const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
      const u_int * p_column0 = sp0.get_column_numbers();
      const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
      const u_int * p_column1 = sp1.get_column_numbers();
      for (i = 0;i < sp.n_rows();i ++) {
	for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	  m.add(i, p_column0[j], m0.global_entry(j));
	for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	  m.add(i, n + p_column1[j], m1.global_entry(j));
      }
    };

  /**
   * 将两个稀疏矩阵按照竖直方式连接成为一个新的稀疏矩阵。
   * 
   * @param m0 上边的矩阵
   * @param m1 下边的矩阵
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */
  template <typename number>
    void vCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok)
    {
      const SparsityPattern& sp0 = m0.get_sparsity_pattern();
      const SparsityPattern& sp1 = m1.get_sparsity_pattern();
      if (! is_pattern_ok)
	vCatSparsityPattern(sp0, sp1, sp);
      m.reinit(sp);
      int i, j, n = sp0.n_rows();
      const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
      const u_int * p_column0 = sp0.get_column_numbers();
      const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
      const u_int * p_column1 = sp1.get_column_numbers();
      for (i = 0;i < n;i ++)
	for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	  m.add(i, p_column0[j], m0.global_entry(j));
      for (i = 0;i < sp1.n_rows();i ++)
	for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	  m.add(i + n, p_column1[j], m1.global_entry(j));
    };

  /**
   * 将两个稀疏矩阵按鞍点问题的方式连接成为一个新的稀疏
   * 矩阵，也就是对于给定的稀疏矩阵 \f$ A \f$ 和 \f$ B
   * \f$，我们得到矩阵
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A & B \\
   *      B^T & 0
   *    \end{array}\right)
   * \f]
   *
   * @param m0 矩阵 \f$ A \f$
   * @param m1 矩阵 \f$ B \f$
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */
  template <typename number>
    void gammaCatSparseMatrix(const SparseMatrix<number>& m0,
			      const SparseMatrix<number>& m1,
			      SparsityPattern& sp,
			      SparseMatrix<number>& m,
			      bool is_pattern_ok)
    {
      const SparsityPattern& sp0 = m0.get_sparsity_pattern();
      const SparsityPattern& sp1 = m1.get_sparsity_pattern();
      if (! is_pattern_ok)
	gammaCatSparsityPattern(sp0, sp1, sp);
      m.reinit(sp);

      int i, j, n = sp0.n_rows();
      const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
      const u_int * p_column0 = sp0.get_column_numbers();
      const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
      const u_int * p_column1 = sp1.get_column_numbers();
      for (i = 0;i < n;i ++) {
	for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	  m.add(i, p_column0[j], m0.global_entry(j));
	for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++) {
	  m.add(i, n + p_column1[j], m1.global_entry(j));
	  m.add(n + p_column1[j], i, m1.global_entry(j));
	}
      }
    };

  /**
   * 将两个稀疏矩阵按照对角方式连接成为一个新的稀疏矩阵。
   * 
   * @param m0 左上角矩阵
   * @param m1 右下角矩阵
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则不对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */
  template <typename number>
    void dCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok)
    {
      const SparsityPattern& sp0 = m0.get_sparsity_pattern();
      const SparsityPattern& sp1 = m1.get_sparsity_pattern();
      if (! is_pattern_ok)
	dCatSparsityPattern(sp0, sp1, sp);
      m.reinit(sp);

      int i, j;
      int r0 = sp0.n_rows(), r1 = sp1.n_rows();
      int c0 = sp0.n_cols();
      const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
      const u_int * p_column0 = sp0.get_column_numbers();
      const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
      const u_int * p_column1 = sp1.get_column_numbers();
      for (i = 0;i < r0;i ++)
	for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	  m.add(i, p_column0[j], m0.global_entry(j));
      for (i = 0;i < r1;i ++)
	for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	  m.add(r0 + i, c0 + p_column1[j], m1.global_entry(j));
    };

  /** 
   * 将四个稀疏矩阵完全连接为一个新的稀疏矩阵
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A_{00} & A_{01} \\
   *      A_{01} & A_{11}
   *    \end{array}\right)
   * \f]
   * 
   * @param m00 左上角矩阵
   * @param m01 右上角矩阵
   * @param m10 左下角矩阵
   * @param m11 右下角矩阵
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */  
  template <typename number>
    void fullCatSparseMatrix(const SparseMatrix<number>& m00,
			     const SparseMatrix<number>& m01,
			     const SparseMatrix<number>& m10,
			     const SparseMatrix<number>& m11,
                 SparsityPattern& sp,
			     SparseMatrix<number>& m,
			     bool is_pattern_ok)
    {
      SparsityPattern sp0, sp1;
      SparseMatrix<number> m0, m1;
      hCatSparseMatrix(m00, m01, sp0, m0, false);
      hCatSparseMatrix(m10, m11, sp1, m1, false);
      vCatSparseMatrix(m0, m1, sp, m, is_pattern_ok);
    };

};

/**
 * end of file
 * 
 */
