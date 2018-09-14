/**
 * SparseMatrixTool.cpp
 * 
 */

#include "SparseMatrixTool.templates.h"

namespace SparseMatrixTool {

  /**
   * 将两个稀疏矩阵模板按水平方向连接成为一个新的模板
   * 要求这两个矩阵具有相同的行数。
   * 
   * @param sp0 左边的模板
   * @param sp1 右边的模板
   * @param sp 得到的新模板
   */
  void hCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows();
    u_int n0 = sp0.n_cols();

    /**
     * 检查这个两个模板是否都已经压缩过了。
     * 
     */
    Assert(sp0.is_compressed(), SparsityPattern::ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * 检查这两个模板是否具有相同的行数。
     * 
     */
    Assert (m0 == sp1.n_rows(), ExcDimensionDontMatch(m0, sp1.n_rows()));
  
    int i, j;
    u_int n1 = sp1.n_cols();
    std::vector<u_int> row_length(m0, 0);
    for (i = 0;i < m0;i ++)
      row_length[i] = sp0.row_length(i) + sp1.row_length(i);
    sp.reinit(m0, n0 + n1, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++) {
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(i, n0 + p_column1[j]);
    }
    sp.compress();
  };

  /**
   * 将两个稀疏矩阵模板按竖直方向连接成为一个新的模板
   * 要求这两个矩阵具有相同的列数。
   *
   * @param sp0 上边的模板
   * @param sp1 下边的模板
   * @param sp 得到的新模板
   */
  void vCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows();
    u_int n0 = sp0.n_cols();

    /**
     * 检查这个两个模板是否都已经压缩过了。
     * 
     */
    Assert(sp0.is_compressed(), ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * 检查这两个模板是否具有相同的列数。
     * 
     */
    Assert (n0 == sp1.n_cols(), ExcDimensionDontMatch(m0, sp1.n_cols()));
  
    u_int i, j, m1 = sp1.n_rows();
    std::vector<u_int> row_length(m0 + m1, 0);
    for (i = 0;i < m0;i ++)
      row_length[i] = sp0.row_length(i);
    for (i = 0;i < m1;i ++)
      row_length[i + m0] = sp1.row_length(i);
    sp.reinit(m0 + m1, n0, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++)
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
    for (i = 0;i < m1;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(m0 + i, p_column1[j]);
    sp.compress();
  };

  /**
   * 将两个稀疏矩阵模板按鞍点问题的方式连接成为一个新的
   * 模板，也就是对于给定的稀疏矩阵 \f$ A \f$ 和 \f$ B
   * \f$，我们得到矩阵
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A & B \\
   *      B^T & 0
   *    \end{array}\right)
   * \f]
   *
   * 的模板。要求 \f$ A \f$ 和 \f$ B \f$ 具有相同的行数。
   * 
   * @param sp0 矩阵 \f$ A \f$ 的模板
   * @param sp1 矩阵 \f$ B \f$ 的模板
   * @param sp 得到的新模板
   */
  void gammaCatSparsityPattern(const SparsityPattern& sp0,
			       const SparsityPattern& sp1,
			       SparsityPattern& sp)
  {
    u_int m = sp0.n_rows();

    /**
     * 检查这个两个模板是否都已经压缩过了。
     * 
     */
    Assert(sp0.is_compressed(), ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * 检查这两个模板是否具有相同的行数以及第一个模板是否
     * 是一个方阵的模板。
     * 
     */
    Assert (m == sp0.n_cols(), ExcDimensionDontMatch(m, sp0.n_cols()));
    Assert (m == sp1.n_rows(), ExcDimensionDontMatch(m, sp1.n_rows()));

    int i, j;
    u_int n1 = sp1.n_cols();
    std::vector<u_int> row_length(m + n1, 1);
    for (i = 0;i < m;i ++)
      row_length[i] = sp0.row_length(i) + sp1.row_length(i);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	row_length[m + p_column1[j]] += 1;
    sp.reinit(m + n1, m + n1, row_length);
    for (i = 0;i < m;i ++) {
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++) {
	sp.add(i, m + p_column1[j]);
	sp.add(m + p_column1[j], i);
      }
    }
    sp.compress();
  };

  /** 
   * 将两个稀疏矩阵模板按对角方向连接成为一个新的模板
   * 
   * @param sp0 左上角的模板
   * @param sp1 右下角的模板
   * @param sp 得到的新模板
   */
  void dCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows(), m1 = sp1.n_rows();
    u_int n0 = sp0.n_cols(), n1 = sp1.n_cols();

    /**
     * 检查这个两个模板是否都已经压缩过了。
     * 
     */
    Assert(sp0.is_compressed(), SparsityPattern::ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    int i = 0, j;
    if (m0 + m1 == n0 + n1) i = 1;
    std::vector<u_int> row_length(m0 + m1, i);
    for (i = 0;i < m0;i ++)
      row_length[i] += sp0.row_length(i); 
    for (i = 0;i < m1;i ++)
      row_length[i + m0] += sp1.row_length(i);
    sp.reinit(m0 + m1, n0 + n1, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++)
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
    for (i = 0;i < m1;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(m0 + i, n0 + p_column1[j]);
    sp.compress();
  };


  /** 
   * 将四个稀疏矩阵模板完全连接为一个新的模板
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A_{00} & A_{01} \\
   *      A_{01} & A_{11}
   *    \end{array}\right)
   * \f]
   * 
   * @param sp00 左上角矩阵的模板
   * @param sp01 右上角矩阵的模板
   * @param sp10 左下角矩阵的模板
   * @param sp11 右下角矩阵的模板
   * @param sp 得到的新模板
   */  
  void fullCatSparsityPattern(const SparsityPattern& sp00,
			      const SparsityPattern& sp01,
			      const SparsityPattern& sp10,
			      const SparsityPattern& sp11,
			      SparsityPattern& sp)
  {
    SparsityPattern sp0, sp1;
    hCatSparsityPattern(sp00, sp01, sp0);
    hCatSparsityPattern(sp10, sp11, sp1);
    vCatSparsityPattern(sp0, sp1, sp);
    sp.compress();
  };

#define number double
  template void hCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);
  template void vCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);
  template void gammaCatSparseMatrix<number>(const SparseMatrix<number>&,
					     const SparseMatrix<number>&,
					     SparsityPattern&,
					     SparseMatrix<number>&,
					     bool);
  template void dCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);

  template void fullCatSparseMatrix(const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    SparsityPattern&,
				    SparseMatrix<number>&,
				    bool is_pattern_ok);
#undef number
};

/**
 * end of file
 * 
 */
