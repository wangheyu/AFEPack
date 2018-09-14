/**
 * @file   SparseMatrixTool.h
 * @author Robert Lie
 * @date   Sun Nov  4 09:06:11 2007
 * 
 * @brief  
 * 
 * 
 */

#include <base/exceptions.h>
#include <lac/sparse_matrix.h>

#include "Miscellaneous.h"

/**
 * 定义一些对稀疏矩阵和稀疏矩阵模板进行操作的工具函数
 * 
 */
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
			   SparsityPattern& sp); 

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
			   SparsityPattern& sp);

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
   * 的模板。要求 \f$ A \f$ 和 \f$ B \f$ 具有相同的行数
   * 以及 \f$ A \f$ 是一个方阵。
   * 
   * @param sp0 矩阵 \f$ A \f$ 的模板
   * @param sp1 矩阵 \f$ B \f$ 的模板
   * @param sp 得到的新模板
   */
  void gammaCatSparsityPattern(const SparsityPattern& sp0,
			       const SparsityPattern& sp1,
			       SparsityPattern& sp);

  /**
   * 将两个稀疏矩阵模板按照对角方式连接成为一个新的模板。
   * 
   * @param sp0 左上角矩阵的模板
   * @param sp1 右下角矩阵的模板
   * @param sp 得到的新模板
   */
  void dCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp);

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
			      SparsityPattern& sp);

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
			  bool is_pattern_ok = true);

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
			  bool is_pattern_ok = true);

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
			      bool is_pattern_ok = true);

  /**
   * 将两个稀疏矩阵按照对角方式连接成为一个新的稀疏矩阵。
   * 
   * @param m0 左上角矩阵
   * @param m1 右下角矩阵
   * @param sp 得到的新矩阵的模板
   * @param m 得到的新矩阵
   * @param is_pattern_ok 如果为真，则对模板进行构造；否则
   *                      假定模板已经进行了正确的构造；缺
   *                      省值为真。
   */
  template <typename number>
    void dCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok = true);

  /** 
   * 将四个稀疏矩阵完全连接为一个新的稀疏矩阵
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A_{00} & A_{01} \\
   *      A_{10} & A_{11}
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
			     bool is_pattern_ok = true);

  DeclException0(ExcNotCompressed);
  DeclException2(ExcDimensionDontMatch, int, int,
		 << "The dimensions " << arg1 << " and " << arg2
		 << " do not match properly.");
};

/**
 * end of file
 * 
 */
