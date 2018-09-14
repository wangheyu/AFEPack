/**
 * @file   BoundaryCondition.h
 * @author Robert Lie
 * @date   Mon Nov 27 12:17:41 2006
 * 
 * @brief  
 * 
 * 
 */

#ifndef __BoundaryCondition_h__
#define __BoundaryCondition_h__

#include <cstdarg>
#include <vector>
#include <set>
#include <map>

#include <lac/sparse_matrix.h>
#include <lac/vector.h>

#include "Miscellaneous.h"

/**
 * 边界条件的基类。这个类拥有一个边界类型整数 _type，用以表征它描述
 * 的边界条件的类型，目前，这个类型事实上主要是具有象征意义，我们能
 * 够自动处理的还是只有狄氏边界条件。
 *
 * 另外，这个类还管理着一个材料标识的数组表示具有这些材料标识的自由
 * 度都将会使用本边界条件，这个数组可以通过 add_mark 函数来进行管理。
 * 注意：边界条件的材料标识只能加不能减的！
 *
 * 两个虚函数 value 和 gradient 用来给派生类做接口，一个计算边界函数
 * 的函数值，一个计算边界函数的梯度向量。
 *
 */
class BCondition {
 public:
  /// 三种基本的边界类型，作为本类的静态变量
  static const int DIRICHLET;
  static const int NEUMANN;
  static const int ROBIN;

 public:
  /*@\{*/
  /// 构造函数和析构函数
  BCondition() {}
  BCondition(int type) : _type(type) {}
  BCondition(const BCondition& b) :
    _type(b._type), _bmark(b._bmark) {}
  virtual ~BCondition() {}
  /*@\}*/

 public:
  /// 拷贝操作符
  BCondition& operator=(const BCondition& b) {
    _type = b._type;
    _bmark = b._bmark;

    return *this;
  }

  /*@\{*/
  /// 读写边界的类型
  int type() const {return _type;}
  int& type() {return _type;}
  /*@\}*/

  /*@\{*/
  /// 读写材料标识数组
  const std::set<int, std::less<int> >& bound_mark() const {return _bmark;}
  std::set<int, std::less<int> >& bound_mark() {return _bmark;}
  /*@\}*/

  /*@\{*/
  /// 往材料标识数组中加入材料标识
  void add_mark(int bm) {add_one_mark(bm);}
  void add_mark(const std::vector<int>& bm) {
    u_int n = bm.size();
    for (u_int i = 0;i < n;++ i) {
      add_one_mark(bm[i]);
    }
  }
  template <class V>
    void add_mark(u_int n, const V& bm) {
    for (u_int i = 0;i < n;i ++) {
      add_one_mark(bm[i]);
    }
  }
  void add_mark(u_int n, int bm0, int bm1, ...) {
    add_one_mark(bm0);
    add_one_mark(bm1);

    va_list ap;
    va_start(ap, bm1);
    for (u_int i = 2;i < n;i ++) {
      add_one_mark(va_arg(ap, int));
    }
    va_end(ap);
  }
  /*@\}*/

 private:
  void add_one_mark(int bm) {
    if (bm == 0) {
      std::cerr << "警告：材料标识 0 表示区域内部，不能加给边界条件。"
                << std::endl;
      return;
    }
    _bmark.insert(bm);
  }

 public:
  /// 边界条件函数的求值函数
  virtual void value(const void * p, void * v) const {}
  /// 边界条件函数的求梯度函数
  virtual void gradient(const void * p, void * g) const {}

 private:
  int _type; /// 边界的类型
  std::set<int, std::less<int> > _bmark; /// 材料标识数组
};



/**
 * 通过使用函数指针来实现 value 和 gradient 两个函数的边界条件。
 *
 */
template <class P, class V, class G = std::vector<V> >
class BCFunction : public BCondition {
  public:
  typedef void (*val_fun_t)(const P&, V&);
  typedef void (*grad_fun_t)(const P&, G&);

  private:
  static void _default_val(const P& p, V& v) {v = V();}
  static void _default_grad(const P& p, G& g) {g = G();}

  val_fun_t _val_fun; /// 求值函数指针
  grad_fun_t _grad_fun; /// 求梯度函数指针

  public:
  BCFunction(val_fun_t vf = &_default_val, 
             grad_fun_t gf = &_default_grad) :
  _val_fun(vf), 
  _grad_fun(gf) {}

  BCFunction(int type, 
             val_fun_t vf = &_default_val, 
             grad_fun_t gf = &_default_grad) : 
  BCondition(type), 
  _val_fun(vf), 
  _grad_fun(gf) {}

  virtual ~BCFunction() {}

  public:
  /*@\{*/
  /// 读写两个函数指针变量。
  const val_fun_t& value_fun_ptr() const {
    return _val_fun;
  }
  val_fun_t& value_fun_ptr() {
    return _val_fun;
  }

  const grad_fun_t& gradient_fun_ptr() const {
    return _grad_fun;
  }
  grad_fun_t& gradient_fun_ptr() {
    return _grad_fun;
  }
  /*@\}*/

  public:
  /*@\{*/
  /// 通过调用函数指针的变量来进行求值和求梯度。
  virtual void value(const void * p, void * v) const {
    (*_val_fun)(*(const P *)p, *(V *)v);
  }
  virtual void gradient(const void * p, void * g) const {
    (*_grad_fun)(*(const P *)p, *(G *)g);
  }
  /*@\}*/
};


/**
 * 边界条件管理器。这个类管理着一堆的边界条件，然后在具体特定的材料
 * 标识的自由度处使用相应的边界条件。目前，这个类还很粗糙，只能自动
 * 处理一定情况下的狄氏边值。这个一定情况指的是：
 *
 *   - 单元上的基函数在它们的插值点上是正交的；
 *   - 基函数在自己插值点上的值取 1；
 *
 * 在其它情况下，这里的处理方法会给出错误的结果。
 * 
 * 使用本函数中定义的这几个类的典型代码行如下：
 *
 * <pre>
 *
 *    BCFunction<Point<DIM>, double> bc(BCondition::DIRICHLET, u_b);
 *    bc.add_mark(5, 1, 2, 3, 4, 5);
 *    BCAdmin bc_admin;
 *    bc_admin.add(bc);
 *    bc_admin.apply(fem_space, mat, u_h, rhs);
 *
 * </pre>
 *
 * 其中定义了一个狄氏边值，取值方式由函数 u_b 给定，然后加入了 5 个边
 * 界标识，分别为 1, 2, 3, 4, 5。这个边值加入到管理器中，然后应用到矩
 * 阵和向量上。
 *
 */
class BCAdmin {
 public:
  typedef BCondition * bc_ptr_t;

 private:
  std::map<int, bc_ptr_t, std::less<int> > _map;

 public:
  /**
   * 应用狄氏边值条件在一个线性系统上。
   * 
   */
  template <
    class SP,    /// 有限元空间的类型
    class MAT,   /// 稀疏矩阵的类型，现在是设为 SparseMatrix<double>
    class XVEC,  /// 解向量的类型，现在是设为 Vector<double>
    class BVEC   /// 右端向量的类型，现在是设为 Vector<double>
    >
    void apply(const SP& sp,
               MAT& A, 
               XVEC& u, 
               BVEC& f, 
               bool preserve_symmetry = true) 
    {
      u_int n_dof = sp.n_dof();
      const SparsityPattern& spA = A.get_sparsity_pattern();
      const std::size_t * rowstart = spA.get_rowstart_indices();
      const u_int * colnum = spA.get_column_numbers();
      for (u_int i = 0;i < n_dof;++ i) {
        int bm = sp.dofInfo(i).boundary_mark;
        if (bm == 0) continue; /// 0 缺省指区域内部
        const bc_ptr_t bc = this->find(bm);
        /// 如果这个边界条件不由本对象处理，或者不是狄氏边界，则跳过去。
        if (bc == NULL) continue;
        if (bc->type() != BCondition::DIRICHLET) continue;
        bc->value((const void *)(&sp.dofInfo(i).interp_point), 
                  (void *)(&u(i)));

        f(i) = A.diag_element(i)*u(i);
        for (u_int j = rowstart[i] + 1;j < rowstart[i + 1];++ j) {
          A.global_entry(j) -= A.global_entry(j);
        }
        if (preserve_symmetry) {
          for (u_int j = rowstart[i] + 1;j < rowstart[i + 1];++ j) {
            u_int k = colnum[j];
            const u_int * p = std::find(&colnum[rowstart[k] + 1],
                                        &colnum[rowstart[k + 1]], i);
            if (p != &colnum[rowstart[k+1]]) {
              u_int l = p - &colnum[rowstart[0]];
              f(k) -= A.global_entry(l)*f(i)/A.diag_element(i);
              A.global_entry(l) -= A.global_entry(l);
            }
          }
        }
      }
    }

  /**
   * 将向量 f 中由本对象负责处理的自由度相同的指标上的元素清零。
   * 
   */
  template <
    class SP,  /// 有限元空间的类型
    class VEC  /// 向量的类型，现在是设为 Vector<double>
    >
    void clear_entry(const SP& sp,
                     VEC& f)
    {
      u_int n_dof = sp.n_dof();
      for (u_int i = 0;i < n_dof;++ i) {
        int bm = sp.dof_info(i).bound_mark();
        if (bm == 0) continue; /// 0 缺省指区域内部
        const bc_ptr_t bc = find(bm);
        if (bc != NULL) f(i) -= f(i);
      }
    }

  /**
   * 往本对象中加入一个边界条件。一个边界条件被加入到管理器中以后，
   * 如果另外调用了 add_mark 以后，是不会自动起作用的。如果有此要求，
   * 需要重新调用本函数。
   * 
   */
  template <class BC> 
    void add(BC& b) 
    {
      const std::set<int, std::less<int> >& bm = b.bound_mark();
      std::set<int, std::less<int> >::const_iterator 
        the_bm = bm.begin(), end_bm = bm.end();
      for (;the_bm != end_bm;++ the_bm) {
        _map[*the_bm] = &b;
      }
    }

  /**
   * 找出对应于材料标识 bm 的边界条件的指针，如果没有找到，则返回
   * NULL。
   * 
   */
  bc_ptr_t find(int bm) const 
    {
      std::map<int, bc_ptr_t, std::less<int> >::const_iterator
        the_ptr = _map.find(bm);
      if (the_ptr != _map.end()) {
        return the_ptr->second;
      } else {
        return NULL;
      }
    }
};

#endif // __BoundaryCondition_h__

/**
 * end of file
 * 
 */
