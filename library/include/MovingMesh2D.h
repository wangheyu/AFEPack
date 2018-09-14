/**
 * @file   MovingMesh2D.h
 * @author Robert Lie
 * @date   Thu Feb 22 11:48:08 2007
 * 
 * @brief  declaration of class MovingMesh2D
 * 
 * 
 */

#ifndef __MovingMesh2D_h__
#define __MovingMesh2D_h__

#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include "AMGSolver.h"
#include "Geometry.h"
#include "EasyMesh.h"

AFEPACK_OPEN_NAMESPACE

class MovingMesh2D : public EasyMesh
{
 public:
  typedef GeometryBM::bmark_t bound_t;
  struct Vertex : public Point<2> {
    int index;
    int boundary_mark;
  };
  struct Edge {
    int index;
    int vertex[2];
    int boundary_mark;
    double normal[2]; /// 物理区域边上的法向
    double logical_normal[2]; /// 计算区域边上的法向
  };
  struct Domain {
    int n_vertex;
    int n_edge;
    std::vector<Vertex> physical_domain_vertex;
    std::vector<Vertex> logical_domain_vertex;
    std::vector<Edge> edge;
  };
  class Solver { /**< 特殊设计的代数多重网格求解器来求解网格方程 */
  private:
    bool					is_initialized;
    int						n_project;
    std::vector<SparseMatrix<double> *>		project_matrix;
    std::vector<SparseMatrix<double> *>	        project_matrix_r;
    std::vector<SparseMatrix<double> *>	        projected_matrix;
    std::vector<const std::vector<int> *>	boundary_mark;
    u_int					smooth_step;
    u_int					min_order;
    const Domain *                              domain;
  public:
    Solver();
    ~Solver();
  public:
    u_int smoothStep() const {return smooth_step;}
    u_int& smoothStep() {return smooth_step;}
    u_int minimalOrder() const {return min_order;}
    u_int& minimalOrder() {return min_order;}

    void reinit(SparseMatrix<double>&,
		const std::vector<int>&,
		const Domain&);
    void clear();
    void solve(std::vector<Vector<double> >& x, 
	       const std::vector<Vector<double> >& r,
	       u_int steps = 5) const;
  private:
    void Project(const SparseMatrix<double>& M,
		 const std::vector<int>& bm,
		 SparseMatrix<double> *& P, 
		 SparseMatrix<double> *& PMPT,
		 SparseMatrix<double> *& Pt, 
		 std::vector<bound_t> *& pbm);
    void Init(SparseMatrix<double>&,
	      const std::vector<int>&,
	      const Domain&);
    void GaussSidel(const SparseMatrix<double>& M, 
		    std::vector<Vector<double> >& x, 
		    const std::vector<Vector<double> >& r,
		    const std::vector<int>& bm,
		    const u_int& s) const;
    SparseMatrix<double> * getPMPT(const SparseMatrix<double> & P, 
                                   const SparseMatrix<double> & M,
                                   const SparseMatrix<double> & Pt) const;
    void lazyPMPT(const SparseMatrix<double>& P, 
		  const SparseMatrix<double>& M,
		  const SparseMatrix<double>& Pt,
		  SparseMatrix<double>& A) const;
  };

 private:
  Domain domain;
  std::vector<Point<2> > logical_node;
  double move_step_length;
  int n_move_step;
  std::vector<Point<2> > move_direction;
  std::vector<Point<2> > logical_move_direction;
  std::vector<double> mon;

  std::vector<std::vector<int> > mb_node; /// 第 i 条边上的节点的指标
  std::vector<int> boundary_mark; /// 给求解器使用的材料标识

  SparsityPattern spM;
  SparseMatrix<double> M;
  Solver solver;

  static int primes[]; /// 素数表
  u_int solve_step; /// 求解器求解时候使用的迭代步数，缺省为 5
  /// 网格移动的迭代最大迭代步数，缺省为 100；0 为不设限制
  u_int max_step; 
  /// 网格移动结束的相对误差容忍度，缺省为 10%
  double tol; 
  /// 解调和映照后，能量减少的百分比
  double energy_reduction; 

 public:
  MovingMesh2D();
  virtual ~MovingMesh2D();

 public:
  const std::vector<double>& monitor() const {return mon;};
  std::vector<double>& monitor() {return mon;};
  const double& monitor(const int& i) const {return mon[i];};
  double& monitor(const int& i) {return mon[i];}

  const std::vector<Point<2> >& moveDirection() const {return move_direction;};
  std::vector<Point<2> >& moveDirection() {return move_direction;};
  const Point<2>& moveDirection(const int& i) const {return move_direction[i];};
  Point<2>& moveDirection(const int& i) {return move_direction[i];};
	
  std::vector<double> moveDirection(const Point<2>& point, 
                                    const int& element) const;
  std::vector<std::vector<double> > moveDirection(const std::vector<Point<2> >& point, 
                                                  const int& element) const;
  double moveDirectionDivergence(const int& element) const;

  double moveStepLength() const {return move_step_length;};
  double& moveStepLength() {return move_step_length;};

  int moveStep() const {return n_move_step;};
  int& moveStep() {return n_move_step;};

  u_int solveStep() const {return solve_step;}
  u_int& solveStep() {return solve_step;}

  double tolerence() const {return tol;}
  double& tolerence() {return tol;}

  u_int maxStep() const {return max_step;}
  u_int& maxStep() {return max_step;}
	
  void moveMesh();
  void outputPhysicalMesh(const std::string& file);
  void outputLogicalMesh(const std::string& file);
  virtual void getMonitor();
  virtual void smoothMonitor(int step = 1);
  virtual void updateMesh();
  virtual void updateSolution() = 0;
  virtual void outputSolution() = 0;
  virtual void getMoveStepLength();
  void readDomain(const std::string& file);

 private:
  void getLogicalMesh();
  void getMoveDirection();
  void readDummy(std::ifstream& is);
  void parseBoundary();

  void outputMoveDirection(const std::string&); /// 调试使用的函数
};

AFEPACK_CLOSE_NAMESPACE

#endif /**< __MovingMesh2D_h__ */

/**
 * end of file
 * 
 */

