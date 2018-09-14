/**
 * @file   MovingMeshFB3D.h
 * @author Robert Lie
 * @date   Mon May 28 22:25:20 2007
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MovingMeshFB3D_h__
#define __MovingMeshFB3D_h__

#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include "AMGSolver.h"
#include "Geometry.h"

AFEPACK_OPEN_NAMESPACE

/**
 * 三维的边界上的点固定的移动网格方法
 *
 */
class MovingMeshFB3D : public Mesh<3,3>
{
 public:
  typedef GeometryBM::bmark_t bound_t;
 private:
  std::vector<Point<3> > logical_node;
  double move_step_length;
  u_int n_move_step;
  std::vector<Point<3> > move_direction;
  std::vector<Point<3> > logical_move_direction;
  std::vector<float> mon;

  u_int n_interior_node;
  u_int n_boundary_node;
  std::vector<int> index;
  std::vector<u_int> interior_node_index;
  std::vector<u_int> boundary_node_index;

  SparsityPattern spM, spN;
  SparseMatrix<double> M, N;
  AMGSolver solver;

  double tol;

 public:
  MovingMeshFB3D();
  virtual ~MovingMeshFB3D();

 public:
  const std::vector<float>& monitor() const {return mon;};
  std::vector<float>& monitor() {return mon;};
  const float& monitor(const u_int& i) const {return mon[i];};
  float& monitor(const u_int& i) {return mon[i];}

  const std::vector<Point<3> >& moveDirection() const {return move_direction;};
  std::vector<Point<3> >& moveDirection() {return move_direction;};
  const Point<3>& moveDirection(const u_int& i) const {return move_direction[i];};
  Point<3>& moveDirection(const u_int& i) {return move_direction[i];};
	
  std::vector<double> moveDirection(const Point<3>& point, 
                                    const int& element) const;
  std::vector<std::vector<double> > moveDirection(const std::vector<Point<3> >& point, 
                                                  const int& element) const;
  double moveDirectionDivergence(const u_int& element) const;

  double moveStepLength() const {return move_step_length;};
  double& moveStepLength() {return move_step_length;};

  u_int moveStep() const {return n_move_step;};
  u_int& moveStep() {return n_move_step;};

  double tolerence() const {return tol;}
  double& tolerence() {return tol;}
	
  void moveMesh();
  void outputPhysicalMesh(const std::string& file);
  void outputLogicalMesh(const std::string& file);
  void readDomain(const std::string& file);

  virtual void getMonitor();
  virtual void smoothMonitor(u_int step = 1);
  virtual void updateMesh();
  virtual void updateSolution() = 0;
  virtual void outputSolution() = 0;
  virtual void getMoveStepLength();
  virtual void getLogicalMesh();

 private:
  void getMoveDirection();
  void readDummy(std::ifstream& is);
};

AFEPACK_CLOSE_NAMESPACE

#endif /**< __MovingMeshFB3D_h__ */

/**
 * end of file
 * 
 */

