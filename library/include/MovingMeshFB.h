/**
 * @file   MovingMeshFB.h
 * @author Robert Lie
 * @date   Mon May 28 22:20:25 2007
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MovingMeshFB_h__
#define __MovingMeshFB_h__

#include "AMGSolver.h"
#include "Geometry.h"
#include "EasyMesh.h"

AFEPACK_OPEN_NAMESPACE

/**
 * 使用边界上不移动的方式的移动网格方法。
 *
 */
class MovingMeshFB : public EasyMesh
{
	std::vector<Point<2> > logical_node;
	double move_step_length;
	int n_move_step;
	std::vector<Point<2> > move_direction;
	std::vector<Point<2> > logical_move_direction;
	std::vector<float> mon;

	int n_interior_node;
	int n_boundary_node;
	std::vector<int> index;
	std::vector<int> interior_node_index;
	std::vector<int> boundary_node_index;

	SparsityPattern spM;
	SparsityPattern spN;
        SparseMatrix<double> M, N;
        AMGSolver solver;
public:
	MovingMeshFB();
	virtual ~MovingMeshFB();
public:
	const std::vector<float>& monitor() const {return mon;};
	std::vector<float>& monitor() {return mon;};
	const float& monitor(const int& i) const {return mon[i];};
	float& monitor(const int& i) {return mon[i];}
	const std::vector<Point<2> >& moveDirection() const {return move_direction;};
	std::vector<Point<2> >& moveDirection() {return move_direction;};
	const Point<2>& moveDirection(const int& i) const {return move_direction[i];};
	Point<2>& moveDirection(const int& i) {return move_direction[i];};
	
	std::vector<double> moveDirection(const Point<2>& point, const int& element) const;
	std::vector<std::vector<double> > moveDirection(const std::vector<Point<2> >& point, const int& element) const;
	double moveDirectionDivergence(const int& element) const;
	const double& moveStepLength() const {return move_step_length;};
	double& moveStepLength() {return move_step_length;};
	const int& n_moveStep() const {return n_move_step;};
	int& n_moveStep() {return n_move_step;};
	
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
	virtual void getLogicalMesh();
private:
	void getMoveDirection();
	void readDummy(std::ifstream& is);
};

AFEPACK_CLOSE_NAMESPACE

#endif // __MovingMeshFB_h__

/**
 * end of file
 * 
 */



