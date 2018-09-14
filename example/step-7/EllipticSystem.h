

// EllipticSystem.h

#ifndef _ELLIPTICSYSTEM_H_

#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/BilinearOperator.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>

#define DIM 2

class EllipticSystem {
public:
	// this class is used to construct matrix $A^0$ and $A^1$
	class Matrix_A : public StiffMatrix<DIM,double> {
	private:
		double (*a[2][2])(const double *);
		double (*c)(const double *);
	public:
		Matrix_A(double (*a00)(const double *), double (*a01)(const double *),
			double (*a10)(const double *), double (*a11)(const double *),
			double (*c0)(const double *), FEMSpace<double,DIM>& sp) :
				StiffMatrix<DIM,double>(sp) {
			a[0][0] = a00;
			a[0][1] = a01;
			a[1][0] = a10;
			a[1][1] = a11;
			c = c0;
		};
		virtual ~Matrix_A() {};
	public:
		virtual void getElementMatrix(const Element<double,DIM>&,
			const Element<double,DIM>&,
			const ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL);
 	};

	// this class is used to construct matrix $C^0$ and $C^1$
	class Matrix_C : public L2InnerProduct<DIM,double> {
	private:
		double (*c)(const double *);
	public:
		Matrix_C(double (*c1)(const double *), FEMSpace<double,DIM>& sp0, FEMSpace<double,DIM>& sp1) :
				c(c1), L2InnerProduct<DIM,double>(sp0, sp1) {};
		~Matrix_C() {};
	public:
		virtual void getElementMatrix(const Element<double,DIM>&,
  			const Element<double,DIM>&,
			const ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL);
	};

private:
	std::string			mesh_file;			// file name of mesh
	HGeometryTree<DIM>		h_geometry_tree;		// the geometry tree
	IrregularMesh<DIM>		irregular_mesh[DIM];		// the meshes
	FEMSpace<double,DIM>		fem_space[DIM];			// the finite element space
	FEMFunction<double,DIM>		solution[DIM];			// the solution

	Indicator<DIM>			indicator[DIM];			// indicator to adapt the meshes

	// variables used to contruct template elements
	std::vector<TemplateElement<double,DIM,DIM> >	template_element;
	TemplateGeometry<DIM>				triangle_template_geometry;
	CoordTransform<DIM,DIM>				triangle_coord_transform;
	TemplateDOF<DIM>				triangle_1_template_dof;
	BasisFunctionAdmin<double,DIM,DIM>		triangle_1_basis_function;
	TemplateGeometry<DIM>				twin_triangle_template_geometry;
	CoordTransform<DIM,DIM>				twin_triangle_coord_transform;
	TemplateDOF<DIM>				twin_triangle_1_template_dof;
	BasisFunctionAdmin<double,DIM,DIM>		twin_triangle_1_basis_function;
public: // contructors and destructor
	EllipticSystem(const std::string& filename);
	~EllipticSystem() {};
public: //
	// initialize those memeber variables
	void init();
	// the engine to run the program
	void run();
	// adapt mesh
	void adaptMesh();
	// build finite element space
	void buildFEMSpace();
	// construct the linear system and solve it
	void solve();
	// calculate the indicator used to adapt the mesh using Kelly estimator
	void getIndicator();
	// calculate the error of the solution
	void getError();
	// save data
	void saveData();
};

#endif // _ELLIPTICSYSTEM_H_

