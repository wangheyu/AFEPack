////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>

#include <AFEPack/LOBPCG.h>

#define PI (4.0*atan(1.0))
#define DIM 3 


int main(int argc, char * argv[])
{

  //Mesh<DIM> mesh;
  //mesh.readData(argv[1]);

  HGeometryTree<DIM> h_tree;
  IrregularMesh<DIM>* ir_mesh;

  h_tree.readMesh(argv[1]);
  ir_mesh = new IrregularMesh<DIM>(h_tree);
  ir_mesh->reinit(h_tree);
  ir_mesh->globalRefine(1);
  ir_mesh->semiregularize();
  ir_mesh->regularize(false);

  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();

  TemplateGeometry<DIM>	triangle_template_geometry;
  triangle_template_geometry.readData("tetrahedron.tmp_geo");
  CoordTransform<DIM,DIM>	triangle_coord_transform;
  triangle_coord_transform.readData("tetrahedron.crd_trs");
  TemplateDOF<DIM>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("tetrahedron.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("tetrahedron.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function);

  FEMSpace<double,DIM> fem_space(mesh, template_element);
 // DGFEMSpace<double, DIM> fem_space (mesh, template_element, edge_template_element);
  //FEMSpace<double,2> fem_space;
  //fem_space = new FEMSpace<double, 2>(mesh, template_element);
  //fem_space.reinit(mesh, template_element);  
	
  int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  ////////////////////////////////////////////////
  // Begin to construct the eigenvalue problem
  StiffMatrix<DIM,double> A(fem_space);
  MassMatrix<DIM,double> B(fem_space);
  A.algebricAccuracy() = 4;
  A.build();
  B.algebricAccuracy() = 4;
  B.build();


  std::vector<int> boundaryDOFIndex;
  std::vector<int> internalDOFIndex;

  boundaryDOFIndex.clear();
  internalDOFIndex.clear();

  std::set<int> boundaryDOFSet;
  std::set<int> internalDOFSet;
  const int& n_ele = fem_space.n_element();
  for(int i = 0;i < n_ele;++ i){
    Element<double, DIM>& the_ele = fem_space.element(i);
    const int& ele_idx = the_ele.index();

    std::vector<int>& ele_dof = the_ele.dof();
    const int& n_ele_dof = ele_dof.size();

    for(int j = 0;j < n_ele_dof;++ j){
      if(fem_space.dofBoundaryMark(ele_dof[j]) == 0){
        internalDOFSet.insert(ele_dof[j]);
      }
      else{
        boundaryDOFSet.insert(ele_dof[j]);
      }
    }
  }

  /// copy boundaryDOFset to boundaryDOFIndex to realize the parallel OpenMP
  int bndDof_idx = 0;
  int intDof_idx = 0;
  boundaryDOFIndex.resize(boundaryDOFSet.size());
  internalDOFIndex.resize(internalDOFSet.size());
  std::set<int>::iterator it;
  for(it = boundaryDOFSet.begin();it != boundaryDOFSet.end();++ it){
    boundaryDOFIndex[bndDof_idx ++] = *it;
  }

  for(it = internalDOFSet.begin();it != internalDOFSet.end();++ it){
    internalDOFIndex[intDof_idx ++] = *it;
  }

      /// boundary condition

      double bnd_value = 0.;
      const int& n_boundaryDOFIndex = boundaryDOFIndex.size();

      const int& n_dof = fem_space.n_dof();
      for(int i = 0;i < n_boundaryDOFIndex;++ i){
        //        if(fem_space->dofInfo(i).boundary_mark == 0) continue;
        const int& dof_idx = boundaryDOFIndex[i];

        {
	  const SparsityPattern& spA = A.get_sparsity_pattern();
          const std::size_t * row_start
            = spA.get_rowstart_indices();
          const unsigned int * column
            = spA.get_column_numbers();
          for (int j = row_start[dof_idx] + 1;j < row_start[dof_idx + 1];j ++) {
            A.global_entry(j) = 0.0;
            int k = column[j];
            for (int l = row_start[k] + 1;l < row_start[k + 1];l ++) {
              if(column[l] == dof_idx){
                A.global_entry(l) = 0.0;
                break;
              }
            }
          }
        }
        {
	  const SparsityPattern& spB = B.get_sparsity_pattern();
          const std::size_t * row_start
                        = spB.get_rowstart_indices();
          const unsigned int * column
            = spB.get_column_numbers();
          for (int j = row_start[dof_idx] + 1;j < row_start[dof_idx + 1];j ++) {
            B.global_entry(j) = 0.0;
            int k = column[j];
            for (int l = row_start[k] + 1;l < row_start[k + 1];l ++) {
              if(column[l] == dof_idx){
                B.global_entry(l) = 0.0;
                break;
              }
            }
          }
        }
      }

  int n_orbital = 10;
  std::vector<Vector<double>> Orbital(n_orbital);
  std::vector<double> Eigenvalue(n_orbital);
  for(int i = 0;i < n_orbital;++ i){
    Orbital[i].reinit(n_dof);
  }

  for(int i = 0;i < n_orbital;++ i){
    for(int j = 0;j < n_dof;++ j){
      Orbital[i](j) = rand()%100;
    }
  }

  double time4LOBPCG;
  time4LOBPCG = omp_get_wtime();
  LOBPCG lobpcg_solver;
  std::string str1("Residual.txt");
  lobpcg_solver.open_log(str1);
  lobpcg_solver.reinit(A, B, Orbital, Eigenvalue);
  double lobpcg_tol = 1.0e-5;
  int lobpcg_iter = 100;
  //lobpcg_solver.solve(lobpcg_tol, lobpcg_iter, false);
  lobpcg_solver.solve(lobpcg_tol, lobpcg_iter);
  lobpcg_solver.close_log();
  time4LOBPCG = omp_get_wtime() - time4LOBPCG;
  std::ofstream output("Time.txt");
  output << time4LOBPCG << std::endl;
  output.close();
  return 0;
};


//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
