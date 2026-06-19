/**
 * @file   ex26.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Nov 27 12:28:09 2009
 * 
 * @brief  测试使用 PETSc 来求解一个椭圆型方程。
 *
 * 本例中的新内容包括：
 *
 *   - 构造 PETSc 的稀疏矩阵和分布式向量；
 *   - 调用 PETSc 的求解器求解线性系统；
 *
 * 运行语法：
 *
 *   mpiexec -np # ./ex26 mesh_dir refine_round
 * 
 */

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_DOF.h>
#include <AFEPack/mpi/MPI_ULoadBalance.h>
#include <petscksp.h>

#define DIM 2
#define DOW DIM
#define PI M_PIl

double _f_(const double * p) {
  return (1.0 + DIM*4.0*PI*PI)*sin(2*PI*p[0])*cos(2*PI*p[1]);
}

static char help[] = "Solving an elliptic equation.";

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef AFEPack::FEMSpace<double,DIM,DOW> fe_space_t;  
  typedef AFEPack::MPI::DOF::GlobalIndex<forest_t, fe_space_t> global_index_t;

  PetscInitialize(&argc, &argv, (char *)NULL, help);

  forest_t forest(PETSC_COMM_WORLD);
  ir_mesh_t ir_mesh;
  AFEPack::MPI::load_mesh(argv[1], forest, ir_mesh); /// 从一个目录中读入网格数据

  int round = 0;
  if (argc >= 3) round = atoi(argv[2]);

  ir_mesh.globalRefine(round);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  setenv("AFEPACK_TEMPLATE_PATH", "/usr/local/AFEPack/template/triangle", 1);

  TemplateGeometry<DIM> tri;
  tri.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> tri_ct;
  tri_ct.readData("triangle.crd_trs");
  TemplateDOF<DIM> tri_td(tri);
  tri_td.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> tri_bf(tri_td);
  tri_bf.readData("triangle.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > tmp_ele(1);
  tmp_ele[0].reinit(tri, tri_td, tri_ct, tri_bf);

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  fe_space_t fem_space(mesh, tmp_ele);
  u_int n_ele = mesh.n_geometry(DIM);
  fem_space.element().resize(n_ele);
  for (int i = 0;i < n_ele;i ++) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  std::cout << "Building global indices ... " << std::flush;
  global_index_t global_index(forest, fem_space);
  global_index.build();
  std::cout << "OK!" << std::endl;

  std::cout << "Building the linear system ... " << std::flush;
  Mat A;
  Vec x, b;
  MatCreateAIJ(PETSC_COMM_WORLD,
	       global_index.n_primary_dof(), global_index.n_primary_dof(),
	       PETSC_DECIDE, PETSC_DECIDE, 
	       0, PETSC_NULL, 0, PETSC_NULL, &A);
  VecCreateMPI(PETSC_COMM_WORLD, global_index.n_primary_dof(), PETSC_DECIDE, &b);
  fe_space_t::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    double vol = the_ele->templateElement().volume();
    const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(5);
    std::vector<Point<DIM> > q_pnt = the_ele->local_to_global(qi.quadraturePoint());
    int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
    std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_pnt);
    std::vector<std::vector<std::vector<double> > > bas_grad = the_ele->basis_function_gradient(q_pnt);

    const std::vector<int>& ele_dof = the_ele->dof();
    u_int n_ele_dof = ele_dof.size();
    FullMatrix<double> ele_mat(n_ele_dof, n_ele_dof);
    Vector<double> ele_rhs(n_ele_dof);
    for (u_int l = 0;l < n_q_pnt;++ l) {
      double JxW = vol*jac[l]*qi.weight(l);
      double f_val = _f_(q_pnt[l]);
      for (u_int i = 0;i < n_ele_dof;++ i) {
        for (u_int j = 0;j < n_ele_dof;++ j) {
          ele_mat(i, j) += JxW*(bas_val[i][l]*bas_val[j][l] +
                                innerProduct(bas_grad[i][l], bas_grad[j][l]));
        }
        ele_rhs(i) += JxW*f_val*bas_val[i][l];
      }
    }
    /**
     * 此处将单元矩阵和单元载荷先计算好，然后向全局的矩阵和载荷向量上
     * 集中，可以提高效率。
     */

    std::vector<int> indices(n_ele_dof);
    for (u_int i = 0;i < n_ele_dof;++ i) {
      indices[i] = global_index(ele_dof[i]);
    }
    MatSetValues(A, n_ele_dof, &indices[0], n_ele_dof, &indices[0], &ele_mat(0,0), ADD_VALUES);
    VecSetValues(b, n_ele_dof, &indices[0], &ele_rhs(0), ADD_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
  VecDuplicate(b, &x);
  std::cout << "OK!" << std::endl;

  KSP solver;
  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A, SAME_NONZERO_PATTERN);
  KSPSetType(solver, KSPCG);
  KSPSetFromOptions(solver);
  KSPSolve(solver, b, x);

  if (forest.rank() == 0) {
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver,&reason);
    if (reason == KSP_DIVERGED_INDEFINITE_PC) {
      printf("\nDivergence because of indefinite preconditioner;\n");
      printf("Run the executable again but with -pc_ilu_shift option.\n");
    } else if (reason<0) {
      printf("\nOther kind of divergence: this should not happen.\n");
    } else {
      PetscInt its;
      KSPGetIterationNumber(solver,&its);
      printf("\nConvergence in %d iterations.\n",(int)its);
    }
    printf("\n");
  }

  MatDestroy(&A);
  VecDestroy(&b);
  KSPDestroy(&solver);

  /// 准备解函数
  FEMFunction<double,DIM> u_h(fem_space);
  Vec X;
  VecCreateSeqWithArray(PETSC_COMM_SELF, 1, global_index.n_local_dof(), &u_h(0), &X);

  /// 将 PETSc 解出来的向量取出到有限元函数 u_h 中来
  std::vector<int> primary_idx(global_index.n_primary_dof());
  global_index.build_primary_index(&primary_idx[0]);
  IS is;
  ISCreateGeneral(forest.communicator(), global_index.n_local_dof(),
		  &global_index(0), PETSC_USE_POINTER, &is);
  VecScatter scatter;
  VecScatterCreate(x, is, X, PETSC_NULL, &scatter);
  VecScatterBegin(scatter, x, X, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scatter, x, X, INSERT_VALUES, SCATTER_FORWARD);

  /// 清理 PETSc 的变量
  VecDestroy(&x);
  VecDestroy(&X);
  VecScatterDestroy(&scatter);
  ISDestroy(&is);

  char filename[1024];
  sprintf(filename, "u_h%d.dx", forest.rank());
  u_h.writeOpenDXData(filename);

  PetscFinalize();

  return 0;
}

/**
 * end of file
 * 
 */

