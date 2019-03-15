/**
 * @file   ex38.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Apr 30 13:51:43 2010
 * 
 * @brief  测试使用 Hypre 来求解一个椭圆型方程。
 *
 * 本例中的新内容包括：
 *
 *   - 使用 GlobalIndex 来构造 Hypre 的分布式的矩阵和向量；
 *   - 调用 Hypre 的 BoomerAMG 来求解方程组；
 *   - 将 Hypre 的向量散发成为内部的有限元函数；
 * 
 * 运行语法：
 *
 *   mpirun -np # ./ex38 mesh_dir refine_round
 *
 */

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_DOF.h>
#include <AFEPack/mpi/MPI_ULoadBalance.h>

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>

#define DIM 2
#define DOW DIM
#define PI M_PIl

double _f_(const double * p) {
  return (1.0 + DIM*4.0*PI*PI)*sin(2*PI*p[0])*cos(2*PI*p[1]);
}

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef AFEPack::FEMSpace<double,DIM,DOW> fe_space_t;  
  typedef AFEPack::MPI::DOF::GlobalIndex<forest_t, fe_space_t> global_index_t;

  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  ir_mesh_t ir_mesh;
  AFEPack::MPI::load_mesh(argv[1], forest, ir_mesh); /// 从一个目录中读入网格数据

  int round = 0;
  if (argc >= 3) round = atoi(argv[2]);

  ir_mesh.globalRefine(round);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

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

  /// 构造 Hypre 的分布式稀疏矩阵模板
  std::cout << "Build sparse matrix ... " << std::flush;
  HYPRE_IJMatrix A;
  HYPRE_IJMatrixCreate(forest.communicator(),
                       global_index.first_index(), global_index.last_index(),
                       global_index.first_index(), global_index.last_index(),
                       &A);
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(A);

  HYPRE_IJVector b;
  HYPRE_IJVectorCreate(forest.communicator(),
                       global_index.first_index(), global_index.last_index(),
                       &b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  fe_space_t::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    double vol = the_ele->templateElement().volume();
    const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(5);
    std::vector<AFEPack::Point<DIM> > q_pnt = the_ele->local_to_global(qi.quadraturePoint());
    int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
    std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_pnt);
    std::vector<std::vector<std::vector<double> > > bas_grad = the_ele->basis_function_gradient(q_pnt);

    const std::vector<int>& ele_dof = the_ele->dof();
    int n_ele_dof = ele_dof.size();
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
    global_index.translate(ele_dof, indices);
    for (u_int i = 0;i < n_ele_dof;++ i) {
      HYPRE_IJMatrixAddToValues(A, 1, &n_ele_dof, 
                                &indices[i], &indices[0], 
                                &ele_mat(i, 0));
    }
    HYPRE_IJVectorAddToValues(b, n_ele_dof, &indices[0], &ele_rhs(0));
  }
  HYPRE_IJMatrixAssemble(A);
  HYPRE_IJVectorAssemble(b);
  std::cout << "OK!" << std::endl;

  /// 准备解向量。
  HYPRE_IJVector x;
  HYPRE_IJVectorCreate(forest.communicator(), 
                       global_index.first_index(),
                       global_index.last_index(),
                       &x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);
  HYPRE_IJVectorAssemble(x);

  HYPRE_ParCSRMatrix A0;
  HYPRE_IJMatrixGetObject(A, (void **)(&A0));
  HYPRE_ParVector b0;
  HYPRE_IJVectorGetObject(b, (void **)(&b0));
  HYPRE_ParVector x0;
  HYPRE_IJVectorGetObject(x, (void **)(&x0));

  /// 调用 BoomerAMG 求解器。
  std::cout << "Solving the linear system ..." << std::flush;
  HYPRE_Solver solver;
  HYPRE_BoomerAMGCreate(&solver);
  HYPRE_BoomerAMGSetTol(solver, 1.0e-12);
  HYPRE_BoomerAMGSetMaxIter(solver, 500);
  HYPRE_BoomerAMGSetup(solver, A0, b0, x0);
  HYPRE_BoomerAMGSolve(solver, A0, b0, x0);
  HYPRE_BoomerAMGDestroy(solver);
  std::cout << "OK!" << std::endl;

  FEMFunction<double,DIM> u_h(fem_space);
  global_index.scatter_hypre_vector(x, u_h);

  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);

  char filename[1024];
  sprintf(filename, "u_h%d.dx", forest.rank());
  u_h.writeOpenDXData(filename);

  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

