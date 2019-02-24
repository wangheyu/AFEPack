/**
 * @file   ex25.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Nov 30 13:14:11 2009
 * 
 * @brief  测试使用 Trilinos 来求解一个椭圆型方程。
 *
 * 本例中的新内容是处理边界条件，我们加上了一个狄氏边值的处理代码。
 * 
 * 运行语法：
 *
 *   mpirun -np # ./ex24 mesh_dir refine_round
 *
 */

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_DOF.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>
#include <Epetra_LinearProblem.h>
#include <ml_MultiLevelOperator.h>
#include <ml_epetra_preconditioner.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <AztecOO.h>

#define DIM 2
#define DOW DIM
#define PI M_PIl

double _f_(const double * p) {
  return DIM*4.0*PI*PI*sin(2*PI*p[0])*cos(2*PI*p[1]);
}

double _u_b_(const double * p) {
  return sin(2*PI*p[0])*cos(2*PI*p[1]);
}

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef FEMSpace<double,DIM,DOW> fe_space_t;  
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

  Epetra_MpiComm comm(forest.communicator());
  Epetra_Map map((int)global_index.n_global_dof(), global_index.n_primary_dof(), 0, comm);
  global_index.build_epetra_map(map);

  /// 构造 Epetra 的分布式稀疏矩阵模板
  std::cout << "Build sparsity pattern ... " << std::flush;
  Epetra_FECrsGraph G(Copy, map, 10);
  fe_space_t::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    u_int n_ele_dof = ele_dof.size();

    /**
     * 建立从局部自由度数组到全局自由度数组的映射表，这是实现分布式并行
     * 状态下的数据结构的关键一步。
     */
    std::vector<int> indices(n_ele_dof);
    for (u_int i = 0;i < n_ele_dof;++ i) {
      indices[i] = global_index(ele_dof[i]);
    }
    G.InsertGlobalIndices(n_ele_dof, &indices[0], n_ele_dof, &indices[0]);
  }
  G.GlobalAssemble();
  std::cout << "OK!" << std::endl;

  /// 准备构造 Epetra 的分布式稀疏矩阵和计算分布式右端项
  std::cout << "Build sparse matrix ... " << std::flush;
  Epetra_FECrsMatrix A(Copy, G);
  Epetra_FEVector b(map);
  the_ele = fem_space.beginElement();
  for (;the_ele != end_ele;++ the_ele) {
    double vol = the_ele->templateElement().volume();
    const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(5);
    std::vector<AFEPack::Point<DIM> > q_pnt = the_ele->local_to_global(qi.quadraturePoint());
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
          ele_mat(i, j) += JxW*(innerProduct(bas_grad[i][l], bas_grad[j][l]));
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
    A.SumIntoGlobalValues(n_ele_dof, &indices[0], n_ele_dof, &indices[0], &ele_mat(0,0));
    b.SumIntoGlobalValues(n_ele_dof, &indices[0], &ele_rhs(0));
  }
  A.GlobalAssemble();
  b.GlobalAssemble();
  std::cout << "OK!" << std::endl;

  /// 准备解向量。
  Epetra_FEVector x(map);

  /// 加上狄氏边值条件
  u_int n_bnd_dof = 0; /// 首先清点边界上自由度的个数
  for (u_int i = 0;i < fem_space.n_dof();++ i) {
    if (fem_space.dofBoundaryMark(i) > 0) {
      /// 如果不是在主几何体上就不做
      if (! global_index.is_dof_on_minimal_geometry(i)) continue;

      n_bnd_dof += 1;
    }
  }

  /// 准备空间存储边界上全局标号、自变量和右端项
  std::vector<int> bnd_idx(n_bnd_dof);
  std::vector<double> x_entry(n_bnd_dof), rhs_entry(n_bnd_dof);

  /// 对自由度做循环
  for (u_int i = 0, j = 0;i < fem_space.n_dof();++ i) {
    if (fem_space.dofBoundaryMark(i) > 0) { /// 边界上的自由度?
      /// 如果不是在主几何体上就不做
      if (! global_index.is_dof_on_minimal_geometry(i)) continue;

      const int& idx = global_index(i); /// 行的全局标号
      bnd_idx[j] = idx; 

      /// 修改矩阵
      int lrid = A.LRID(idx);
      int row_nnz, *row_idx;
      double *row_entry, row_diag;
      A.ExtractMyRowView(lrid, row_nnz, row_entry, row_idx); /// 取出矩阵的行
      for (int k = 0;k < row_nnz;++ k) { /// 对矩阵的行进行修改
        if (A.LCID(row_idx[k]) != lrid) {   /// 如果不是对角元
          row_entry[k] = 0.0;      /// 则将矩阵元素清零
        } else {                   /// 而对角元保持不变
          row_diag = row_entry[k]; /// 并记录下对角元
        }
      }

      /// 计算并记下自变量和右端项，假设自由度值为插值量
      double u_b_val = _u_b_(fem_space.dofInfo(i).interp_point);
      x_entry[j] = u_b_val;
      rhs_entry[j] = row_diag*u_b_val;

      j += 1;
    }
  } 
  std::cout << "# DOF on the boundary: " << n_bnd_dof << std::endl;

  /// 修改解变量和右端项
  x.ReplaceGlobalValues(n_bnd_dof, &bnd_idx[0], &x_entry[0]);
  b.ReplaceGlobalValues(n_bnd_dof, &bnd_idx[0], &rhs_entry[0]);

  /// 调用 AztecOO 的求解器。
  std::cout << "Solving the linear system ..." << std::flush;
  Epetra_LinearProblem problem(&A, &x, &b);
  AztecOO solver(problem);
  ML_Epetra::MultiLevelPreconditioner precond(A, true);
  solver.SetPrecOperator(&precond);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 100);
  solver.Iterate(5000, 1.0e-12);
  std::cout << "OK!" << std::endl;

  Epetra_Map fe_map(-1, global_index.n_local_dof(), &global_index(0), 0, comm);
  FEMFunction<double,DIM> u_h(fem_space);
  Epetra_Import importer(fe_map, map);
  Epetra_Vector X(View, fe_map, &u_h(0));
  X.Import(x, importer, Add);

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

