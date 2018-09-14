/**
 * @file   AMGSolver.cpp
 * @author Ruo Li
 * @date   Sun Mar 13 16:45:52 2002
 * 
 * @brief  algebraic multigrid solver for positive defined system.
 * 
 * 
 */

#include "AMGSolver.h"
#include "assert.h"

AFEPACK_OPEN_NAMESPACE

typedef AMGSolver::Matrix Matrix;

void AMGSolver::Project(const Matrix& M, 
			Matrix *& P, 
			Matrix *& PMPt)
{
  Assert(M.m() == M.n(), ExcDimensionMismatch(M.m(), M.n()));
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();
  /**
   * choose the core points and influenced points. the procedure is divided into
   * two steps: the first step is to choose the dependent and influenced points.
   * 
   */
  u_int n_core_point = 0;
  std::vector<int> flag(M.m(), 0); /**
				    *  0: candidate to be C-point
				    *  1: F-point
				    * -1: C-point
				    */
  std::vector<double> beta(M.m(), 0.0);
  for (u_int i = 0;i < M.m();i ++) {
    /**
     * calculate the threshhold to construct the influence set
     * 
     */
    beta[i] = 0.0;
    for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++)
      beta[i] = std::max(fabs(M.global_entry(j)), beta[i]);
    beta[i] *= alpha;
    if (flag[i] != 0) continue; /** if this point is not a candidate */
    /**
     * label the F-points around this C-point
     * 
     */
    for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++) {
      const u_int& k = M_colnums[j];
      if (flag[k] == -1) continue;
      if (fabs(M.global_entry(j)) < beta[i]) continue; //!
      flag[k] = 1;
    }
    /**
     * label this point as C-point
     * 
     */
    flag[i] = -1;
    n_core_point ++;
  }
  /**
   * the second step is to check if there are influenced points without common
   * dependent point. if such case occured, we add one of the point to be a
   * dependent point.
   * 
   */
  for (u_int i = 0;i < M.m();i ++) {
    if (flag[i] != 1) continue; /** if this point is not a F-point */
    /**
     * calculate the threshhold to construct the influence set
     * and collect the C-points around this point
     * 
     */
    std::list<u_int> dep_pnt;
    for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++) {
      const u_int& k = M_colnums[j];
      if (fabs(M.global_entry(j)) < beta[i]) continue; //!
      if (flag[k] == -1) dep_pnt.push_back(k);
    }
    /**
     * check its connected points if there are F-F connection
     * 
     */
    bool is_FF = false;
    for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++) {
      const u_int& k = M_colnums[j];
      if (   (flag[k] == -1)                     /** this is a C-point */
	  || (k < i)                             /** this case should have been handled */
	  || (fabs(M.global_entry(j)) < beta[i]) /** this is a weak connection */
	 ) continue;
      /**
       * the point selected is now strongly connected to this point and
       * it is a F-point. we check if there are common entries between
       * the C-points around this point and the point selected. if found, 
       * this is not a F-F connection; otherwise, it is.
       * 
       */
      is_FF = true;
      for (u_int j1 = M_rowstart[k] + 1;j1 < M_rowstart[k + 1];j1 ++) {
	const u_int& k1 = M_colnums[j1];
	if (flag[k1] != -1) continue;
	std::list<u_int>::iterator result = 
	  std::find(dep_pnt.begin(), dep_pnt.end(), k1);
	if (result != dep_pnt.end()) {
	  is_FF = false;
	  break;
	}
      }
      if (is_FF) break;
    }
    /**
     * if there is a F-F connection, we label this point to be C-point
     * 
     */
    if (!is_FF) continue;
    flag[i] = -1;
    n_core_point ++;
  }
  /**
   * then collect the information to construct the project matrix
   * 
   */
  std::vector<u_int> core_point(n_core_point);
  std::vector<int> core_at(M.m(), -1); //! 
  std::vector<u_int> row_length(n_core_point, 0);
  std::vector<u_int> col_length(M.m(), 0);
  std::vector<double> sigma(M.m(), 0.0);
  for (u_int i = 0, n = 0;i < M.m();i ++) {
    if (flag[i] != -1) continue;
    Assert(flag[i] == -1, ExcInternalError());
    core_point[n] = i;
    core_at[i] = n; //! 
    for (u_int j = M_rowstart[i];j < M_rowstart[i+1];j ++) {
      const u_int& k = M_colnums[j];
      row_length[n] ++;
      col_length[k] ++;
      sigma[k] += fabs(M.global_entry(j));
    }
    n ++;
  }
  /**
   * construct the sparsity pattern for the project matrix
   * 
   */
  SparsityPattern& spP = *(new SparsityPattern(n_core_point, M.m(), row_length));
  SparsityPattern spPt = SparsityPattern(M.m(), n_core_point, col_length);
  for (u_int i = 0;i < n_core_point;i ++) {
    const u_int& k = core_point[i];
    spP.add(i, k);
    spPt.add(k, i);
    for (u_int j = M_rowstart[k] + 1;j < M_rowstart[k+1];j ++) {
      const u_int& l = M_colnums[j];
      if (flag[l] == -1) continue;
      spP.add(i, l);
      spPt.add(l, i);
    }
  }
  spP.compress();
  spPt.compress();
  /**
   * fill the entries for the project matrix
   * 
   */
  const std::size_t * P_rowstart = spP.get_rowstart_indices();
  const u_int * P_colnums = spP.get_column_numbers();
  const std::size_t * Pt_rowstart = spPt.get_rowstart_indices();
  const u_int * Pt_colnums = spPt.get_column_numbers();
  P = new Matrix(spP);
  Matrix Pt(spPt);
  for (u_int i = 0;i < M.m(); i ++) {
    const int &k = core_at[i];
    if( k >= 0 ) {
      Pt.add(i, k, 1.0);
      P->add(k, i, 1.0);
    }
    else {
      double aii = M.diag_element(i);
      for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++) {
	if (fabs(M.global_entry(j)) < beta[i])
	  aii += M.global_entry(j);
      }
      const u_int * dep_begin = &Pt_colnums[Pt_rowstart[i]];
      const u_int * dep_end = &Pt_colnums[Pt_rowstart[i + 1]];

      for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++) {
	const u_int &l = M_colnums[j];
	if (fabs(M.global_entry(j)) < beta[i]) 
	  continue;	
	if (core_at[l] >= 0)  {
	  Pt.add(i, core_at[l], fabs(M.global_entry(j)/aii) );
	  P->add(core_at[l], i, fabs(M.global_entry(j)/aii) );
	}
	else {
	  double d = 0.0;
	  for (u_int i1 = Pt_rowstart[l];i1 < Pt_rowstart[l + 1];i1 ++) {
	    const u_int &c1 = Pt_colnums[i1];
	    if (std::find(dep_begin, dep_end, c1) != dep_end) {
	      d += M.el(l, core_point[c1]);
	    }
	  }
	  if(fabs(d*aii) < 1e-16) {
	    //std::cerr << " F-F assert failure!" << std::endl;
	    continue;  //! Maybe error, F-F has no common dependence
	  }
	  for (u_int i1 = Pt_rowstart[l];i1 < Pt_rowstart[l + 1];i1 ++) {
	    const u_int &c1 = Pt_colnums[i1];
	    if (std::find(dep_begin, dep_end, c1) != dep_end) {
	      Pt.add(i, c1, fabs(M.el(i, l)*M.el(l, core_point[c1])/(d*aii)));
	      P->add(c1, i, fabs(M.el(i, l)*M.el(l, core_point[c1])/(d*aii)));
	    }
	  }
	}
      }
    }
  }

  /**
   * prepare data for constructing the projected matrix
   * 
   */
  PMPt = getPMPt(*P, M, Pt);
}



/**
 * This is the Project subroutine for fast project to solve without high
 * accuracy requirement.
 * 
 */
void AMGSolver::lazyProject(const Matrix& M, 
			    Matrix *& P,
			    Matrix *& PMPt,
			    Matrix *& Pt)
{
  std::vector<u_int> counter(M.m());
  std::vector<u_int> core_point(M.m());
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();

  u_int n_core_point = 0;
  const std::size_t * rowstart_indices = spM.get_rowstart_indices();
  for (u_int i = 0;i < M.m();i ++) {
    if (counter[i] > 0) continue; /**< this row is handled already */
    for (u_int j = M_rowstart[i] + 1;j < M_rowstart[i + 1];j ++)
      counter[M_colnums[j]] ++;
    core_point[n_core_point ++] = i;
  }

  SparsityPattern& spP = *new SparsityPattern(n_core_point, M.n(), spM.max_entries_per_row());
  SparsityPattern& spPt = *new SparsityPattern(M.n(), n_core_point, spM.max_entries_per_row());
  for (u_int i = 0;i < n_core_point;i ++) {
    const u_int& k = core_point[i];
    spP.add(i, k);
    spPt.add(k, i);
    for (u_int j = M_rowstart[k] + 1;j < M_rowstart[k + 1];j ++) {
      const u_int& l = M_colnums[j];
      spP.add(i, l);
      spPt.add(l, i);
    }
  }
  spP.compress();
  spPt.compress();

  P = new Matrix(spP);
  Pt = new Matrix(spPt);
  for (u_int i = 0;i < n_core_point;i ++) {
    const u_int& k = core_point[i];
    P->add(i, k, 1.0);
    Pt->add(k, i, 1.0);
    for (u_int j = M_rowstart[k] + 1;j < M_rowstart[k + 1];j ++) {
      const u_int& l = M_colnums[j];
      P->add(i, l, 1.0/counter[l]);
      Pt->add(l, i, 1.0/counter[l]);
    }
  }

  PMPt = getPMPt(*P, M, *Pt);
}



Matrix * AMGSolver::getPMPt(const Matrix& P, 
			    const Matrix& M,
			    const Matrix& Pt) const
{
  const SparsityPattern& spP = P.get_sparsity_pattern();
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const SparsityPattern& spPt = Pt.get_sparsity_pattern();
  const std::size_t * P_rowstart = spP.get_rowstart_indices();
  const u_int * P_colnums = spP.get_column_numbers();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();
  const std::size_t * Pt_rowstart = spPt.get_rowstart_indices();
  const u_int * Pt_colnums = spPt.get_column_numbers();

  std::vector<u_int> row_length(P.m(), 0);
  std::vector<bool> flag(P.m(), true);
  std::vector<u_int> index(P.m());
  std::vector<std::vector<u_int> > col_index(P.m());
  for (u_int i = 0;i < P.m();i ++) {
    row_length[i] = 1; /**< add the diagonal entry at first */
    flag[i] = false;
    index[0] = i;
    for (u_int j = P_rowstart[i];j < P_rowstart[i+1];j ++) {
      const u_int& a = P_colnums[j];
      for (u_int k = M_rowstart[a];k < M_rowstart[a+1];k ++) {
	const u_int& b = M_colnums[k];
	for (u_int l = Pt_rowstart[b];l < Pt_rowstart[b+1];l ++) {
	  const u_int&  c = Pt_colnums[l];
	  if (flag[c]) {
	    index[row_length[i] ++] = c;
	    flag[c] = false;
	  }
	}
      }
    }
    col_index[i].resize(row_length[i]);
    for (u_int j = 0;j < row_length[i];j ++) {
      col_index[i][j] = index[j];
      flag[index[j]] = true;
    }
  }

  SparsityPattern& spA = *(new SparsityPattern(P.m(), row_length));
  for (u_int i = 0;i < P.m();i ++) {
    for (u_int j = 0;j < row_length[i];j ++) {
      spA.add(i, col_index[i][j]);
    }
  }
  spA.compress();

  Matrix *A = new Matrix(spA);
  std::vector<double> row_entry(P.m(), 0.0);
  for (u_int i = 0;i < P.m();i ++) {
    for (u_int j = P_rowstart[i];j < P_rowstart[i+1];j ++) {
      const u_int& a = P_colnums[j];
      for (u_int k = M_rowstart[a];k < M_rowstart[a+1];k ++) {
	const u_int& b = M_colnums[k];
	for (u_int l = Pt_rowstart[b];l < Pt_rowstart[b+1];l ++) {
	  const u_int& c = Pt_colnums[l];
	  row_entry[c] += P.global_entry(j)*M.global_entry(k)*Pt.global_entry(l);
	}
      }
    }
    for (u_int j = 0;j < row_length[i];j ++) {
      A->add(i, col_index[i][j], row_entry[col_index[i][j]]);
      row_entry[col_index[i][j]] = 0.0;
    }
  }

  return A;
}



void AMGSolver::GaussSidel(const Matrix& M, 
			   Vector<double>& x, 
			   const Vector<double>& r,
			   const int& s) const
{
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * rowstart = spM.get_rowstart_indices();
  const u_int * colnums = spM.get_column_numbers();
  for (u_int i = 0;i < s;i ++) {
    for (u_int j = 0;j < M.m();j ++) {
      double r0 = r(j);
      for (u_int k = rowstart[j] + 1;k < rowstart[j+1];k ++) {
	r0 -= x(colnums[k]) * M.global_entry(k);
      }
      x(j) = r0 / M.global_entry(rowstart[j]);
    }
  }
}


void AMGSolver::reinit(const Matrix& M)
{
  clear();
	
  std::cerr << "AMGSolver initializing ..." << std::flush;
  const SparsityPattern& spM = M.get_sparsity_pattern();
  Assert(spM.is_compressed(), SparsityPattern::ExcNotCompressed());
  Assert(M.m() == M.n(), ExcDimensionMismatch(M.m(), M.n()));
  u_int k, m;
  k = m = M.m();
  float sparse_degree = M.n_nonzero_elements()/(float(m)*m);
  projected_matrix.push_back(&M);
  k = m = M.m();
  n_project = 0;
  is_most_project_full = true;
  while(m > int(n_most_project_order) || sparse_degree < n_most_project_sparse_degree) {
    Matrix *PMPt, *P;
    Project(*projected_matrix[n_project], P, PMPt);
    m = PMPt->m();
    if (k == m) {
      is_most_project_full = false;
      const SparsityPattern& spP = P->get_sparsity_pattern();
      delete P; delete (&spP);
      const SparsityPattern& spM = PMPt->get_sparsity_pattern();
      delete PMPt; delete (&spM);
      break;
    }
    k = m;
    project_matrix.push_back(P);
    projected_matrix.push_back(PMPt);
    sparse_degree = PMPt->n_nonzero_elements()/(float(m)*m);
    n_project ++;
#ifdef DEBUG
    std::cerr << "\tlevel " << n_project - 1
	      << " with order " << m
	      << " and sparse degree " << sparse_degree
	      << std::endl;
#endif // DEGUG
  };
  if (is_most_project_full) {
    M_n.reinit(projected_matrix[n_project]->m(), projected_matrix[n_project]->n());
    const std::size_t * row_start = projected_matrix[n_project]
      ->get_sparsity_pattern().get_rowstart_indices();
    const u_int * col_index = projected_matrix[n_project]
      ->get_sparsity_pattern().get_column_numbers();
    for (k = 0;k < projected_matrix[n_project]->m();k ++) {
      for (m = row_start[k];m < row_start[k + 1];m ++) {
	M_n(k, col_index[m]) = projected_matrix[n_project]->global_entry(m);
      }
    }
    M_n.gauss_jordan();
  }
  is_initialized = true;
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



void AMGSolver::lazyInit(const Matrix& M)
{
  std::cerr << "AMGSolver initializing in lazy mode ..." << std::flush;
  const SparsityPattern& spM = M.get_sparsity_pattern();
  Assert(spM.is_compressed(), SparsityPattern::ExcNotCompressed());
  Assert(M.m() == M.n(), ExcDimensionMismatch(M.m(), M.n()));
  u_int k, m;
  k = m = M.m();
  float sparse_degree = M.n_nonzero_elements()/(float(m)*m);
  projected_matrix.push_back(&M);
  n_project = 0;
  is_most_project_full = true;
  while(m > n_most_project_order || sparse_degree < n_most_project_sparse_degree) {
    Matrix *PMPt, *P, *Pt;
    lazyProject(*projected_matrix[n_project], P, PMPt, Pt);
    m = PMPt->m();
    if (k == m) {
      is_most_project_full = false;
      const SparsityPattern& spP = P->get_sparsity_pattern();
      delete P; delete (&spP);
      const SparsityPattern& spPt = Pt->get_sparsity_pattern();
      delete Pt; delete (&spPt);
      const SparsityPattern& spM = PMPt->get_sparsity_pattern();
      delete PMPt; delete (&spM);
      break;
    }
    k = m;
    project_matrix.push_back(P);
    project_matrix_r.push_back(Pt);
    projected_matrix.push_back(PMPt);
    sparse_degree = PMPt->n_nonzero_elements()/(float(m)*m);
    n_project ++;
#ifdef DEBUG
    std::cerr << "\tlevel " << n_project - 1
	      << " with order " << m
	      << " and sparse degree " << sparse_degree
	      << std::endl;
#endif // DEGUG
  };
  if (is_most_project_full) {
    M_n.reinit(projected_matrix[n_project]->m(), projected_matrix[n_project]->n());
    const std::size_t * row_start = projected_matrix[n_project]
      ->get_sparsity_pattern().get_rowstart_indices();
    const u_int * col_index = projected_matrix[n_project]
      ->get_sparsity_pattern().get_column_numbers();
    for (k = 0;k < projected_matrix[n_project]->m();k ++) {
      for (m = row_start[k];m < row_start[k + 1];m ++) {
	M_n(k, col_index[m]) = projected_matrix[n_project]->global_entry(m);
      }
    }
    M_n.gauss_jordan();
  }
  is_initialized = true;
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



/**
 * to solve a linear system with the same sparsity pattern as already initialized
 * matrix
 * 
 */
void AMGSolver::lazyReinit(const Matrix& M)
{
  const SparsityPattern& spM = M.get_sparsity_pattern();
  if (is_initialized) {
    if (&(projected_matrix[0]->get_sparsity_pattern()) != &spM) {
      std::cerr << "The solver is not initialized with the same sparsity pattern." 
		<< std::endl;
      abort();
    }
  }
  else {
    lazyInit(M);
    return;
  }

  int i, j, k, m;
  std::cerr << "AMGSolver reinitializing in lazy mode ..." << std::flush;
  projected_matrix[0] = &M;
  for (j = 0;j < n_project;j ++) {
    const SparsityPattern& spP = projected_matrix[j + 1]->get_sparsity_pattern();
    delete projected_matrix[j + 1];
    delete (&spP);
    projected_matrix[j + 1] = getPMPt(*project_matrix[j], 
				      *projected_matrix[j], 
				      *project_matrix_r[j]);
  };
  if (is_solve_most_project_exactly && is_most_project_full) {
    M_n.reinit(projected_matrix[n_project]->m(), projected_matrix[n_project]->n());
    const std::size_t * row_start = projected_matrix[n_project]
      ->get_sparsity_pattern().get_rowstart_indices();
    const u_int * col_index = projected_matrix[n_project]
      ->get_sparsity_pattern().get_column_numbers();
    for (i = 0;i < projected_matrix[n_project]->m();i ++) {
      for (j = row_start[i];j < row_start[i + 1];j ++) {
	M_n(i, col_index[j]) = projected_matrix[n_project]->global_entry(j);
      }
    }
    M_n.gauss_jordan();
  }
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



void AMGSolver::solve(Vector<double>& x, 
		      const Vector<double>& r, 
		      double tol, 
		      u_int s1,
		      int mode) const
{
  Assert(is_initialized == true, ExcNotInitialized());
  Assert(projected_matrix[0]->m() == r.size(), 
         ExcDimensionMismatch(projected_matrix[0]->m(), r.size()));
  Assert(x.size() == r.size(), 
	 ExcDimensionMismatch(x.size(), r.size()));
  if (tol == 0.0) tol = tolerence();
  Vector<double> r1(r), v0, v1, v2(r.size());
  projected_matrix[0]->vmult(v2, x);
  v2.add(-1.0, r);
  double residual = v2.l1_norm();
  double init_residual = residual;
  if ((mode&0x0001) == 0) { // mode == 0 is the solve mode, mode == 1 is the precondition mode
    std::cerr << "AMGSolver begin with initial residual " << residual << " ..." << std::endl;
  }
  u_int iteration_step = 0;
  std::vector<Vector<double> *> projected_r(n_project + 1);
  std::vector<Vector<double> *> projected_x(n_project + 1);
  projected_x[0] = &x;
  projected_r[0] = &r1;
  for (u_int i = 1;i <= n_project;i ++) {
    projected_r[i] = new Vector<double>(projected_matrix[i]->m());
    projected_x[i] = new Vector<double>(projected_matrix[i]->m());
  }
  while(residual >= tol*init_residual) {
    for (int i = 0;i < n_project;i ++) {
      GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
      v0.reinit(projected_x[i]->size());
      projected_matrix[i]->vmult(v0, *projected_x[i]);
      v1 = *projected_r[i];
      v1 -= v0;
      project_matrix[i]->vmult(*projected_r[i+1], v1);
      (*projected_x[i+1]) = 0;
    }
		
    if (is_solve_most_project_exactly) {
      if (is_most_project_full) {
        M_n.vmult(*projected_x[n_project], *projected_r[n_project]);
      }
      else {
        for (u_int j = 0;j < projected_matrix[n_project]->m();j ++) {
          (*projected_x[n_project])(j) = (*projected_r[n_project])(j)
            / projected_matrix[n_project]->diag_element(j);
        }
      }
    } else {
      GaussSidel(*projected_matrix[n_project],
                 *projected_x[n_project], 
                 *projected_r[n_project],
                 smooth_step);
    }
		
    for (int i = n_project-1;i >= 0;i --) {
      v0.reinit(projected_x[i]->size());
      project_matrix[i]->Tvmult(v0, *projected_x[i+1]);
      (*projected_x[i]) += v0;
      GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], smooth_step);
    }
		
    projected_matrix[0]->vmult(v2, x);
    v2.add(-1.0, r);
    residual = v2.l1_norm();
    iteration_step ++;
    if (iteration_step > s1) break;
    if ((mode&0x0002) != 0) {
      std::cerr << "\r\tresidual = " << residual << std::flush;
    }
  };
  for (u_int i = 1;i <= n_project;i ++) {
    delete projected_r[i];
    delete projected_x[i];
  }
  if ((mode&0x0001) == 0) {
    if (residual < tol*init_residual) {
      std::cerr << "\r\tconverge with residual " << residual
		<< " at step " << iteration_step << "." << std::endl;
    }
    else {
      std::cerr << "\r\tfailed to converge with residual " << residual
		<< " at step " << iteration_step << "." << std::endl;
    }
  }
}



void AMGSolver::clear()
{	
  if (!is_initialized) return;
  for (u_int i = 0;i < n_project;i ++) {
    const SparsityPattern& spP = project_matrix[i]->get_sparsity_pattern();
    delete project_matrix[i];
    delete (&spP);
    const SparsityPattern& spM = projected_matrix[i+1]->get_sparsity_pattern();
    delete projected_matrix[i+1];
    delete (&spM);
  }
  project_matrix.clear();
  projected_matrix.clear();

  if (project_matrix_r.size() > 0) {
    for (u_int i = 0;i < n_project;i ++) {
      const SparsityPattern& spP = project_matrix_r[i]->get_sparsity_pattern();
      delete project_matrix_r[i];
      delete (&spP);
    }
    project_matrix_r.clear();
  }
  is_initialized = false;
}



AMGSolver::AMGSolver() : 
  is_initialized(false),
  toler(1.e-12),
  n_project(0),
  smooth_step(3),
  n_most_project_order(50),
  n_most_project_sparse_degree(0.382),
  alpha(0.25),
  is_solve_most_project_exactly(true)
{}



AMGSolver::AMGSolver(const Matrix& m, 
		     double tol, 
		     u_int s,
		     u_int nmpo,
		     double nmpsd,
		     double alp) :
  is_initialized(false),
  toler(tol),
  n_project(0),
  smooth_step(s),
  n_most_project_order(nmpo),
  n_most_project_sparse_degree(nmpsd),
  alpha(alp),
  is_solve_most_project_exactly(true)
{
  reinit(m);
}



AMGSolver::~AMGSolver()
{
  clear ();
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

