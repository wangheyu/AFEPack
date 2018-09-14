/**
 * FEMSpace.vector_value.templates.h : by R.Lie
 * 
 */

#ifndef _FEMSpace_vector_value_templates_h_
#define _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type        nVector<vector_length, double>
#define vector_zero       value_type(0.0)
#define Number            double

template <>
value_type FEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const Point<DOW>& p, 
							      const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, k;
  value_type val;
  const std::vector<int>& dof = e.dof();
  std::vector<value_type> basis_value = e.basis_function_value(p);
  j = dof.size();
  for (i = 0, val = vector_zero;i < j;i ++) {
    for (k = 0;k < vector_length;k ++) {
      val[k] += (*this)(dof[i])*basis_value[i][k];
    }
  }
  return val;
}

template <>
std::vector<value_type> 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const Point<DOW>& p,
						      const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, k, l;
  std::vector<value_type> val(DOW, vector_zero);
  const std::vector<int>& dof = e.dof();
  std::vector<std::vector<value_type> > basis_gradient = e.basis_function_gradient(p);
  j = dof.size();
  for (i = 0;i < j;i ++) {
    for (k = 0;k < DOW;k ++) {
      for (l = 0;l < vector_length;l ++) {
	val[k][l] += (*this)(dof[i])*basis_gradient[i][k][l];
      }
    }
  }
  return val;
}

template <>
std::vector<value_type> 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const std::vector<Point<DOW> >& p,
						   const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, i1, j1, k;
  i = p.size();
  std::vector<value_type> val(i, vector_zero);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  std::vector<std::vector<value_type> > basis_value = e.basis_function_value(p);
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++) {
      for (k = 0;k < vector_length;k ++) {
	val[i1][k] += (*this)(dof[j1])*basis_value[j1][i1][k];
      }
    }
  }
  return val;
}

template <>
std::vector<std::vector<value_type> > 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const std::vector<Point<DOW> >& p,
						      const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, i1, j1, l, m;
  i = p.size();
  std::vector<std::vector<value_type> > val(i, std::vector<value_type>(DOW, vector_zero));
  const std::vector<int>& dof = e.dof();
  std::vector<std::vector<std::vector<value_type> > > basis_gradient = e.basis_function_gradient(p);
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++) {
      for (m = 0;m < DOW;m ++) {
	for (l = 0;l < vector_length;l ++) {
	  val[i1][m][l] += (*this)(dof[j1])*basis_gradient[j1][i1][m][l];
	}
      }
    }
  }
  return val;
}


/**
 * for optimization. Assume the basis function value and gradient are known,
 * the following functions will speedup the calculation to get the value and
 * gradient of the finite element function.
 */
template <>
value_type 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const std::vector<value_type>& basis_value, 
						   const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, l;
  value_type val;
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i = 0, val = vector_zero;i < j;i ++) {
    for (l = 0;l < vector_length;l ++) {
      val[l] += (*this)(dof[i])*basis_value[i][l];
    }
  }
  return val;
}

template <>
std::vector<value_type> 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const std::vector<std::vector<value_type> >& basis_gradient,
						      const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, k, l;
  std::vector<value_type> val(DOW, vector_zero);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i = 0;i < j;i ++) {
    for (k = 0;k < DOW;k ++) {
      for (l = 0;l < vector_length;l ++) {
	val[k][l] += (*this)(dof[i])*basis_gradient[i][k][l];
      }
    }
  }
  return val;
}

template <>
std::vector<value_type> 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const std::vector<std::vector<value_type> >& basis_value,
						   const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, i1, j1, l;
  i = basis_value[0].size();
  std::vector<value_type> val(i, vector_zero);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++) {
      for (l = 0;l < vector_length;l ++) {
	val[i1][l] += (*this)(dof[j1])*basis_value[j1][i1][l];
      }
    }
  }
  return val;
}

template <>
std::vector<std::vector<value_type> > 
FEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const std::vector<std::vector<std::vector<value_type> > >& basis_gradient,
						      const Element<value_type,DIM,DOW,TDIM>& e) const
{
  int i, j, i1, j1, l, m;
  i = basis_gradient[0].size();
  std::vector<std::vector<value_type> > val(i, std::vector<value_type>(DIM, vector_zero));
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++) {
      for (m = 0;m < DOW;m ++) {
	for (l = 0;l < vector_length;l ++) {
	  val[i1][m][l] += (*this)(dof[j1])*basis_gradient[j1][i1][m][l];
	}
      }
    }
  }
  return val;
}

#undef value_type
#undef vector_zero
#undef Number

AFEPACK_CLOSE_NAMESPACE

#endif // _FEMSpace_vector_value_templates_h_

/**
 * end of file
 * 
 */
