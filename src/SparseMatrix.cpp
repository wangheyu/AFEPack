#include <AFEPack/SparseMatrix.h>

AFEPACK_OPEN_NAMESPACE

template <typename T>
SparseMatrix<T>::SparseMatrix(SparsityPattern& _sp)
{
  initialized = 1;
  sparsity_pattern = &_sp;

  std::vector<T>::resize(sparsity_pattern->get_n_nonzeros(), 0.0);
}

template <typename T>
SparseMatrix<T>::SparseMatrix(const SparsityPattern& _sp)
{
  initialized = 1;
  //sparsity_pattern = const_cast<SparsityPattern*>(&_sp);
  sparsity_pattern = &_sp;
  std::vector<T>::resize(sparsity_pattern->get_n_nonzeros(), 0.0);
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator = (const T s)
{
  std::vector<T>::resize(sparsity_pattern->get_n_nonzeros(), s);
  return *this;
}

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator = (const SparseMatrix<T>& src)
{
  copy_from(src);
  return *this;
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::copy_from (const SparseMatrix<T> &src)
{
  assert(("The size is not match", n_nonzero_elements() == src.n_nonzero_elements()));
  initialized = 1;
  std::copy (src.get_data(), &(src.get_data()[n_nonzero_elements()-1]), get_data_impl());
  return *this;
}


template <typename T>
void SparseMatrix<T>::reinit(SparsityPattern& _sp)
{
  sparsity_pattern = &_sp;

  std::vector<T>::resize(sparsity_pattern->get_n_nonzeros(), 0.0);
}

template <typename T>
void SparseMatrix<T>::reinit(const SparsityPattern& _sp) 
{
  sparsity_pattern = &_sp;

  std::vector<T>::resize(sparsity_pattern->get_n_nonzeros(), 0.0);
}

template <typename T>
void SparseMatrix<T>::clear_impl()
{
  is_square = false;
  sparsity_pattern = NULL;
  initialized = 0;

  std::vector<T>::clear();
}

template <typename T>
//SparsityPattern* SparseMatrix<T>::get_sparsity_pattern() // original vecma version;
const SparsityPattern* SparseMatrix<T>::get_sparsity_pattern_ptr() const 
{
  return sparsity_pattern;
}

template <typename T>
const SparsityPattern& SparseMatrix<T>::get_sparsity_pattern() const
{
  return *sparsity_pattern;
}

template <typename T>
void SparseMatrix<T>::set_entry_impl(const unsigned int& row_idx,
			  const unsigned int& col_idx,
			  const T& entry)
{
  const std::size_t * row = sparsity_pattern->get_row();
  const unsigned int * col = sparsity_pattern->get_col();

  for(std::size_t i = row[row_idx];i < row[row_idx + 1];++ i){
    if(col[i] == col_idx){
      std::vector<T>::operator[](i) = entry;
      return;
    }
  }

  throw "Something wrong with SparseMatrix::add function";
}

template <typename T>
void SparseMatrix<T>::add_entry_impl(const unsigned int& row_idx,
			  const unsigned int& col_idx,
			  const T& entry)
{
  const std::size_t * row = sparsity_pattern->get_row();
  const unsigned int * col = sparsity_pattern->get_col();

  for(unsigned int i = row[row_idx];i < row[row_idx + 1];++ i){
    if(col[i] == col_idx){
      std::vector<T>::operator[](i) += entry;
      return;
    }
  }

  throw "Something wrong with SparseMatrix::add function";
}


template <typename T>
void SparseMatrix<T>::print_sparseMatrix()
{
  sparsity_pattern->print();

  const unsigned int n_nonzeros = sparsity_pattern->get_n_nonzeros();
  for(unsigned int i = 0;i < n_nonzeros;++ i){
    std::cout << std::vector<T>::operator[](i) << " ";
  }
  std::cout << std::endl;
}

template <typename T>
void SparseMatrix<T>::print_impl()
{
  unsigned int m = sparsity_pattern->get_m();
  unsigned int n = sparsity_pattern->get_n();
  const std::size_t * row = sparsity_pattern->get_row();
  const unsigned int * col = sparsity_pattern->get_col();  

  std::vector<T> row_vector(n, 0.0);    
  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = row[i];j < row[i + 1];++ j){
      row_vector[col[j]] = std::vector<T>::operator[](j);
    }

    std::cout << "| ";
    for(unsigned int j = 0;j < n;++ j){
      std::cout << row_vector[j] << " ";
    }
    std::cout << " |" <<std::endl;

    for(unsigned int j = row[i];j < row[i + 1];++ j){
      row_vector[col[j]] = 0.;
    }
    
  }
}

template <typename T>
T& SparseMatrix<T>::global_entry(const unsigned int& _n)
{
  return (*this)[_n];
}

template <typename T>
const T& SparseMatrix<T>::global_entry(const unsigned int& _n) const
{
  return (*this)[_n];
}

template <typename T>
T SparseMatrix<T>::el(const unsigned int& _m, const unsigned int& _n) const
{
  const SparsityPattern& sp = get_sparsity_pattern();
  const std::size_t * row_idx = sp.get_rowstart_indices();
  const unsigned int * col_idx = sp.get_column_numbers();

  for(unsigned int i  = row_idx[_m];i < row_idx[_m + 1];++ i){
    if(col_idx[i] == _n){
      return (*this)[i];
    }
  }

  return 0;
}

template <typename T>
T& SparseMatrix<T>::diag_element_impl(const unsigned int& _n)
{
  /// is_square should be true.

  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  return (*this)[row_idx[_n]];
}

template <typename T>
const T& SparseMatrix<T>::diag_element_impl(const unsigned int& _n) const
{
  /// is_square should be true.

  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  return (*this)[row_idx[_n]];
}

/// old version (src, des); M*src = des;
/// new version (des, src)
template <typename T>
//void SparseMatrix<T>::mult_v_impl(Vector<T>&src, Vector<T>&des)
void SparseMatrix<T>::mult_v_impl(Vector<T>&des, Vector<T>&src) const
{
  unsigned int m = sparsity_pattern->get_m();
  unsigned int n = sparsity_pattern->get_n();  
  assert(n == src.size());
  //assert(src.size() == des.size());
  des.reinit(m);
  assert(m == des.size());
  
  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  T * des_data = des.get_data();
  T * src_data = src.get_data();

  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = row_idx[i];j < row_idx[i + 1];++ j){
      des[i] += std::vector<T>::operator[](j) * src[col_idx[j]];
    }
  }
}

template <typename T>
void SparseMatrix<T>::mult_v_impl(Vector<T>&des,const Vector<T>&src) const
{
  unsigned int m = sparsity_pattern->get_m();
  unsigned int n = sparsity_pattern->get_n();  
  assert(n == src.size());
  des.reinit(m);
  assert(m == des.size());
  
  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  T * des_data = des.get_data();
  const T * src_data = src.get_data();

  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = row_idx[i];j < row_idx[i + 1];++ j){
      des[i] += std::vector<T>::operator[](j) * src[col_idx[j]];
    }
  }
}

/// des := alpha * A * src + beta * y,
template <typename T>
void SparseMatrix<T>::mult_v_impl(const T& alpha, Vector<T>& src, const T& beta, Vector<T>& des) const
{
  unsigned int m = sparsity_pattern->get_m();
  unsigned int n = sparsity_pattern->get_n();  
  assert(n == src.size());
  assert(src.size() == des.size());
  
  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  T * des_data = des.get_data();
  const T * src_data = src.get_data();

  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = row_idx[i];j < row_idx[i + 1];++ j){
      des[i] += alpha * std::vector<T>::operator[](j) * src[col_idx[j]] + beta * des[i];
    }
  }  
}

template <typename T>
void SparseMatrix<T>::Tmult_v_impl(Vector<T>&des, Vector<T>&src)
{
  unsigned int m = sparsity_pattern->get_m();
  unsigned int n = sparsity_pattern->get_n();  
  assert(m == src.size());
  //assert(src.size() == des.size());
  des.reinit(n);
  assert(n == des.size());
  
  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();

  T * des_data = des.get_data();
  T * src_data = src.get_data();
  
  for (std::size_t i = 0; i < m ; i++)
    {
      for (unsigned int j=row_idx[i]; j<row_idx[i+1] ; j++)
        {
          const std::size_t p = col_idx[j];
          des[p] += std::vector<T>::operator[](j) * src(i);
        }
    }

}

/// des = alpha * srcA * srcB + beta * des
template <typename T>
void SparseMatrix<T>::mult_m_impl(const T& alpha,
			     SparseMatrix<T>& src,
			     const T& beta,
			     DenseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  assert(get_m_impl() == des.get_m_impl());
  assert(src.get_n_impl() == des.get_n_impl());

  unsigned int m = get_m_impl();
  unsigned int n = src.get_n_impl();
  const std::size_t * row_idx = sparsity_pattern->get_row();
  const unsigned int * col_idx = sparsity_pattern->get_col();
  T* data = std::vector<T>::data();
  T* src_data = src.get_data();
  T* des_data = des.get_data();


  
  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int k = 0;k < n;++ k){    
      for(unsigned int j = row_idx[i];j < row_idx[i + 1];++ j){
	des.add_entry_impl(i, k, alpha * data[j] * src_data[col_idx[j] * n + k] + beta * des_data[col_idx[j] * n + k]);
      }
    }
  }
}

/// The unique case is des = alpha * srcA * srcB, say we need to build
/// a brand-new des. That means it does not matter what beta is.
template <typename T>
void SparseMatrix<T>::mult_m_impl(const T& alpha,
			     SparseMatrix<T>& src,
			     const T& beta,
			     SparseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  //  assert(srcA.get_m() == des.get_m());
  //  assert(srcB.get_n() == des.get_n());

  unsigned int m = get_m_impl();
  unsigned int n = get_n_impl();
  unsigned int m_src = src.get_m_impl();
  unsigned int n_src = src.get_n_impl();  
  
  /// to prepare info
  const std::size_t * rowIdx = sparsity_pattern->get_row();
  const unsigned int * colIdx = sparsity_pattern->get_col();
  T * data = get_data_impl();
  //SparsityPattern * sparsity_pattern_src = src.get_sparsity_pattern_ptr();
  const SparsityPattern * sparsity_pattern_src = src.get_sparsity_pattern_ptr();
  const std::size_t * rowIdx_src = sparsity_pattern_src->get_row();
  const unsigned int * colIdx_src = sparsity_pattern_src->get_col();
  T * data_src = src.get_data_impl();
  
  //const SparsityPattern * des_sp = des.get_sparsity_pattern_ptr();
  

  /// to generate sparsity pattern of des
  std::vector<std::set<unsigned int> > row_nonzeros(m);
  std::vector<unsigned int> n_row_nonzeros(m);  
  
  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = rowIdx[i];j < rowIdx[i + 1];++ j){
      unsigned int col_idx = colIdx[j];
      for(unsigned int k = rowIdx_src[col_idx];k < rowIdx_src[col_idx + 1];++ k){
	row_nonzeros[i].insert(colIdx_src[k]);
      }
    }
    n_row_nonzeros[i] = row_nonzeros[i].size();
  }

  SparsityPattern * des_sp;
  des_sp = new SparsityPattern(m, n_src, n_row_nonzeros);

  std::set<unsigned int>::iterator it;
  for(std::size_t i = 0;i < m;++ i){
    for(it = row_nonzeros[i].begin();it != row_nonzeros[i].end();++ it){
      des_sp->add(i, *it);
    }
  }

  des_sp->compress();

  /// fill entries of res.
  des.reinit(*des_sp);

  for(std::size_t i = 0;i < m;++ i){
    for(unsigned int j = rowIdx[i];j < rowIdx[i + 1];++ j){
      unsigned int col_idx = colIdx[j];
      for(unsigned int k = rowIdx_src[col_idx];k < rowIdx_src[col_idx + 1];++ k){
	des.add_entry_impl(i, colIdx_src[k], data[j] * data_src[k]);
      }
    }
  }   
}


template <typename T>
T* SparseMatrix<T>::get_data_impl()
{
  return std::vector<T>::data();
}

template <typename T>
const T* SparseMatrix<T>::get_data_impl() const
{
  return std::vector<T>::data();
}

template <typename T>
const unsigned int& SparseMatrix<T>::get_m_impl() const
{
  return sparsity_pattern->get_m();
}

template <typename T>
const unsigned int& SparseMatrix<T>::m() const
{
  return sparsity_pattern->get_m();
}

template <typename T>
const unsigned int& SparseMatrix<T>::get_n_impl() const
{
  return sparsity_pattern->get_n();
}

template <typename T>
const unsigned int& SparseMatrix<T>::n() const
{
  return sparsity_pattern->get_n();
}

template <typename T>
//const unsigned int& SparseMatrix<T>::n_nonzero_elements() const
const unsigned int SparseMatrix<T>::n_nonzero_elements() const
{
  return sparsity_pattern->get_n_nonzeros();
}

template class SparseMatrix<double>;
AFEPACK_CLOSE_NAMESPACE
