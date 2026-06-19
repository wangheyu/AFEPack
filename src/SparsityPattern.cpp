#include <AFEPack/SparsityPattern.h>

AFEPACK_OPEN_NAMESPACE

SparsityPattern::SparsityPattern()
{
  ;
}

SparsityPattern::SparsityPattern(const u_int& _m,
				 const u_int& _n,
				 const u_int& _max_entries_per_row)
{
  m = _m;
  n = _n;
  //max_entries_per_row = _max_entries_per_row;
  max_row_length = _max_entries_per_row;

  row.resize(m + 1);
  row[0] = 0;
  for(u_int i = 1;i < m;++ i){
    //row[i] = row[i - 1] + max_entries_per_row;
    row[i] = row[i - 1] + max_row_length;
  }
  //row[m] = m * max_entries_per_row;
  //col.resize(m * max_entries_per_row, -1);
  row[m] = m * max_row_length;
  col.resize(m * max_row_length, -1);
  compressed = false;
  //n_nonzeros = m * max_entries_per_row;
  n_nonzeros = m * max_row_length;
}

SparsityPattern::SparsityPattern(const u_int& _m,
				 const u_int& _n,
				 const std::vector<u_int>& _entries_per_row)
{
  m = _m;
  n = _n;
  //max_entries_per_row = _max_entries_per_row;

  row.resize(m + 1);
  row[0] = 0;
  n_nonzeros = 0;
  for(u_int i = 1;i < m;++ i){
    row[i] = row[i - 1] + _entries_per_row[i - 1];
    n_nonzeros += _entries_per_row[i - 1];
  }
  n_nonzeros += _entries_per_row[m - 1];
  row[m] = n_nonzeros;
  col.resize(n_nonzeros, -1);

  compressed = false;
}

SparsityPattern::SparsityPattern(const u_int& _m,
				 const std::vector<u_int>& _entries_per_row)
{
  reinit(_m, _m, _entries_per_row);
}

SparsityPattern::~SparsityPattern()
{
  ;	  
}

void SparsityPattern::reinit(const u_int& _m,
			     const u_int& _n)
{
  m = _m;
  n = _n;
}

void SparsityPattern::reinit(const u_int& _m,
			     const u_int& _n,
			     const u_int& _max_entries_per_row)
{
  m = _m;
  n = _n;
  //max_entries_per_row = _max_entries_per_row;
  max_row_length = _max_entries_per_row;
  row.resize(m + 1);
  row[0] = 0;
  for(u_int i = 1;i < m;++ i){
    //row[i] = row[i - 1] + max_entries_per_row;
    row[i] = row[i - 1] +max_row_length;
  }
  //row[m] = m * max_entries_per_row;
  row[m] = m * max_row_length;
    //col.resize(m * max_entries_per_row, -1);
  col.resize(m * max_row_length, -1);

  compressed = false;
  // n_nonzeros = m * max_entries_per_row;
  n_nonzeros = m * max_row_length;
}

void SparsityPattern::reinit(const u_int& _m,
			     const u_int& _n,
			     const std::vector<u_int>& _entries_per_row)
{
  m = _m;
  n = _n;
  //max_entries_per_row = _max_entries_per_row;

  row.resize(m + 1);
  row[0] = 0;
  n_nonzeros = 0;
  for(u_int i = 1;i < m;++ i){
    row[i] = row[i - 1] + _entries_per_row[i - 1];
    n_nonzeros += _entries_per_row[i - 1];
  }
  n_nonzeros += _entries_per_row[m - 1];
  row[m] = n_nonzeros;
  col.resize(n_nonzeros, -1);

  compressed = false;
}

const u_int& SparsityPattern::get_m() const
{
  return m;
}

const u_int& SparsityPattern::get_n() const
{
  return n;
}

const u_int SparsityPattern::get_n_nonzeros() const
{
  return n_nonzeros;
}

void SparsityPattern::set_max_entries_per_row(const u_int& _max_entries_per_row)
{
  // max_entries_per_row = _max_entries_per_row;
  max_row_length =  _max_entries_per_row;
}

/// In solving PDEs, the diagnol element always play an important role
/// in GS iteration, etc. Hence, we put the diagnol element at the
/// first poisition in each row, to improve the efficiency

void SparsityPattern::put_diagEle_first_in_each_row()
{
  for(u_int i = 0;i < m;++ i){
    if(i == col[row[i]]) continue;
    
    for(u_int j = row[i] + 1;j < row[i + 1];++ j){
      if(col[j] == i) std::swap(col[row[i]], col[j]);
    }
  }
}


/// Initially, the vector col is filled with -1. During the
/// formation of the sparsity pattern, the number -1 would be
/// replaced by column numbers appeared in the sparse matrix. So in
/// compress(), all element with entry -1 will be removed from the
/// vector col, with the following three steps.
/// 1: using std::remove() to "delete" all elements with entry -1;
/// 2: using vec.erase() to change the size of vector
/// 3: using vec.shrink_to_fit() to change the capacity of the vector
/// Complexity:...
void SparsityPattern::compress()
{
  assert(!compressed);

  /// revise row
  //std::vector<u_int> new_rowIdx(row);
  std::vector<size_t> new_rowIdx(row); // updated version by Liu C.
  u_int idx = 0;
  for(u_int i = 0;i < m;++ i){
    for(u_int j = row[i];j < row[i + 1];++ j){
      if(col[j] == -1){
	new_rowIdx[i + 1] = idx;
	break;
      }
      else if (j == row[i + 1] - 1){
	idx ++;
	new_rowIdx[i + 1] = idx;
      }
      else {
	idx ++;	
      }
    }
  }
  row = new_rowIdx;
  
  /// revise col
  col.erase(std::remove(col.begin(), col.end(), -1), col.end());
  col.shrink_to_fit();

  /// finally, set compressed be true;
  compressed = true;
  n_nonzeros = col.size();

  /// if the matrix is a square matrix, we call
  /// put_diagEle_first_in_each_row() to defaulty set the diag element
  /// be the first one in each row
  if(m == n) put_diagEle_first_in_each_row();
}

/// When compressed is false, the vector col would have
/// max_entries_per_row elements in each row. add(x,x) would check
/// if the given position is available in the current status. If
/// yes, just return; otherwise, add the given position in the
/// pattern.
void SparsityPattern::add(const u_int& row_idx, const u_int& col_idx)
{
  assert(!compressed);
    
  for(u_int i = row[row_idx];i < row[row_idx + 1];++ i){
    if(col[i] == col_idx) return;
    if(col[i] == -1) {
      col[i] = col_idx;
      return;
    }
  }

  throw "Pretty much like that max_entries_per_row is too small...";
    
}

std::size_t * SparsityPattern::get_row()
{
  return row.data();
}
const std::size_t * SparsityPattern::get_row() const
{
  return row.data();
}

std::size_t * SparsityPattern::get_rowstart_indices()
{
  return row.data();
}
const std::size_t * SparsityPattern::get_rowstart_indices() const
{
  return row.data();
}

u_int * SparsityPattern::get_col()
{
  return col.data();
}
const u_int * SparsityPattern::get_col() const
{
  return col.data();
}


u_int * SparsityPattern::get_column_numbers()
{
  return col.data();
}
const u_int * SparsityPattern::get_column_numbers() const
{
  return col.data();
}

const void SparsityPattern::print() const
{
  std::cout << "The row vector is: " << std::endl;
  const u_int& row_size = row.size();
  for(u_int i = 0;i < row_size;++ i){
    std::cout << row[i] << " ";
  }
  std::cout << std::endl;

  const u_int& col_size = col.size();
  for(u_int i = 0;i < col_size;++ i){
    std::cout << col[i] << " ";
  }
  std::cout << std::endl;
}

u_int SparsityPattern::max_entries_per_row() const
{
  if (!compressed)
    return max_row_length;

  // if compress() was called, we
  // use a better algorithm which
  // gives us a sharp bound
  size_t m = 0;
  for (size_t i=1; i<row.size(); ++i)
    m = std::max (m, static_cast<size_t>(row[i]-row[i-1]));

  return m;
}

bool SparsityPattern::is_compressed()
{
  return compressed;
}

const bool SparsityPattern::is_compressed() const
{
  return compressed;
}

AFEPACK_CLOSE_NAMESPACE
