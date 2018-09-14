/**
 * @file   PropertyTable.cpp
 * @author Robert Lie
 * @date   Tue Nov 14 13:25:13 2006
 * 
 * @brief  
 * 
 * 
 */

#include "PropertyTable.h"

#define TEMPLATE 

#ifndef INLINE
#define INLINE
#endif

#define THIS details::PropertyTableBase

/// 全局的性质表数据表
THIS::prop_tbl_t THIS::_prop_tbl;

TEMPLATE
THIS::~PropertyTableBase()
{
  clear_property();
}

TEMPLATE
INLINE
void THIS::clear_property() {
  prop_iter_iter_t
    the_iter = _prop_iter.begin(),
    end_iter = _prop_iter.end();
  for (;the_iter != end_iter;++ the_iter) {
    free_property(the_iter);
  }
}

TEMPLATE
INLINE
void * THIS::new_property(u_int i)
{
  THIS::prop_tbl_entry_t& entry = _prop_tbl.find(i)->second;
  THIS::prop_tbl_list_t& lst = entry._list;
  void * p_data = entry._alloc->allocate(); /// 分配内存
  std::pair<u_int, THIS::prop_iter_entry_t> pit;
  pit.first = i;
  pit.second._global_entry = &entry;
  pit.second._iter = lst.insert(lst.end(), THIS::property_entry_t(this, p_data));
  _prop_iter.insert(pit);
  return p_data;
}

TEMPLATE
INLINE
void THIS::free_property(prop_iter_iter_t& iter)
{
  THIS::prop_tbl_entry_t& entry = *(iter->second._global_entry);
  THIS::prop_tbl_list_t& lst = entry._list;
  const THIS::prop_tbl_list_iterator_t& the_prop = iter->second._iter;
  entry._alloc->deallocate(the_prop->data_ptr());
  lst.erase(the_prop);
  _prop_iter.erase(iter);
}

#undef THIS
#undef INLINE
#undef TEMPLATE

/**
 * end of file
 * 
 */
