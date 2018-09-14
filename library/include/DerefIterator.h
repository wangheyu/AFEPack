/**
 * @file   DerefIterator.h
 * @author Ruo Li <rli@aztec>
 * @date   Mon Sep 28 12:46:18 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __DerefIterator_h__
#define __DerefIterator_h__

template <class IT, class VT>
  struct _Deref_iterator : public IT {
  typedef _Deref_iterator<IT,VT> _Self;
  typedef IT _Base;

  typedef VT value_type;
  typedef value_type* pointer;
  typedef value_type& reference;

  _Deref_iterator() : _Base() {}
  _Deref_iterator(const _Base& __x) : _Base(__x) {}
  _Deref_iterator(const _Self& __x) : _Base(__x) {}

  template <class OIT>
    _Self& operator=(const OIT& __x) {
    _Base::operator=(__x);
    return *this;
  }
  
  reference deref() const { return dynamic_cast<reference>(*(_Base::operator*())); }

  reference operator*() const { return deref(); }
  pointer operator->() const { return &deref(); }
};

#endif // __DerefIterator_h__

/**
 * end of file
 * 
 */
