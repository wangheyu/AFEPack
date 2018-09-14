/**
 * @file   PropertyTable.h
 * @author Robert Lie
 * @date   Mon Nov 13 20:52:26 2006
 * 
 * @brief  
 * 
 * 
 */

#ifndef __PropertyTable_h__
#define __PropertyTable_h__

#include <iostream>
#include <vector>
#include <list>
#include <map>

#include "Miscellaneous.h"

template <class T> class property_id_t; /// 前置声明
class PropertyTable;

namespace details {
  /**
   * 资源 ID 的基类，其虚函数用来使得派生类中的内存分配和释放函数能够
   * 被正确的调用。
   */
  class property_id_allocator_t
  {
  public:
    virtual ~property_id_allocator_t() {}

    virtual void * allocate() const {return NULL;}
    virtual void deallocate(void *) const {}
  };

  struct null_type {}; /// 用一个不占内存的类型做缺省类型

  class PropertyTableBase;

  /** 
   * 申请一个资源 ID
   * 
   * @param id 申请到的资源 ID 会写在其中
   * 
   * @return 是否申请到
   *
   */
  template <class T> 
    void _new_property_id(property_id_t<T>& id);

  /** 
   * 释放一个资源 ID
   * 
   * @param id 被释放的资源 ID
   *
   */
  template <class T> 
    void _free_property_id(property_id_t<T>& id);

};


/**
 * 资源 ID 类，其模板参数 T 是该资源存储的数据的数据类型。
 * 
 */
template <class T = details::null_type>
  class property_id_t : public details::property_id_allocator_t {
 public:
 typedef T value_type;
 private:
 typedef property_id_t<value_type> _Self;
 public:
 property_id_t() : _id(0xffff) {}
 property_id_t(const _Self& id) : _id(id.id()) {}
 virtual ~property_id_t() {
   if (_id != 0xffff) {
     details::_free_property_id<value_type>(*this);
   }
 }

 u_int id() const {return _id;}

 public:
 virtual void * allocate() const {
   return (void *)(new value_type());
 }
 virtual void deallocate(void * p_t) const {
   delete (value_type *)(p_t);
 }

 private:
 friend void details::_new_property_id<>(_Self&);
 friend void details::_free_property_id<>(_Self&);

 private:
 u_int& id() {return _id;}
 u_int _id;
};

namespace details {

  class PropertyTableBase {
  public:
    /// 全局性质表中存储的数据的类型
    class property_entry_t {
    public:
    property_entry_t() : p_obj(NULL), p_data(NULL) {}
    property_entry_t(PropertyTableBase * po, void * pd) :
      p_obj(po), p_data(pd) {}
    property_entry_t(const property_entry_t& pe) :
      p_obj(pe.p_obj), p_data(pe.p_data) {}

    public:
      property_entry_t& operator=(const property_entry_t& pe) {
        p_obj = pe.p_obj;
        p_data = pe.p_data;

        return *this;
      }

      PropertyTableBase * object_ptr() const {return p_obj;}
      PropertyTableBase *& object_ptr() {return p_obj;}

      void * data_ptr() const {return p_data;}
      void *& data_ptr() {return p_data;}

    private:
      PropertyTableBase * p_obj;
      void * p_data;
    };

    /// 每一个性质使用一个全局的 std::list 类型来存储及其遍历器类型
    typedef std::list<property_entry_t> prop_tbl_list_t;
    typedef prop_tbl_list_t::iterator prop_tbl_list_iterator_t;

    /// 全局的所有性质表列表都放在一个 map 中，map 的值类型中包含内存分配器的
    /// 指针和性质数据的列表。
    struct prop_tbl_entry_t {
      property_id_allocator_t * _alloc;
      prop_tbl_list_t _list;
    };
    typedef std::map<u_int, prop_tbl_entry_t> prop_tbl_t;

    /// 每个对象的性质表记录其在全局中的应用
    struct prop_iter_entry_t {
      prop_tbl_entry_t * _global_entry;
      prop_tbl_list_iterator_t _iter;
    };
    typedef std::map<u_int, prop_iter_entry_t> prop_iter_t;
    typedef prop_iter_t::iterator prop_iter_iter_t;
    typedef prop_iter_t::const_iterator prop_iter_const_iter_t;

  private:
    //@{
    /// 构造函数和析构函数使用私有保护，避免被直接使用
    PropertyTableBase() {}
    PropertyTableBase(const PropertyTableBase& pt) {}
    virtual ~PropertyTableBase();
    //@}

    /// 拷贝算子
    PropertyTableBase& operator=(const PropertyTableBase& pt) { return *this; }

    template <class T>
      bool has_property(const property_id_t<T>& i, 
                        prop_iter_const_iter_t& iter) const
      {
        iter = _prop_iter.find(i.id());
        return (iter != _prop_iter.end());
      }
    template <class T>
      bool has_property(const property_id_t<T>& i, 
                        prop_iter_iter_t& iter)
      {
        iter = _prop_iter.find(i.id());
        return (iter != _prop_iter.end());
      }

  public:
    /** 
     * 返回此对象是否具有该性质资源
     * 
     * @param i 性质资源 ID
     * 
     * @return 逻辑值
     */  
    template <class T>
      bool has_property(const property_id_t<T>& i) const
      {
        return (_prop_iter.find(i.id()) != _prop_iter.end());
      }

    /** 
     * 分配对象上的资源空间
     * 
     * @param id 资源 ID
     * 
     * @return 资源的指针
     */
    template <class T> 
      T * new_property(const property_id_t<T>& i)
      {
        assert(! has_property(i)); /// 确认没有此性质
        return (T *)new_property(i.id());
      }

    /** 
     * 获取对象上的资源指针
     * 
     * @param id 资源 ID
     * 
     * @return 资源的指针
     */
    template <class T> 
      T * get_property(const property_id_t<T>& i) const
      {
        prop_iter_const_iter_t iter;
        if (! has_property(i, iter)) {
          return (T *)NULL;
        }
        else {
          return (T *)(iter->second._iter->data_ptr());
        }
      }

    /**
     * 释放对象上的资源指针
     *
     * @param id 资源 ID
     * 
     */
    template <class T> 
      void free_property(const property_id_t<T>& i)
      {
        prop_iter_iter_t iter;
        if (has_property(i, iter)) { /// 确认有此性质
          free_property(iter);
        }
      }

    /**
     * 释放对象上的所有资源指针。
     * 
     */
    void clear_property();

  public:
    template <class T> friend void _new_property_id(property_id_t<T>& i);
    template <class T> friend void _free_property_id(property_id_t<T>& i);
    friend class ::PropertyTable;

  private:
    //@{
    /// 根据 id 号来分配和释放资源
    void * new_property(u_int i); /// 根据 id 号分配资源
    void free_property(prop_iter_iter_t& it); /// 根据 id 号释放资源
    //@}

    prop_iter_t _prop_iter; /**< 本对象的性质在全局表中的位置 */

    static prop_tbl_t _prop_tbl; /**< 全局的性质表 */
  };
};

/*!

  本类用来将性质表始终做成可写的。对于几何体这样的对象来说，如果直接从
  性质表派生，则当我们得到一个有 const 修饰符这样的对象的时候，就不能改
  变其性质表，我们通过本类来做一个包装，使得被 const 修饰符修饰的对象也
  能获得可写的性质表。

*/
class PropertyTable {
 private:
  typedef details::PropertyTableBase Base;
  void * p_tbl;
 public:
  //@{
  /// 构造函数和析构函数
  PropertyTable() {
    p_tbl = new Base();
  }
  /// 注意：拷贝构造函数实际上不做任何拷贝！
  PropertyTable(const PropertyTable& ptp) {
    p_tbl = new Base();
  }
  ~PropertyTable() {
    Base * p = (Base *)p_tbl;
    delete p;
  }
  //@}

  /// 注意：拷贝算子啥都不做的！
  PropertyTable& operator=(const PropertyTable& ptp) {
    return *this;
  }

 public:
  /** 
   * 检查对象上的资源空间
   * 
   * @param id 资源 ID
   * 
   * @return 是否有相应的资源
   */
  template <class T> 
    bool has_property(const property_id_t<T>& id) const {
    Base * p = (Base *)p_tbl;
    return p->has_property(id);
  }

  /** 
   * 分配对象上的资源空间
   * 
   * @param id 资源 ID
   * 
   * @return 资源的指针
   */
  template <class T> 
    T * new_property(const property_id_t<T>& id) const {
    Base * p = (Base *)p_tbl;
    return p->new_property(id);
  }

  /** 
   * 获取对象上的资源指针
   * 
   * @param id 资源 ID
   * 
   * @return 资源的指针
   */
  template <class T> 
    T * get_property(const property_id_t<T>& id) const {
    Base * p = (Base *)p_tbl;
    return p->get_property(id);
  }

  /**
   * 释放对象上的资源指针
   *
   * @param id 资源 ID
   * 
   */
  template <class T> 
    void free_property(const property_id_t<T>& id) const {
    Base * p = (Base *)p_tbl;
    p->free_property(id);
  }

  void clear_property() {
    Base * p = (Base *)p_tbl;
    p->clear_property();
  }
};

namespace details {
  template <class T> 
    void _new_property_id(property_id_t<T>& i) {
    typedef details::PropertyTableBase PTB;
    /// 搜索空闲 ID
    i.id() = 0;
    while (PTB::_prop_tbl.find(i.id()) != PTB::_prop_tbl.end()) {
      ++ i.id();
    }

    std::pair<u_int, typename PTB::prop_tbl_entry_t> entry;
    entry.first = i.id();
    entry.second._alloc = &i;
    PTB::_prop_tbl.insert(entry);
  }

  template <class T> 
    void _free_property_id(property_id_t<T>& i) 
    {
      typedef details::PropertyTableBase PTB;
      typename PTB::prop_tbl_t::iterator entry = PTB::_prop_tbl.find(i.id());
      if (entry == PTB::_prop_tbl.end()) return; /// 此资源已经释放过了

      typename PTB::prop_tbl_list_t& lst = entry->second._list;
      typename PTB::prop_tbl_list_t::iterator the_ptr = lst.begin();
      typename PTB::prop_tbl_list_t::iterator end_ptr = lst.end();
      for (;the_ptr != end_ptr;) {
        typename PTB::prop_tbl_list_iterator_t this_ptr = the_ptr;
        ++ the_ptr;
        this_ptr->object_ptr()->free_property(i);
      }
      PTB::_prop_tbl.erase(entry);
      i.id() = 0xffff;
    }
}

template <class T> void 
new_property_id(property_id_t<T>& pid) {
  details::_new_property_id(pid);
}

template <class T> void 
free_property_id(property_id_t<T>& pid) {
  details::_free_property_id(pid);
}


#endif // __PropertyTable_h__

/**
 * end of file
 * 
 */

