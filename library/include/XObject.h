/**
 * @file   XObject.h
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Fri Nov 15 15:54:50 2013
 * 
 * @brief  可以利用性质表进行随意构建的类的基类
 * 
 * 
 */

/*<

我们可以模仿下面的例子来构造具有灵活数据的类：

- 类的定义：
  struct A : XObject<A> {
    static const MemberProperty<A,int> member;
  };

  const MemberProperty<A,int> A::member;

- 类的使用：

  A a;
  a.alloc(a.member);
  int& member = a.get(a.member);
  member = 0;
  a.dealloc(a.member);

 */

#ifndef __XObject_h__
#define __XObject_h__

#include <string>
#include <map>

#include <AFEPack/PropertyTable.h>

template <class O, class T> class MemberProperty;

template <class DERIVED>
class XObject : public PropertyTable {
 public:
  template <class O, class T>
    T& alloc(const MemberProperty<O,T>& mp) const {
    return *new_property(mp.pid());
  }
  template <class O, class T>
    T& get(const MemberProperty<O,T>& mp) const {
    return *get_property(mp.pid());
  }
  template <class O, class T>
    void dealloc(const MemberProperty<O,T>& mp) const {
    free_property(mp.pid());
  }
};

template <class O, class T>
  class MemberProperty {
 public:
  typedef T value_t;
  typedef O object_t;

 public:
  MemberProperty<O,T>() {}

  const property_id_t<value_t>& pid() const {
    return _impl._pid;
  }

 private:
  struct _Impl {
    property_id_t<value_t> _pid;
    _Impl() {
      new_property_id<value_t>(_pid);
    }
  };
  const _Impl _impl;
};

#endif // __XObject_h__

/**
 * end of file
 * 
 */
