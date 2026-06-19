/**
 * @file   MPI_SyncProp.h
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Fri Dec 31 09:13:08 2010
 * 
 * @brief  在分区间对几何体上的一个性质进行同步。
 * 
 */

/*<

  对表达成为几何体上的性质的数据进行分区间的同步。此类接受三个模板参数：

    - FOREST：几何遗传树的类型；
    - PROPOUT：发送出去的性质的数据类型；
    - PROPIN：接收进来的性质的数据类型；

  此类使用的方式为以下步骤：

  1. 用户程序定义 property_id_t<PROPOUT> 和 property_id_t<PROPIN> 两个
     性质ID，并为部分几何体设置PROPOUT性质；

  2. 声明一个本类的对象，并使用几何遗传数和两个性质ID对其进行初始化；

  3. 调用本类对象的函数 sync<dim>() 对dim维的共享几何体上的数据进行交换；

  类型 PROPOUT 和 PROPIN 需要提供下面的流操作符：

     AFEPack::ostream<>& operator<<(AFEPack::ostream<>& os,
                                    const PROPOUT& prop_out)

     AFEPack::istream<>& operator>>(AFEPack::istream<>& is,
                                    PROPIN& prop_in)

 */

#ifndef __MPI_SyncProp_h__
#define __MPI_SyncProp_h__

#include "MPI_HGeometry.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <class FOREST, class PROPOUT, class PROPIN = PROPOUT>
    class PropSyncer {
  private:
  typedef PropSyncer<FOREST,PROPOUT,PROPIN> this_t;
  const FOREST * _p_forest;
  const property_id_t<PROPOUT> * _p_pid_out;
  const property_id_t<PROPIN> * _p_pid_in;

  public:
  PropSyncer() : _p_forest(NULL), _p_pid_out(NULL), _p_pid_in(NULL) {}
  PropSyncer(const FOREST& forest,
             const property_id_t<PROPOUT>& pid_out)
  : _p_forest(&forest), _p_pid_out(&pid_out), _p_pid_in(&pid_out) {}
  PropSyncer(const FOREST& forest,
             const property_id_t<PROPOUT>& pid_out,
             const property_id_t<PROPIN>& pid_in)
  : _p_forest(&forest), _p_pid_out(&pid_out), _p_pid_in(&pid_in) {}

  public:
  void set_forest(const FOREST& forest) {
    _p_forest = &forest;
  }

  void set_pid(const property_id_t<PROPOUT>& pid) {
    _p_pid_out = &pid;
    _p_pid_in = &pid;
  }
  void set_pid(const property_id_t<PROPOUT>& pid_out,
               const property_id_t<PROPIN>& pid_in) {
    _p_pid_out = &pid_out;
    _p_pid_in = &pid_in;
  }

  void set_pid_out(const property_id_t<PROPOUT>& pid) {
    _p_pid_out = &pid;
  }
  void set_pid_in(const property_id_t<PROPIN>& pid) {
    _p_pid_in = &pid;
  }

  template <int D>
  void sync() {
    sync_data(_p_forest->communicator(),
              _p_forest->template get_shared_list<D>(),
              *this,
              &this_t::template pack<D>,
              &this_t::template unpack<D>,
              &this_t::template is_pack_info<D>);
  }
  void sync(int D) {
    switch (D) {
    case 0: sync<0>(); break;
    case 1: sync<1>(); break;
    case 2: sync<2>(); break;
    case 3: sync<3>(); break;
    }
  }

  public:
  template <int D> bool
  is_pack_info(HGeometry<D,FOREST::dow> * p_geo) {
    return (p_geo->get_property(*_p_pid_out) != NULL);
  }
  template <int D>
  void pack(HGeometry<D,FOREST::dow> * p_geo,
            int remote_rank,
            Migration::ostream<>& os) {
    const PROPOUT * p_buf = p_geo->get_property(*_p_pid_out);
    assert (p_buf != NULL);
    os << *p_buf;
  }
  template <int D>
  void unpack(HGeometry<D,FOREST::dow> * p_geo,
              int remote_rank,
              Migration::istream<>& is) {
    PROPIN * p_buf = p_geo->get_property(*_p_pid_in);
    if (p_buf == NULL) {
      p_buf = p_geo->new_property(*_p_pid_in);
    }
    is >> *p_buf;
  }

  private:
  const void * p_op; /// 操作符的指针
  public:
  template <int D, class OP>
  void sync(const OP& op) {
    p_op = &op;
    sync_data(_p_forest->communicator(),
              _p_forest->template get_shared_list<D>(),
              *this,
              &this_t::template pack<D>,
              &this_t::template varunpack<D,OP>,
              &this_t::template is_pack_info<D>);
  }
  template <class OP>
  void sync(int D, const OP& op) {
    switch (D) {
    case 0: sync<0,OP>(op); break;
    case 1: sync<1,OP>(op); break;
    case 2: sync<2,OP>(op); break;
    case 3: sync<3,OP>(op); break;
    }
  }

  public:
  template <int D, class OP>
  void varunpack(HGeometry<D,FOREST::dow> * p_geo,
                 int remote_rank,
                 Migration::istream<>& is) {
    PROPIN * p_buf = p_geo->get_property(*_p_pid_in);
    if (p_buf == NULL) {
      p_buf = p_geo->new_property(*_p_pid_in);
      is >> *p_buf;
    } else {
      const OP& op = *((const OP *)p_op);
      PROPIN buf_in;
      is >> buf_in;
      *p_buf = op(*p_buf, buf_in);
    }
  }

  };
}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_SyncProp_h__

/**
 * end of file
 * 
 */
