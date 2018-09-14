/**
 * @file   MPI.h
 * @author Ruo Li <rli@aztec>
 * @date   Fri Sep 25 14:58:29 2009
 * 
 * @brief  对MPI的接口
 * 
 * 我们对 MPI 提供一个简单而一致的接口，所有的数据传送均通过此方式进行，
 * 整个数据发送和接收的实现原理是：
 *
 * 1. 构造一个共享对象列表，表示一个实际的数据事实上在不同的进程均有其
 *    备份，如果两个进程上有共享对象，那么我们判定这两个进程之间需要交
 *    换数据；
 *
 * 2. 通过这个共享对象列表，可以构建所谓数据发送映射图表；
 *
 * 3. 在此数据映射发送图表上，我们完成需要交换的数据的打包、交换和解包
 *    操作。
 *
 * 在实际的实现和使用上：
 *
 * 1. 共享对象列表有两个实现，一个是 Shared_list，另一个是
 *    Shared_ptr_list。共享的对象对我们的应用来说，事实上就是几何遗传
 *    树中的几何体，这样的列表的构建一般都由几何遗传树自己完成。作为用
 *    户，可以在几何遗传树提供的列表上进行筛选，获得子列表以提升效率；
 *
 * 2. 数据发送映射图表的实现为类 Transmit_map，其构建方式提供多个接口
 *    函数，均命名为 build。其中可以支持对共享对象进行过滤，以获得效率
 *    上的提升。
 *
 * 3. 我们称一个类为数据打包器（DATA_PACKER），如果其提供了一个打包函数，
 *    一个解包函数，选择性地，其还可以提供一个共享对象过滤函数（参考
 *    sync_data 函数的说明）。通过一个 Transmit_map 和一个数据打包器就
 *    能够最终完成数据的通讯。
 * 
 * 作为一个典型的使用方式，可以参考 MPI::FaceData 名空间中的类的实现方
 * 式，其中将上述的实现进一步包装成为给最终用户使用的一个非常友好的接
 * 口，使得最终用户仅仅需要指定哪些数据需要传送就够了。当然，这也带来
 * 了一些限制。
 *
 * 其它的使用例子，可以参考 HGeometryMatcher 中在自适应的时候对几何体
 * 进行匹配的部分的实现，以及 HLoadBalance 中的两处应用。
 *
 */

#ifndef __MPI_h__
#define __MPI_h__

#include <string>
#include <vector>
#include <list>
#include <map>
#include <mpi.h>

#include "../../include/DerefIterator.h"
#include "../../include/BinaryBuffer.h"
#include "../../include/Migration.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  /**
   * 取一个大于 RAND_MAX/2 的随机数做标签，用来避免大部分的数据冲突。
   */
  int get_comm_tag(MPI_Comm comm);
  int get_comm_int(MPI_Comm comm);

  /**
   * 进程间发送和接收数据。此函数根据区域分解的基本方式来完成数据的传
   * 送，我们假设：对两个进程 i 和 j，要么互相发送数据，要么根本没有关
   * 系。发送和接收的数据都被打包成为二进制形式存储在一个使用
   * BinaryBuffer 实现的数组中，整个发送和接收的操作分为两步进行：第一
   * 步是告知对方自己要发送的数据的大小；第二步是发送要发送的数据。
   * 
   */
  template <class DataIterator, class TargetIterator>
    void sendrecv_data(MPI_Comm comm, 
                       int n,
                       DataIterator start_send_data,
                       DataIterator start_recv_data,
                       TargetIterator start_target) {
    int tag = get_comm_tag(comm);
    if (n == 0) return;

    MPI_Request request[2*n];
    MPI_Status status[2*n];

    /// 交换发送数据的大小
    int send_data_size[n], recv_data_size[n];
    DataIterator the_send_data = start_send_data;
    for (int i = 0;i < n;++ i, ++ the_send_data) {
      send_data_size[i] = the_send_data->size();
    }

    TargetIterator the_target = start_target;
    for (int i = 0;i < n;++ i, ++ the_target) {
      MPI_Isend(&send_data_size[i], 1, MPI_INT,
                *the_target, tag, comm, &request[i]);
      MPI_Irecv(&recv_data_size[i], 1, MPI_INT,
                *the_target, tag, comm, &request[i + n]);
    }
    MPI_Waitall(2*n, request, status);

    /// 接收数据缓冲区准备空间
    DataIterator the_recv_data = start_recv_data;
    for (int i = 0;i < n;++ i, ++ the_recv_data) {
      the_recv_data->resize(recv_data_size[i]);
    }

    /// 发送和接收数据
    int n_request = 0;
    the_target = start_target;
    the_send_data = start_send_data;
    the_recv_data = start_recv_data;
    for (int i = 0;i < n;++ i) {
      if (send_data_size[i] > 0) {
        MPI_Isend(the_send_data->start_address(), send_data_size[i], MPI_CHAR,
                  *the_target, tag, comm, &request[n_request ++]);
      }
      if (recv_data_size[i] > 0) {
        MPI_Irecv(the_recv_data->start_address(), recv_data_size[i], MPI_CHAR,
                  *the_target, tag, comm, &request[n_request ++]);
      }
      ++ the_target, ++ the_send_data, ++ the_recv_data;
    }
    if (n_request > 0) MPI_Waitall(n_request, request, status);
  }

  /**
   * 类型为 T 的远程指针。
   */
  template <class T>
    struct Remote_pointer {
      int type; /// 指针共享类型，约定 0 为缺省值，表示两个对象完全一致
      T * ptr;  /// 指针的远程内存地址
    Remote_pointer() : type(0), ptr(NULL) {}
    Remote_pointer(int _type, T * _ptr) :
      type(_type), ptr(_ptr) {}
    Remote_pointer(const Remote_pointer<T>& rp) :
      type(rp.type), ptr(rp.ptr) {}
      Remote_pointer<T>& operator=(const Remote_pointer<T>& rp) {
        type = rp.type;
        ptr = rp.ptr;
      }
    };

  namespace Shared_type_filter {
    struct all {
      bool operator()(int type) const {
        return true;
      }
    };
    template <int D0, int D1>
      struct between {
        bool operator()(int type) const {
          return (type >= D0)&&(type < D1);
        }
      };
    template <int D>
      struct only {
        bool operator()(int type) const {
          return (type == D);
        }
      };
    template <int D>
      struct except {
        bool operator()(int type) const {
          return (type != D);
        }
      };
    template <int D>
      struct greater_than {
        bool operator()(int type) const {
          return (type > D);
        }
      };
    template <int D>
      struct less_than {
        bool operator()(int type) const {
          return (type < D);
        }
      };
    template <class FILTER>
      struct negate {
        FILTER filter;

        negate() {}
        negate(const FILTER& _filter) : filter(_filter) {}

        bool operator()(int type) const {
          if (filter(type)) return false;
          else return true;
        }
      };
  }

  /**
   * 一个类型为 T 的对象的拷贝的远程指针数组。这个信息完整的描述一个对
   * 象的共享情况。本地对象的指针另外存储，不在此数组中。其元素是
   *
   *           <rank, remote_pointer>
   *
   * 的对。
   */
  template <class T>
    struct Shared_object : public std::multimap<int,Remote_pointer<T> > {
    typedef Remote_pointer<T> pointer_t;
    typedef std::pair<int,pointer_t> pair_t;
    typedef std::multimap<int,pointer_t> _Base;

    T * _ptr; /// 本地对象的指针

    Shared_object() {}
    Shared_object(T& t) : _ptr(&t) {}
    bool add_clone(int rank, T* ptr) {
      return this->add_clone(rank, pointer_t(0, ptr));
    }
    bool add_clone(int rank, int type, T* ptr) {
      return this->add_clone(rank, pointer_t(type, ptr));
    }
    bool add_clone(int rank, const pointer_t& ptr) {
      bool result = false;
      if (! is_duplicate_entry(rank, ptr)) {
        this->insert(pair_t(rank, ptr));
        result = true;
      } 
      return result;
    }

    T *& local_pointer() { return _ptr; }
    T * local_pointer() const { return _ptr; }
    T& local_object() const { return *_ptr; }

    /**
     * 检查是否已经有此条目。
     */
    bool is_duplicate_entry(int rank, 
                            const pointer_t& ptr) const {
      typedef typename _Base::const_iterator it_t;
      std::pair<it_t,it_t> range = _Base::equal_range(rank);
      it_t the_ptr = range.first, end_ptr = range.second;
      for (;the_ptr != end_ptr;++ the_ptr) {
        if (the_ptr->second.ptr == ptr.ptr) {
          return true;
        }
      }
      return false;
    }

    //@{
    /**
     * 对一个对象来说，其所有的拷贝所在的秩中最小的那个秩称为"首秩"。由
     * 于不是完全一致的共享对象(共享类型不为 0)在一个进程上可以有多个拷
     * 贝，一个对象在首秩上的拷贝可以有多个。
     */

    /**
     * 返回对象的首秩。其中rank必须是本进程的秩。
     */
    int primary_rank(int rank) const {
      return std::min(_Base::begin()->first, rank);
    }

    /**
     * 返回对象是否在首秩上。由于此处没有存储秩，因此需要传入秩作为参数。
     * rank 必须是本进程的秩，确保在 rank 上有对象的拷贝。由于multimap
     * 是有序的容器，实现时只需要比较第一个元素的秩即可。
     */
    bool is_on_primary_rank(int rank) const {
      return (_Base::begin()->first >= rank);
    }
    //@}

    /**
     * 对一个在首秩上的对象，其可能还有多个拷贝，其中那个内存地址最小
     * 的对象称为"主对象"。
     *
     * 判断此对象本身是否是主对象。基于上面同样的原因，也传入本进程的
     * 秩作为参数。
     * 
     */
    bool is_primary_object(int rank) const {
      int first_rank = _Base::begin()->first;
      bool result = true;
      if (first_rank < rank) {
        result = false; /// 不在首秩上
      } else if (first_rank == rank) { /**
                                        * 在此情况下，本地有多个拷贝，
                                        * 下面在几个本地拷贝中进行比较。
                                        */
        typedef typename _Base::const_iterator it_t;
        it_t the_ptr = _Base::begin();
        it_t end_ptr = _Base::upper_bound(rank);
        for (;the_ptr != end_ptr;++ the_ptr) {
          assert (the_ptr->second.ptr != _ptr);
          if (the_ptr->second.ptr < _ptr) {
            result = false;
            break;
          }
        }
      } else { /// (first_rank > rank)
        result = true; /// 是本进程上的唯一拷贝
      }
      return result;
    }
  };

  /**
   * 一系列类型为 T 的对象的远程指针数组的列表。我们将此信息收集起来好
   * 能够构造发送信息的数据结构。
   */
  template <class T, 
    template <class C, typename ALLOC = std::allocator<C> > class CNT = std::list>
    struct Shared_list : public CNT<Shared_object<T> > {};

  /**
   * 也可以使用指针作为列表的元素，通过 _Deref_iterator 使得遍历器的值
   * 仍然返回 Shared_object<T> 类型。
   */
  template <class T,
    template <class C, typename ALLOC = std::allocator<C> > class CNT = std::list>
    struct Shared_ptr_list : public CNT<Shared_object<T> *> {
    typedef CNT<Shared_object<T> *> base_t;
    typedef _Deref_iterator<typename base_t::iterator, Shared_object<T> > iterator;
    typedef _Deref_iterator<typename base_t::const_iterator, const Shared_object<T> > const_iterator;
    iterator begin() { return base_t::begin(); }
    iterator end() { return base_t::end(); }
    const_iterator begin() const { return base_t::begin(); }
    const_iterator end() const { return base_t::end(); }
    typename base_t::iterator begin_ptr() { return base_t::begin(); }
    typename base_t::iterator end_ptr() { return base_t::end(); }
    typename base_t::const_iterator begin_ptr() const { return base_t::begin(); }
    typename base_t::const_iterator end_ptr() const { return base_t::end(); }
  };

  /**
   * 数据发送映射图表。其中数据的意义为：
   *
   * std::map<秩, std::pair<信息条数, std::list<std::pair<本地指针，远程指针> > > >
   *
   */
  template <class T, class SHARED_TYPE_FILTER=Shared_type_filter::all>
    struct Transmit_map : 
    public std::map<int, std::pair<int, std::list<std::pair<T*, T*> > > > {
    typedef std::list<std::pair<T*, T*> > value_t;
    typedef std::pair<int, value_t> pair_t;
    typedef std::map<int, pair_t> _Base;
    typedef SHARED_TYPE_FILTER type_filter_t;

    type_filter_t type_filter;

    /**
     * 基于一个共享数据表创建数据发送接收映像图，要求 CONTAINER 的遍历
     * 器的值类型为 Shared_object<T>。
     */
    template <class CONTAINER>
      void build(const CONTAINER& shlist) {
      _Base::clear();

      typename CONTAINER::const_iterator 
        the_obj = shlist.begin(),
        end_obj = shlist.end();
      for (;the_obj != end_obj;++ the_obj) {
        this->add_object(*the_obj);
      }
    }

    template <class CONTAINER>
      void build(const CONTAINER& shlist,
                 bool (*filter)(T *)) {
      _Base::clear();

      typename CONTAINER::const_iterator 
        the_obj = shlist.begin(),
        end_obj = shlist.end();
      for (;the_obj != end_obj;++ the_obj) {
        this->add_object(*the_obj, 
                         (*filter)(the_obj->local_pointer()));
      }
    }

    template <class CONTAINER, class DATA_PACKER>
      void build(const CONTAINER& shlist,
                 DATA_PACKER& data_packer,
                 bool (DATA_PACKER::*filter)(T *)) {
      _Base::clear();

      typename CONTAINER::const_iterator 
        the_obj = shlist.begin(),
        end_obj = shlist.end();
      for (;the_obj != end_obj;++ the_obj) {
        this->add_object(*the_obj, 
                         (data_packer.*filter)(the_obj->local_pointer()));
      }
    }

    template <class CONTAINER, class DATA_PACKER>
      void build(const CONTAINER& shlist,
                 const DATA_PACKER& data_packer,
                 bool (DATA_PACKER::*filter)(T *) const) {
      _Base::clear();

      typename CONTAINER::const_iterator 
        the_obj = shlist.begin(),
        end_obj = shlist.end();
      for (;the_obj != end_obj;++ the_obj) {
        this->add_object(*the_obj, 
                         (data_packer.*filter)(the_obj->local_pointer()));
      }
    }

    /**
     * 基于一个遍历器创建数据发送接收映像图，要求 ITERATOR 的值类型为
     * Shared_object<T>。
     */
    template <class ITERATOR>
      void build(ITERATOR& begin, ITERATOR& end) {
      _Base::clear();

      ITERATOR the_obj(begin);
      for (;the_obj != end;++ the_obj) {
        this->add_object(*the_obj);
      }
    }

    /**
     * 加入一个共享对象。主要作私有用途。
     */
    void add_object(const Shared_object<T>& obj, 
                    bool is_add_entry = true) {
      typename Shared_object<T>::const_iterator
        the_ptr = obj.begin(),
        end_ptr = obj.end();
      if (! is_add_entry) {
        /**
         * 如果不加入条目，rank一定需要加上去，否则会导致通讯死锁。
         */
        for (;the_ptr != end_ptr;++ the_ptr) {
          int rank = the_ptr->first;
          if (this->find(rank) == this->end()) {
            (*this)[rank] = pair_t(0, value_t());
          }
        }
      } else {
        T * obj_ptr = obj.local_pointer();
        for (;the_ptr != end_ptr;++ the_ptr) {
          const int& rank = the_ptr->first;
          if (this->find(rank) == this->end()) {
            (*this)[rank] = pair_t(0, value_t());
          }

          if (type_filter(the_ptr->second.type)) {
            pair_t& pair = (*this)[rank];
            pair.first += 1;
            pair.second.push_back(std::pair<T*,T*>(obj_ptr, 
                                                   the_ptr->second.ptr));
          }
        }
      }
    }
  };

  /**
   * 在进程之间传递数据的函数的使用方式将会成对出现的，发送数据的时候，
   * 由远程指针中的 rank 指定发送目标，数据包由一系列的数据段组成。
   *
   *    typedef void (*data_pack_t)(T * local_obj,
   *                                int remote_rank,
   *                                Migration::ostream<>& os);
   *    typedef void (*data_unpack_t)(T * local_obj,
   *                                  int remote_rank,
   *                                  Migration::istream<>& is);
   *
   *
   * 此类有一个模板参数 DATA_PACKER，这个类需要提供两个成员函数类型如
   * 上面的 data_pack_t 和 data_unpack_t，来进行数据的编码和解码工作。
   * 
   * @param comm 通讯器
   * @param data_packer 可以完成数据编码和解码的类的对象
   * @param pack 进行编码的函数指针
   * @param unpack 进行解码的函数指针
   */
  template <class T, class DATA_PACKER, class SHARED_TYPE_FILTER>
    void sync_data(MPI_Comm comm,
                   Transmit_map<T,SHARED_TYPE_FILTER>& map,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&)) {
    typedef Transmit_map<T,SHARED_TYPE_FILTER> map_t;
    typedef typename map_t::value_t value_t;

    std::list<int> target_list;
    std::list<BinaryBuffer<> > data_buffer_in, data_buffer_out;

    int n = 0;
    typename map_t::iterator
      the_pair = map.begin(),
      end_pair = map.end();
    for (;the_pair != end_pair;++ the_pair, ++ n) {
      int rank = the_pair->first;
      target_list.push_back(rank);

      data_buffer_in.push_back(BinaryBuffer<>());
      data_buffer_out.push_back(BinaryBuffer<>());

      Migration::ostream<> os(data_buffer_out.back());
      int n_item = the_pair->second.first;
      if (n_item == 0) continue;

      value_t& lst = the_pair->second.second;
      os << n_item; /// 输出信息条目数
      typename value_t::iterator
        the_ptr = lst.begin(), end_ptr = lst.end();
      for (;the_ptr != end_ptr;++ the_ptr) {
        T *& local_obj = the_ptr->first;
        T *& remote_obj = the_ptr->second;
        os << remote_obj;
        (data_packer.*pack)(local_obj, rank, os);
      }
    }

    sendrecv_data(comm, n, data_buffer_out.begin(), data_buffer_in.begin(),
                  target_list.begin());

    typename std::list<BinaryBuffer<> >::iterator 
      the_buf = data_buffer_in.begin();
    the_pair = map.begin();
    for (;the_pair != end_pair;++ the_pair, ++ the_buf) {
      if (the_buf->size() == 0) continue;

      int rank = the_pair->first;
      Migration::istream<> is(*the_buf);
      int n_item;
      T * local_obj;
      is >> n_item; /// 信息条目数
      for (int i = 0;i < n_item;++ i) {
        is >> local_obj;
        (data_packer.*unpack)(local_obj, rank, is);
      }
    }
  }

  /** 
   * 使用 SHARED_LIST 中的信息进行数据同步化，同步化操作来自于
   * DATA_PACKER 及其提供的方法。
   * 
   * @param comm 通信器
   * @param shlist 共享对象列表
   * @param data_packer 可以完成数据编码和解码的类的对象
   * @param pack 进行编码的函数指针
   * @param unpack 进行解码的函数指针
   */
  template <class T, class SHARED_LIST, class DATA_PACKER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&)) {
    Transmit_map<T> map;
    map.build(shlist);
    sync_data(comm, map, data_packer, pack, unpack);
  }

  template <class T, class SHARED_LIST, class DATA_PACKER, class SHARED_TYPE_FILTER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&),
                   const SHARED_TYPE_FILTER& stf) {
    Transmit_map<T,SHARED_TYPE_FILTER> map;
    map.build(shlist);
    sync_data(comm, map, data_packer, pack, unpack);
  }

  /**
   * 使用成员函数做过滤器的版本。
   */
  template <class T, class SHARED_LIST, class DATA_PACKER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&),
                   bool (DATA_PACKER::*filter)(T *)) {
    Transmit_map<T> map;
    map.build(shlist, data_packer, filter);
    sync_data(comm, map, data_packer, pack, unpack);
  }

  template <class T, class SHARED_LIST, class DATA_PACKER, class SHARED_TYPE_FILTER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&),
                   bool (DATA_PACKER::*filter)(T *),
                   const SHARED_TYPE_FILTER& stf) {
    Transmit_map<T,SHARED_TYPE_FILTER> map;
    map.build(shlist, data_packer, filter);
    sync_data(comm, map, data_packer, pack, unpack);
  }

  /**
   * 使用常型成员函数做过滤器的版本。
   */
  template <class T, class SHARED_LIST, class DATA_PACKER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&),
                   bool (DATA_PACKER::*filter)(T *) const) {
    Transmit_map<T> map;
    map.build(shlist, data_packer, filter);
    sync_data(comm, map, data_packer, pack, unpack);
  }

  template <class T, class SHARED_LIST, class DATA_PACKER, class SHARED_TYPE_FILTER>
    void sync_data(MPI_Comm comm,
                   SHARED_LIST& shlist,
                   DATA_PACKER& data_packer,
                   void (DATA_PACKER::*pack)(T *,int,Migration::ostream<>&),
                   void (DATA_PACKER::*unpack)(T *,int,Migration::istream<>&),
                   bool (DATA_PACKER::*filter)(T *) const,
                   const SHARED_TYPE_FILTER& stf) {
    Transmit_map<T,SHARED_TYPE_FILTER> map;
    map.build(shlist, data_packer, filter);
    sync_data(comm, map, data_packer, pack, unpack);
  }

} // namespace MPI

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_h__

/**
 * end of file
 * 
 */
