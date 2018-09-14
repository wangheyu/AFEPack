/**
 * @file   Migration.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Sep 24 11:51:29 2009
 * 
 * @brief  数据迁移
 * 
 * 
 */

#ifndef __Migration_h__
#define __Migration_h__

#include <map>
#include "BinaryBuffer.h"

AFEPACK_OPEN_NAMESPACE

/**
 * 在进行整个环境的序列化以前，我们先构造序列化的环境。我们设定基本前提
 * 为：假定进程上仅仅有一棵几何遗传树，在此树上多个非正则网格，在每个非
 * 正则网格上建立了一个正则网格，在正则网格上则建立了多个有限元空间，每
 * 个有限元空间上有多个有限元函数。
 *
 * 由于在进行网格的拆分和合并的时候，正则网格需要重建，从而有限元空间也
 * 需要进行重建，为了能够对有限元函数的信息进行迁移，我们在猜分以前，将
 * 有限元函数的信息存储到相应的几何遗传树的几何体上去。几何遗传树中的几
 * 何体继承了 HGeomeryBuffer 这个类，从而我们可以将信息存储在这个类的缓
 * 冲区中。比如从有限元函数的出发，每个自由度都能够找到自身所依附的正则
 * 网格中的几何体，而每个正则网格中的几何体则需要能够找寻到在正则化的过
 * 程中相应的几何遗传树中的几何体，我们对正则网格进行了改进，使得这个检
 * 索成为可能。这个是通过在正则网格中加入了h_geometry_ptr数组来完成的。
 *
 * 在此基础上，我们假设在每一个几何体上附着了一组数据，这些数据按照数
 * 据 ID 被分成组，比如一个有限元函数的数据我们可以分为一组。一个几何
 * 体上，属于同一个数据 ID 的数据则必须由用户自己负责分析，也就是说，
 * 这组数据本身必须在上下文中是可以自解释的。具体参考 export_fe_func
 * 和 import_fe_func 这对函数的实现。
 *
 */
namespace Migration {

  template <typename BUFFER>
    class stream_base {
  protected:
    typedef typename BUFFER::char_t char_t;
    typedef typename BUFFER::const_iterator const_iterator;
  public:
    /// 空构造函数
  stream_base() : _buf(NULL) {}
    /// 提供一个数据缓存的构造函数 
  stream_base(BUFFER& buf) : _buf(&buf) {}
    /// 析构函数
    ~stream_base() {}

  public:
    /** 
     * 提供数据缓存对象
     * 
     * @param buf 
     */
    void set_buffer(BUFFER& buf) { _buf = &buf; }

    /** 
     * 查询得到数据缓存对象
     * @return  数据缓存
     */
    const BUFFER& get_buffer() const { return *_buf; }

  protected:
    /// 被操作的缓存对象
    BUFFER * _buf;
  };

  /// 支持数据写入缓存。
  template <typename BUFFER = BinaryBuffer<> >
    class ostream : public stream_base<BUFFER> {
  private:
  typedef stream_base<BUFFER> _Base;
  typedef typename _Base::char_t char_t;
  using _Base::_buf;
  public:
  ostream() : _Base() {}
  ostream(BUFFER& buf) : _Base(buf) {}
  ~ostream() {}

  /** 
   * 将类型为T数据t 放到当前缓存向量的末尾注意这里的T只支持基本数据类型的写
   * 入，其他任何非c++基本类型的写入需要用户自己给定写入格式。
   * 
   * @param t 待写入缓存的数据
   */
  template <class T>
  void encode(const T& t) {
    int n = sizeof(t);
    const char_t * ptr = (const char_t *) (&t);
    _buf->insert(_buf->end(), ptr, ptr + n);
  }

  void encode_binary(void * data, int n) {
    const char_t * ptr = (const char_t *)data;
    _buf->insert(_buf->end(), ptr, ptr + n); 
  }

  };

  /// 支持数据从缓存中读取。此对象构造完成后，目前只能从头到尾遍历读取一次。若要多次读取，需要重新构造流对象。
  template <typename BUFFER = BinaryBuffer<> >
    class istream : public stream_base<const BUFFER>{
  private:
  typedef stream_base<const BUFFER> _Base;
  typedef typename _Base::char_t char_t;
  typedef typename _Base::const_iterator iterator;
  using _Base::_buf;
  public:
  istream() : _Base() {}
  istream(const BUFFER& buf) : _Base(buf) { reset(); }
  ~istream() {}

  /// 设定作为数据源的缓存对象
  void set_buffer(const BUFFER& buf) {
    _Base::set_buffer(buf);
    reset();
  }

  /** 
   *  从当前向量读取数据到变量t中。这里同样只支持c++基本数据类型，若有需要处
   *  理其他类型的数据，用户需要提供带有特定数据类型参数的同名函数。
   * 
   * @param t 被读取的数据，同时被转换成类型T
   */
  template <class T>
  void decode(T& t) {
    std::size_t n = sizeof(T);
    std::copy(_pos, _pos + n, (char_t *)(&t));
    _pos += n;
  }
  void decode_binary(void * data, int n) {
    std::copy(_pos, _pos + n, (char_t *)data);
    _pos += n;
  }

  private:
  /** 
   * 将指示子放到数据缓存的开始位置，准备读取缓存数据
   * 
   */
  void reset() {
    if (_buf != NULL) {
      _pos = _buf->begin();
    }
  }

  private:
  /// 读取缓存是逐步进行的，因此需要一个迭代子指示当前的位置
  iterator _pos;
  }; 

  template <class OSTREAM, class T>
    void encode(OSTREAM& os, const T& t) {
    os.encode(t);
  }
  template <class ISTREAM, class T>
    void decode(ISTREAM& is, T& t) {
    is.decode(t);
  }

  template <class OSTREAM, class T>
    void encode_binary(OSTREAM& os, const T& t) {
    os.encode_binary(t);
  }
  template <class ISTREAM, class T>
    void decode_binary(ISTREAM& is, T& t) {
    is.decode_binary(t);
  }

#define BASIC_BINARY_IO(TYPE)                                   \
  template <class BUF>                                          \
    ostream<BUF>& operator<<(ostream<BUF>& os, const TYPE& t) { \
    os.encode(t);                                               \
    return os;                                                  \
  }                                                             \
  template <class BUF>                                          \
  istream<BUF>& operator>>(istream<BUF>& is, TYPE& t) {         \
    is.decode(t);                                               \
    return is;                                                  \
  }

  BASIC_BINARY_IO(bool)
  BASIC_BINARY_IO(char)
  BASIC_BINARY_IO(int)
  BASIC_BINARY_IO(unsigned int)
  BASIC_BINARY_IO(long)
  BASIC_BINARY_IO(unsigned long)
  BASIC_BINARY_IO(double)
  BASIC_BINARY_IO(long long)
  BASIC_BINARY_IO(unsigned long long)
  BASIC_BINARY_IO(float)
  BASIC_BINARY_IO(void *)

#undef BASIC_BINARY_IO

  template <class BUF, class T>
    ostream<BUF>& operator<<(ostream<BUF>& os, const T*& t) {
    os.encode(t);
    return os;
  }
  template <class BUF, class T>
    ostream<BUF>& operator<<(ostream<BUF>& os, T*& t) {
    os.encode(t);
    return os;
  }
  template <class BUF, class T>
    istream<BUF>& operator>>(istream<BUF>& is, const T*& t) {
    is.decode(t);
    return is;
  }
  template <class BUF, class T>
    istream<BUF>& operator>>(istream<BUF>& is, T*& t) {
    is.decode(t);
    return is;
  }

  template <class BUF, class T>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<T>& t) {
    u_int n = t.size();
    os << n;
    for (u_int i = 0;i < n;++ i) {
      os << t[i];
    }
    return os;
  }
  template <class BUF, class T>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<T>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    for (u_int i = 0;i < n;++ i) {
      is >> t[i];
    }
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<double>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(double)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<double>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(double)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<int>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(int)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<int>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(int)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const std::vector<u_int>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&t[0]), sizeof(u_int)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, std::vector<u_int>& t) {
    u_int n;
    is >> n;
    t.resize(n);
    is.decode_binary(&t[0], sizeof(u_int)*n);
    return is;
  }

  template <class BUF>
    ostream<BUF>& operator<<(ostream<BUF>& os, const Vector<double>& t) {
    u_int n = t.size();
    os << n;
    os.encode_binary((void *)(&(*t.begin())), sizeof(double)*n);
    return os;
  }
  template <class BUF>
    istream<BUF>& operator>>(istream<BUF>& is, Vector<double>& t) {
    u_int n;
    is >> n;
    t.reinit(n);
    is.decode_binary(&t(0), sizeof(double)*n);
    return is;
  }

  template <class BUF, typename CHAR>
    ostream<BUF>& operator<<(ostream<BUF>& os, const BinaryBuffer<CHAR>& buf) {
    const std::vector<CHAR>& v(buf);
    os << v;
    return os;
  }
  template <class BUF, typename CHAR>
    istream<BUF>& operator>>(istream<BUF>& is, BinaryBuffer<CHAR>& buf) {
    std::vector<CHAR>& v(buf);
    is >> v;
    return is;
  }

  typedef int data_id_t; /// 数据 ID 的类型：为有符号整型
  typedef std::string data_name_t; /// 数据名称类型：字符串
  typedef BinaryBuffer<> data_buffer_t; /// 数据缓冲区类型
  /**
   * 数据 ID 到缓冲区映射表，这个是用在最终存储数据位置的类型。
   */
  typedef std::map<data_id_t, data_buffer_t> buffer_t; 

  data_id_t name_to_id(const data_name_t& dn); /// 从数据名称到ID的转换
  data_id_t register_data_name(const data_name_t& dn, bool); /// 登记数据名称
  void initialize(); /// 初始化数据迁移环境
  bool is_valid(const data_id_t&); /// 判断一个数据ID是否合法

  /**
   * 用来进行数据迁移的基本缓冲区对象。
   */
  struct HBuffer {
    buffer_t buffer;
    virtual ~HBuffer() {}
  };

  namespace details {
    /**
     * 获取树结构中几何体上的数据缓冲区。对于存储的情形，如果相应的缓冲
     * 区不存在，则新建立改缓冲区。
     * 
     * @param geo 树结构中的几何体 
     * @param is_save 存储还是读入
     * 
     * @return 相应的缓冲区
     */
    template <class HGEO> BinaryBuffer<>& 
      get_buffer(const HGEO& geo,
                 const data_id_t& data_id,
                 bool is_save) {
      /**
       * 强行将const指针转换成了非const的指针，似乎没法子呢...
       */
      buffer_t& buffer = *(buffer_t *)(&geo.buffer);
      if (is_save) {
        return buffer[data_id];
      } else {
        buffer_t::iterator it = buffer.find(data_id);
        if (it == buffer.end()) {
          return *((BinaryBuffer<> *)NULL);
        }
        return it->second;
      }
    }

    /**
     * 获取几何体上的数据缓冲区。对于存储的情形，如果相应的缓冲区不存在，
     * 则新建立改缓冲区。
     * 
     * @param mesh 正则网格
     * @param data_id 数据 ID
     * @param dimension 几何体维数
     * @param geo_idx 几何体序号
     * @param is_save 存储还是读入
     * 
     * @return 相应的缓冲区
     */
    template <class MESH> BinaryBuffer<>& 
      get_buffer(MESH& mesh,
                 const data_id_t& data_id,
                 u_int dimension,
                 u_int geo_idx,
                 bool is_save) {
      buffer_t& buffer = mesh.h_geometry(dimension, geo_idx)->buffer;
      if (is_save) {
        return buffer[data_id];
      } else {
        buffer_t::iterator it = buffer.find(data_id);
        if (it == buffer.end()) {
          return *((BinaryBuffer<> *)NULL);
        }
        return it->second;
      }
    }

  }

  /** 
   * 获取树结构中的几何体上的一个输出流，用来输出数据。
   * 
   * @param geo 树结构中的几何体
   * @param data_id 数据 ID
   * @param os 获得的输出流
   */
  template <class HGEO, class STREAM> void
    get_export_stream(HGEO& geo,
                      const data_id_t& data_id,
                      STREAM& os) {
    os.set_buffer(details::get_buffer(geo, data_id, true));
  }
  
} // end of namespace Migration

AFEPACK_CLOSE_NAMESPACE

#endif // __Migration_h__

/**
 * end of file
 * 
 */
 
