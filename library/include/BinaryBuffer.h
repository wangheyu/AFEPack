/**
 * @file   BinaryBuffer.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Sep 24 09:46:48 2009
 * 
 * @brief  二进制缓冲区及其上构造的流。
 * 
 * 
 */

#ifndef __BinaryBuffer_h__
#define __BinaryBuffer_h__

#include <vector>
#include "Miscellaneous.h"

AFEPACK_OPEN_NAMESPACE

template <typename CHAR = char>
  class BinaryBuffer : public std::vector<CHAR> {
 public:
 typedef CHAR char_t;
 typedef std::vector<CHAR> _Base;
 typedef typename _Base::iterator iterator;
 typedef typename _Base::const_iterator const_iterator;

 BinaryBuffer() : _Base() {}
 ~BinaryBuffer() {}

 /** 
  * 返回当前缓存的起始地址，供MPI传递函数使用
  */
 void * start_address() const { return (void *)&(*this)[0]; }

};

AFEPACK_CLOSE_NAMESPACE

#endif // __BinaryBuffer_h__

/**
 * End of File
 * 
 */
