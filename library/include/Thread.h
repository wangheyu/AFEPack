/**
 * @file   Thread.h
 * @author Ruo Li
 * @date   Sun Feb 13 21:57:00 2005
 * 
 * @brief  for convenient using of POSIX thread in C++
 * 
 * 
 */

#ifdef MULTITHREAD

#ifndef __Thread_h__
#define __Thread_h__

#include <pthread.h>
#include <iostream>
#include <list>

struct NullThread {};

#define ARG 0
#define CLASSNAME Thread0
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 1
#define CLASSNAME Thread1
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 2
#define CLASSNAME Thread2
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 3
#define CLASSNAME Thread3
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 4
#define CLASSNAME Thread4
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 5
#define CLASSNAME Thread5
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 6
#define CLASSNAME Thread6
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 7
#define CLASSNAME Thread7
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 8
#define CLASSNAME Thread8
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 9
#define CLASSNAME Thread9
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 10
#define CLASSNAME Thread10
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

/** 
 * 设置内在的线程个数，当大于1时，库的有些部分会自动进行多
 * 线程化。
 * 
 * @param int 线程个数
 */
void setThread(unsigned int);
/** 
 * 取得内在的线程个数
 * 
 * 
 * @return 线程个数
 */
unsigned int getThread();

/*!
\class ThreadManager
\brief 进行线程管理

进行线程的创建以及管理的工作，事实上我只是写了一个join函数。

*/
class ThreadManager
{
 private:
  std::list<pthread_t> threads;	/**< 线程的句柄 */
  std::list<void *>    args;	/**< argments of the thread */
  bool                 is_joined; /**< 是否已经join过了 */
 public:
  ThreadManager() : 
    is_joined(false) {};
  ~ThreadManager() {
    if (!is_joined && !threads.empty()) {
      std::cerr << "Thread manager is not joined before destory."
		<< std::endl;
      abort();
    }
  };
  /** 
   * 启动一个线程
   * 
   * @param fun_data 传递给线程入口程序的参数表
   * @param attr 线程属性
   * 
   * @return 成功返回0，否则-1
   */
  template <class T>
  int start(T fun_data,
	    pthread_attr_t * attr = NULL) {
    pthread_t th;
    T * fun_data_copy = new T(fun_data); /** 
					  * it looks that it is more safe 
					  * if we make a copy here.
					  */
    int error_number = pthread_create(&th, 
				      attr,
				      fun_data.thread_entry, 
				      (void *)fun_data_copy);
    if (error_number == 0) {
      threads.push_back(th);
      args.push_back(fun_data_copy);
#ifdef DEBUG
      //std::cerr << "new thread created successfully. " 
      //	<< std::endl;
#endif
      return 0;
    }
    else {
      std::cout << "thread creating failure with error_number "
		<< error_number << std::endl;
      exit(-1);
    }
  };
  /** 
   * 等待线程执行结束
   * 
   * 
   * @return 成功返回0，否则-1
   */
  template <class T>
  int join(T arg_collector) {
    typedef typename T::FunDataType FunData;
    std::list<pthread_t>::iterator the_th = threads.begin();
    std::list<pthread_t>::iterator end_th = threads.end();
    std::list<void *>::iterator the_arg = args.begin();
    for (;the_th != end_th;the_th ++, the_arg ++) {
      int error_number = pthread_join(*the_th, NULL);
      if (error_number) {
	std::cout << "thread join error with error_number "
		  << error_number << std::endl;
	exit(-1);
      }
      delete (FunData *)(*the_arg);
    }
    clean();
    return 0;
  };
 private:
  void clean() {
    threads.clear();
    args.clear();
    is_joined = false;
  };
};

#endif // __Thread_h__

#endif 

/**
 * end of file
 * 
 */
