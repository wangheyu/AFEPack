/**
 * @file   ThreadBase.h
 * @author Ruo Li
 * @date   Sun Feb 13 20:08:19 2005
 * 
 * @brief  for convenient using of POSIX thread in C++
 * 
 * 
 */

/*!
\class Thread?
\brief 这个类的目标是为了能够比较方便的在AFEPack中使用POSIX线程提供一个接口。

这个类的实现比较复杂，我主要是借鉴了deal.II中关于多线程支持的方法，如果你非常
想了解其中的技术，可以去参考那里的文档。学习一下可以提高C++方面的水平，建议
有空的时候去看一眼。我这个代码写得比较难读，建议不要读我的这个代码，浪费时间。
如果你读懂了，说明你的水平比我高，;-)。另外，我不提供任何关于线程安全性方面
的保证，全依赖你自己的代码去实现了。

deal.II中的多线程依赖于ACE，我在这里
直接自己写了一个比较简单的线程管理器。下面主要说明一下使用的方法，这对于大多
数用户来说是已经足够了。这个多线程的支持，可以帮助你方便的产生一组线程，协同
完成一个比较大的任务。线程的入口程序可以是类成员函数，也可以不是类成员函数。
我们先来介绍入口程序不是类成员函数的情况的用法：

\verbatim
／／ 这是线程的入口函数，计算第i个Chebyshev多项式在点0.64处的值
void fun(int i, const double& d)
{
  d = cos(i*acos(0.64));
};

int main(int argc, char * argv[])
{
  double val[5];
  ThreadManager th_man;        ／／ 创建线程管理器对象
  for (int i = 0;i < 5;i ++) { ／／ 我们创建5个线程来计算
    ／／ 创建线程，事实上所有的秘密都隐藏在这一行中
    th_man.start(encapsulate(&fun).collectArgs(i, val[i])); 
  }
  th_man.join(encapsulate(&fun)); ／／ 等待线程结束
  ／／ 上面这一行实现得非常难看，可是我实在是没有办法了
  ／／ 如果你有什么高招，请无论如何告诉我
  for (int i = 0;i < 5;i ++) { ／／ 输出计算结果
    std::cout << "Chebyshev(" << i 
	      << ", 0.64) = " << val[i]
	      << std::endl;
  }
  return 0;
};
\endverbatim

如果是类成员函数做入口函数，很显然的一件事情就是一定要把对象传
进去，请仿照下面的例子：

\verbatim
class Chebyshev {
 private:
  double x;
 public: ／／ 线程的入口函数必须是公有的
  Chebyshev(const double& _x) : x(_x) {};
  void fun(int i, const double& d) {
    d = cos(i*acos(x));
  };
};

int main(int argc, char * argv[])
{
  double val[5];
  Chebyshev ch(0.64);          ／／ 创建类的对象
  ThreadManager th_man;        ／／ 创建线程管理器对象
  for (int i = 0;i < 5;i ++) { ／／ 我们创建5个线程来计算
    ／／ 创建线程，我们现在要把对象传进去
    th_man.start(encapsulate(&Chebyshev::fun).collectArgs(&ch, i, val[i])); 
  }
  th_man.join(encapsulate(&Chebyshev::fun)); ／／ 等待线程结束
  for (int i = 0;i < 5;i ++) { ／／ 输出计算结果
    std::cout << "Chebyshev(" << i 
	      << ", 0.64) = " << val[i]
	      << std::endl;
  }
  return 0;
};
\endverbatim

这样的实现，唯一的目标就是为了使用起来方便。

如果在进行有限元计算的时候，你实现调用一行
\verbatim
setThread(n); ／／ n>1
\endverbatim
那么在进行有限元空间的构造的时候的有一部分代码会使用多线程进行，以后
我会将更多的部分进行多线程化，但是需要时间，请耐心等待哦！

我们使用AFEPack写的有限元的程序，常常有下面的函数
\verbatim
void fun()
{
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    ... ... ／／我们会在这个单元上进行一个比较复杂的数值积分的操作
  }
};

／／在另外的地方，我们调用上面这个函数
fun();
\endverbatim
对于这样的一段程序，我们就可以使用下面的方式来进行多线程化了：
\verbatim
void fun(int rank, int n_thread)
{
  int n_element = fem_space.n_elemen();
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
  int stride = n_element/n_thread;
  the_ele += stride*rank;
  if (rank + 1 < n_thread)
    end_ele = the_ele + stride;
  for (;the_ele != end_ele;++ the_ele) { ／／现在只是在一部分单元上了
    ... ... ／／还是在这个单元上进行一个比较复杂的数值积分的操作
  }
};

／／当然这时候在其他的地方调用这个程序就稍微复杂一点了
int n_thread = getThread();
ThreadManager th_man;
for (int rank = 1;rank < n_thread;rank ++) { ／／启动一堆线程
  th_man.start(encapsulate(&fun).collectArgs(rank, n_thread));
}
fun(0, n_thread); ／／自己也不要闲着
th_man.join(encapsulate(&fun)); ／／等着工作完成
\endverbatim
至少我自己觉得这样的使用方式还算是比较方便的，;-)

*/
template <typename Return
#if ARG>=1 
,typename Arg1 
#endif
#if ARG>=2 
,typename Arg2 
#endif
#if ARG>=3 
,typename Arg3 
#endif
#if ARG>=4 
,typename Arg4 
#endif
#if ARG>=5 
,typename Arg5 
#endif
#if ARG>=6 
,typename Arg6 
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
,typename Class=NullThread>
struct CLASSNAME 
{
  public:

typedef Return (*FunPtr)(
#if ARG>=1 
			 Arg1 
#endif
#if ARG>=2 
			 ,Arg2 
#endif
#if ARG>=3 
			 ,Arg3 
#endif
#if ARG>=4 
			 ,Arg4 
#endif
#if ARG>=5 
			 ,Arg5 
#endif
#if ARG>=6 
			 ,Arg6 
#endif
#if ARG>=7 
			 ,Arg7
#endif
#if ARG>=8 
			 ,Arg8 
#endif
#if ARG>=9 
			 ,Arg9 
#endif
#if ARG>=10 
			 ,Arg10 
#endif
			 );

  struct FunData 
  {
    public:
    FunPtr     fun_ptr;
#if ARG>=1 
    Arg1 arg1;
#endif
#if ARG>=2 
    Arg2 arg2;
#endif
#if ARG>=3 
    Arg3 arg3;
#endif
#if ARG>=4 
    Arg4 arg4;
#endif
#if ARG>=5 
    Arg5 arg5;
#endif
#if ARG>=6 
    Arg6 arg6;
#endif
#if ARG>=7 
    Arg7 arg7;
#endif
#if ARG>=8 
    Arg8 arg8;
#endif
#if ARG>=9 
    Arg9 arg9;
#endif
#if ARG>=10 
    Arg10 arg10;
#endif
    FunData(FunPtr _fun_ptr
#if ARG>=1 
	    ,Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) :
      fun_ptr(_fun_ptr)
#if ARG>=1 
	 ,arg1(_arg1)
#endif
#if ARG>=2 
	 ,arg2(_arg2)
#endif
#if ARG>=3 
	 ,arg3(_arg3)
#endif
#if ARG>=4 
	 ,arg4(_arg4)
#endif
#if ARG>=5 
	 ,arg5(_arg5)
#endif
#if ARG>=6 
	 ,arg6(_arg6)
#endif
#if ARG>=7 
	 ,arg7(_arg7)
#endif
#if ARG>=8 
	 ,arg8(_arg8)
#endif
#if ARG>=9 
	 ,arg9(_arg9)
#endif
#if ARG>=10 
	 ,arg10(_arg10)
#endif
    {};
    FunData(const FunData& _fd) :
      fun_ptr(_fd.fun_ptr)
#if ARG>=1 
	 ,arg1(_fd.arg1)
#endif
#if ARG>=2 
	 ,arg2(_fd.arg2)
#endif
#if ARG>=3 
	 ,arg3(_fd.arg3)
#endif
#if ARG>=4 
	 ,arg4(_fd.arg4)
#endif
#if ARG>=5 
	 ,arg5(_fd.arg5)
#endif
#if ARG>=6 
	 ,arg6(_fd.arg6)
#endif
#if ARG>=7 
	 ,arg7(_fd.arg7)
#endif
#if ARG>=8 
	 ,arg8(_fd.arg8)
#endif
#if ARG>=9 
	 ,arg9(_fd.arg9)
#endif
#if ARG>=10 
	 ,arg10(_fd.arg10)
#endif
    {};
    FunData& operator=(const FunData& _fd) 
{
fun_ptr = _fd.fun_ptr;
#if ARG>=1 
 arg1 = _fd.arg1;
#endif
#if ARG>=2 
 arg2 = _fd.arg2;
#endif
#if ARG>=3 
 arg3 = _fd.arg3;
#endif
#if ARG>=4 
 arg4 = _fd.arg4;
#endif
#if ARG>=5 
 arg5 = _fd.arg5;
#endif
#if ARG>=6 
 arg6 = _fd.arg6;
#endif
#if ARG>=7 
 arg7 = _fd.arg7;
#endif
#if ARG>=8 
 arg8 = _fd.arg8;
#endif
#if ARG>=9 
 arg9 = _fd.arg9;
#endif
#if ARG>=10 
 arg10 = _fd.arg10;
#endif
 return *this;
};
    static void * thread_entry(void *arg_ptr) {
      FunData * mfd = reinterpret_cast<FunData *>(arg_ptr);
      (*(mfd->fun_ptr))(
#if ARG>=1 
			    mfd->arg1
#endif
#if ARG>=2 
			    , mfd->arg2
#endif
#if ARG>=3 
			    , mfd->arg3
#endif
#if ARG>=4 
			    , mfd->arg4
#endif
#if ARG>=5 
			    , mfd->arg5
#endif
#if ARG>=6 
			    , mfd->arg6
#endif
#if ARG>=7 
			    , mfd->arg7
#endif
#if ARG>=8 
			    , mfd->arg8
#endif
#if ARG>=9 
			    , mfd->arg9
#endif
#if ARG>=10 
			    , mfd->arg10
#endif
			  );
    return NULL;
  };

  };

    class ArgCollector
  {
  public:
    typedef FunData FunDataType;
    FunPtr fun_ptr;

    ArgCollector (FunPtr _fun_ptr_) {
      fun_ptr = _fun_ptr_;
    };
        
    FunData collectArgs (
#if ARG>=1 
	    Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) {
      return FunData(fun_ptr
#if ARG>=1 
		     , _arg1
#endif
#if ARG>=2 
		     , _arg2
#endif
#if ARG>=3 
		     , _arg3
#endif
#if ARG>=4 
		     , _arg4
#endif
#if ARG>=5 
		     , _arg5
#endif
#if ARG>=6 
		     , _arg6
#endif
#if ARG>=7 
		     , _arg7
#endif
#if ARG>=8 
		     , _arg8
#endif
#if ARG>=9 
		     , _arg9
#endif
#if ARG>=10 
		     , _arg10
#endif
		      );
    };

  };

typedef Return (Class::*MemFunPtr)(
#if ARG>=1 
				   Arg1 
#endif
#if ARG>=2 
				   ,Arg2 
#endif
#if ARG>=3 
				   ,Arg3 
#endif
#if ARG>=4 
				   ,Arg4 
#endif
#if ARG>=5 
				   ,Arg5 
#endif
#if ARG>=6 
				   ,Arg6 
#endif
#if ARG>=7 
				   ,Arg7
#endif
#if ARG>=8 
				   ,Arg8 
#endif
#if ARG>=9 
				   ,Arg9 
#endif
#if ARG>=10 
				   ,Arg10 
#endif
				   );

  struct MemFunData 
  {
    public:
    MemFunPtr     mem_fun_ptr;
    Class *       object;
#if ARG>=1 
    Arg1 arg1;
#endif
#if ARG>=2 
    Arg2 arg2;
#endif
#if ARG>=3 
    Arg3 arg3;
#endif
#if ARG>=4 
    Arg4 arg4;
#endif
#if ARG>=5 
    Arg5 arg5;
#endif
#if ARG>=6 
    Arg6 arg6;
#endif
#if ARG>=7 
    Arg7 arg7;
#endif
#if ARG>=8 
    Arg8 arg8;
#endif
#if ARG>=9 
    Arg9 arg9;
#endif
#if ARG>=10 
    Arg10 arg10;
#endif
    MemFunData(MemFunPtr _mem_fun_ptr
	       ,Class * _object
#if ARG>=1 
	       ,Arg1 _arg1
#endif
#if ARG>=2 
	       ,Arg2 _arg2
#endif
#if ARG>=3 
	       ,Arg3 _arg3
#endif
#if ARG>=4 
	       ,Arg4 _arg4
#endif
#if ARG>=5 
	       ,Arg5 _arg5
#endif
#if ARG>=6 
	       ,Arg6 _arg6
#endif
#if ARG>=7 
	       ,Arg7 _arg7
#endif
#if ARG>=8 
	       ,Arg8 _arg8
#endif
#if ARG>=9 
	       ,Arg9 _arg9
#endif
#if ARG>=10 
	       ,Arg10 _arg10
#endif
	       ) :
      mem_fun_ptr(_mem_fun_ptr)
	 ,object(_object)
#if ARG>=1 
	 ,arg1(_arg1)
#endif
#if ARG>=2 
	 ,arg2(_arg2)
#endif
#if ARG>=3 
	 ,arg3(_arg3)
#endif
#if ARG>=4 
	 ,arg4(_arg4)
#endif
#if ARG>=5 
	 ,arg5(_arg5)
#endif
#if ARG>=6 
	 ,arg6(_arg6)
#endif
#if ARG>=7 
	 ,arg7(_arg7)
#endif
#if ARG>=8 
	 ,arg8(_arg8)
#endif
#if ARG>=9 
	 ,arg9(_arg9)
#endif
#if ARG>=10 
	 ,arg10(_arg10)
#endif
    {};
    MemFunData(const MemFunData& _mfd) :
      mem_fun_ptr(_mfd.mem_fun_ptr)
	 ,object(_mfd.object)
#if ARG>=1 
	 ,arg1(_mfd.arg1)
#endif
#if ARG>=2 
	 ,arg2(_mfd.arg2)
#endif
#if ARG>=3 
	 ,arg3(_mfd.arg3)
#endif
#if ARG>=4 
	 ,arg4(_mfd.arg4)
#endif
#if ARG>=5 
	 ,arg5(_mfd.arg5)
#endif
#if ARG>=6 
	 ,arg6(_mfd.arg6)
#endif
#if ARG>=7 
	 ,arg7(_mfd.arg7)
#endif
#if ARG>=8 
	 ,arg8(_mfd.arg8)
#endif
#if ARG>=9 
	 ,arg9(_mfd.arg9)
#endif
#if ARG>=10 
	 ,arg10(_mfd.arg10)
#endif
    {};
    MemFunData& operator=(const MemFunData& _mfd)
    {
      mem_fun_ptr = _mfd.mem_fun_ptr;
      object = _mfd.object;
#if ARG>=1 
      arg1 = _mfd.arg1;
#endif
#if ARG>=2 
      arg2 = _mfd.arg2;
#endif
#if ARG>=3 
      arg3 = _mfd.arg3;
#endif
#if ARG>=4 
      arg4 = _mfd.arg4;
#endif
#if ARG>=5 
      arg5 = _mfd.arg5;
#endif
#if ARG>=6 
      arg6 = _mfd.arg6;
#endif
#if ARG>=7 
      arg7 = _mfd.arg7;
#endif
#if ARG>=8 
      arg8 = _mfd.arg8;
#endif
#if ARG>=9 
      arg9 = _mfd.arg9;
#endif
#if ARG>=10 
      arg10 = _mfd.arg10;
#endif
      return *this;
    };
    static void * thread_entry(void *arg_ptr) {
      MemFunData * mfd = reinterpret_cast<MemFunData *>(arg_ptr);
      ((mfd->object)->*(mfd->mem_fun_ptr))(
#if ARG>=1 
					   mfd->arg1
#endif
#if ARG>=2 
					   , mfd->arg2
#endif
#if ARG>=3 
					   , mfd->arg3
#endif
#if ARG>=4 
					   , mfd->arg4
#endif
#if ARG>=5 
					   , mfd->arg5
#endif
#if ARG>=6 
					   , mfd->arg6
#endif
#if ARG>=7 
					   , mfd->arg7
#endif
#if ARG>=8 
					   , mfd->arg8
#endif
#if ARG>=9 
					   , mfd->arg9
#endif
#if ARG>=10 
					   , mfd->arg10
#endif
					   );
      return NULL;
    };

  };

    class MemArgCollector
  {
  public:
    typedef MemFunData FunDataType;
    MemFunPtr mem_fun_ptr;

    MemArgCollector (MemFunPtr _mem_fun_ptr_) {
      mem_fun_ptr = _mem_fun_ptr_;
    };
        
    MemFunData collectArgs (Class * object
#if ARG>=1 
	    ,Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) {
      return MemFunData(mem_fun_ptr
			,object
#if ARG>=1 
			, _arg1
#endif
#if ARG>=2 
			, _arg2
#endif
#if ARG>=3 
			, _arg3
#endif
#if ARG>=4 
			, _arg4
#endif
#if ARG>=5 
			, _arg5
#endif
#if ARG>=6 
			, _arg6
#endif
#if ARG>=7 
			, _arg7
#endif
#if ARG>=8 
			, _arg8
#endif
#if ARG>=9 
			, _arg9
#endif
#if ARG>=10 
			, _arg10
#endif
			);
    };
  };
};

template <typename Return
#if ARG>=1
,typename Arg1
#endif
#if ARG>=2
,typename Arg2
#endif
#if ARG>=3
,typename Arg3
#endif
#if ARG>=4
,typename Arg4
#endif
#if ARG>=5
,typename Arg5
#endif
#if ARG>=6
,typename Arg6
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
>
typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
>::ArgCollector encapsulate(Return (*_fun_ptr)(
#if ARG>=1
Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
)) {
    return typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
>::ArgCollector(_fun_ptr);
};

template <typename Return
#if ARG>=1
,typename Arg1
#endif
#if ARG>=2
,typename Arg2
#endif
#if ARG>=3
,typename Arg3
#endif
#if ARG>=4
,typename Arg4
#endif
#if ARG>=5
,typename Arg5
#endif
#if ARG>=6
,typename Arg6
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
,typename Class>
typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
,Class>::MemArgCollector encapsulate(Return (Class::*_mem_fun_ptr)(
#if ARG>=1
Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
)) {
  return typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
,Class>::MemArgCollector(_mem_fun_ptr);
};


/**
 * end of file
 * 
 */
