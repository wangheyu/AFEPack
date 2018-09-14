/**
 * @file   MPI_PeriodHandler.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Mar 11 16:02:25 2010
 * 
 * @brief  处理周期区域的工具
 *
 * 面对周期区域，我们需要处理的事情包括如下几项：
 *
 * 1. 在进行自适应的时候，几何体的匹配方式需要修改，具有周期对应性质的
 *    几何体要进行匹配；
 *
 * 2. 在进行了负载平衡操作以后，匹配了的几何体因为只存储了一个备份，因
 *    此几何体的坐标可能出错，我们需要对其做一个校正；
 *
 * 3. 在进行自由度匹配的时候，也要对具有周期对应性质的自由度进行归并；
 *
 * 在以上几项均完成的情况下，整个程序仍然可能出错：如果两个具有周期对应
 * 性质的几何体在负载平衡的时候被放到了同一个进程上，那么在数据展开时它
 * 将失去一些不应该归并的拷贝，使得整个数据结构被破坏。这样的情形如何避
 * 免尚未考虑好。
 *
 * 对于前面进行处理的三个步骤，我们可以采用如下的方式完成：
 *
 * a. 在1和3中，进行匹配的程序使用计算点与点之间的距离的模式。我们在这
 *    些地方进行扩展，使用一个型为 match_point(p0, p1, h) 的函数来代替
 *    过去的 distance(p0, p1) < h 的语句，使得当 p0 和 p1 是具有周期对
 *    应性质的点进行判断的时候返回真值；
 *
 * b. 在2的处理中，我们对周期对应的几何体存储几个备份。为了做到这一点，
 *    我们在进行负载平衡时，对所有的周期匹配方式建立的共享关系都忽略。
 * 
 * 我们来详细分析一下做负载平衡时所出现的问题：这个问题原则上是由于共享
 * 的几何体事实上有两种不同的共享方式，但是我们却没有分开考虑造成的。一
 * 种方式是两个共享的几何体本质上就是同一个几何体，其所有的数据都完全一
 * 致；另一种方式是两个共享的几何体从物理上来说是同一个几何体，但是具有
 * 不一样的数据(在周期边界条件下，表现为坐标值不相同)。这个问题事实上在
 * 并非分布式并行的时候就已经存在了，即所谓Physical Entity和Geometric
 * Entity的一对多和多对一的关系处理问题。要彻底解决这个问题，需要从网格
 * 数据结构部分开始修改。所以签于其复杂性，我们暂时先不做此考虑。
 * 
 * 目前，我们将远程共享指针中加入一个整型变量表征所谓"对象共享类型"。
 * 如果两个对象是完全相同的，那么对象必须在不同的进程上，共享类型为 0。
 * 对于共享类型不是 0 的情形，我们认为两个对象不是完全的复制品，并允许
 * 其出现在同一个进程上。
 * 
 * 在做负载平衡的时候，对所有非 0 的共享关系均忽略。因此我们的处理方式
 * 相当于在负载平衡的时候先拆散所有的周期匹配关系，然后在负载平衡完成
 * 以后，将所有周期匹配关系重新建立起来。这将会浪费一些运行时间，但是
 * 可以节省大量的程序开发上的工作量。
 *
 * 注：尽管此处处理的是周期情形，但是对于用户定义的特殊匹配器，其实可以
 * 做非常不同方式的匹配。参考 ex37 中的写法。
 *
 */

#ifndef __MPI_PeriodHandler_h__
#define __MPI_PeriodHandler_h__

AFEPACK_OPEN_NAMESPACE

namespace MPI {
  namespace Periodic {

    /**
     * 对周期区域下的点进行匹配的匹配器。约定匹配函数非负的返回值是匹配
     * 上，负值是没有匹配上。在用户自身写匹配器的情况请注意这一约定。
     */
    template <int DOW>
      struct PointDistance {
        int bmark;
        std::vector<double> period;
        /** 
         * 对点 p0 和 p1 进行坐标匹配，匹配的误差忍量为 tol。
         */
        int value(const Point<DOW>& p0,
                  const Point<DOW>& p1,
                  double tol) const {
          int type = 0;
          for (int i = 0;i < DOW && (type >= 0);++ i) {
            double a = fabs(p0[i] - p1[i]);
            if (a > tol) {
              a = fabs(a - period[i]);
              type = (a < tol)?1:-1;
            }
          }
          return type;
        }
        /** 
         * 读入周期描述的配置文件。
         *
         * 我们定义配置文件的格式如下：
         *
         * bnd_mark # 一个整数，表示需要进行处理的几何体的边界标识，我
         *          # 们要求需要做周期处理的几何体具有相同的边界标识。
         *
         * T1 T2 ... TD # D个实数，表示D个方向上每个方向的周期。
         *
         * 文件中任何位置出现字符 #，那么该字符直到行尾都被视为注释
         *
         */
         void readConfigFile(const std::string& file) {
          filtering_istream is;
          OpenFilteredStream(file, is);
          is >> bmark;
          period.resize(DOW);
          for (int i = 0;i < DOW;++ i) {
            is >> period[i];
          }
        }
      };

    namespace details {

      /**
       * 对几何体geo及其低维几何体中的周期几何体进行收集。
       */
      template <int DIM, int DOW>
        void collectGeometry(HGeometry<DIM,DOW>& geo,
                             int bnd_mark,
                             property_id_t<bool>& pid,
                             std::vector<std::vector<void *> >& bnd_geo) {
        if (geo.get_property(pid) == NULL) {
          if (geo.bmark == bnd_mark) {
            geo.new_property(pid);
            bnd_geo[DIM].push_back((void *)(&geo));
            if (geo.isRefined()) { /// 对孩子进行匹配
              for (u_int i = 0;i < geo.n_child;++ i) {
                collectGeometry(*(geo.child[i]), bnd_mark, pid, bnd_geo);
              }
            }
          }
          if (DIM > 1) { /// 对边界进行匹配
            for (u_int i = 0;i < geo.n_boundary;++ i) {
              collectGeometry(*(geo.boundary[i]), bnd_mark, pid, bnd_geo);
            }
          } else if (DIM == 1) { /// 对顶点进行匹配
            for (u_int i = 0;i < geo.n_vertex;++ i) {
              collectGeometry(*(geo.vertex[i]), bnd_mark, pid, bnd_geo);
            }
          }
        }
      }

      /**
       * 计算一个几何体做匹配时候的参考坐标点。
       */
      template <class GEO>
        void geometry_reference_point(const GEO& geo, 
                                      Point<GEO::dow>& pnt) {
        int n_vtx = geo.n_vertex;
        for (int n = 0;n < GEO::dow;++ n) {
          pnt[n] = 0.0;
          for (int i = 0;i < n_vtx;++ i) {
            pnt[n] += (*geo.vertex[i])[n];
          }
          pnt[n] /= n_vtx;
        }
      }
      template <int DOW>
        void geometry_reference_point(const HGeometry<0,DOW>& geo, 
                                      Point<DOW>& pnt) {
        pnt = geo;
      }

      struct matchTolerence {
        double operator()() const {
          return 1.0e-08; /// 此参量如何传递进来？
        }
      };

      /** 
       * 对本地几何体和远程几何体进行匹配。参数中
       *
       * forest 是几何遗传树；
       * local_* 分别是本地的进程号，几何体的指针数组和参考坐标点。
       * remote_* 分别是远程的进程号，几何体的指针数组和参考坐标点。
       *
       * 需要注意的是远程进程也可能就是本地进程自身，此处是统一处理的。
       * 当匹配的时候，如果匹配上，但是类型为 0 的话，则两个几何体是本
       * 质相同的，应该已经匹配过了，所以此处也屏蔽掉这样的情况。约定
       * 匹配函数非负的返回值是匹配上，负值是没有匹配上。在用户自身写
       * 匹配器的情况请注意这一约定。
       * 
       */
      template <int DIM, class FOREST>
        void matchRankGeometryBuffer(FOREST& forest,
                                     int local_rank,
                                     std::vector<void *>& local_geo,
                                     std::vector<Point<FOREST::dow> >& local_grp,
                                     int remote_rank,
                                     std::vector<void *>& remote_geo,
                                     std::vector<Point<FOREST::dow> >& remote_grp) {
        typedef HGeometry<DIM,FOREST::dow> geometry_t;
        matchTolerence h;

        u_int n_remote_geo = remote_geo.size();
        u_int n_local_geo = local_geo.size();
        for (u_int i = 0;i < n_local_geo;++ i) {
          geometry_t * p_geo_i = (geometry_t *)(local_geo[i]);
          Point<FOREST::dow>& pi = local_grp[i];
          for (u_int j = 0;j < n_remote_geo;++ j) {
            geometry_t * p_geo_j = (geometry_t *)(remote_geo[j]);
            Point<FOREST::dow>& pj = remote_grp[j];

            int type = forest.matcher().value(pi, pj, h());
            if (type > 0) { /// 匹配上但是不是同一个几何体的情形
              Shared_object<geometry_t> * 
                p_info = forest.get_shared_info(*p_geo_i);
              if (p_info == NULL) p_info = forest.new_shared_info(*p_geo_i);
              p_info->add_clone(remote_rank, type, p_geo_j);
            }
            /**
             * 注意如果上面匹配上了以后也还要继续进行匹配，因为一个几
             * 何体在同一个进程上可能有多个不完全一致的拷贝。
             */
          }
        }
      }

      /** 
       * 对第DIM维的周期几何体进行匹配。此处仅仅是发送和接收数据，然后
       * 调用 matchRankGeometryBuffer 完成具体的匹配操作。
       */
      template <int DIM, class FOREST>
        void matchRankGeometry(FOREST& forest,
                               int shift,
                               std::vector<void *>& geo,
                               std::vector<Point<FOREST::dow> >& grp,
                               BinaryBuffer<>& buf, 
                               int tag) {
        if (shift == 0) {
          matchRankGeometryBuffer<DIM,FOREST>(forest, shift, geo, grp, forest.rank(), geo, grp);
        } else {
          BinaryBuffer<> ibuf;
          int rank = forest.rank(), n_rank = forest.n_rank();
          int src = rank - shift;
          int dst = rank + shift;
          if (src < 0) src += n_rank;
          if (dst >= n_rank) dst -= n_rank;
          
          MPI_Comm comm = forest.communicator();
          int send_size = buf.size(), recv_size;
          MPI_Status status;
          MPI_Sendrecv(&send_size, 1, MPI_INT, dst, tag, 
                       &recv_size, 1, MPI_INT, src, tag, 
                       comm, &status);

          ibuf.resize(recv_size);
          MPI_Sendrecv(buf.start_address(), send_size, MPI_CHAR, dst, tag,
                       ibuf.start_address(), recv_size, MPI_CHAR, src, tag, 
                       comm, &status);

          Migration::istream<BinaryBuffer<> > is(ibuf);
          u_int n_remote_geo;
          is >> n_remote_geo;
          if (n_remote_geo > 0) {
            std::vector<void *> remote_geo(n_remote_geo);
            std::vector<Point<FOREST::dow> > remote_grp(n_remote_geo);
            for (u_int i = 0;i < n_remote_geo;++ i) {
              is >> remote_geo[i] >> remote_grp[i];
            }
            matchRankGeometryBuffer<DIM,FOREST>(forest, rank, geo, grp, 
                                                src, remote_geo, remote_grp);
          }
        }
      }

      /**
       * 对DIM维的周期几何体进行匹配。
       */
      template <int DIM, class FOREST>
        void matchGeometry(FOREST& forest,
                           std::vector<void *>& geo) {
        typedef HGeometry<DIM,FOREST::dow> geometry_t;

        /**
         * 准备发送数据的缓冲区。
         */
        BinaryBuffer<> buf;
        Migration::ostream<BinaryBuffer<> > os(buf);
        u_int n = geo.size();
        std::vector<Point<FOREST::dow> > grp(n);
        os << n;
        for (u_int i = 0;i < n;++ i) {
          geometry_t * p_geo = (geometry_t *)(geo[i]);
          geometry_reference_point(*p_geo, grp[i]);
          os << geo[i] << grp[i];
        }

        int tag = 9996;

        int n_rank = forest.n_rank();
        for (int i = 0; i < n_rank;++ i) {
          matchRankGeometry<DIM,FOREST>(forest, i, geo, grp, buf, tag);
        }
      }

    }

    /**
     * 几何遗传树读入匹配器中的周期的值。
     */
    template <class FOREST>
      void readConfigFile(FOREST& forest,
                          const std::string& file) {
      forest.matcher().readConfigFile(file);
    }

    /**
     * 分析出几何遗传树中所有几何体的周期型共享入口。实现的方法为：
     *
     *   1. 将每一维上的所有边界标识为 bnd_mark 的几何体搜集起来；
     *   2. 进程间交换数据；
     *   3. 对每一维的几何体从低维到高维进行比较，对满足周期条件的几何
     *      体进行匹配；
     *
     * 此处的计算量是每一维标识的几何体个数的平方，是比较费时间的算法。
     * 整个匹配的计算量如下：假定总共有P个进程，在秩n进程上，有A_{D,n}
     * 个D-维几何体边界标识为 bnd_mark，则在秩n进程上进行匹配的计算量为
     *
     *   \sum_{D=0}^DIM \sum_{m=0}^P A_{D,n}*A_{D,m}
     *
     * 秩n进程接收和发送的整个数据通讯的量分别为
     *
     *   \sum_{D=0}^DIM \sum_{m=0}^P A_{D,m}
     *
     *   \sum_{D=1}^DIM P*A_{D,n}
     *
     */
    template <class FOREST>
      void periodizeForest(FOREST& forest) {
      int bnd_mark = forest.matcher().bmark;

      /// 对背景网格中的周期几何体进行收集
      std::vector<std::vector<void *> > bnd_geo(FOREST::dim + 1);
      property_id_t<bool> pid;
      new_property_id(pid);
      typename FOREST::RootIterator 
        the_geo = forest.beginRootElement(),
        end_geo = forest.endRootElement();
      for (;the_geo != end_geo;++ the_geo) {
        details::collectGeometry(*the_geo, bnd_mark,
                                 pid, bnd_geo);
      }
      free_property_id(pid);

      /// 对周期几何体进行匹配
      for (int i = 0;i <= FOREST::dim;++ i) {
        switch (i) {
        case 0: details::matchGeometry<0,FOREST>(forest, bnd_geo[i]); break;
        case 1: details::matchGeometry<1,FOREST>(forest, bnd_geo[i]); break;
        case 2: details::matchGeometry<2,FOREST>(forest, bnd_geo[i]); break;
        case 3: details::matchGeometry<3,FOREST>(forest, bnd_geo[i]); break;
        }
      }
    }

    /**
     * 对内存中的几何遗传树，根据说明网格的周期结构的文本型配置文件，并
     * 分析和建立所有周期型共享入口。
     */
    template <class FOREST>
      void periodizeForest(FOREST& forest,
                           const std::string& file) {
      /// 读入配置文件
      readConfigFile(forest, file);
      periodizeForest(forest);
    }

  }
}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_PeriodHandler_h__
/**
 * end of file
 * 
 */
