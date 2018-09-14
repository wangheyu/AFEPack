/**
 * @file   MPI_LoadBalance.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu May  6 08:27:21 2010
 * 
 * @brief  
 * 
 * 
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <mpi.h>

#include "../include/MPI.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  namespace lb_details {
    struct sync_file_buf {
      enum {SEND, RECV} type;
      std::string filename;
      long filesize;
      void * buffer;
      MPI_Request * request;
    };

    void send_file_begin(MPI_Comm comm,
                         int target,
                         int tag,
                         const char * filename,
                         long filesize,
                         sync_file_buf& sfb) {
      size_t dummy;
      sfb.type = sync_file_buf::SEND;
      sfb.buffer = malloc(filesize);

      FILE * fp = fopen(filename, "rb");
      dummy = fread(sfb.buffer, filesize, 1, fp);
      fclose(fp);

      MPI_Isend(sfb.buffer, filesize, MPI_CHAR, target, tag, comm, sfb.request);
    }
    void send_file_end(sync_file_buf& sfb) {
      int dummy;
      free(sfb.buffer);
      char command[2048];
      sprintf(command, "rm -f %s", sfb.filename.c_str());
      dummy = system(command);
    }

    void recv_file_begin(MPI_Comm comm,
                         int source,
                         int tag,
                         const char * filename,
                         long filesize,
                         sync_file_buf& sfb) {
      sfb.type = sync_file_buf::RECV;
      sfb.buffer = malloc(filesize);
      sfb.filename = filename;
      sfb.filesize = filesize;

      MPI_Irecv(sfb.buffer, filesize, MPI_CHAR, source, tag, comm, sfb.request);
    }
    void recv_file_end(sync_file_buf& sfb) {
      size_t dummy;
      FILE * fp = fopen(sfb.filename.c_str(), "wb");
      dummy = fwrite(sfb.buffer, sfb.filesize, 1, fp);
      fclose(fp);
      free(sfb.buffer);
    }

    long get_file_size(const char * filename) {
      struct stat st;
      long length = 0;
      if (stat(filename, &st) == 0) {
        length = st.st_size;
      }
      return length;
    }

    bool is_file_exist(const char * filename) {
      struct stat st;
      return (stat(filename, &st) == 0);
    }

    void bcast_small_file(MPI_Comm comm,
                          int root,
                          int tag,
                          const char * filename) {
      int rank, dummy;
      long filesize;
      MPI_Comm_rank(comm, &rank);
      if (rank == root) {
        filesize = get_file_size(filename);
      }
      MPI_Bcast(&filesize, 1, MPI_LONG, root, comm);

      void * buffer = malloc(filesize);
      if (rank == root) {
        FILE * fp = fopen(filename, "rb");
        dummy = fread(buffer, filesize, 1, fp);
        fclose(fp);
      }
      MPI_Bcast(buffer, filesize, MPI_CHAR, root, comm);

      if (rank != root) {
        if (! is_file_exist(filename)) {
          FILE * fp = fopen(filename, "wb");
          dummy = fwrite(buffer, filesize, 1, fp);
          fclose(fp);
        }
      }
      free(buffer);
    }

  } /// namespace lb_details

  /**
   * 假定在数据存储的时候，所有数据都存储在每个计算节点的局部存储上。本
   * 函数将所有这些文件都收集到网络文件系统上的同一个目录中。此实现中直
   * 接使用shell命令来完成，数据传输由网络文件系统来完成，事实上不可能
   * 做到实时，而且可能有很久的延时。
   */
  void lb_collect_local_data_dir(MPI_Comm comm,
                                 const std::string& src_dir,
                                 const std::string& dst_dir) {
    int rank, dummy;
    MPI_Comm_rank(comm, &rank);
    char command[2048];
    if (rank == 0) {
      sprintf(command, "mkdir -p %s", dst_dir.c_str());
      dummy = system(command);
      sprintf(command, "cp -f %s/.config %s/", src_dir.c_str(), dst_dir.c_str());
      dummy = system(command);
      sprintf(command, "cp -f %s/.migration.cfg %s/", src_dir.c_str(), dst_dir.c_str());
      dummy = system(command);
    }
    MPI_Barrier(comm);
    sprintf(command, 
            "for i in `ls %s`; \
               do mkdir -p %s/$i ; \
               file=\"%s/$i/%d.dat\" ; \
               if [ -f $file ]; then \
                 cp $file %s/$i/ ; \
               fi \
             done", 
            src_dir.c_str(), 
            dst_dir.c_str(), 
            src_dir.c_str(), rank, 
            dst_dir.c_str());
    dummy = system(command);
    MPI_Barrier(comm);
  }

  /**
   * 在做负载平衡的时候，假定数据被存储在每个计算节点的局部存储上。我
   * 们使用下面的函数将进程之间的文件进行交换，使得每个进程都拿到其需
   * 要载入的数据文件。
   */
  void lb_sync_local_data_dir(MPI_Comm comm,
                              const std::string& src_dir,
                              const std::string& dst_dir) {
    int rank, n_rank, dummy;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &n_rank);

    std::vector<long> send_size(n_rank);
    std::vector<long> recv_size(n_rank);

    std::vector<MPI_Request> request(2*n_rank);
    std::vector<MPI_Status> status(2*n_rank);

    char filename[2048], command[2048];

    int tag = get_comm_tag(comm);
    for (int i = 0;i < n_rank;++ i) {
      sprintf(filename, "%s/%d/%d.dat", src_dir.c_str(), i, rank);
      send_size[i] = lb_details::get_file_size(filename);
    }
    MPI_Alltoall(&send_size[0], 1, MPI_LONG,
                 &recv_size[0], 1, MPI_LONG,
                 comm);

    /**
     * 先将目录建立好并进行同步，使得后面可以节省数据传输的开销。
     */
    sprintf(command, "mkdir -p %s/%d", dst_dir.c_str(), rank);
    dummy = system(command);
    MPI_Barrier(comm);

    int n_req = 0;
    std::vector<lb_details::sync_file_buf> sfb(2*n_rank);
    for (int i = 0;i < n_rank;++ i) {     
      char src_file[2048], dst_file[2048];
      sprintf(src_file, "%s/%d/%d.dat", src_dir.c_str(), i, rank);
      sprintf(dst_file, "%s/%d/%d.dat", dst_dir.c_str(), rank, i);

      if (send_size[i] > 0) {
	char dst_rank_dir[2048];
        sprintf(dst_rank_dir, "%s/%d", dst_dir.c_str(), i);
        /**
         * 如果目标目录存在，说明 i 和 rank 共用同样的存储空间，此时本
         * 进程直接将旧文件移动成为新文件即可。
         */
        if (lb_details::is_file_exist(dst_rank_dir)) {
          sprintf(command, "mv %s %s/", src_file, dst_rank_dir);
          dummy = system(command);
        } else { /// 否则通过网络发送数据
          sfb[n_req].request = &request[n_req];
          lb_details::send_file_begin(comm, i, tag, src_file, send_size[i], sfb[n_req ++]);
        }
      }
      if (recv_size[i] > 0) {
	char src_rank_dir[2048];
        sprintf(src_rank_dir, "%s/%d", dst_dir.c_str(), i);
        /**
         * 如果源目录存在，则说明 i 和 rank 共用同样的存储空间，此时进
         * 程 i 将会负责将旧文件移动成为新文件，这里什么都不必做。
         */
        if (! lb_details::is_file_exist(src_rank_dir)) {
          sfb[n_req].request = &request[n_req];
          lb_details::recv_file_begin(comm, i, tag, dst_file, recv_size[i], sfb[n_req ++]);
        }
      }
    }

    MPI_Waitall(n_req, &request[0], &status[0]);
    for (int i = 0;i < n_req;++ i) {
      switch (sfb[i].type) {
      case lb_details::sync_file_buf::SEND:
        lb_details::send_file_end(sfb[i]); break;
      case lb_details::sync_file_buf::RECV:
        lb_details::recv_file_end(sfb[i]); break;
      default:
        ; /// 不会运行到这里的
      }
    }

    /**
     * 最后对目录中的两个配置文件进行同步。目前，配置文件的名称为
     * .config 和 .migration.cfg。
     */
    if (rank == 0) {
      sprintf(command, "cp %s/.config %s/.config", src_dir.c_str(), dst_dir.c_str());
      dummy = system(command);
      sprintf(command, "cp %s/.migration.cfg %s/.migration.cfg", src_dir.c_str(), 
              dst_dir.c_str());
      dummy = system(command);
    }
    lb_details::bcast_small_file(comm, 0, tag, (dst_dir + "/.config").c_str());
    lb_details::bcast_small_file(comm, 0, tag, (dst_dir + "/.migration.cfg").c_str());
    dummy = system("sync");
  }

}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
