/**
 * @file   MPI.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Sat Oct  3 20:48:26 2009
 * 
 * @brief  
 * 
 * 
 */

#include <cstdlib>
#include <mpi.h>

#include "../include/MPI.h"

#ifndef MPI_TAG_UB
#define MPI_TAG_UB 4096
#endif

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  int get_comm_tag(MPI_Comm comm) {
    int rank, tag;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      tag = abs(rand());
      tag %= MPI_TAG_UB; /// 模掉MPI标签值的上界
    }
    MPI_Bcast(&tag, 1, MPI_INT, 0, comm);
    
    return tag;
  }

  int get_comm_int(MPI_Comm comm) {
    int rank, tag;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      tag = abs(rand());
    }
    MPI_Bcast(&tag, 1, MPI_INT, 0, comm);
    
    return tag;
  }

}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
