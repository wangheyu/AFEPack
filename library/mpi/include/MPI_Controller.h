/**
 * @file   MPI_Controller.h
 * @author Ruo Li <rli@aztec>
 * @date   Tue Nov  3 07:29:05 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MPI_Controller_h__
#define __MPI_Controller_h__

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <mpi.h>

#include "../../include/Miscellaneous.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  namespace details {
    struct controller_packer {
      std::string name;
      std::string description;
      virtual void call(const char * args) const = 0;
      virtual ~controller_packer() {}
    };
    template <class T>
      struct controller_packer_impl : public controller_packer {
      T * obj_ptr;
      void (T::*fun_ptr)(const char *);
      virtual void call(const char * args) const {
        (obj_ptr->*fun_ptr)(args);
      }
      virtual ~controller_packer_impl() {}
    };
    template <class T>
      struct controller_packer_impl<const T> : public controller_packer {
      const T * obj_ptr;
      void (T::*fun_ptr)(const char *) const;
      virtual void call(const char * args) const {
        (obj_ptr->*fun_ptr)(args);
      }
      virtual ~controller_packer_impl() {}
    };
    void setControllerEntry(controller_packer *);
  }

  void registerController(const char * name,
                          void (*function)(const char *),
                          const char * desc = NULL);
  void controllerScript(const char * filename, MPI_Comm = MPI_COMM_WORLD);
  void getControl();

  template <class T>
    void registerController(const char * name,
                            T * obj,
                            void (T::*function)(const char *),
                            const char * desc = NULL)
    {
      details::controller_packer_impl<T> * 
        cnt_ptr = new details::controller_packer_impl<T>;
      cnt_ptr->name = name;
      cnt_ptr->obj_ptr = obj;
      cnt_ptr->fun_ptr = function;
      if (desc != NULL) {
        cnt_ptr->description = desc;
      }
      details::setControllerEntry(cnt_ptr);
    }

  template <class T>
    void registerController(const char * name,
                            const T * obj,
                            void (T::*function)(const char *) const,
                            const char * desc = NULL)
    {
      details::controller_packer_impl<const T> * 
        cnt_ptr = new details::controller_packer_impl<const T>;
      cnt_ptr->name = name;
      cnt_ptr->obj_ptr = obj;
      cnt_ptr->fun_ptr = function;
      if (desc != NULL) {
        cnt_ptr->description = desc;
      }
      details::setControllerEntry(cnt_ptr);
    }
}

AFEPACK_CLOSE_NAMESPACE

#endif /* __MPI_Controller_h__ */

/**
 * end of file
 * 
 */
