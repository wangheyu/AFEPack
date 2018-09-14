/**
 * @file   MPI_Controller.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Nov  3 07:30:07 2009
 * 
 * @brief  
 * 
 * 
 */

#include "../include/MPI_Controller.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  namespace __global_controller_environment {
    struct Controller {
      MPI_Comm comm;
      int rank;
      int n_entry;
      details::controller_packer * entry[32];
    };

    struct Controller controller;
  }

  namespace details {
    template <>
    struct controller_packer_impl<void> : public controller_packer {
      void (*fun_ptr)(const char *);
      virtual void call(const char * args) const {
        (*fun_ptr)(args);
      }
    };
  }

  using __global_controller_environment::Controller;

  void details::setControllerEntry(details::controller_packer * cnt_ptr) {
    Controller& controller = __global_controller_environment::controller;
    controller.entry[controller.n_entry] = cnt_ptr;
    controller.n_entry ++;
  }

  void init_controller(MPI_Comm comm)
  {
    Controller& controller = __global_controller_environment::controller;
    controller.comm = comm;
    MPI_Comm_rank(comm, &controller.rank);
  }

  void registerController(const char * name, 
                          void (*function)(const char *), 
                          const char * desc)
  {
    details::controller_packer_impl<void> * 
      cnt_ptr = new details::controller_packer_impl<void>;
    cnt_ptr->name = name;
    cnt_ptr->fun_ptr = function;
    if (desc != NULL) {
      cnt_ptr->description = desc;
    }
    details::setControllerEntry(cnt_ptr);
  }

  int getOneControl(const int i)
  {
    Controller& controller = __global_controller_environment::controller;
    int control, dummy, ret;
    size_t n;
    char * argments = NULL;
    if (controller.rank == 0) {
      char file[128];
      FILE * fp;
      sprintf(file, ".control.%d", i);
      if ((fp = fopen(file, "r")) == NULL) {
        ret = 0;
      } else {
        dummy = fscanf(fp, "%d", &control);
        dummy = getline(&argments, &n, fp);
        fclose(fp);
        sprintf(file, "rm -f .control.%d", i);
        dummy = system(file);
        ret = 1;
      }
    }
    MPI_Bcast(&ret, 1, MPI_INT, 0, controller.comm);
    if (ret == 0) return 0;

    MPI_Bcast(&control, 1, MPI_INT, 0, controller.comm);
    MPI_Bcast(&n, sizeof(size_t), MPI_CHAR, 0, controller.comm);
    if (controller.rank != 0) argments = (char *)malloc(n);
    MPI_Bcast(argments, n, MPI_CHAR, 0, controller.comm);

    controller.entry[control]->call(argments);
    if (argments) free(argments);
    return 1;
  }

  void getControl()
  {
    Controller& controller = __global_controller_environment::controller;
    int i = 1, dummy;
    do {
      if (!getOneControl(i ++)) {
        break;
      }
    } while (1);
  }

  void controllerScript(const char * file,
                        MPI_Comm comm) {
    init_controller(comm);

    Controller& controller = __global_controller_environment::controller;
    if (controller.rank != 0) return;

    int i, dummy;
    FILE * fp;
    char buffer[256];
    fp = fopen(file, "w");
    fprintf(fp, "#! /bin/bash\n");
    fprintf(fp, "\n");
    fprintf(fp, "RETVAL=0\n");
    fprintf(fp, "for ((i = 1;i > 0;i = i+1)); do \n");
    fprintf(fp, "  if ! [ -f .control.$i ]; then\n");
    fprintf(fp, "    file=\".control.$i\" ;\n");
    fprintf(fp, "    echo \"This is the No.\" $i \" unmanaged requests.\" ;\n");
    fprintf(fp, "    i=-1 ;\n");
    fprintf(fp, "  fi\n");
    fprintf(fp, "done\n\n");
    fprintf(fp, "case \"$1\" in\n");
    for (i = 0;i < controller.n_entry;i ++) {
      fprintf(fp, "  %s)\n", controller.entry[i]->name.c_str());
      fprintf(fp, "    echo \"%d $2 $3 $4 $5\" >$file\n", i);
      fprintf(fp, "    ;;\n");
    }
    fprintf(fp, "  *)\n");
    fprintf(fp, "    echo \"Unmatched argument, valid arguments are:\"\n");
    for (i = 0;i < controller.n_entry;i ++) {
      fprintf(fp, "    echo \"\t%s: %s\"\n", controller.entry[i]->name.c_str(), 
              controller.entry[i]->description.c_str());
    }
    fprintf(fp, "    RETVAL=1\n");
    fprintf(fp, "esac\n");
    fprintf(fp, "\n");
    fprintf(fp, "exit $RETVAL\n");
    fclose(fp);
    sprintf(buffer, "chmod 755 %s", file);
    dummy = system(buffer);
    printf("\nUsing %s to manage the job.\n", file);
  }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
