/**
 * @file   main.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 23:22:55 2009
 * 
 * @brief  
 * 
 * 
 */

#include "SCL.h"

extern void init_data(SCL& scl, int argc, char * argv[]);
SCL app;

void save_data(const char * dummy) {
  app.save_data();
}

void dump_data(const char * dummy) {
  int n_rank;
  std::string dir;
  std::istringstream is(dummy);
  is >> n_rank >> dir;
  app.dump(n_rank, dir);
}

void load_data(const char * dummy) {
  std::string dir;
  std::istringstream is(dummy);
  is >> dir;
  app.load(dir);
}

int main(int argc, char * argv[])
{
  try {
    MPI_Init(&argc, &argv);

    AFEPack::MPI::registerController("save", &save_data, "save solution.");
    AFEPack::MPI::registerController("dump", &dump_data, "dump the runtime environment.");
    AFEPack::MPI::registerController("load", &load_data, "load the runtime environment.");
    AFEPack::MPI::controllerScript("run.sh", MPI_COMM_WORLD);

    init_data(app, argc, argv);
    app.run();

    MPI_Finalize();
  }
  catch(std::exception& e) {
    std::cerr << "Exception caughted:" << std::endl
      << e.what ()
      << std::endl;
  }
  catch(...) {
    std::cerr << "Exception caughted:" << std::endl
      << "unknown exception caughted."
      << std::endl;
  }

  return 0;
};

/**
 * end of file
 * 
 */

