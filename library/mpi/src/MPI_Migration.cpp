/**
 * @file   MPI_Migration.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Oct 19 16:03:55 2009
 * 
 * @brief  
 * 
 * 
 */

#include "../include/MPI_Migration.h"

AFEPACK_OPEN_NAMESPACE

namespace Migration {
  namespace details {

    class _mpi_access {
      static MPI_Comm comm;
    public:
      static void initialize(MPI_Comm comm);
      static data_id_t register_data_name(const data_name_t& dn);
      static void load_config(const std::string& filename);
      static void save_config(const std::string& filename); 
      static void ensured_open_fstream(const std::string& filename,
                                       std::ifstream& is);
      static void ensured_open_filtering_stream(const std::string& filename,
                                                filtering_istream& is);
   };

    MPI_Comm _mpi_access::comm;

    void _mpi_access::initialize(MPI_Comm _comm) {
      comm = _comm;
      _global_environment::initialize();
    }

    data_id_t _mpi_access::register_data_name(const data_name_t& dn) {
      std::map<data_name_t, data_id_t>::const_iterator 
        it = _global_environment::_data_name2id_table.find(dn);

      int n_rank;
      MPI_Comm_size(comm, &n_rank);

      int is_rank_found = (it != _global_environment::_data_name2id_table.end())?1:0;
      int is_found = 0;
      MPI_Allreduce(&is_rank_found, &is_found, 1, MPI_INT, MPI_SUM, comm);
      assert ((is_found == 0 || is_found == n_rank));
      if (is_found == n_rank) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
          std::cerr << "Warning: data name \"" << dn
                    << "\" is already registered." << std::endl;
        }
        return it->second;
      } else {
        data_id_t max_data_id;
        MPI_Allreduce(&_global_environment::_data_id, &max_data_id, 1, 
                      MPI_INT, MPI_MAX, comm);
        _global_environment::_data_id = max_data_id;

        _global_environment::_data_name2id_table[dn] = _global_environment::_data_id;
        _global_environment::_data_id2name_table[_global_environment::_data_id] = dn;
        it = _global_environment::_data_name2id_table.find(dn);
        _global_environment::_data_id += 1;
        return it->second;
      }
    }

    void _mpi_access::ensured_open_fstream(const std::string& filename,
                                           std::ifstream& is) {
      int dummy;
      do {
        is.open(filename.c_str());
        if (! is) {
          is.clear();
          dummy = system("sync");
        } else break;
      } while (true);
    }

    void _mpi_access::ensured_open_filtering_stream(const std::string& filename,
                                                    filtering_istream& is) {
      int dummy;
      do {
        OpenFilteredStream(filename, is);
        if (! is) {
          is.reset();
          dummy = system("sync");
        } else break;
      } while (true);
    }

    void _mpi_access::load_config(const std::string& dirname) {
      char filename[1024];
      sprintf(filename, "%s/.migration.cfg", dirname.c_str());
      filtering_istream is;
      ensured_open_filtering_stream(filename, is);
      
      int n_item;
      is >> n_item;

      _global_environment::initialize();
      data_id_t max_data_id = 0;
      for (int i = 0;i < n_item;++ i) {
        data_id_t data_id;
        std::string data_name;
        is >> data_id >> data_name;
        _global_environment::_data_name2id_table[data_name] = data_id;
        _global_environment::_data_id2name_table[data_id] = data_name;
        if (max_data_id < data_id) {
          max_data_id = data_id;
        }
      }
      _global_environment::_data_id = max_data_id + 1;
    }

    void _mpi_access::save_config(const std::string& dirname) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      if (rank == 0) {
        char filename[1024];
        sprintf(filename, "%s/.migration.cfg", dirname.c_str());
        std::ofstream os(filename);
        int n_item = _global_environment::_data_name2id_table.size();
        os << n_item << " \t# number of entries\n";
        std::map<data_name_t, data_id_t>::const_iterator 
          it = _global_environment::_data_name2id_table.begin();
        for (int i = 0;i < n_item;++ i, ++ it) {
          const data_id_t& id = it->second;
          const data_name_t& name = it->first;
          os << id << " \t" << name << " \t# data ID/NAME pair\n" ;
        }
        os.close();
      }
    }

  }

  void initialize(MPI_Comm comm) {
    details::_mpi_access::initialize(comm);
  }

  data_id_t register_data_name(const data_name_t& dn) {
    return details::_mpi_access::register_data_name(dn);
  }

  void load_config(const std::string& filename) {
    details::_mpi_access::load_config(filename);
  }

  void save_config(const std::string& filename) {
    details::_mpi_access::save_config(filename);
  }

  void ensured_open_fstream(const std::string& filename,
                            std::ifstream& is) {
    details::_mpi_access::ensured_open_fstream(filename, is);
  }

  void ensured_open_filtering_stream(const std::string& filename,
                                     filtering_istream& is) {
    details::_mpi_access::ensured_open_filtering_stream(filename, is);
  }

}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
