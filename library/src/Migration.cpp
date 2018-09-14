/**
 * @file   Migration.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Sep 24 14:38:11 2009
 * 
 * @brief  
 * 
 * 
 */

#include <map>
#include "HGeometry.h"

#include "Migration.details.h"

AFEPACK_OPEN_NAMESPACE

namespace Migration {
  namespace details {

    data_id_t _global_environment::_data_id;
    std::map<data_id_t, data_name_t> _global_environment::_data_id2name_table;
    std::map<data_name_t, data_id_t> _global_environment::_data_name2id_table;

    data_id_t _global_environment::name_to_id(const data_name_t& dn) {
      std::map<data_name_t, data_id_t>::const_iterator 
        it = _data_name2id_table.find(dn);
      if (it == _data_name2id_table.end()) {
        std::cerr << "Warning: data name \"" << dn
                  << "\" is not registered." << std::endl;
        return -1;
      }
      return it->second;
    }

    data_id_t _global_environment::register_data_name(const data_name_t& dn, 
                                                      bool flag) {
      std::map<data_name_t, data_id_t>::const_iterator 
        it = _data_name2id_table.find(dn);
      if (it != _data_name2id_table.end()) {
        std::cerr << "Data name \"" << dn
                  << "\" is already registered." << std::endl;
        return it->second;
      } else {
        _data_name2id_table[dn] = _data_id;
        _data_id2name_table[_data_id ++] = dn;
        it = _data_name2id_table.find(dn);
        return it->second;
      }
    }

    void _global_environment::initialize() {
      _data_id = 0;
      _data_name2id_table.clear();
      _data_id2name_table.clear();
    }

  }

  data_id_t name_to_id(const data_name_t& dn) {
    return details::_global_environment::name_to_id(dn);
  }
  data_id_t register_data_name(const data_name_t& dn, bool flag) {
    return details::_global_environment::register_data_name(dn, flag);
  }
  void initialize() {
    details::_global_environment::initialize();
  }
  bool is_valid(const data_id_t& id) {
    return (id >= 0);
  }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
