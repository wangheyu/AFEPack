/**
 * @file   Migration.details.h
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct 20 07:47:02 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __Migration_details_h__
#define __Migration_details_h__

#include "Migration.h"

AFEPACK_OPEN_NAMESPACE

namespace Migration {
  typedef AFEPack::Migration::data_id_t data_id_t;
  typedef AFEPack::Migration::data_name_t data_name_t;

  namespace details {

    class _mpi_access;

    struct _global_environment {
    private:
      static data_id_t _data_id;
      static std::map<data_id_t, data_name_t> _data_id2name_table;
      static std::map<data_name_t, data_id_t> _data_name2id_table;

    public:
      /**
       * 检索一个数据名称的 ID。
       */
      static data_id_t name_to_id(const data_name_t& dn);
      /**
       * 注册一个数据名称的 ID。
       */
      static data_id_t register_data_name(const data_name_t& dn, bool);
      /**
       * 数据迁移环境进行初始化。
       */
      static void initialize();

      friend class _mpi_access;
    };

  }
}

AFEPACK_CLOSE_NAMESPACE

#endif // __Migration_details_h__

/**
 * end of file
 * 
 */
