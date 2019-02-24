/**
 * @file   MPI_UGeometry_archive.h
 * @author Ruo Li <rli@aztec>
 * @date   Sun Nov  1 17:42:32 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MPI_UGeometry_archive_h__
#define __MPI_UGeometry_archive_h__

#include <sstream>

#include <boost/version.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_array.hpp>

#define BOOST_ARCHIVE_SOURCE
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/binary_object.hpp>

#if BOOST_VERSION >= 103801
#if BOOST_VERSION < 106501
#include <boost/serialization/pfto.hpp>
#endif
#include <boost/archive/detail/register_archive.hpp>
#include <boost/archive/impl/archive_serializer_map.ipp>
#else
#include <boost/pfto.hpp>
#include <boost/archive/impl/archive_pointer_oserializer.ipp>
#include <boost/archive/impl/archive_pointer_iserializer.ipp>
#endif 

#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>
#include <boost/archive/impl/basic_binary_oarchive.ipp>
#include <boost/archive/impl/basic_binary_iarchive.ipp>

#include "../../include/Serialization.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <class T>
    struct _is_HGeometry_type {
      enum { value = false };
    };
  template <int DIM, int DOW>
    struct _is_HGeometry_type<AFEPack::HGeometry<DIM,DOW> > {
    enum { value = true };
  };
  template <int DIM, int DOW>
    struct _is_HGeometry_type<AFEPack::HGeometry<DIM,DOW> *> {
    enum { value = true };
  };
  template <int DIM, int DOW>
    struct _is_HGeometry_type<const AFEPack::HGeometry<DIM,DOW> > {
    enum { value = true };
  };
  template <int DIM, int DOW>
    struct _is_HGeometry_type<AFEPack::HGeometry<DIM,DOW>* const> {
    enum { value = true };
  };

#if BOOST_VERSION>105800
#define COMMA_0
#else
#define COMMA_0 , 0
#endif
  
#if BOOST_VERSION>=103400
  template <class LB, 
    class Elem = std::ostream::char_type, 
    class Tr = std::ostream::traits_type>
    class HGeometry_oarchive : 
  public boost::archive::binary_oarchive_impl<HGeometry_oarchive<LB,Elem,Tr>, Elem, Tr> {
    typedef HGeometry_oarchive<LB,Elem,Tr> derived_t;
    typedef boost::archive::binary_oarchive_impl<derived_t,Elem,Tr> base_t;
#else
    template <class LB>
    class HGeometry_oarchive : 
    public boost::archive::binary_oarchive_impl<HGeometry_oarchive<LB> > {
      typedef HGeometry_oarchive<LB> derived_t;
      typedef boost::archive::binary_oarchive_impl<derived_t> base_t;
#endif
#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    public:
#else
      friend class boost::archive::detail::interface_oarchive<derived_t>;
      friend class boost::archive::basic_binary_oarchive<derived_t>;
      friend class boost::archive::basic_binary_oprimitive<derived_t, std::ostream>;
      friend class boost::archive::save_access;
#endif
      template <class T>
        void save_override(T& t, int = 0) {
        base_t::save_override(t COMMA_0);

        BOOST_STATIC_ASSERT(! (_is_HGeometry_type<T>::value) );
      }
      template <int DIM, int DOW>
        void save_override(const HGeometry<DIM,DOW>& geo, int = 0) {
        base_t::save_override(geo.buffer COMMA_0);
        base_t::save_override(geo.index COMMA_0);
        base_t::save_override(geo.bmark COMMA_0);

        this->save_override(geo.parent, 0);
        for (u_int i = 0;i < geo.n_vertex;++ i) {
          this->save_override(geo.vertex[i], 0);
        }
        for (u_int i = 0;i < geo.n_boundary;++ i) {
          this->save_override(geo.boundary[i], 0);
        }
        for (u_int i = 0;i < geo.n_child;++ i) {
          this->save_override(geo.child[i], 0);
        }

        /// 调试用代码
        if (geo.isRefined()) {
          if (lb().is_on_this_new_rank(*geo.child[0], new_rank())) {
            for (u_int i = 1;i < geo.n_child;++ i) {
              assert (lb().is_on_this_new_rank(*geo.child[i], new_rank()));
            }
          } else {
            for (u_int i = 1;i < geo.n_child;++ i) {
              assert ((! lb().is_on_this_new_rank(*geo.child[i], new_rank())));
            }
          }
        }

        lb().set_is_saved(geo);
      }
      template <int DOW>
        void save_override(const HGeometry<0,DOW>& geo, int = 0) {
        base_t::save_override(geo.buffer COMMA_0);
        base_t::save_override(boost::serialization::base_object<Point<DOW> >(geo) COMMA_0);
        base_t::save_override(geo.index COMMA_0);
        base_t::save_override(geo.bmark COMMA_0);

        lb().set_is_saved(geo);
      }
      template <int DIM, int DOW>
        void save_override(HGeometry<DIM,DOW>* const& p_geo, int = 0) {
        int type = 0;
        if (p_geo == NULL) { /// 空指针
          type = 0;
          base_t::save_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
        } else {
          if (! (lb().is_on_this_new_rank(*p_geo, new_rank()))) {
            type = 0;
            base_t::save_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
          } else if (lb().is_save_on_this_rank(*p_geo, new_rank())) {
            unsigned long * global_idx = lb().get_global_idx(*p_geo);
            if (global_idx == NULL) { /// 普通指针
              type = (1<<0);
              base_t::save_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
              base_t::save_override(p_geo COMMA_0);
            } else { /// 本进程上存储的共享型指针
              type = (1<<1);
              base_t::save_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
              base_t::save_override(p_geo COMMA_0);
              base_t::save_override(boost::serialization::make_binary_object(global_idx, sizeof(unsigned long)) COMMA_0);

              typename LB::rank_map_t * p_map = lb().get_rank_map(*p_geo);
              int n_share = p_map->size() - 1;
              base_t::save_override(boost::serialization::make_binary_object(&n_share, sizeof(int)) COMMA_0);
              if (n_share > 0) {
                typename LB::rank_map_t::iterator
                  the_pair = p_map->begin(),
                  end_pair = p_map->end();
                for (;the_pair != end_pair;++ the_pair) {
                  int new_rank = the_pair->first;
                  if (new_rank == this->new_rank()) continue;
                  base_t::save_override(boost::serialization::make_binary_object(&new_rank, sizeof(int)) COMMA_0);
                }
              }
            }
          } else { /// 其他进程上存储的指针
            type = (1<<2);
            base_t::save_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
            unsigned long * global_idx = lb().get_global_idx(*p_geo);
            assert (global_idx != NULL);
            base_t::save_override(boost::serialization::make_binary_object(global_idx, sizeof(unsigned long)) COMMA_0);
          }
        }
      }
    public:
    HGeometry_oarchive(LB& lb, int new_rank, std::ostream& os, unsigned flags = 0) : 
      base_t(os, flags), _lb(&lb), _new_rank(new_rank) {}
    public:
      LB& lb() const { return *_lb; }
      int new_rank() const { return _new_rank; }
    private:
      LB * _lb;
      int _new_rank;
    };


#if BOOST_VERSION>=103400
    template <class LB,
    class Elem = std::istream::char_type, 
    class Tr = std::istream::traits_type>
    class HGeometry_iarchive :
    public boost::archive::binary_iarchive_impl<HGeometry_iarchive<LB,Elem,Tr>,Elem,Tr> {
      typedef HGeometry_iarchive<LB,Elem,Tr> derived_t;
      typedef boost::archive::binary_iarchive_impl<derived_t,Elem,Tr> base_t;
#else
      template <class LB>
      class HGeometry_iarchive :
      public boost::archive::binary_iarchive_impl<HGeometry_iarchive<LB> > {
        typedef HGeometry_iarchive<LB> derived_t;
        typedef boost::archive::binary_iarchive_impl<derived_t> base_t;
#endif
#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
      public:
#else
        friend class boost::archive::detail::interface_iarchive<derived_t>;
        friend class boost::archive::basic_binary_iarchive<derived_t>;
        friend class boost::archive::basic_binary_iprimitive<derived_t, std::istream>;
        friend class boost::archive::load_access;
#endif
        template<class T>
          void load_override(T& t, int = 0){
          base_t::load_override(t COMMA_0);

          BOOST_STATIC_ASSERT(! (_is_HGeometry_type<T>::value) );
        }
        template<int DIM, int DOW>
          void load_override(HGeometry<DIM,DOW> &geo, int = 0){
          base_t::load_override(geo.buffer COMMA_0);
          base_t::load_override(geo.index COMMA_0);
          base_t::load_override(geo.bmark COMMA_0);

          this->load_override(geo.parent, 0);
          for (u_int i = 0;i < geo.n_vertex;++ i) {
            this->load_override(geo.vertex[i], 0);
          }
          for (u_int i = 0;i < geo.n_boundary;++ i) {
            this->load_override(geo.boundary[i], 0);
          }
          for (u_int i = 0;i < geo.n_child;++ i) {
            this->load_override(geo.child[i], 0);
          }
        }
        template <int DOW>
          void load_override(HGeometry<0,DOW>& geo, int = 0) {
          base_t::load_override(geo.buffer COMMA_0);
          base_t::load_override(boost::serialization::base_object<Point<DOW> >(geo) COMMA_0);
          base_t::load_override(geo.index COMMA_0);
          base_t::load_override(geo.bmark COMMA_0);
        }
        template <int DIM, int DOW>
          void load_override(HGeometry<DIM,DOW> *& p_geo, int = 0) {
          int type;
          unsigned long global_idx;
          base_t::load_override(boost::serialization::make_binary_object(&type, sizeof(int)) COMMA_0);
          if (type == 0) { /// 空指针
            p_geo = NULL;
          } else if (type == (1<<0)) { /// 普通指针
            base_t::load_override(p_geo COMMA_0);
          } else if (type == (1<<1)) { /// 共享型指针
            base_t::load_override(p_geo COMMA_0);

            base_t::load_override(boost::serialization::make_binary_object(&global_idx, sizeof(unsigned long)) COMMA_0);
            lb().merge_global_pointer(type, global_idx, p_geo);

            int n_share;
            base_t::load_override(boost::serialization::make_binary_object(&n_share, sizeof(int)) COMMA_0);
            if (n_share > 0) {
              for (int i = 0;i < n_share;++ i) {
                int rank;
                base_t::load_override(boost::serialization::make_binary_object(&rank, sizeof(int)) COMMA_0);
                assert (rank != this->new_rank());
                lb().share_global_pointer(rank, global_idx, p_geo);
              }
            }
          } else { /// 其它进程存储的指针
            assert (type == (1<<2));
            base_t::load_override(boost::serialization::make_binary_object(&global_idx, sizeof(unsigned long)) COMMA_0);
            lb().merge_global_pointer(type, global_idx, p_geo);
          }
        }
      public:
      HGeometry_iarchive(LB& lb, int new_rank, std::istream & is, unsigned flags = 0) :
        base_t(is, flags), _lb(&lb), _new_rank(new_rank) {}
      public:
        LB& lb() const { return *_lb; }
        int new_rank() const { return _new_rank; }
      private:
        LB * _lb;
        int _new_rank;
      };

} // MPI

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_UGeometry_archive_h__

/**
 * end of file
 * 
 */
