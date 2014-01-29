/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2008  Todd Wease <->

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Based on a an original implementation kindly released by Todd
    Wease.
  */
#pragma once

#include <magnet/containers/iterator_pair.hpp>

namespace magnet {
  namespace containers {
    /*! \brief This container roughly approximates a multimap
      implementation but uses a vector data structure.
      
      The key restriction for this container is that the number of
      keys is known, the keys start at zero, and are sequential. 
    */
    template <typename InnerSet>
    class Vector_Multimap {
      std::vector<InnerSet> _data;
    public:
      typedef typename InnerSet::const_iterator const_iterator;

      void erase(size_t cell, size_t particle) {
        _data[cell].erase(particle);
      }

      void insert(size_t cell, size_t particle) {
	_data[cell].insert(particle);
      }

      typedef magnet::containers::IteratorPairRange<const_iterator> RangeType;
      RangeType getKeyContents(const size_t key) const {
#ifdef MAGNET_DEBUG
	if (key >= _data.size()) M_throw() << "Access out of range (key="<<key << ", size=" << _data.size();
#endif
	return RangeType(_data[key].begin(), _data[key].end());
      }

      void resize(uint32_t keycount) { _data.resize(keycount); }

      void clear() { _data.clear(); }
    };


    /*! \brief This container roughly approximates a multimap
      implementation but uses a set data structure.
      
      The restriction for this container are that the type of the keys
      and mapped data are uint32_t types. You must provide a set
      container which stores a uint64_t type as the first template
      argument. The number of keys is known, the keys start at zero,
      and are sequential.
    */
    template <typename Set>
    class Set_Multimap {
      static_assert(std::is_same<typename Set::value_type, uint64_t>::value, "The Set used must store a uint64_t type");
      static uint64_t toKey(uint32_t cell, uint32_t particle) { return (uint64_t(cell) << 32) | particle; }
      
      Set _data;
    public:
      struct const_iterator : public Set::const_iterator
      {
	const_iterator(const typename Set::const_iterator& iter): Set::const_iterator(iter) {}
	typename Set::const_iterator::value_type operator*() const { return Set::const_iterator::operator*() & (~uint32_t(0)); }
	typename Set::const_iterator::value_type const * operator->() const = delete;
      };

      void erase(uint32_t cell, uint32_t particle) {
        _data.erase(toKey(cell, particle));
      }

      void insert(uint32_t cell, uint32_t particle) {
	_data.insert(toKey(cell, particle));
      }

      typedef magnet::containers::IteratorPairRange<const_iterator> RangeType;
      RangeType getKeyContents(const uint32_t cellID) const {
	return RangeType(_data.lower_bound(toKey(cellID, 0)), _data.lower_bound(toKey(cellID + 1, 0)));
      }

      void resize(uint32_t cellcount) {}
      size_t size() const { return _data.size(); }
      void clear() { _data.clear(); }
    };

  }
}
