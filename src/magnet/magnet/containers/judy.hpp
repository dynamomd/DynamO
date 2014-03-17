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

#include <Judy.h>
#include <magnet/exception.hpp>

namespace magnet {
  namespace containers {
    namespace detail {
      /*! \brief A class used to provide read-only iterators for 
          Judy array sets and maps.
       */
      template<typename Container, typename value_type>
      class ConstJudyIterator : public std::iterator<std::bidirectional_iterator_tag, value_type>
      {
      public:
	ConstJudyIterator(const Container& container, value_type value):
	  _container(container), _value(value) {}
	ConstJudyIterator& operator++() { return *this = _container.next(*this); }
	ConstJudyIterator& operator--() { return *this = _container.prev(*this); }
	const value_type& operator*() const { return _value; }
	value_type const * operator->() const { return &_value; }
	bool operator!=(const ConstJudyIterator& o) const { return !(*this == o); }
	bool operator==(const ConstJudyIterator& o) const { 
#ifdef MAGNET_DEBUG
	  if (&_container != &o._container) M_throw() << "Cannot compare iterators from different containers";
#endif
	  return _value == o._value; 
	}
      protected:
	//Private as we can't allow iterator assignment as we have a
	//reference to the container
	ConstJudyIterator& operator=(const ConstJudyIterator& o) {
#ifdef MAGNET_DEBUG
	  if (&_container != &o._container) M_throw() << "Cannot assign iterators from different containers";
#endif
	  _value = o._value;
	  return *this;
	}
	const Container& _container;
	value_type _value;
      };
    }

    /*! \brief A Judy array set which contains types which are the
      size of a machine word (32/64 bit). 
      
      The stored type is reinterpreted as an unsigned integer and is
      sorted in increasing order.
     */
    template<typename KeyType = Word_t, Word_t endIndex = ~Word_t(0)>
    class JudySet
    {
    public:
      JudySet(): _array((Pvoid_t) NULL), _count(0) {}

      typedef KeyType key_type;
      typedef key_type value_type;
      static_assert(sizeof(KeyType) == sizeof(Word_t), "Must use a machine word as KeyType");

    private:
      Pvoid_t _array;
      size_t _count;

    public:
      typedef detail::ConstJudyIterator<JudySet, value_type> const_iterator;

      const_iterator find(const key_type& key) const { 
	if (Judy1Test(_array, key, NULL))
	  return const_iterator(*this, key);
	else
	  return end();
      }

      const_iterator lower_bound(key_type key) const {
	if (Judy1First(_array, reinterpret_cast<Word_t*>(&key), NULL) != 0)
	  return const_iterator(*this, key);
	else 
	  return end();
      }

      const_iterator upper_bound(key_type key) const { return next(key); }

      const_iterator next(const const_iterator& it) const {
	key_type key = *it;
	if (Judy1Next(_array, reinterpret_cast<Word_t*>(&key), NULL) != 0)
	  return const_iterator(*this, key);
	else 
	  return end();
      }

      const_iterator prev(const const_iterator& it) const {
	key_type key = *it;
	if (Judy1Prev(_array, reinterpret_cast<Word_t*>(&key), NULL) != 0)
	  return const_iterator(*this, key);
	else
	  return end();
      }

      std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const {
	const_iterator it1 = find(key);
	if (it1 == end()) 
	  return std::pair<const_iterator, const_iterator>(end(), end());
	const_iterator it2 = it1;
	return std::pair<const_iterator, const_iterator>(it1, ++it2);
      }

      const_iterator findNth(const size_t i) {
	key_type key;
	if (Judy1ByCount(_array, i+1, reinterpret_cast<Word_t*>(&key), NULL))
	  return const_iterator(*this, key);
	else
	  return end();
      }
      
      const_iterator begin() const { return lower_bound(0); }
      const_iterator end() const { return const_iterator(*this, endIndex); }
      
    public:
      ~JudySet() { clear(); }
      void clear() { Judy1FreeArray(&_array, NULL); _count = 0; }
      size_t size() const { return _count; }
      bool empty() const { return _count == 0; }
      void insert(const key_type key) { _count += Judy1Set(&_array, key, NULL); }
      void erase(const key_type key) { _count -= Judy1Unset(&_array, key, NULL); }
      bool count(const key_type key) { return find(key) != end(); }
    };

    /*! \brief A Judy array map which between types which are the size
      of a machine word (32/64 bit).
      
      The key and stored type are is reinterpreted as an unsigned
      integer and is sorted in increasing order.
     */
    template<typename KeyType = Word_t, typename MappedType = Word_t, Word_t endIndex = ~Word_t(0)>
    class JudyMap
    {
    public:
      static_assert(sizeof(KeyType) == sizeof(Word_t), "Must use a machine word as KeyType");
      static_assert(sizeof(MappedType) == sizeof(Word_t), "Must use a machine word as MappedType");
      typedef KeyType key_type;
      typedef MappedType mapped_type;
      typedef std::pair<key_type, mapped_type> value_type;

      JudyMap(): _array((Pvoid_t) NULL), _count(0) {}
    private:
      Pvoid_t _array;
      size_t _count;

    public:
      typedef detail::ConstJudyIterator<JudyMap, value_type> const_iterator;
      
      const_iterator find(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLGet(_array, idx, NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator lower_bound(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLFirst(_array, reinterpret_cast<Word_t*>(&idx), NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator upper_bound(key_type idx) const {
	return next(idx);
      }

      const_iterator next(const const_iterator& it) const { return next(it->first); }
      const_iterator next(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLNext(_array, reinterpret_cast<Word_t*>(&idx), NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator prev(const const_iterator& it) const { return prev(it->first); }
      const_iterator prev(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLPrev(_array, reinterpret_cast<Word_t*>(&idx), NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      void insert(const value_type value) {
	(*this)[value.first] = value.second;
      }

      void erase(const key_type key) {
	_count -= JudyLDel(&_array, key, NULL);
      }

      const_iterator begin() const { return lower_bound(0); }
      const_iterator end() const { return const_iterator(*this, value_type(endIndex, 0)); }

      std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const {
	const_iterator it1 = find(key);
	if (it1 == end())
	  return std::pair<const_iterator, const_iterator>(end(), end());
	const_iterator it2 = it1;
	return std::pair<const_iterator, const_iterator>(it1, ++it2);
      }

      ~JudyMap() { clear(); }
      void clear() { JudyLFreeArray(&_array, NULL); _count = 0; }
      size_t size() const { return _count; }
      size_t empty() const { return _count == 0; }
      bool count(key_type key) const { return find(key) != end(); }

      const_iterator findNth(const size_t i) const {
	key_type index;
	mapped_type* PValue = (mapped_type*)JudyLByCount(_array, i+1, reinterpret_cast<Word_t*>(&index), NULL);
	if (PValue != NULL)
	  return const_iterator(*this, value_type(index, *PValue));
	else
	  return end();
      }

      /*!\brief This proxy is used to double check if an assignment of
	zero is done, and delete the entry if it is. */
      struct EntryProxy {
      public:
	EntryProxy(JudyMap& container, const key_type& key):_container(container), _key(key) {}
	operator const mapped_type() const {
	  return *_container.getPtr(_key);
	}
	EntryProxy& operator=(mapped_type newval) {
	  *(_container.getPtr(_key)) = newval;
	  return *this;
	}
      private:
	JudyMap& _container;
	const key_type _key;
      };

      /*! \brief This non-const array access operator uses EntryProxy
	to check if any values assigned are zero so they may be deleted. */
      EntryProxy operator[](key_type key) {
	return EntryProxy(*this, key); 
      }

    protected:
      mapped_type* getPtr(key_type key) {
	//Try finding it first
	mapped_type* PValue = (mapped_type*)JudyLGet(_array, key, NULL);
	if (PValue == NULL)
	  {//Key does not exist
	    ++_count;
	    PValue = (mapped_type*)JudyLIns(&_array, key, NULL);
	  }
	return PValue;
      }
    };
  }
}
