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
	ConstJudyIterator& operator++() { return *this = _container.findNext(_value); }
	ConstJudyIterator& operator--() { return *this = _container.findPrev(_value); }
	const value_type& operator*() const { return _value; }
	value_type const * operator->() const { return &_value; }
	bool operator!=(const ConstJudyIterator& o) const { return !(*this == o); }
	bool operator==(const ConstJudyIterator& o) const { return (&_container == &o._container) && (_value == o._value); }
      private:
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
      typedef KeyType key_type;
      typedef key_type value_type;
      static_assert(sizeof(KeyType) == sizeof(Word_t), "Must use a machine word as KeyType");

    private:
      Pvoid_t _array = (Pvoid_t) NULL;
      size_t _count = 0;

    public:
      typedef detail::ConstJudyIterator<JudySet, value_type> const_iterator;

      const_iterator find(key_type key) const { 
	if (Judy1Test(_array, key, NULL))
	  return const_iterator(*this, key);
	else
	  return end();
      }
      const_iterator findFirst(key_type key) const {
	if (Judy1First(_array, &key, NULL) != 0)
	  return const_iterator(*this, key);
	else 
	  return end();
      }

      const_iterator findNext(key_type key) const {
	if (Judy1Next(_array, &key, NULL) != 0)
	  return const_iterator(*this, key);
	else 
	  return end();
      }

      const_iterator findPrev(key_type key) const {
	if (Judy1Prev(_array, &key, NULL) != 0)
	  return const_iterator(*this, key);
	else
	  return end();
      }

      const_iterator findNth(const size_t i) {
	key_type key;
	Judy1ByCount(_array, i+1, &key, NULL);
	return const_iterator(*this, key);
      }
      
      const_iterator begin() const { return findFirst(0); }
      const_iterator end() const { return const_iterator(*this, endIndex); }
      
    public:
      ~JudySet() { clear(); }
      void clear() { Judy1FreeArray(&_array, NULL); _count = 0; }
      inline size_t size() const { return _count; }
      inline void insert(const key_type key) { _count += Judy1Set(&_array, key, NULL); }
      inline void erase(const key_type key) { _count -= Judy1Unset(&_array, key, NULL); }
      inline bool count(const key_type key) { return find(key) != end(); }
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

    private:
      Pvoid_t _array = (Pvoid_t) NULL;
      size_t _count = 0;

    public:
      typedef detail::ConstJudyIterator<JudyMap, value_type> const_iterator;
      
      const_iterator find(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLGet(_array, idx, NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator findFirst(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLFirst(_array, (Word_t*)&idx, NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator findNext(value_type val) const { return findNext(val.first); }
      const_iterator findNext(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLNext(_array, (Word_t*)&idx, NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      const_iterator findPrev(value_type val) const { return findNext(val.first); }
      const_iterator findPrev(key_type idx) const {
	mapped_type* ptr = (mapped_type*)JudyLPrev(_array, (Word_t*)&idx, NULL);
	if (ptr != NULL)
	  return const_iterator(*this, value_type(idx, *ptr));
	else
	  return end();
      }

      inline void insert(const value_type value) {
	(*this)[value.first] = value.second;
      }

      inline void erase(const key_type key) {
	_count -= JudyLDel(&_array, key, NULL);
      }

      const_iterator begin() const { return findFirst(0); }
      const_iterator end() const { return const_iterator(*this, value_type(endIndex, 0)); }

      ~JudyMap() { clear(); }
      void clear() { JudyLFreeArray(&_array, NULL); _count = 0; }
      size_t size() const { return _count; }
      bool count(key_type key) const { return find(key) != end(); }

      const_iterator findNth(const size_t i) const {
	key_type index;
	mapped_type* PValue = (mapped_type*)JudyLByCount(_array, i+1, &index, NULL);
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
      inline EntryProxy operator[](key_type key) {
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
