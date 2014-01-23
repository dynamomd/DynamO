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
    class JudyPairSet
    {
    public:
      typedef std::pair<size_t, size_t> key_type;
      typedef key_type value_type;

    private:
      static const Word_t endIndex = ~Word_t(0);
      static const size_t half_shift = sizeof(Word_t) * 4;
      static const size_t mask = (size_t(1) << half_shift) - 1;

      static inline Word_t keyToID(key_type key) {
	return Word_t((std::min(key.first, key.second) << half_shift) | std::max(key.first, key.second)); 
      }
      
      static inline key_type IDtoKey(Word_t key) {
	return key_type((size_t(key) >> half_shift) & mask,  size_t(key) & mask);
      }

    public:
      class const_iterator : public std::iterator<std::bidirectional_iterator_tag, value_type> {
      public:
	const_iterator(const JudyPairSet& container, Word_t idx = 0): _container(container), _idx(idx) {
	  
	  int Rc = Judy1First(_container._array, &_idx, NULL);
	  if (Rc == 0) _idx = endIndex;
	}

	const_iterator& operator++() { 
	  int Rc = Judy1Next(_container._array, &_idx, NULL);
	  if (Rc == 0) _idx = endIndex;
	  return *this;
	}

	const_iterator& operator--() {
	  int Rc = Judy1Prev(_container._array, &_idx, NULL);
	  if (Rc == 0) _idx = endIndex;
	  return *this;
	}

	value_type operator*() const { return IDtoKey(_idx); }
	bool operator!=(const const_iterator& o) const { return !(*this == o); }
	bool operator==(const const_iterator& o) const { return _idx == o._idx;}

      private:
	const JudyPairSet& _container;
	Word_t _idx;
      };

      const_iterator begin() const { return const_iterator(*this); }
      const_iterator end() const { return const_iterator(*this, endIndex); }
      
      Pvoid_t _array = (Pvoid_t) NULL;
      size_t _count = 0;
    public:
      ~JudyPairSet() { clear(); }

      void clear() { Judy1FreeArray(&_array, NULL); _count = 0; }

      inline size_t size() const { return _count; }
      
      inline void insert(const key_type key) {
	_count += Judy1Set(&_array, keyToID(key), NULL);
      }
      
      inline void erase(const key_type key) {
	_count -= Judy1Unset(&_array, keyToID(key), NULL);
      }
      
      inline bool count(const key_type key) {
	return Judy1Test(_array, keyToID(key), NULL);
      }
      
      inline key_type operator[](size_t i) {
	size_t ID;
	Judy1ByCount(_array, i+1, &ID, NULL);
	return IDtoKey(ID);
      }
    };

    /*!\brief This is a container that stores a single size_t
       identified by a pair of size_t values.
       
       To efficiently store this state, a Judy array is used and
       entries are only stored if the state is non-zero.  To
       facilitate the storage only if non-zero behaviour, the array
       access operator is overloaded to automatically return a size_t
       0 for any entry which is missing. It also returns a proxy which
       deletes entries when they are set to 0.
    */
    class JudyPairMap
    {
    public:
      typedef std::pair<size_t, size_t> key_type;
      typedef Word_t mapped_type;
      typedef std::pair<key_type, mapped_type> value_type;
      static const Word_t endIndex = ~Word_t(0);

    private:
      static const size_t half_shift = sizeof(Word_t) * 4;
      static const size_t mask = (size_t(1) << half_shift) - 1;
      
      Pvoid_t _array = (Pvoid_t) NULL;
      size_t _count = 0;
      
      static inline Word_t keyToID(key_type key) {
	return Word_t((std::min(key.first, key.second) << half_shift) | std::max(key.first, key.second)); 
      }
  
      static inline key_type IDtoKey(Word_t key) {
	return key_type((size_t(key) >> half_shift) & mask,  size_t(key) & mask);
      }

    public:
      class const_iterator : public std::iterator<std::bidirectional_iterator_tag, value_type> {
      public:
	const_iterator(const JudyPairMap& container, size_t idx = 0):_container(container), _idx(idx) {
	  Word_t* ptr = (Word_t*)JudyLFirst(_container._array, &_idx, NULL);
	  if (ptr == NULL) 
	    {_idx = endIndex; _value = 0;}
	  else
	    _value = *ptr;
	}

	const_iterator& operator++() { 
	  Word_t* ptr = (Word_t*)JudyLNext(_container._array, &_idx, NULL);
	  if (ptr == NULL) 
	    {_idx = endIndex; _value = 0;}
	  else
	    _value = *ptr;
	  return *this;
	}

	const_iterator& operator--() { 
	  Word_t* ptr = (Word_t*)JudyLPrev(_container._array, &_idx, NULL);
	  if (ptr == NULL) 
	    {_idx = endIndex; _value = 0;}
	  else
	    _value = *ptr;
	  return *this;
	}

	value_type operator*() const { return value_type(IDtoKey(_idx), _value); }
	bool operator!=(const const_iterator& o) const { return !(*this == o); }
	bool operator==(const const_iterator& o) const { return _idx == o._idx;}

      private:
	const JudyPairMap& _container;
	Word_t _idx;
	Word_t _value;
      };

      const_iterator begin() const { return const_iterator(*this); }
      const_iterator end() const { return const_iterator(*this, endIndex); }

      ~JudyPairMap() { clear(); }

      void clear() { JudyLFreeArray(&_array, NULL); _count = 0; }

      /*!\brief This proxy is used to double check if an assignment of
	zero is done, and delete the entry if it is. */
      struct EntryProxy {
      public:
	EntryProxy(JudyPairMap& container, const key_type& key):
	  _container(container), _ID(keyToID(key)) {}

	operator const mapped_type() const {
	  mapped_type* PValue = NULL;
	  PValue = (mapped_type*)JudyLGet(_container._array, _ID, NULL);
	  return (PValue == NULL) ? 0 : *PValue;
	}
	
	EntryProxy& operator=(mapped_type newval) {
	  if (newval == 0)
	    _container.erase(_ID);
	  else
	    {
	      mapped_type& Value = *(mapped_type*)JudyLIns(&(_container._array), _ID, NULL);
	      if (Value == 0) ++_container._count;
	      Value = newval;
	    }
	  return *this;
	}
	
      private:
	JudyPairMap& _container;
	const Word_t _ID;
      };

      size_t size() const { return _count; }
      
      value_type operator[](const size_t i) const {
	Word_t ID;
	mapped_type val = *(mapped_type*)JudyLByCount(_array, i+1, &ID, NULL);
	return value_type(IDtoKey(ID), val);
      }
      
      /*! \brief This non-const array access operator uses EntryProxy
	to check if any values assigned are zero so they may be deleted. */
      EntryProxy operator[](const key_type& key) {
	return EntryProxy(*this, key); 
      }

      mapped_type operator[](const key_type& key) const {
	mapped_type* PValue = NULL;
	PValue = (mapped_type*)JudyLGet(_array, keyToID(key), NULL);
	return (PValue == NULL) ? 0 : *PValue;
      }

      inline void insert(const value_type value) {
	mapped_type& PValue = *(mapped_type*)JudyLIns(&_array, keyToID(value.first), NULL);
	if (PValue == 0) ++_count;
	PValue = value.second;
      }

      inline void erase(const key_type key) {
	erase(keyToID(key));
      }

    protected:
      inline void erase(const Word_t key) {
	_count -= JudyLDel(&_array, key, NULL);
      }
    };
  }
}
