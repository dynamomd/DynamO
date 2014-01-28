/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <magnet/containers/iterator_pair.hpp>
#include <magnet/math/morton_number.hpp>
#include <array>

namespace magnet {
  namespace containers {
    template <size_t NDim>
    class RowMajorOrdering {
    public:
      typedef std::array<size_t, NDim> ArrayType;

      size_t toIndex(const ArrayType& loc) const  {
	size_t index = loc[NDim - 1];
	for (size_t i(1); i < NDim; ++i)
	  index = index * _dimensions[NDim - i - 1] + (loc[NDim - i - 1] % _dimensions[NDim - i - 1]);
	return index;
      }

      ArrayType toCoord(size_t index) const {
	ArrayType coord;
	for (size_t i(0); i < NDim; ++i)
	  {
	    coord[i] = index % _dimensions[i];
	    index /= _dimensions[i];
	  }
	return coord;
      }

      class const_iterator : public std::iterator<std::forward_iterator_tag, size_t>
      {
      public:
	const_iterator(const RowMajorOrdering& container, const ArrayType& start, const ArrayType& distance):
	  _container(container), _start(start), _distance(distance), _pos({{0,0,0}}) {}

	const_iterator(const RowMajorOrdering& container, const ArrayType& start, const ArrayType& distance, const ArrayType& pos):
	  _container(container), _start(start), _distance(distance), _pos(pos) {}

	const_iterator& operator++() {
	  ++_pos[0];
	  for (size_t i(1); (i < NDim) && (_pos[i-1] == _distance[i-1]); ++i)
	    {
	      //Increment the next coordinate
	      ++_pos[i];
	      //reset the last coordinate
	      _pos[i - 1] = 0;
	    }
	  return *this;
	}

	value_type operator*() const { 
	  ArrayType loc;
	  for (size_t i(0); i < NDim; ++i)
	    loc[i] = _start[i] + _pos[i];
	  return _container.toIndex(loc); 
	}
	bool operator!=(const const_iterator& o) const { return !(*this == o); }
	bool operator==(const const_iterator& o) const {
#ifdef MAGNET_DEBUG
	  if (&_container != &o._container) M_throw() << "Cannot compare iterators from different containers";
#endif

	  return _pos == o._pos;
	}
	//protected:
	const RowMajorOrdering& _container;
	const ArrayType _start;
	const ArrayType _distance;
	ArrayType _pos;
	size_t _index;
      };

      RowMajorOrdering(const ArrayType& dimensions): _dimensions(dimensions), _end()
      {
	_end[NDim - 1] = _dimensions[NDim - 1];
      }

      template<typename... Args>
      RowMajorOrdering(Args... args): _dimensions(std::forward<Args>(args)...), _end()
      {
	_end[NDim - 1] += _dimensions[NDim - 1] + 1;
      }
  
      size_t length() const {
	size_t length = 1;
	for (size_t i(0); i < NDim; ++i)    
	  length *= _dimensions[i];
	return length;
      }
  
      const_iterator begin() const { 
	return const_iterator(*this, {{0,0,0}}, _dimensions); 
      }

      const_iterator end() const { 
	return const_iterator(*this, {{0,0,0}}, _dimensions, _end); 
      }
  
      typedef magnet::containers::IteratorPairRange<const_iterator> IndexRange;
      IndexRange getIndices(const ArrayType& start, const ArrayType& distance) const {
	ArrayType end = ArrayType();
	end[NDim - 1] = distance[NDim - 1];
	return IndexRange(const_iterator(*this, start, distance), const_iterator(*this, start, distance, end));
      }

      const ArrayType& getDimensions() const { return _dimensions; }
      
    private:
      const ArrayType _dimensions;
      ArrayType _end;
    };

    template <size_t NDim>
    class MortonOrdering {
    public:
      typedef std::array<size_t, NDim> ArrayType;
      
      size_t toIndex(const ArrayType& loc) const  {
	size_t index = 0;
	for (size_t i(0); i < NDim; ++i)
	  index += magnet::math::DilatedInteger<NDim>(loc[i] % _dimensions[i]).getDilatedValue() << i;
	return index;
      }

      ArrayType toCoord(const magnet::math::MortonNumber<NDim>& index) const {
	ArrayType coord;
	for (size_t i(0); i < NDim; ++i)
	  coord[i] = index[i].getRealValue();
	return coord;
      }

      ArrayType newLocation(ArrayType loc, size_t dim, int direction) const  {
	loc[dim] += _dimensions[dim] + direction;
	return index;
      }

      class const_iterator : public std::iterator<std::forward_iterator_tag, size_t>
      {
      public:
	const_iterator(const MortonOrdering& container, const ArrayType& start, const ArrayType& distance, const ArrayType& pos = ArrayType()):
	  _container(container), _start(start), _distance(distance), _pos(pos) {}

	const_iterator& operator++() {
	  ++_pos[0];
	  for (size_t i(1); (i < NDim) && (_pos[i-1] == _distance[i-1]); ++i)
	    {
	      //Increment the next coordinate
	      ++_pos[i];
	      //reset the last coordinate
	      _pos[i - 1] = 0;
	    }
	  return *this;
	}

	value_type operator*() const { 
	  ArrayType loc;
	  for (size_t i(0); i < NDim; ++i)
	    loc[i] = _start[i] + _pos[i];
	  return _container.toIndex(loc); 
	}
	bool operator!=(const const_iterator& o) const { return !(*this == o); }
	bool operator==(const const_iterator& o) const {
#ifdef MAGNET_DEBUG
	  if (&_container != &o._container) M_throw() << "Cannot compare iterators from different containers";
#endif

	  return _pos == o._pos;
	}
	//protected:
	const MortonOrdering& _container;
	const ArrayType _start;
	const ArrayType _distance;
	ArrayType _pos;
	size_t _index;
      };

      MortonOrdering(const ArrayType& dimensions): _dimensions(dimensions), _end()
      {
	_end[NDim - 1] = _dimensions[NDim - 1];
      }

      template<typename... Args>
      MortonOrdering(Args... args): _dimensions(std::forward<Args>(args)...), _end()
      {
	_end[NDim - 1] += _dimensions[NDim - 1] + 1;
      }
  
      size_t length() const {
	size_t length = 0;
	for (size_t i(0); i < NDim; ++i)
	  length += magnet::math::DilatedInteger<NDim>(_dimensions[i]).getDilatedValue() << i;
	return length;
      }
  
      const_iterator begin() const { 
	return const_iterator(*this, ArrayType(), _dimensions);
      }

      const_iterator end() const { 
	return const_iterator(*this, ArrayType(), _dimensions, _end);
      }
  
      typedef magnet::containers::IteratorPairRange<const_iterator> IndexRange;

      IndexRange getIndices(const ArrayType& start, const ArrayType& distance) const {
	ArrayType end = ArrayType();
	end[NDim - 1] = distance[NDim - 1];
	return IndexRange(const_iterator(*this, start, distance), const_iterator(*this, start, distance, end));
      }

      IndexRange getSurroundingIndices(const ArrayType& center, const ArrayType& distance) const {
	ArrayType zero;
	ArrayType range;
	for (size_t i(0); i < NDim; ++i)
	  {
	    zero[i] = (center[i] + _dimensions[i] - distance[i]) % _dimensions[i];
	    range[i] = 2 * distance[i] + 1;
	  }
	
	return getIndices(zero, range);
      }

      const ArrayType& getDimensions() const { return _dimensions; }
      
    private:
      ArrayType _dimensions;
      ArrayType _end;
    };
  }
}
