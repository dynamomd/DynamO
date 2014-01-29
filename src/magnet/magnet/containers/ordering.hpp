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
#include <magnet/math/dilated_int.hpp>
#include <array>

namespace magnet {
  namespace containers {
    namespace detail {
      /*! \brief Base class for any derived Ordering classes.

	Ordering classes help with storing multidimensional arrays in
	memory (which is linear).
	
	\tparam NDim How many dimensions the array has.

	\tparam Derived This class uses the Curiously Recurring Template Pattern
	(CRTP) to determine its Derived class type, so that it can
	access the specialised functions it needs.
       */
      template<size_t NDim, typename Derived>
      class OrderingBase {
      public:
	typedef std::array<size_t, NDim> ArrayType;

	class const_iterator : public std::iterator<std::forward_iterator_tag, size_t>
	{
	public:
	  const_iterator(const Derived& container, const ArrayType& start, const ArrayType& distance):
	    _container(container), _start(start), _distance(distance), _pos({{0,0,0}}) {}
	  
	  const_iterator(const Derived& container, const ArrayType& start, const ArrayType& distance, const ArrayType& pos):
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
	  const Derived& _container;
	  const ArrayType _start;
	  const ArrayType _distance;
	  ArrayType _pos;
	  size_t _index;
	};
	
	OrderingBase(const ArrayType& dimensions): _dimensions(dimensions), _end() {
	  _end[NDim - 1] = _dimensions[NDim - 1];
	}

	template<typename... Args>
	OrderingBase(Args... args): _dimensions(std::forward<Args>(args)...), _end() {
	  _end[NDim - 1] += _dimensions[NDim - 1] + 1;
	}

	/*! \brief How many unique coordinates exist in the system. */
	size_t size() const {
	  size_t val = 1;
	  for (size_t i(0); i < NDim; ++i)
	    val *= _dimensions[i];
	  return val;
	}

	const ArrayType& getDimensions() const { return _dimensions; }

	const_iterator begin() const { 
	  return const_iterator(static_cast<const Derived&>(*this), {{0,0,0}}, _dimensions); 
	}
	
	const_iterator end() const { 
	  return const_iterator(static_cast<const Derived&>(*this), {{0,0,0}}, _dimensions, _end); 
	}
	
	typedef magnet::containers::IteratorPairRange<const_iterator> IndexRange;
	IndexRange getIndices(const ArrayType& start, const ArrayType& distance) const {
	  ArrayType end = ArrayType();
	  end[NDim - 1] = distance[NDim - 1];
	  return IndexRange(const_iterator(static_cast<const Derived&>(*this), start, distance), 
			    const_iterator(static_cast<const Derived&>(*this), start, distance, end));
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

      protected:
	ArrayType _dimensions;
	ArrayType _end;
      };
    }

    /*! \brief Standard Row-Major ordering of elements in memory.
      
      \tparam NDim The dimensionality of the array.
    */
    template <size_t NDim>
    class RowMajorOrdering : public detail::OrderingBase<NDim, RowMajorOrdering<NDim> > {
      typedef typename detail::OrderingBase<NDim, RowMajorOrdering<NDim> > Base;
    public:
      using typename Base::ArrayType;

      RowMajorOrdering(const ArrayType& dimensions): Base(dimensions) { }

      template<typename... Args> RowMajorOrdering(Args... args): Base(std::forward<Args>(args)...) {}

      size_t toIndex(const ArrayType& loc) const  {
	size_t index = 0;
	for (size_t i(0); i < NDim; ++i)
	  {
	    index = index * Base::_dimensions[NDim - 1 - i] + (loc[NDim - i - 1] % Base::_dimensions[NDim - 1 - i]);
	  }
	return index;
      }

      ArrayType toCoord(size_t index) const {
	ArrayType coord;
	for (size_t i(0); i < NDim; ++i)
	  {
	    coord[i] = index % Base::_dimensions[i];
	    index /= Base::_dimensions[i];
	  }
	return coord;
      }
  
      /*! \brief How many elements are needed to store the array. */
      size_t length() const { return Base::size(); }
    };

    /*! \brief Morton ordering of elements in memory.
      
      Morton ordering is a type of ordering that ensures that elements
      close to each other in space are more likely to be closer
      together in (linear) memory spaces.

      \tparam NDim The dimensionality of the array.
    */
    template <size_t NDim>
    class MortonOrdering : public detail::OrderingBase<NDim, MortonOrdering<NDim> > {
      typedef typename detail::OrderingBase<NDim, MortonOrdering<NDim> > Base;
    public:
      using typename Base::ArrayType;

      MortonOrdering(const ArrayType& dimensions): Base(dimensions) { }

      template<typename... Args> MortonOrdering(Args... args): Base(std::forward<Args>(args)...) {}
      
      size_t toIndex(const ArrayType& loc) const  {
	size_t index = 0;
	for (size_t i(0); i < NDim; ++i)
	  {
	    index += magnet::math::DilatedInteger<NDim>(loc[i] % Base::_dimensions[i]).getDilatedValue() << i;
	  }
	return index;
      }

      ArrayType toCoord(const size_t index) const {
	ArrayType coord;
	for (size_t i(0); i < NDim; ++i)
	  {
	    magnet::math::DilatedInteger<NDim> dilatedint;
	    dilatedint.setDilatedValue(index >> i);
	    coord[i] = dilatedint.getRealValue();
	  }
	return coord;
      }

      /*! \brief How many elements are needed to store the array. */
      size_t length() const {
	size_t length = 0;
	for (size_t i(0); i < NDim; ++i)
	  length += magnet::math::DilatedInteger<NDim>(Base::_dimensions[i]).getDilatedValue() << i;
	return length;
      }
    };
  }
}
