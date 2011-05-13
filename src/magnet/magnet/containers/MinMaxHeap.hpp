/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#include <magnet/exception.hpp>
#include <boost/array.hpp>

namespace magnet {
  namespace containers {    
    template <typename T, std::size_t N> class MinMaxHeap;
    
    //The unsupported sizes of the MinMax heap
    template <typename T> class MinMaxHeap<T,0> {};
    template <typename T> class MinMaxHeap<T,1> {};

    //! \brief A fixed size MinMax heap implementation. 
    //!
    //! This type of
    //! container is optimized for access to the minimum and maximum
    //! elements of the heap.
    //! \tparam T Stored type in the heap.
    //! \tparam N Size of the heap.
    template <typename T, std::size_t N>
    class MinMaxHeap
    {
      typedef typename boost::array<T, N> Container;
      typedef typename Container::value_type Comparable;
      //! The heap array
      Container _array; 
      //! Number of elements currently in heap
      size_t _currentSize;

    public:
      MinMaxHeap(): _currentSize(0) {}

      typedef typename Container::iterator iterator;
      typedef typename Container::const_iterator const_iterator;

      inline iterator begin() { return _array.begin(); }
      inline const_iterator begin() const { return _array.begin(); }
      inline iterator end() { return begin() + _currentSize; }
      inline const_iterator end() const { return begin() + _currentSize; }

      inline void pop()
      {
#ifdef MAGNET_DEBUG
	if( empty() )
	  M_throw() << "*** DeleteMin failed: Heap is empty ***";
#endif
    
	_array[0] = _array[--_currentSize];
	percolateDown(1);
      }

      // FindMin
      // Finds and returns the minimum item in the heap
      // Returns: minimum item in the heap
      inline const Comparable & top() const
      {
#ifdef MAGNET_DEBUG
	if (empty())
	  M_throw() << "*** FindMin failed: Heap is empty ***";
#endif
	return _array[0];
      }


      // FindMax
      // Finds and returns the maximum item in the heap
      // Returns: maximum item in the heap
      inline const Comparable & bottom() const
      {
#ifdef MAGNET_DEBUG
	if ( empty() )
	  M_throw() << "*** FindMax failed: Heap is empty ***";
#endif
	if (_currentSize <= 1)
	  return _array[0];
	else if (_currentSize == 2)
	  return _array[1];
	else if (N > 2) //This is just to stop a compiler warning for N=2 containers
	  return _array[1] > _array[2] ? _array[1] : _array[2];
	else 
	  M_throw() << "Bad _currentSize when accessing heap bottom";
      }


      inline Comparable & unsafe_bottom()
      {
#ifdef MAGNET_DEBUG
	if ( empty() )
	  M_throw() << "*** FindMax failed: Heap is empty ***";
#endif
	if (_currentSize == 1)
	  return _array[0];
	else if ( _currentSize == 2 )
	  return _array[1];
	else if (N > 2) //This is just to stop a compiler warning for N=2 containers
	  return _array[1] > _array[2] ? _array[1] : _array[2];
	else 
	  M_throw() << "Bad _currentSize when accessing heap bottom in N=2 MinMaxHeap";
      }

      // Insert
      // Inserts item into the min-max heap, maintaining heap order
      //   Duplicates are allowed.
      // Parameter x: item to be inserted
      // Side Effects: throws Overflow if heap is full
      inline void insert( const Comparable & x )
      {
#ifdef MAGNET_DEBUG
	if (full())
	  M_throw() << "*** Insert failed: Heap is full ***";
#endif
	_array[_currentSize] = x;    
	percolateUp(++_currentSize);
      }


      // DeleteMin
      // Remove the smallest item from the min-max heap and
      //   place it in minItem.
      // Parameter minItem: reference from calling function
      //   to put the smallest item in.
      // Side Effects: throws Underflow if heap is empty
      inline void deleteMin( Comparable & minItem )
      {
#ifdef MAGNET_DEBUG
	if( empty() )
	  M_throw() << "*** DeleteMin failed: Heap is empty ***";
#endif
    
	minItem = _array[ 1 ];
	pop();
      }

      // DeleteMax
      // Remove the largest item from the min-max heap and
      //   place it in maxItem.
      // Parameter maxItem: reference from calling function
      //   to put the largest item in.
      // Side Effects: throws Underflow if heap is empty
      inline void deleteMax( Comparable & maxItem)
      {
	size_t maxIndex;
    
#ifdef MAGNET_DEBUG
	if ( empty() )
	  M_throw() << "*** DeleteMax failed: Heap is empty ***";
#endif
    
	if (_currentSize ==1)
	  maxIndex = 0;
	else if ( _currentSize == 2 )
	  maxIndex = 1;
	else
	  maxIndex = _array[1] > _array[2] ? 1 : 2;
    
	maxItem = _array[maxIndex];
	_array[maxIndex] = _array[--_currentSize];
	percolateDown(maxIndex+1);   
      }


      inline void replaceMax(const Comparable & newMaxItem)
      {
	size_t maxIndex;

#ifdef MAGNET_DEBUG
	if ( empty() )
	  M_throw() << "*** DeleteMax failed: Heap is empty ***";
#endif

	if (_currentSize == 1)
	  maxIndex = 0;
	else if (_currentSize == 2)
	  maxIndex = 1;
	else
	  maxIndex = _array[1] > _array[2] ? 1 : 2;

  
	_array[ maxIndex ] = _array[--_currentSize];
	percolateDown(maxIndex + 1);
	insert(newMaxItem);
      }

      inline void clear() { _currentSize = 0; }

      inline size_t size() const { return _currentSize; }

      // IsEmpty
      // Checks to see if heap is logically empty
      // Returns: true if empty, false if not
      inline bool empty() const { return _currentSize == 0; }

      // IsFull
      // Checks to see if heap is logically full
      // Returns: true if full, false if not
      inline bool full() const { return _currentSize == _array.size(); }

      inline void swap(MinMaxHeap<T, N>& rhs)
      {
	std::swap(_currentSize, rhs._currentSize);
	_array.swap(rhs._array);
      }

    private:

      // PercolateUp
      // Used to maintain min-max heap order after insertion.
      //   Determines whether current heap level is a min level
      //   or a max level and calls percolateUpMin or percolateUpMax.
      // Parameter hole: index in array where item was inserted
      inline void percolateUp( size_t hole )
      {
	size_t parent = hole / 2;
	size_t level = 0;
    
	for ( size_t i = 2; i <= hole; i *= 2 )
	  level++;
    
	// min level
	if ( level % 2 == 0 ) {
      
	  if ( parent > 0 && _array[hole - 1] > _array[parent-1] ) {
	    swapElements(hole, parent);
	    percolateUpMax(parent);
	  }
	  else
	    percolateUpMin(hole);
	}
	// max level
	else {
      
	  if ( parent > 0 && _array[hole-1] < _array[parent-1] ) {
	    swapElements( hole, parent );
	    percolateUpMin(parent);
	  }
	  else
	    percolateUpMax(hole);
	}
      }

      // PercolateUpMin
      // Called by percolateUp to maintain min-max heap order on 
      //   the min levels.
      // Parameter hole: index in array of item to be percolated
      inline void percolateUpMin( size_t hole )
      {
	size_t grandparent = hole / 4;

	if ( grandparent > 0 && _array[hole-1] < _array[grandparent-1] ) {
	  swapElements(hole, grandparent);
	  percolateUpMin(grandparent);
	}
      }

      // PercolateUpMax
      // Called by percolateUp to maintain min-max heap order on
      //   the max levels.
      // Parameter hole: index in array of item to be percolated
      inline void percolateUpMax(size_t hole)
      {
	size_t grandparent = hole / 4;

	if ( grandparent > 0 && _array[hole -1] > _array[grandparent-1] ) {
	  swapElements( hole, grandparent );
	  percolateUpMax( grandparent );
	}
      }

      // PercolateDown
      // Used to maintain min-max heap order after deletion.
      //   Determines whether current heap level is a min level
      //   or a max level and calls percolateDownMin or percolateDownMax.
      // Parameter hole: index in array where item was deleted
      inline void percolateDown( size_t hole )
      {
	size_t level = 0;
    
	for (size_t i = 2; i <= hole; i *= 2, ++level);
    
	if ( level % 2 == 0 )
	  percolateDownMin(hole);
	else
	  percolateDownMax(hole);
      }

      // PercolateDownMin
      // Called by percolateDown to maintain min-max heap order on
      //   the min levels.
      // Parameter hole: index in array of item to be percolated
      inline void percolateDownMin( size_t hole )
      {
	// find minimum value of children and grandchildren
	// hole * 2 = index of first child if it exists
	// hole * 4 = index of first grandchild if it exists
	size_t minIndex = findMinDescendent( hole * 2, hole * 4 );

	// at least one descendent
	if ( minIndex > 0 ) {

	  // min descendent is a grandchild
	  if ( minIndex >= hole * 4 ) {

	    // if less than grandparent, i.e value at hole, swap
	    if ( _array[minIndex-1] < _array[hole-1] ) {
	      swapElements( hole, minIndex );

	      // if greater than parent, swap
	      if ( _array[minIndex-1] > _array[minIndex / 2 -1] )
		swapElements( minIndex, minIndex / 2 );

	      percolateDownMin( minIndex );
	    }
	  }
	  // min descendent is a child
	  else {

	    // if less than parent, i.e value at hole, swap
	    if ( _array[minIndex-1] < _array[hole-1] )
	      swapElements(hole, minIndex);
	  }
	}
      }


      // PercolateDownMax
      // Called by percolateDown to maintain min-max heap order on
      //   the max levels.
      // Parameter hole: index in array of item to be percolated
      inline void percolateDownMax(size_t hole)
      {
	size_t maxIndex;

	// find maximum value of children and grandchildren
	// hole * 2 = index of first child if it exists
	// hole * 4 = index of first grandchild if it exists
	maxIndex = findMaxDescendent( hole * 2, hole * 4 );

	// at least one descendent
	if ( maxIndex > 0 ) {

	  // max descendent is a grandchild
	  if ( maxIndex >= hole * 4 ) {

	    // if greater than grandparent, i.e value at hole, swap
	    if ( _array[maxIndex-1] > _array[hole-1] ) {
	      swapElements( hole, maxIndex );

	      // if less than parent, swap
	      if ( _array[maxIndex-1] < _array[maxIndex / 2 - 1] )
		swapElements( maxIndex, maxIndex / 2 );

	      percolateDownMax( maxIndex );
	    }
	  }
	  // max descendent is a child
	  else {

	    // if greater than parent, i.e value at hole, swap
	    if ( _array[maxIndex-1] > _array[hole-1] )
	      swapElements( hole, maxIndex );
	  }
	}
      }
 

      // FindMinDescendent
      // Helper function of percolateDownMin that finds the min
      //   of the children and grandchildren of the item to be
      //   percolated down
      // Parameter child: array index of first child
      // Parameter grandchild: array index of first grandchild
      inline size_t findMinDescendent( size_t child, size_t grandchild )
      {
	size_t minIndex = 0;
	size_t minChild = child;
	size_t minGrandchild = grandchild;

	if ( child <= _currentSize ) {

	  if ( child != _currentSize && _array[child] < _array[minChild-1] )
	    minChild = child + 1;

	  minIndex = minChild;

	  if (grandchild <= _currentSize) {

	    for ( size_t i = 1; grandchild < _currentSize && i < 4; i++, grandchild++ ) {

	      if ( _array[grandchild] < _array[minGrandchild-1] )
		minGrandchild = grandchild + 1;
	    }

	    if ( _array[minGrandchild-1] < _array[minChild-1] )
	      minIndex = minGrandchild;
	  }
	}

	return minIndex;
      }

      // FindMaxDescendent
      // Helper function of percolateDownMax that finds the max
      //   of the children and grandchildren of the item to be
      //   percolated down
      // Parameter child: array index of first child
      // Parameter grandchild: array index of first grandchild
      inline size_t findMaxDescendent( size_t child, size_t grandchild )
      {
	size_t maxIndex = 0;
	size_t maxChild = child;
	size_t maxGrandchild = grandchild;

	if ( child <= _currentSize ) {

	  if ( child != _currentSize && _array[child] > _array[maxChild-1] )
	    maxChild = child + 1;

	  maxIndex = maxChild;

	  if ( grandchild <= _currentSize ) {

	    for ( size_t i = 1; grandchild < _currentSize && i < 4; i++, grandchild++ ) {

	      if ( _array[grandchild] > _array[maxGrandchild-1] )
		maxGrandchild = grandchild + 1;
	    }

	    if ( _array[maxGrandchild-1] > _array[maxChild-1] )
	      maxIndex = maxGrandchild;
	  } 
	}

	return maxIndex;
      }

      // Swap
      // Swaps elements of two indices
      // Parameter indexOne: first index of array of item to be swapped
      // Parameter indexTwo: second index of array of item to be swapped
      inline void swapElements(size_t indexOne, size_t indexTwo)
      {std::swap(_array[indexOne-1], _array[indexTwo-1]);}
    };
  }
}

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<typename T, std::size_t N> inline
  void swap(magnet::containers::MinMaxHeap<T, N>& lhs, 
	    magnet::containers::MinMaxHeap<T, N>& rhs)
  { lhs.swap(rhs); }
}
