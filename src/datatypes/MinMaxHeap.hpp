//------------------------------------------------------------------------
// File:     MinMaxHeap.H
// Project:  CMSC 341  Project 5  Fall 2003
// Author:   Todd Wease, Section 0301
// Created:  11/24/03
// Current:  12/05/03
// 
// The min-max heap is an extension of the binary heap that has both a
// min order and max order thus supporting both deleteMin and deleteMax.
// Min order of the heap is maintained on even levels of the heap, the
// root being at level 0, and max order is maintained on odd levels of the
// heap.  Both orders are intermeshed such that no child of a min level
// is ever less than the parent and no child of a max level is ever greater
// than the parent.
// 
//-------------------------------------------------------------------------
// I have read and I understand the course policy on cheating. By submitting
// the following program, I am stating that the program was produced by my
// individual effort. 
//-------------------------------------------------------------------------

#ifndef MINMAXHEAP_H
#define MINMAXHEAP_H

#include <iostream>
#include <vector>
#include <boost/array.hpp>

using namespace std;

template <class Comparable>
class MinMaxHeap
{
  size_t _currentSize;            // Number of elements currently in heap
  std::vector<Comparable> _array; // The heap array

public:
  MinMaxHeap(const size_t heapSize):_array(heapSize) {}

  typedef typename std::vector<Comparable>::iterator iterator;
  typedef typename std::vector<Comparable>::const_iterator const_iterator;

  inline iterator begin() { return _array.begin() + 1; }
  inline const_iterator begin() const { return _array.begin() + 1; }
  inline iterator end() { return begin() + _currentSize; }
  inline const_iterator end() const { return begin() + _currentSize; }

  inline void pop();

  // FindMin
  // Finds and returns the minimum item in the heap
  // Returns: minimum item in the heap
  inline const Comparable & top() const;

  // FindMax
  // Finds and returns the maximum item in the heap
  // Returns: maximum item in the heap
  const Comparable & bottom() const;

  Comparable & unsafe_bottom();

  // Insert
  // Inserts item into the min-max heap, maintaining heap order
  //   Duplicates are allowed.
  // Parameter x: item to be inserted
  // Side Effects: throws Overflow if heap is full
  void insert( const Comparable & x );

  // DeleteMin
  // Remove the smallest item from the min-max heap and
  //   place it in minItem.
  // Parameter minItem: reference from calling function
  //   to put the smallest item in.
  // Side Effects: throws Underflow if heap is empty
  void deleteMin( Comparable & minItem );

  // DeleteMax
  // Remove the largest item from the min-max heap and
  //   place it in maxItem.
  // Parameter maxItem: reference from calling function
  //   to put the largest item in.
  // Side Effects: throws Underflow if heap is empty
  void deleteMax( Comparable & maxItem);

  void replaceMax(const Comparable & maxItem);

  inline void clear() { _currentSize = 0; }

  inline size_t size() const { return _currentSize; }

  // IsEmpty
  // Checks to see if heap is logically empty
  // Returns: true if empty, false if not
  inline bool empty() const { return _currentSize == 0; }

  // IsFull
  // Checks to see if heap is logically full
  // Returns: true if full, false if not
  inline bool full() const { return _currentSize == _array.size() - 1; }

  inline void swap(MinMaxHeap<Comparable>& rhs)
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
  void percolateUp( size_t hole );

  // PercolateUpMin
  // Called by percolateUp to maintain min-max heap order on 
  //   the min levels.
  // Parameter hole: index in array of item to be percolated
  void percolateUpMin( size_t hole );

  // PercolateUpMax
  // Called by percolateUp to maintain min-max heap order on
  //   the max levels.
  // Parameter hole: index in array of item to be percolated
  void percolateUpMax( size_t hole );

  // PercolateDown
  // Used to maintain min-max heap order after deletion.
  //   Determines whether current heap level is a min level
  //   or a max level and calls percolateDownMin or percolateDownMax.
  // Parameter hole: index in array where item was deleted
  void percolateDown( size_t hole );

  // PercolateDownMin
  // Called by percolateDown to maintain min-max heap order on
  //   the min levels.
  // Parameter hole: index in array of item to be percolated
  void percolateDownMin( size_t hole );

  // PercolateDownMax
  // Called by percolateDown to maintain min-max heap order on
  //   the max levels.
  // Parameter hole: index in array of item to be percolated
  void percolateDownMax( size_t hole );

  // FindMinDescendent
  // Helper function of percolateDownMin that finds the min
  //   of the children and grandchildren of the item to be
  //   percolated down
  // Parameter child: array index of first child
  // Parameter grandchild: array index of first grandchild
  size_t findMinDescendent( size_t child, size_t grandchild );

  // FindMaxDescendent
  // Helper function of percolateDownMax that finds the max
  //   of the children and grandchildren of the item to be
  //   percolated down
  // Parameter child: array index of first child
  // Parameter grandchild: array index of first grandchild
  size_t findMaxDescendent( size_t child, size_t grandchild );

  // Swap
  // Swaps elements of two indices
  // Parameter indexOne: first index of array of item to be swapped
  // Parameter indexTwo: second index of array of item to be swapped
  inline void swapElements( size_t indexOne, size_t indexTwo )
  {std::swap(_array[ indexOne ], _array[ indexTwo]);}
};

namespace std
{
  /*! \brief Template specialisation of the std::swap function for pList*/
  template<typename T> inline
  void swap(MinMaxHeap<T>& lhs, MinMaxHeap<T>& rhs)
  {
    lhs.swap(rhs);
  }
}


#include "MinMaxHeap.cpp"

#endif

