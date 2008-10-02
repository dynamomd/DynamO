/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSSHeap
#define CSSHeap
#include <vector>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>
#include "../../base/is_exception.hpp"

template<typename dataType, typename idType>
struct Node_type
{  
  Node_type(const idType npID, const dataType& nCC):
    pID(npID), ptrData(&nCC) {}
  
  inline bool operator>(const Node_type& node2) const
  {
    return *ptrData > *node2.ptrData;
  }

  inline bool operator<(const Node_type& node2) const
  {
    return *ptrData < *node2.ptrData;
  }

  idType pID;
  const dataType* ptrData;
};

template<typename T, typename idType = unsigned long,
	 class Compare = std::greater<Node_type<T, idType> > >
class CSHeap 
{ 
  typedef typename std::vector<Node_type<T, idType> >::iterator CRanIt;
  typedef typename std::iterator_traits<CRanIt>::difference_type CDiff;
  typedef Node_type<T, idType> Node;

  Compare comp;
  std::vector<Node> local_heap;
  std::vector<T>    data_stack;
  std::vector<CDiff> idMap;
  
public:
  CSHeap() {}
  ~CSHeap() {}

  void resize(size_t size) 
  { 
    clear();
    data_stack.resize(size);
    idMap.resize(size);
  }

  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  inline iterator begin() { return data_stack.begin(); }
  inline const_iterator begin() const { return data_stack.begin(); }
  inline iterator end() { return data_stack.end(); }
  inline const_iterator end() const { return data_stack.end(); }
  inline size_t size() const { return data_stack.size(); }
  inline bool empty() const { return data_stack.empty(); }

  inline const T& operator[](idType a) const 
  {
#ifdef DYNAMO_DEBUG 
    if (data_stack.empty())
      I_throw() << "Heap not yet sized";
#endif
    
    return data_stack[a]; 
  }
  
  inline T& operator[](idType a) 
  {
#ifdef DYNAMO_DEBUG 
    if (data_stack.empty())
      I_throw() << "Heap not yet sized";
#endif

    return data_stack[a]; 
  }

  void init()
  {
#ifdef DYNAMO_DEBUG 
    if (data_stack.empty())
      I_throw() << "Heap not yet sized";
#endif

    idType ID(0);
    BOOST_FOREACH(T& val, data_stack)
      local_heap.push_back(Node(ID++, val));
    
    make_heap(local_heap.begin(), local_heap.end(),comp);
    
    //Iterate over the heap looking at ID numbers and placing the
    //correct map location
    for (CRanIt iPtr = local_heap.begin();
	 iPtr != local_heap.end(); iPtr++)
      idMap[iPtr->pID] = iPtr - local_heap.begin();
  }

  void update(idType ID)
  {
    CRanIt first = local_heap.begin();

#ifdef DYNAMO_DEBUG 
    if (data_stack.empty())
      I_throw() << "Heap not yet sized";

    if (first == local_heap.end())
      I_throw() << "Update called on an empty heap";
 
    if (idMap.empty())
      I_throw() << "Heap not yet initialised";
#endif
    
    unsigned long index = idMap[ID];
    CRanIt pos = first + index; 
    
    if(index > 0 && comp(*(first + ((index - 1)/2)), *pos))
      up_heap(pos);
    else 
      down_heap(pos);
  }

  inline const T& next_Data() const
  {
#ifdef DYNAMO_DEBUG 
    if (idMap.empty())
      I_throw() << "Heap not yet initialised";
#endif
    return *local_heap.front().ptrData;
  }

  inline T& next_Data()
  {
#ifdef DYNAMO_DEBUG 
    if (idMap.empty())
      I_throw() << "Heap not yet initialised";
#endif
    return data_stack[local_heap.front().pID];
  }

  inline const idType next_ID() const
  {
#ifdef DYNAMO_DEBUG 
    if (idMap.empty())
      I_throw() << "Heap not yet initialised";
#endif

    return local_heap.front().pID;
  }

  void clear()
  {
    idMap.clear();
    local_heap.clear();
    data_stack.clear();
  }

protected:

  inline void up_heap(CRanIt pos)
  {    
    CRanIt first = local_heap.begin();

    CDiff parent = (pos - first - 1) / 2;
    CDiff index = pos - first;

    Node mov = *pos;
    
    while(index > 0 && comp(*(first + parent),mov) )
      {
	*(first + index) = *(first + parent);
	//Performed a move, update the ID map
	idMap[(first + index)->pID] = index;

	index = parent;
	parent = (parent - 1) / 2;
      }
    
    if(index != pos - first)
      {
	*(first + index) = mov;
	//Performed a move, update the ID map
	idMap[(first + index)->pID] = index;
      }
  }
  

  inline void down_heap(CRanIt pos)
  {     
    CRanIt first = local_heap.begin();
    CRanIt last = local_heap.end();

    Node mov = *pos;
    CDiff index = pos - first;
    CDiff left = 2 * index + 1;
    CDiff right = 2 * index + 2;
    CDiff len = last - first;
    CDiff largest;
    
    while(left < len)
      {
	largest = ( right < len && comp(*(first + left), *(first + right)) ) ? right : left;

	if( comp(mov, *(first + largest)) )
	  {
	    *(first + index) = *(first + largest); 
	    //Performed a move, update the ID map
	    idMap[(first + index)->pID] = index;
	    index = largest;
	    left = 2 * index + 1;
	    right = 2 * index + 2;
	  } 
	else 
	  break;
      }
    
    if(index != pos - first)
      {
	*(first + index) = mov;
	//Performed the move, update the ID map
	idMap[mov.pID] = index;
      }
  }

};

#endif
