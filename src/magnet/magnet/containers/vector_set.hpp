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
#include <vector>
#include <algorithm>

namespace magnet {
  namespace containers {
    /*! \brief An unorderd Set container based on the std::vector
        class.
	
	This is a simple std::set like implementation using a
	std::vector. This means that certain operations (such as find)
	are O(N), but for small sets this container will greatly
	outperform other approaches as it should fit in the cache.

	The set is unordered to allow a fast erase operation.
     */
    template <typename T>
    class VectorSet : public std::vector<T> {
    public:
      void erase(size_t particle) {
	typename std::vector<T>::iterator pit = find(particle);
#ifdef DYNAMO_DEBUG
	if (pit == end())
	  M_throw() << "Removing a particle " << particle << " which is not in this cell";
#endif
	//Base::erase(pit); 

	//If we erase as above, this forces a shuffle of the entries,
	//as they are unordered we can speed this up by shuffling only
	//one entry. This may result in an extra copy if the set has
	//only one item.
	std::swap(*pit, std::vector<T>::back());
	std::vector<T>::pop_back();
      }

      void insert(size_t particle) { std::vector<T>::push_back(particle); }
      
      typename std::vector<T>::const_iterator find(const T& val) const {
	return std::find(std::vector<T>::begin(), std::vector<T>::end(), val);
      }

      typename std::vector<T>::iterator find(const T& val) {
	return std::find(std::vector<T>::begin(), std::vector<T>::end(), val);
      }

      typename std::vector<T>::const_iterator count(const T& val) {
	return std::find(std::vector<T>::begin(), std::vector<T>::end(), val) != std::vector<T>::end();
      }
    };

  }
}
