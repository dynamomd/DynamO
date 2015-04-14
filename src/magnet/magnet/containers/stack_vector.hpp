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

#include <magnet/exception.hpp>
#include <array>

namespace magnet {
  namespace containers {
    /*! \brief Stack allocated std::vector.

      This class is a fast and lightweight std::vector-like
      container. It is useful for returning small numbers of items
      from functions without allocating memory on the heap.
     */
    template<class T, size_t Nmax>
    class StackVector: public std::array<T, Nmax> {
      typedef std::array<T, Nmax> Base;
    public:
      template<size_t Nmax2>
      StackVector(const StackVector<T, Nmax2>& vec):
	Base()
      {
	static_assert(Nmax2 <= Nmax, "Can only convert to larger StackVector containers");
	_size = vec.size();
	std::copy(vec.begin(), vec.end(), Base::begin());
      }

      StackVector(): Base(), _size(0) {}
      
      StackVector(std::initializer_list<T> _list):
	Base(),
	_size(0)
      {
	auto it = _list.begin();
	for (size_t i(0); (i < Nmax) && (it != _list.end()); ++i, ++it)
	  push_back(*it);
      }
            
      constexpr typename Base::size_type size() const { return _size; }
      constexpr bool empty() const { return size() == 0; }

      typename Base::iterator end() { return typename Base::iterator(Base::data() + _size); }
      typename Base::const_iterator end() const { return typename Base::iterator(Base::data() + _size); }
      typename Base::const_iterator cend() const { return typename Base::iterator(Base::data() + _size); }

      typename Base::reverse_iterator rbegin() { return typename Base::reverse_iterator(Base::end()); }
      typename Base::const_reverse_iterator rbegin() const { return typename Base::const_reverse_iterator(Base::end()); }
      typename Base::const_reverse_iterator crbegin() const { return typename Base::const_reverse_iterator(Base::end()); }
      typename Base::reverse_iterator rend() { return typename Base::reverse_iterator(Base::begin()); }
      typename Base::const_reverse_iterator rend() const { return typename Base::const_reverse_iterator(Base::begin()); }
      typename Base::const_reverse_iterator crend() const { return typename Base::const_reverse_iterator(Base::begin()); }
      
      typename Base::reference back() { return _size ? *(Base::end() - 1) : *Base::end(); }
      typename Base::const_reference back() const { return _size ? *(Base::end() - 1) : *Base::end(); }

      void push_back(const T& val) {
#ifdef MAGNET_DEBUG
	if (_size+1 > Nmax)
	  M_throw() << "Cannot push elements to a filled StackVector " << *this;
#endif
	Base::operator[](_size) = val;
	++_size;
      }

      T pop_back() {
#ifdef MAGNET_DEBUG
	if (empty())
	  M_throw() << "Cannot pop elements from an emptry StackVector " << *this;
#endif
	return Base::operator[](--_size);
      }

      template<size_t Nmax2>
      void extend(const StackVector<T,Nmax2>& ovec) {
	for (const T& a: ovec)
	  push_back(a);
      }
      
    private:
      size_t _size;
    };

    template<class T, size_t Nmax>
    std::ostream& operator<<(std::ostream& os, const StackVector<T,Nmax>&s) {
      os << "StackVector{ ";
      for (const auto& val : s)
	os << val << " ";
      os << "}";
      return os;
    }

    template<class T1, class T2, size_t Nmax>
    std::ostream& operator<<(std::ostream& os, const StackVector<std::pair<T1,T2>,Nmax>&s) {
      os << "StackVector{ ";
      for (const auto& val : s)
	os << "[" << val.first << ", " << val.second << "] ";
      os << "}";
      return os;
    }

    template<std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    tuple_print(const std::tuple<Tp...>& t, std::ostream& os)
    { }
    
    template<std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    tuple_print(const std::tuple<Tp...>& t, std::ostream& os)
    {
      os << std::get<I>(t) << " ";
      tuple_print<I + 1, Tp...>(t, os);
    }
      
    template<size_t Nmax, typename... Tp>
    std::ostream& operator<<(std::ostream& os, const StackVector<std::tuple<Tp...>,Nmax>&s) 
    {
      os << "StackVector{ ";
      for (const auto& val : s) {
	os << "[";
	tuple_print(val, os);
	os << "] ";
      }
      os << "}";
      return os;
    }
  }
}
    
