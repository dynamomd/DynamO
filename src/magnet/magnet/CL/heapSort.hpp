/*    DYNAMO:- Event driven molecular dynamics simulator 
 *    http://www.marcusbannerman.co.uk/dynamo
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <magnet/CL/detail/common.hpp>

namespace magnet {
  namespace CL {
    template<class T>
    class heapSort : public detail::functor<heapSort<T> >
    {
      cl::Kernel _sortKernel, _dataSortKernel;
    
    public:
      void build(cl::CommandQueue queue, cl::Context context)
      {
	detail::functor<heapSort<T> >::build(queue, context, "");

	// set up kernel
	_sortKernel = cl::Kernel(detail::functor<heapSort<T> >::_program, "heapSort");

	_dataSortKernel = cl::Kernel(detail::functor<heapSort<T> >::_program, "heapSortData");	
      }

      void operator()(cl::Buffer input, cl_uint ascending = true)
      {
	const cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);

	cl::KernelFunctor clsort = _sortKernel.bind(detail::functor<heapSort<T> >::_queue, cl::NDRange(1), 
						     cl::NDRange(1));      
	clsort(input, size);
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput)
      {
	const cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	if (size != dataInput.getInfo<CL_MEM_SIZE>() / sizeof(cl_uint))
	  M_throw() << "Data-key buffer size mismatch";
	

	cl::KernelFunctor clsort = _dataSortKernel.bind(detail::functor<heapSort<T> >::_queue, cl::NDRange(1), 
							cl::NDRange(1));      
	clsort(keyInput, dataInput, size);

      }
      static inline std::string kernelSource();
    };
  }
}
#include <magnet/CL/detail/kernels/heapSort.clh>
