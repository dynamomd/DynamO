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
      cl::KernelFunctor _clsort, _cldatasort;
    
    public:
      heapSort(cl::CommandQueue queue, cl::Context context):
	detail::functor<heapSort<T> >(queue, context, "")
      {
	// set up kernel
	_sortKernel = cl::Kernel(detail::functor<heapSort<T> >::_program, "heapSort");

	_dataSortKernel = cl::Kernel(detail::functor<heapSort<T> >::_program, "heapDataSort");
	
	_clsort = _sortKernel.bind(detail::functor<heapSort<T> >::_queue, cl::NDRange(1), 
				   cl::NDRange(1));

	_cldatasort = _dataSortKernel.bind(detail::functor<heapSort<T> >::_queue, cl::NDRange(1), 
				   cl::NDRange(1));
      }

      void operator()(cl::Buffer input, cl_uint ascending = true)
      {
	const cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);

      
	//Small sort on blocks of up to 512
	clsmallsort(input, size);
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput)
      {

      }
      static inline std::string kernelSource();
    };
  }
}
#include <magnet/CL/detail/kernels/heapSort.clh>
