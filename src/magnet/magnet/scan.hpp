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

#include <magnet/detail/common.hpp>

namespace magnet {
  
  template<class T>
  class scan : public detail::functor<scan<T> > 
  {
    cl::Kernel _prescanKernel, _uniformAddKernel;
    
  public:
    scan(cl::CommandQueue queue, cl::Context context):
      detail::functor<scan<T> >(queue, context, "")
    {
      _prescanKernel = cl::Kernel(detail::functor<scan<T> >::_program, "prescan");
      _uniformAddKernel = cl::Kernel(detail::functor<scan<T> >::_program, "uniformAdd");
    }

    void operator()(cl::Buffer input, cl::Buffer output, cl_uint size)
    {
      cl_uint nGroups = (((size / 2) + 256 - 1) / 256);
      cl::KernelFunctor prescan
	= _prescanKernel.bind(detail::functor<scan<T> >::_queue, 
					 cl::NDRange(256 * nGroups), 
					 cl::NDRange(256));
      
      cl::Context context = (detail::functor<scan<T> >::_queue);

      cl::Buffer partialSums(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * (nGroups));
      
      prescan(input, output, partialSums, size);
      
      //Recurse if we've got to scan the sub buffers
      if (nGroups > 1)
	{
	  operator()(partialSums, partialSums, nGroups);
	  
	  cl::KernelFunctor uniformAdd
	    = _uniformAddKernel.bind(detail::functor<scan<T> >::_queue, 
				     cl::NDRange( 256 * nGroups), cl::NDRange(256));
	  
	  uniformAdd(output, output, partialSums, size);
	}
    }

    static inline std::string kernelSource();
  };
};

#include "magnet/scan.clh"

//void scan(cl::CommandQueue& queue, cl::Context& context, cl::Buffer& bufferIn, cl_uint size)
//{
//  cl_uint nGroups = (((size / 2) + 256 - 1) / 256);
//  cl::KernelFunctor prescan
//    = prescanKernel.bind(queue, cl::NDRange( 256 * nGroups), cl::NDRange(256));
//  
//  cl::Buffer partialSums(context, CL_MEM_WRITE_ONLY, sizeof(cl_uint) * (nGroups));
//  
//  events.push_back(prescan(bufferIn, bufferIn, partialSums, size));
//
//  //Recurse if we've got to scan the sub buffers
//  if (nGroups > 1)
//    {
//      scan(queue, context, partialSums, nGroups);
//
//      cl::KernelFunctor uniformAdd
//	= uniformAddKernel.bind(queue, cl::NDRange( 256 * nGroups), cl::NDRange(256));
//
//      events.push_back(uniformAdd(bufferIn, bufferIn, partialSums, size));
//    }
//}
