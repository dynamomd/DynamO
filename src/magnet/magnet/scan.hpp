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
      //Workgroups of 256 work-items process 512 elements
      cl_uint nGroups = (((size / 2) + 256 - 1) / 256);
      
      cl::Buffer partialSums(detail::functor<scan<T> >::_context,
			     CL_MEM_READ_WRITE, sizeof(cl_uint) * (nGroups));
      
      _prescanKernel.bind(detail::functor<scan<T> >::_queue, 
			  cl::NDRange(256 * nGroups), 
			  cl::NDRange(256))
	(input, output, partialSums, size);
      
      //Recurse if we've got more than one workgroup
      if (nGroups > 1)
	{
	  operator()(partialSums, partialSums, nGroups);
	  
	  _uniformAddKernel.bind(detail::functor<scan<T> >::_queue, 
				 cl::NDRange(nGroups * 256), 
				 cl::NDRange(256))
	    (output, output, partialSums, size);
	}
    }

    static inline std::string kernelSource();
  };
};

#include "magnet/scan.clh"
