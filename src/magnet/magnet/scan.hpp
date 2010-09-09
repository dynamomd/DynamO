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
    

    std::vector<cl::Buffer> _partialSumBufferStack;
    cl_uint _lastSize;

  public:
    scan(cl::CommandQueue queue, cl::Context context):
      detail::functor<scan<T> >(queue, context, ""),
      _lastSize(0)
    {
      _prescanKernel = cl::Kernel(detail::functor<scan<T> >::_program, "prescan");
      _uniformAddKernel = cl::Kernel(detail::functor<scan<T> >::_program, "uniformAdd");
    }

    inline void operator()(cl::Buffer input, cl::Buffer output)
    {
      cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);

      //Workgroups of 256 work-items process 512 elements
      cl_uint nGroups = (((size / 2) + 256 - 1) / 256);

      //Check if the size has changed from last call,
      //We cache the partial sum buffers as reallocating buffers is slow
      if (size != _lastSize) 
	{
	  //Rebuild the stack of buffers
	  _partialSumBufferStack.clear();

	  for (cl_uint stageSize = nGroups; stageSize > 1; 
	       stageSize = (stageSize + 511) / 512)
	    _partialSumBufferStack.push_back(cl::Buffer(detail::functor<scan<T> >::_context,
							CL_MEM_READ_WRITE, sizeof(cl_uint) * stageSize));	  

	  _partialSumBufferStack.push_back(cl::Buffer(detail::functor<scan<T> >::_context,
						      CL_MEM_READ_WRITE, sizeof(cl_uint)));	  
	}

      //Start the recursion
      recursionFunction(input, output, size, 0);

      _lastSize = size;
    }

    static inline std::string kernelSource();

  protected:
    inline void recursionFunction(cl::Buffer input, cl::Buffer output, cl_uint size, cl_uint stage)
    {
      cl_uint nGroups = (((size / 2) + 256 - 1) / 256);

      _prescanKernel.bind(detail::functor<scan<T> >::_queue, 
			  cl::NDRange(256 * nGroups), 
			  cl::NDRange(256))
	(input, output, _partialSumBufferStack[stage], size);
      
      
      //Recurse if we've got more than one workgroup
      if (nGroups > 1)
	{
	  recursionFunction(_partialSumBufferStack[stage], _partialSumBufferStack[stage], 
			    nGroups, stage+1);
	  
	  _uniformAddKernel.bind(detail::functor<scan<T> >::_queue,
				 cl::NDRange(nGroups * 256), 
				 cl::NDRange(256))
	    (output, output, _partialSumBufferStack[stage], size);
	}
    }
  };
};

#include "magnet/detail/kernels/scan.clh"
