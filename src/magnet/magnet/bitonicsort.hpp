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
  class bitonicSort : public detail::functor<bitonicSort<T> >
  {
    cl::Kernel _sortKernel, _smallSortKernel, _subSortKernel;
    
  public:
    bitonicSort(cl::CommandQueue queue, cl::Context context):
      detail::functor<bitonicSort<T> >(queue, context, "")
    {
      // set up kernel
      _sortKernel = cl::Kernel(detail::functor<bitonicSort<T> >::_program, "bitonicSort");
      _smallSortKernel = cl::Kernel(detail::functor<bitonicSort<T> >::_program, 
				    "bitonicLocalSortKernel");
      _subSortKernel = cl::Kernel(detail::functor<bitonicSort<T> >::_program, 
				  "bitonicSubStageSort");
    }

    void operator()(cl::Buffer input, cl_uint ascending = true)
    {
      const cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);
      const cl_uint groupSize = 256;

      size_t numStages = 0;
      for (size_t temp = size; temp > 1; temp >>= 1) ++numStages;

      cl_uint powerOfTwo = 1 << numStages;
      
      if (size != powerOfTwo)
	M_throw() << "This bitonic sort only works on power of two sized arrays, size =" 
		  << size;
      
      cl::KernelFunctor clsort
	= _sortKernel.bind(detail::functor<bitonicSort<T> >::_queue, cl::NDRange(powerOfTwo), 
			   cl::NDRange(groupSize));
      
      cl::KernelFunctor clsmallsort
	= _smallSortKernel.bind(detail::functor<bitonicSort<T> >::_queue, 
				cl::NDRange(powerOfTwo / 2), cl::NDRange(groupSize));
      
      cl::KernelFunctor clsubsort
	= _subSortKernel.bind(detail::functor<bitonicSort<T> >::_queue, 
			      cl::NDRange(powerOfTwo / 2), cl::NDRange(groupSize));
      
      //All stages except the last one use the reverse intended direction!
      cl_uint initial_direction = 1-ascending;
      
      std::vector<cl::Event> events;

      ///////////////////////////////////////////////////////////////////////
      //Small sort on blocks of up to 512
      clsmallsort(input, initial_direction);
      
      //Now for the full sort
      for (cl_uint stage = 9; stage < numStages-1; ++stage)
	{
	  //Do the first (stage - 8) passes using the slow kernel
	  for (cl_uint stagePass = 0; stagePass < stage-8/*stage+1*/; ++stagePass)
	    clsort(input, stage, stagePass, size, initial_direction);
	  
	  //The final passes can be done using the small sort kernel
	  clsubsort(input, size, initial_direction, stage);
	}

      {
	cl_uint stage = numStages-1;	
	for (cl_uint stagePass = 0; stagePass < stage+1; ++stagePass)
	  clsort(input, stage, stagePass, size, ascending);
      }
    }

    static inline std::string kernelSource();
  };
};

#include <magnet/detail/kernels/bitonicSort.clh>
