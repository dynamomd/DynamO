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
      _smallSortKernel = cl::Kernel(detail::functor<bitonicSort<T> >::_program, "bitonicLocalSort");
      _subSortKernel = cl::Kernel(detail::functor<bitonicSort<T> >::_program, "bitonicSubStageSort");
    }

    void operator()(cl::Buffer input, cl::Buffer output, bool ascending = true)
    {
      cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);
      
      size_t numStages = 0;
      for (size_t temp = size; temp > 1; temp >>= 1) ++numStages;

      cl_uint powerOfTwo = 1 << numStages;
      
      if (size != powerOfTwo)
	M_throw() << "Cannot use this bitonic sort on non-power of two sized arrays, size = " << size;

      cl::KernelFunctor clsort
	= _sortKernel.bind(detail::functor<bitonicSort<T> >::_queue, cl::NDRange(powerOfTwo), 
			   cl::NDRange(groupSize));
      
      cl::KernelFunctor clsmallsort
	= _smallSortKernel.bind(detail::functor<bitonicSort<T> >::_queue, cl::NDRange(powerOfTwo / 2), 
				cl::NDRange(groupSize));
      
      cl::KernelFunctor clsubsort
	= _subSortKernel.bind(detail::functor<bitonicSort<T> >::_queue, cl::NDRange(powerOfTwo / 2), 
			      cl::NDRange(groupSize));
      
      //All stages except the last one use the reverse intended direction!
      cl_uint initial_direction = 1-ascending;
      
      std::vector<cl::Event> events;

      ///////////////////////////////////////////////////////////////////////
      //Small sort on blocks of up to 512
      clsmallsort(bufferIn, size, initial_direction, 
		  cl::__local(2 * sizeof(cl_uint ) * groupSize));
      
      //Now for the full sort
      for (cl_uint stage = 9; stage < numStages-1; ++stage)
	{
	  //Do the first (stage - 8) passes using the slow kernel
	  for (cl_uint stagePass = 0; stagePass < stage-8/*stage+1*/; ++stagePass)
	    clsort(bufferIn, stage, stagePass, size, initial_direction);
	  
	  //The final passes can be done using the small sort kernel
	  clsubsort(bufferIn, size, initial_direction, 
		    cl::__local(2 * sizeof(cl_uint ) * groupSize),
		    stage);
	}

      {
	cl_uint stage = numStages-1;	
	for (cl_uint stagePass = 0; stagePass < stage+1; ++stagePass)
	  clsort(bufferIn, stage, stagePass, size, ascending);
      }

//
//      //Workgroups of 256 work-items process 512 elements
//      cl_uint nGroups = (((size / 2) + 256 - 1) / 256);
//      
//      cl::Buffer partialSums(detail::functor<bitonicSort<T> >::_context,
//			     CL_MEM_READ_WRITE, sizeof(cl_uint) * (nGroups));
//      
//      _prescanKernel.bind(detail::functor<bitonicSort<T> >::_queue, 
//			  cl::NDRange(256 * nGroups), 
//			  cl::NDRange(256))
//	(input, output, partialSums, size);
//      
//      //Recurse if we've got more than one workgroup
//      if (nGroups > 1)
//	{
//	  operator()(partialSums, partialSums);
//	  
//	  _uniformAddKernel.bind(detail::functor<bitonicSort<T> >::_queue, 
//				 cl::NDRange(nGroups * 256), 
//				 cl::NDRange(256))
//	    (output, output, partialSums, size);
//	}
    }

    static inline std::string kernelSource();
  };
};

#include "magnet/detail/kernels/bitonicSort.clh"
