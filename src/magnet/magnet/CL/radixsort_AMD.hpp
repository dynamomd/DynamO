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
#include <magnet/CL/scan.hpp>

#include <sstream>
#include <fstream>

namespace magnet {
  namespace CL {
    template<class T>
    class radixSortAMD : public detail::functor<radixSortAMD<T> >
    {
      cl::Kernel _histogramKernel, 
	_permuteKernel,
	_dataPermuteKernel;

      scan<cl_uint> _scanFunctor;
    
      cl::Buffer _buckets, _offsets, _doubleBuffer, _dataDoubleBuffer;
    
      cl_uint _lastSize, _lastDataSize;

    public:
      radixSortAMD():
	_lastSize(0),
	_lastDataSize(0)
      {}

      void build(cl::CommandQueue queue, cl::Context context)
      {
	//Build the scan functor
	_scanFunctor.build(queue, context);

	//And build this functor
	detail::functor<radixSortAMD<T> >::build(queue, context, "");

	_histogramKernel 
	  = cl::Kernel(detail::functor<radixSortAMD<T> >::_program, "histogram");
	_permuteKernel 
	  = cl::Kernel(detail::functor<radixSortAMD<T> >::_program, "permute");
	_dataPermuteKernel 
	  = cl::Kernel(detail::functor<radixSortAMD<T> >::_program, "datapermute");
      }

      void operator()(cl::Buffer keyInput, cl::Buffer keyOutput, cl_uint bits_to_sort = 0)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 64;
	cl_uint bitsPerPass = 4;
	cl_uint maxRadixDigit = 1 << bitsPerPass;
	cl_uint keysPerWorkitem = 256;
	cl_uint nWorkGroups = size / (groupSize * keysPerWorkitem);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % (groupSize * keysPerWorkitem))
	  M_throw() << "Radix sort works on whole multiples of " 
		    << groupSize * keysPerWorkitem << " elements only, please pad your data";

	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(detail::functor<radixSortAMD<T> >::_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * groupSize * maxRadixDigit);
	    
	    _doubleBuffer = cl::Buffer(detail::functor<radixSortAMD<T> >::_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }

	cl::KernelFunctor _histogramFunc
	  = _histogramKernel.bind(detail::functor<radixSortAMD<T> >::_queue, 
				  cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
      
	cl::KernelFunctor _permuteFunc
	  = _permuteKernel.bind(detail::functor<radixSortAMD<T> >::_queue, 
				cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
       
	try {
	  detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(keyInput, keyOutput, 0, 0, 
									 sizeof(T) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If keyInput and keyOutput are equal no need to copy
	  }

	static int number = 0;
	++number;
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _histogramFunc(keyOutput, 
			   _buckets, 
			   startBit,
			   cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			   size, keysPerWorkitem, bitsPerPass
			   );

	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	    

	    _permuteFunc(keyOutput, _buckets, startBit, 
			 cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			 _doubleBuffer,
			 size, keysPerWorkitem, bitsPerPass
			 );
	    
	    detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(_doubleBuffer, keyOutput, 0, 0, 
									   sizeof(T) * size);
	  }

	_lastSize = size;
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput,
		      cl::Buffer keyOutput, cl::Buffer dataOutput, 
		      cl_uint bits_to_sort = 0)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 64;
	cl_uint bitsPerPass = 4;
	cl_uint maxRadixDigit = 1 << bitsPerPass;
	cl_uint keysPerWorkitem = 256;
	cl_uint nWorkGroups = size / (groupSize * keysPerWorkitem);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % (groupSize * keysPerWorkitem))
	  M_throw() << "Radix sort works on whole multiples of " 
		    << groupSize * keysPerWorkitem << " elements only, please pad your data";

	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(detail::functor<radixSortAMD<T> >::_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * groupSize * maxRadixDigit);
	    
	    _doubleBuffer = cl::Buffer(detail::functor<radixSortAMD<T> >::_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }

	if (_lastDataSize != size)
	  {
	    _dataDoubleBuffer = cl::Buffer(detail::functor<radixSortAMD<T> >::_context, 
					   CL_MEM_READ_WRITE, sizeof(cl_uint) * size);
	  }

	cl::KernelFunctor _histogramFunc
	  = _histogramKernel.bind(detail::functor<radixSortAMD<T> >::_queue, 
				  cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
      
	cl::KernelFunctor _dataPermuteFunc
	  = _dataPermuteKernel.bind(detail::functor<radixSortAMD<T> >::_queue, 
				    cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
       
	try {
	  detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(keyInput, keyOutput, 0, 0, 
									 sizeof(T) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If keyInput and keyOutput are equal no need to copy
	  }

	try {
	  detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(dataInput, dataOutput, 0, 0, 
									 sizeof(cl_uint) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If dataInput and dataOutput are equal no need to copy
	  }

	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _histogramFunc(keyOutput, 
			   _buckets, 
			   startBit,
			   cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			   size, keysPerWorkitem, bitsPerPass
			   );

	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	    

	    _dataPermuteFunc(keyOutput, _buckets, dataOutput, startBit, 
			     cl::__local(sizeof(cl_uint) * maxRadixDigit * groupSize),
			     _doubleBuffer, _dataDoubleBuffer,
			     size, keysPerWorkitem, bitsPerPass
			     );
	    
	    detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(_doubleBuffer, keyOutput, 0, 0, 
									   sizeof(T) * size);
	    detail::functor<radixSortAMD<T> >::_queue.enqueueCopyBuffer(_dataDoubleBuffer, dataOutput, 0, 0, 
									   sizeof(cl_uint) * size);
	  }

	_lastDataSize = _lastSize = size;
      }

      static inline std::string kernelSource();
    };
  }
}

#include "magnet/CL/detail/kernels/radixsort_AMD.clh"
