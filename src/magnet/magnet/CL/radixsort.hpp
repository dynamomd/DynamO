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

namespace magnet {
  namespace CL {
    template<class T>
    class radixSort : public detail::functor<radixSort<T> >
    {
      cl::Kernel _radixSortKernel, 
	_findRadixOffsetsKernel, 
	_reorderKeysKernel,
	_radixSortDataKernel,
	_reorderKeysDataKernel;

      scan<cl_uint> _scanFunctor;
    
      cl::Buffer _buckets, _offsets, _doubleBuffer, _dataDoubleBuffer;
    
      cl_uint _lastSize, _lastDataSize;

      //cl::KernelFunctor _dataSort;

    public:
      radixSort():
	_lastSize(0),
	_lastDataSize(0)
      {}

      void build (cl::CommandQueue queue, cl::Context context)
      {
	//Build the scan functor
	_scanFunctor.build(queue, context);

	//And build this functor
	detail::functor<radixSort<T> >::build(queue, context, "");

	_radixSortKernel 
	  = cl::Kernel(detail::functor<radixSort<T> >::_program, "radixBlockSortKernel");
	_findRadixOffsetsKernel 
	  = cl::Kernel(detail::functor<radixSort<T> >::_program, "findRadixOffsetsKernel");
	_reorderKeysKernel 
	  = cl::Kernel(detail::functor<radixSort<T> >::_program, "reorderKeys");
	_radixSortDataKernel 
	  = cl::Kernel(detail::functor<radixSort<T> >::_program, "radixBlockSortDataKernel");
	_reorderKeysDataKernel 
	  = cl::Kernel(detail::functor<radixSort<T> >::_program, "reorderKeysData");
      }

      void operator()(cl::Buffer keyInput, cl::Buffer keyOutput, cl_uint bits_to_sort = 0, cl_uint bitsPerPass = 4)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 256;
	cl_uint nWorkGroups = ((size / 4) + groupSize - 1) / groupSize;
	cl_uint maxRadixDigit = 2 << (bitsPerPass - 1);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % 1024)				
	  M_throw() << "Radix sort works on whole multiples of 1024 elements only, please pad your data";

	cl::KernelFunctor _blockSortFunc
	  = _radixSortKernel.bind(detail::functor<radixSort<T> >::_queue, 
				  cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _findRadixOffsetsFunc
	  = _findRadixOffsetsKernel.bind(detail::functor<radixSort<T> >::_queue, 
					 cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _ReorderFunc
	  = _reorderKeysKernel.bind(detail::functor<radixSort<T> >::_queue, 
				    cl::NDRange(size/4), cl::NDRange(groupSize));
       
	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(detail::functor<radixSort<T> >::_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _offsets = cl::Buffer(detail::functor<radixSort<T> >::_context, 
				  CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _doubleBuffer = cl::Buffer(detail::functor<radixSort<T> >::_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }
      
      
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _blockSortFunc(keyInput, _doubleBuffer, size, startBit, bitsPerPass);
	  
	    _findRadixOffsetsFunc(_doubleBuffer, _buckets, _offsets, 
				  size, startBit, bitsPerPass,
				  cl::__local(sizeof(cl_uint) * maxRadixDigit));
	  
	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	  
	    _ReorderFunc(_doubleBuffer, keyOutput, _buckets, _offsets, 
			 size, startBit, bitsPerPass,
			 cl::__local(sizeof(cl_uint) * maxRadixDigit),
			 cl::__local(sizeof(cl_uint) * maxRadixDigit)
			 );
	  }

	_lastSize = size;
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput,
		      cl::Buffer keyOutput, cl::Buffer dataOutput, 
		      cl_uint bits_to_sort = 0, cl_uint bitsPerPass = 4)
      {
	if (bits_to_sort == 0)
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);

	if ((dataInput.getInfo<CL_MEM_SIZE>() / sizeof(cl_uint)) != size)
	  M_throw() << "Key and data set size mismatch";

	cl_uint groupSize = 256;
	cl_uint nWorkGroups = ((size / 4) + groupSize - 1) / groupSize;
	cl_uint maxRadixDigit = 2 << (bitsPerPass - 1);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % 1024)				
	  M_throw() << "Radix sort works on whole multiples of 1024 elements only, please pad your data";

	cl::KernelFunctor _blockSortDataFunc
	  = _radixSortDataKernel.bind(detail::functor<radixSort<T> >::_queue, 
				      cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _findRadixOffsetsFunc
	  = _findRadixOffsetsKernel.bind(detail::functor<radixSort<T> >::_queue, 
					 cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _reorderDataFunc
	  = _reorderKeysDataKernel.bind(detail::functor<radixSort<T> >::_queue, 
					cl::NDRange(size/4), cl::NDRange(groupSize));
      
	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(detail::functor<radixSort<T> >::_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _offsets = cl::Buffer(detail::functor<radixSort<T> >::_context, 
				  CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _doubleBuffer = cl::Buffer(detail::functor<radixSort<T> >::_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }
 

	if (_lastDataSize != size)
	  {
	    _dataDoubleBuffer = cl::Buffer(detail::functor<radixSort<T> >::_context, 
					   CL_MEM_READ_WRITE, sizeof(cl_uint) * size);
	  }
     
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _blockSortDataFunc(keyInput, dataInput, _doubleBuffer, _dataDoubleBuffer, size, 
			       startBit, bitsPerPass)
	      ;
	  
	    _findRadixOffsetsFunc(_doubleBuffer, _buckets, _offsets, 
				  size, startBit, bitsPerPass,
				  cl::__local(sizeof(cl_uint) * maxRadixDigit))
	      ;
	  
	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	  
	    _reorderDataFunc(_doubleBuffer, _dataDoubleBuffer, keyOutput, dataOutput, 
			     _buckets, _offsets, 
			     size, startBit, bitsPerPass,
			     cl::__local(sizeof(cl_uint) * maxRadixDigit),
			     cl::__local(sizeof(cl_uint) * maxRadixDigit)
			     )
	      ;
	  }

	_lastDataSize = _lastSize = size;
      }

      static inline std::string kernelSource();
    };
  }
}

#include "magnet/CL/detail/kernels/radixsort.clh"
