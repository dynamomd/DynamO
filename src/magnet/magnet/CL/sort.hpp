/*    dynamo:- Event driven molecular dynamics simulator 
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

#include <magnet/CL/radixsort_NVIDIA.hpp>
#include <magnet/CL/radixsort_AMD.hpp>
#include <magnet/CL/heapSort.hpp>

namespace magnet {
  namespace CL {
    template<class T>
    class sort
    {
    public:
      enum modeType {CPU, NVIDIA, AMD, UNSET };

    private:
      radixSortNVIDIA<T> _nvSorter;
      radixSortAMD<T> _amdSorter;
      heapSort<T> _cpuSorter;

      modeType _mode;

    public:
      sort():_mode(UNSET) {}

      modeType getMode() { return _mode; }

      size_t padding()
      {
	switch(_mode)
	  {
	  case CPU:
	    return 1;
	  case NVIDIA:
	    return 1024;
	  case AMD:
	    return 64 * 256;
	  default:
	    M_throw() << "Functor has not yet been built";
	  }
      }


      void build (cl::CommandQueue queue, cl::Context context)
      {
	cl::Device dev = queue.getInfo<CL_QUEUE_DEVICE>();

	if (dev.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU)
	  _mode = CPU;
	else
	  if (dev.getInfo<CL_DEVICE_VENDOR>().find("Advanced Micro Devices") != std::string::npos)
	    _mode = AMD;
	  else
	    _mode = NVIDIA;

	switch(_mode)
	  {
	  case CPU:
	    _cpuSorter.build(queue, context);
	    break;
	  case NVIDIA:
	    _nvSorter.build(queue, context);
	    break;
	  case AMD:
	    _amdSorter.build(queue, context);
	    break;
	  default:
	    M_throw() << "Could not determine which sorting algorithm to use";
	  }
      }

      void operator()(cl::Buffer input)
      {
	switch(_mode)
	  {
	  case CPU:
	    _cpuSorter(input);
	    break;
	  case NVIDIA:
	    _nvSorter(input, input);
	    break;
	  case AMD:
	    _amdSorter(input, input);
	    break;
	  default:
	    M_throw() << "Functor has not yet been built";
	  }
      }
	
      void operator()(cl::Buffer keyInput, cl::Buffer dataInput)
      {
	switch(_mode)
	  {
	  case CPU:
	    _cpuSorter(keyInput, dataInput);
	    break;
	  case NVIDIA:
	    _nvSorter(keyInput, dataInput, keyInput, dataInput);
	    break;
	  case AMD:
	    _amdSorter(keyInput, dataInput, keyInput, dataInput);
	    break;
	  default:
	    M_throw() << "Functor has not yet been built";
	  }
      }

    };
  }
}
