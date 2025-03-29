/*    dynamo:- Event driven molecular dynamics simulator
 *    http://www.dynamomd.org
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
#define STRINGIFY(A) #A

namespace magnet {
namespace CL {
template <class T> class bitonicSort : public detail::Program {
  cl::Kernel _sortKernel, _smallSortKernel, _subSortKernel;

public:
  void build(cl::CommandQueue queue, cl::Context context) {
    Program::build(queue, context, "");
    // set up kernel
    _sortKernel = cl::Kernel(_program, "bitonicSort");
    _smallSortKernel = cl::Kernel(_program, "bitonicLocalSortKernel");
    _subSortKernel = cl::Kernel(_program, "bitonicSubStageSort");
  }

  void operator()(cl::Buffer input, cl_uint ascending = true) {
    const cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);
    const cl_uint groupSize = 256;

    size_t numStages = 0;
    for (size_t temp = size; temp > 1; temp >>= 1)
      ++numStages;

    cl_uint powerOfTwo = 1 << numStages;

    if (size != powerOfTwo)
      M_throw()
          << "This bitonic sort only works on power of two sized arrays, size ="
          << size;

    cl::KernelFunctor clsort = _sortKernel.bind(_queue, cl::NDRange(powerOfTwo),
                                                cl::NDRange(groupSize));

    cl::KernelFunctor clsmallsort = _smallSortKernel.bind(
        _queue, cl::NDRange(powerOfTwo / 2), cl::NDRange(groupSize));

    cl::KernelFunctor clsubsort = _subSortKernel.bind(
        _queue, cl::NDRange(powerOfTwo / 2), cl::NDRange(groupSize));

    // All stages except the last one use the reverse intended direction!
    cl_uint initial_direction = 1 - ascending;

    ///////////////////////////////////////////////////////////////////////
    // Small sort on blocks of up to 512
    clsmallsort(input, initial_direction);

    // Now for the full sort
    for (cl_uint stage = 9; stage < numStages - 1; ++stage) {
      // Do the first (stage - 8) passes using the slow kernel
      for (cl_uint stagePass = 0; stagePass < stage - 8 /*stage+1*/;
           ++stagePass)
        clsort(input, stage, stagePass, size, initial_direction);

      // The final passes can be done using the small sort kernel
      clsubsort(input, size, initial_direction, stage);
    }

    {
      cl_uint stage = numStages - 1;
      for (cl_uint stagePass = 0; stagePass < stage + 1; ++stagePass)
        clsort(input, stage, stagePass, size, ascending);
    }
  }

  virtual std::string initKernelSrc() {
    return std::string("#define keyType ") + detail::traits<T>::kernel_type() +
           "\n" STRINGIFY(
               __kernel void bitonicSort(
                   __global keyType * theArray, const uint stage,
                   const uint passOfStage, const uint realsize,
                   const uint direction) {
                 uint threadId = get_global_id(0);

                 uint pairDistance = 1 << (stage - passOfStage);
                 uint blockWidth = 2 * pairDistance;

                 uint leftId = (threadId % pairDistance) +
                               (threadId / pairDistance) * blockWidth;

                 uint rightId = leftId + pairDistance;

                 if ((leftId >= realsize) || (rightId >= realsize))
                   return;

                 uint leftElement = theArray[leftId];
                 uint rightElement = theArray[rightId];

                 uint sameDirectionBlockWidth = 1 << stage;

                 uint sortIncreasing = direction;
                 sortIncreasing = (threadId / sameDirectionBlockWidth) % 2 == 1
                                      ? 1 - sortIncreasing
                                      : sortIncreasing;

                 keyType greater =
                     leftElement > rightElement ? leftElement : rightElement;

                 keyType lesser =
                     leftElement > rightElement ? rightElement : leftElement;

                 theArray[leftId] = sortIncreasing ? lesser : greater;
                 theArray[rightId] = sortIncreasing ? greater : lesser;
               }

               // This sorts small blocks of power of two elements that fit into
               // shared memory
               void bitonicLocalSort(__local keyType * cache, uint direction,
                                     uint size) {
                 for (uint stageStride = 2; stageStride < size;
                      stageStride <<= 1) {
                   uint blockDirection =
                       (get_local_id(0) & (stageStride / 2)) == direction;

                   for (uint passStride = stageStride / 2; passStride > 0;
                        passStride >>= 1) {
                     barrier(CLK_LOCAL_MEM_FENCE);
                     uint pos = 2 * get_local_id(0) -
                                (get_local_id(0) & (passStride - 1));
                     if ((cache[pos] < cache[pos + passStride]) ==
                         blockDirection) {
                       keyType tmp = cache[pos];
                       cache[pos] = cache[pos + passStride];
                       cache[pos + passStride] = tmp;
                     }
                   }
                 }

                 uint blockDirection = (get_group_id(0) & 0x1) == direction;
                 for (uint passStride = size / 2; passStride > 0;
                      passStride >>= 1) {
                   barrier(CLK_LOCAL_MEM_FENCE);
                   uint pos = 2 * get_local_id(0) -
                              (get_local_id(0) & (passStride - 1));
                   if ((cache[pos] < cache[pos + passStride]) ==
                       blockDirection) {
                     keyType tmp = cache[pos];
                     cache[pos] = cache[pos + passStride];
                     cache[pos + passStride] = tmp;
                   }
                 }

                 barrier(CLK_LOCAL_MEM_FENCE);
               }

               __kernel void bitonicLocalSortKernel(
                   __global keyType * inputArray, uint direction) {
                 __local keyType cache[512];

                 int offset =
                     2 * get_local_size(0) * get_group_id(0) + get_local_id(0);

                 inputArray += offset;

                 // Load
                 cache[get_local_id(0)] = inputArray[0];
                 cache[get_local_id(0) + get_local_size(0)] =
                     inputArray[get_local_size(0)];

                 // Sort
                 bitonicLocalSort(cache, direction, 2 * get_local_size(0));

                 // Store
                 inputArray[0] = cache[get_local_id(0)];
                 inputArray[get_local_size(0)] =
                     cache[get_local_id(0) + get_local_size(0)];
               }

               __kernel __attribute__((reqd_work_group_size(
                   256, 1,
                   1))) void bitonicSubStageSort(__global keyType * inputArray,
                                                 uint realsize, uint direction,
                                                 uint stage) {
                 //__local block must be 2*get_local_size(0) big
                 // stage must be >= 8
                 __local keyType cache[512];

                 // offset the data array pointer
                 int offset =
                     2 * get_local_size(0) * get_group_id(0) + get_local_id(0);
                 inputArray += offset;

                 if (offset < realsize)
                   cache[get_local_id(0)] = inputArray[0];

                 if ((offset + get_local_size(0)) < realsize)
                   cache[get_local_id(0) + get_local_size(0)] =
                       inputArray[get_local_size(0)];

                 uint stageStride = 2 << stage;
                 {
                   uint blockDirection =
                       (get_local_id(0) & ((stageStride / 2) << stage)) != 0;

                   for (uint passStride = get_local_size(0); passStride > 0;
                        passStride >>= 1) {
                     barrier(CLK_LOCAL_MEM_FENCE);
                     uint pos = 2 * get_local_id(0) -
                                (get_local_id(0) & (passStride - 1));
                     if ((offset + pos + passStride - get_local_id(0)) <
                         realsize)
                       if ((cache[pos] < cache[pos + passStride]) ==
                           blockDirection) {
                         uint tmp = cache[pos];
                         cache[pos] = cache[pos + passStride];
                         cache[pos + passStride] = tmp;
                       }
                   }
                 }

                 uint blockDirection =
                     !(get_group_id(0) & (0x1 << (stage - 8)));
                 for (uint passStride = get_local_size(0); passStride > 0;
                      passStride >>= 1) {
                   barrier(CLK_LOCAL_MEM_FENCE);
                   uint pos = 2 * get_local_id(0) -
                              (get_local_id(0) & (passStride - 1));
                   if ((offset + pos + passStride - get_local_id(0)) < realsize)
                     if ((cache[pos] < cache[pos + passStride]) ==
                         blockDirection) {
                       uint tmp = cache[pos];
                       cache[pos] = cache[pos + passStride];
                       cache[pos + passStride] = tmp;
                     }
                 }

                 barrier(CLK_LOCAL_MEM_FENCE);
                 if (offset < realsize)
                   inputArray[0] = cache[get_local_id(0)];

                 if ((offset + get_local_size(0)) < realsize)
                   inputArray[get_local_size(0)] =
                       cache[get_local_id(0) + get_local_size(0)];
               });
  }
};
} // namespace CL
} // namespace magnet
#undef STRINGIFY
