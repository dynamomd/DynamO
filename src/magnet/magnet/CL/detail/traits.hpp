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

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

// Simple macro to convert a token to a string
#define STRINGIFY(A) #A

// This is for defining scalar types in CL_TYPE_FACTORY
#define SCALAR_TYPE(F, kernel_type)                                            \
  F(kernel_type, cl_##kernel_type, 1, cl_##kernel_type)

// This is for defining vector types in CL_TYPE_FACTORY
#define VEC_TYPE(F, kernel_type, kernel_bitshift_type)                         \
  F(kernel_type, cl_##kernel_type, 1, cl_##kernel_type, cl_##kernel_type##2,   \
    cl_##kernel_type##4, cl_##kernel_type##8, cl_##kernel_type##16,            \
    cl_##kernel_bitshift_type)                                                 \
  F(kernel_type##2, cl_##kernel_type##2, 2, cl_##kernel_type,                  \
    cl_##kernel_type##4, cl_##kernel_type##8, cl_##kernel_type##16, void,      \
    cl_##kernel_bitshift_type##2)                                              \
  F(kernel_type##4, cl_##kernel_type##4, 4, cl_##kernel_type,                  \
    cl_##kernel_type##8, cl_##kernel_type##16, void, void,                     \
    cl_##kernel_bitshift_type##4)                                              \
  F(kernel_type##8, cl_##kernel_type##8, 8, cl_##kernel_type,                  \
    cl_##kernel_type##16, void, void, void, cl_##kernel_bitshift_type##8)      \
  F(kernel_type##16, cl_##kernel_type##16, 16, cl_##kernel_type, void, void,   \
    void, void, cl_##kernel_bitshift_type##16)

// This macro factory lets us generate type traits in an easy manner
// The format for the passed macro function F is
// F(kernel_type,host_type,tensororder,base_type,vec2_type,vec4_type,vec8_type,vec16_type,bitshift_type)

#define BASE_CL_TYPE_FACTORY(F)                                                \
  VEC_TYPE(F, char, char)                                                      \
  VEC_TYPE(F, uchar, uchar)                                                    \
  VEC_TYPE(F, short, short)                                                    \
  VEC_TYPE(F, ushort, ushort)                                                  \
  VEC_TYPE(F, int, int)                                                        \
  VEC_TYPE(F, uint, uint)                                                      \
  VEC_TYPE(F, long, long)                                                      \
  VEC_TYPE(F, ulong, ulong)                                                    \
  VEC_TYPE(F, float, int)

// Need a way to check for the double data type, and the half data type!
#define CL_TYPE_FACTORY(F)                                                     \
  BASE_CL_TYPE_FACTORY(F)                                                      \
  VEC_TYPE(F, double, ulong)

namespace magnet {
namespace CL {
namespace detail {
template <class T> struct traits {
  static const bool isCLType = false;
};

#ifndef DOXYGEN_SHOULD_IGNORE_THIS
// We use the CL_TYPE_FACTORY to generate our type traits
#define TRAIT_FACTORY(cl_type, hosttype, tensororder, basetype, vec2type,      \
                      vec4type, vec8type, vec16type, bitshifttype)             \
  template <> struct traits<hosttype> {                                        \
    static const bool is_CL_type = true;                                       \
    static const int tensor_order = tensororder;                               \
    static inline std::string kernel_type() { return STRINGIFY(cl_type); }     \
    typedef vec2type vec2_type;                                                \
    typedef vec4type vec4_type;                                                \
    typedef vec8type vec8_type;                                                \
    typedef vec16type vec16_type;                                              \
    typedef bitshifttype bitshiftable_type;                                    \
  };

// Call the factory
CL_TYPE_FACTORY(TRAIT_FACTORY)

#undef TRAIT_FACTORY

#endif
} // namespace detail
} // namespace CL
} // namespace magnet

#undef CL_TYPE_FACTORY
#undef STRINGIFY
