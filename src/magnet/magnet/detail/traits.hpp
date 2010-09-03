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

#include <CL/cl.hpp>

#define STRINGIFY(A) #A

//This macro factory lets us generate type traits in a quick manner
//The format is (kernel_type, host_type, tensor_order, basetype)
#define CL_TYPE_FACTORY(F)					   \
  F(char,  cl_char,  1, char)					   \
  F(char2, cl_char2, 2, char)					   \
  F(char4, cl_char4, 4, char)					   \
  F(char8, cl_char8, 8, char)					   \
  F(char16,cl_char16,16,char)					   \
  F(uchar,cl_uchar,1,uchar)					   \
  F(uchar2,cl_uchar2,2,uchar)					   \
  F(uchar4,cl_uchar4,4,uchar)					   \
  F(uchar8,cl_uchar8,8,uchar)					   \
  F(uchar16,cl_uchar16,16,uchar)				   \
  F(short,cl_short,1,short)					   \
  F(short2,cl_short2,2,short)					   \
  F(short4,cl_short4,4,short)					   \
  F(short8,cl_short8,8,short)					   \
  F(short16,cl_short16,16,short)				   \
  F(ushort,cl_ushort,1,ushort)			    		   \
  F(ushort2,cl_ushort2,2,ushort)				   \
  F(ushort4,cl_ushort4,4,ushort)				   \
  F(ushort8,cl_ushort8,8,ushort)				   \
  F(ushort16,cl_ushort16,16,ushort)				   \
  F(int,cl_int,1,int)						   \
  F(int2,cl_int2,2,int)						   \
  F(int4,cl_int4,4,int)						   \
  F(int8,cl_int8,8,int)						   \
  F(int16,cl_int16,16,int)					   \
  F(uint,cl_uint,1,uint)					   \
  F(uint2,cl_uint2,2,uint)					   \
  F(uint4,cl_uint4,4,uint)					   \
  F(uint8,cl_uint8,8,uint)					   \
  F(uint16,cl_uint16,16,uint)					   \
  F(long,cl_long,1,long)					   \
  F(long2,cl_long2,2,long)					   \
  F(long4,cl_long4,4,long)					   \
  F(long8,cl_long8,8,long)					   \
  F(long16,cl_long16,16,long)					   \
  F(ulong,cl_ulong,1,ulong)					   \
  F(ulong2,cl_ulong2,2,ulong)					   \
  F(ulong4,cl_ulong4,4,ulong)					   \
  F(ulong8,cl_ulong8,8,ulong)					   \
  F(ulong16,cl_ulong16,16,ulong)				   \
  F(float,cl_float,1,float)					   \
  F(float2,cl_float2,2,float)					   \
  F(float4,cl_float4,4,float)					   \
  F(float8,cl_float8,8,float)					   \
  F(float16,cl_float16,16,float)				   \
  F(double,cl_double,1,double)					   \
  F(double2,cl_double2,2,double)				   \
  F(double4,cl_double4,4,double)				   \
  F(double8,cl_double8,8,double)				   \
  F(double16,cl_double16,16,double)				   \

//  F(half,cl_half,1,half)					   
//  F(half2,cl_half2,2,half)					   
//  F(half4,cl_half4,4,half)					   
//  F(half8,cl_half8,8,half)					   
//  F(half16,cl_half16,16,half)					   

//  F(bool,magnet::_traits_detail::no_host_type,1,bool)
//  F(size_t,magnet::_traits_detail::no_host_type,1,size_t)
//  F(ptrdiff_t,magnet::_traits_detail::no_host_type,1,ptrdiff_t)
//  F(intptr_t,magnet::_traits_detail::no_host_type,1,intptr_t)
//  F(uintptr_t,magnet::_traits_detail::no_host_type,1,uintptr_t)


namespace magnet {
  namespace _traits_detail {
    
    struct no_host_type {};    
    template<typename T>
    struct isvoid {
      static const bool eval = false;
    };
    
    template<>
    struct isvoid<void> {
      static const bool eval = true;
    };
  }
  
  template<class T>
  struct traits
  {
    static const bool isCLType = false;
  };
   
  //Now we first create typedefs of the cl types, as some
  //implementations just use macro definitions of the types, but we
  //want traits!
#define TYPE_FACTORY(cl_type,host_type,tensor_order,base_type) \
  typedef host_type CL_##cl_type;

  //Call the factory
  CL_TYPE_FACTORY(TYPE_FACTORY)

#undef TYPE_FACTORY

  //Now that we have our own CL_* types, we can use traits!
#define TRAIT_FACTORY(cl_type,hosttype,tensororder,basetype) \
  template<> struct traits<CL_##cl_type>			\
  {								\
    static const bool is_CL_type = true;			\
    static const int  tensor_order = tensororder;		\
    static const std::string kernel_type;			\
    typedef CL_##cl_type base_type;				\
  };
  
  //Call the factory
  CL_TYPE_FACTORY(TRAIT_FACTORY)

#undef TRAIT_FACTORY
}
