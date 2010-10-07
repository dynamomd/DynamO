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

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include <magnet/exception.hpp>
#include <magnet/CL/detail/extension_wrangler.hpp>

namespace magnet {
  namespace CL {
    namespace detail {
    
      /* This is a CRTP base class that builds kernels into functors
       *
       * It requires that the type that inherits it, specifies its own
       * type in the template parameter (T) and defines a static member
       * function T::kernelSrc().
       */
      template<class T>
      class functor {
      protected:
	inline void build(::cl::CommandQueue queue,cl::Context context, std::string buildFlags)
	{
	  _queue = queue;
	  _context = context;

	 cl::Program::Sources sources;

	  std::string extensions;

	  //Load the double extension if it is present
	  if (magnet::CL::detail::detectExtension(context,"cl_khr_fp64"))
	    extensions += "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
	  else if (magnet::CL::detail::detectExtension(context,"cl_amd_fp64"))
	    extensions += "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n";

	  sources.push_back(std::make_pair(extensions.c_str(),
					   extensions.size()));

	  //Now load the kernel source
	  std::string kernelSrc = format_code(T::kernelSource());	
	  sources.push_back(std::make_pair(kernelSrc.c_str(),
					   kernelSrc.size()));
	
	  //Attempt to build the source
	  _program =cl::Program(_context, sources);

	  std::vector<cl::Device> devices = _context.getInfo<CL_CONTEXT_DEVICES>();

	  try {
	    _program.build(devices, buildFlags.c_str());
	  } catch(::cl::Error& err) {
	    for(size_t dev = 0; dev < devices.size(); ++dev) {
	      std::string msg 
		= format_code(_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[dev]));
	    
	      if (!msg.size())
		continue;
	    
	      M_throw() << "Compilation failed for device " 
			<< devices[dev].getInfo<CL_DEVICE_NAME>()
			<< "\nBuild Log:\n" << msg;
	    }
	  }
	}

	::cl::Program _program;
	::cl::CommandQueue _queue;
	::cl::Context _context;

	/*! \brief Can search and replace elements in a std::string. */
	inline std::string 
	format_code(std::string in)
	{
	  return search_replace(in,";",";\n");
	}

	inline std::string 
	search_replace(std::string in, const std::string& from, const std::string& to)
	{
	  if (!in.empty())
	    {
	      std::string::size_type toLen = to.length();
	      std::string::size_type frLen = from.length();
	      std::string::size_type loc = 0;
	    
	      while (std::string::npos != (loc = in.find(from, loc)))
		{
		  in.replace(loc, frLen, to);
		  loc += toLen;
		
		  if (loc >= in.length())
		    break;
		}
	    }
	  return in;
	}

      };
    }
  }
}
