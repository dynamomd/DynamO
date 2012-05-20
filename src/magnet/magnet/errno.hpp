/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <errno.h>
#include <magnet/exception.hpp>

#define ERRNO_ENUM_FACTORY(F)			\
F(E2BIG, "Argument list too long.") \
F(EACCES, "Permission denied.") \
F(EADDRINUSE, "Address in use.") \
F(EADDRNOTAVAIL, "Address not available.") \
F(EAFNOSUPPORT, "Address family not supported.") \
F(EAGAIN, "Resource unavailable, try again") \
F(EALREADY, "Connection already in progress.") \
F(EBADF, "Bad file descriptor.") \
F(EBADMSG, "Bad message.") \
F(EBUSY, "Device or resource busy.") \
F(ECANCELED, "Operation canceled.") \
F(ECHILD, "No child processes.") \
F(ECONNABORTED, "Connection aborted.") \
F(ECONNREFUSED, "Connection refused.") \
F(ECONNRESET, "Connection reset.") \
F(EDEADLK, "Resource deadlock would occur.") \
F(EDESTADDRREQ, "Destination address required.") \
F(EDOM, "Mathematics argument out of domain of function.") \
F(EDQUOT, "Reserved.") \
F(EEXIST, "File exists.") \
F(EFAULT, "Bad address.") \
F(EFBIG, "File too large.") \
F(EHOSTUNREACH, "Host is unreachable.") \
F(EIDRM, "Identifier removed.") \
F(EILSEQ, "Illegal byte sequence.") \
F(EINPROGRESS, "Operation in progress.") \
F(EINTR, "Interrupted function.") \
F(EINVAL, "Invalid argument.") \
F(EIO, "I/O error.") \
F(EISCONN, "Socket is connected.") \
F(EISDIR, "Is a directory.") \
F(ELOOP, "Too many levels of symbolic links.") \
F(EMFILE, "File descriptor value too large.") \
F(EMLINK, "Too many links.") \
F(EMSGSIZE, "Message too large.") \
F(EMULTIHOP, "Reserved.") \
F(ENAMETOOLONG, "Filename too long.") \
F(ENETDOWN, "Network is down.") \
F(ENETRESET, "Connection aborted by network.") \
F(ENETUNREACH, "Network unreachable.") \
F(ENFILE, "Too many files open in system.") \
F(ENOBUFS, "No buffer space available.") \
F(ENODATA, " No message is available on the STREAM head read queue.") \
F(ENODEV, "No such device.") \
F(ENOENT, "No such file or directory.") \
F(ENOEXEC, "Executable file format error.") \
F(ENOLCK, "No locks available.") \
F(ENOLINK, "Reserved.") \
F(ENOMEM, "Not enough space.") \
F(ENOMSG, "No message of the desired type.") \
F(ENOPROTOOPT, "Protocol not available.") \
F(ENOSPC, "No space left on device.") \
F(ENOSR, " No STREAM resources.") \
F(ENOSTR, " Not a STREAM.") \
F(ENOSYS, "Function not supported.") \
F(ENOTCONN, "The socket is not connected.") \
F(ENOTDIR, "Not a directory.") \
F(ENOTEMPTY, "Directory not empty.") \
F(ENOTRECOVERABLE, "State not recoverable.") \
F(ENOTSOCK, "Not a socket.") \
F(ENOTSUP, "Not supported.") \
F(ENOTTY, "Inappropriate I/O control operation.") \
F(ENXIO, "No such device or address.") \
F(EOVERFLOW, "Value too large to be stored in data type.") \
F(EOWNERDEAD, "Previous owner died.") \
F(EPERM, "Operation not permitted.") \
F(EPIPE, "Broken pipe.") \
F(EPROTO, "Protocol error.") \
F(EPROTONOSUPPORT, "Protocol not supported.") \
F(EPROTOTYPE, "Protocol wrong type for socket.") \
F(ERANGE, "Result too large.") \
F(EROFS, "Read-only file system.") \
F(ESPIPE, "Invalid seek.") \
F(ESRCH, "No such process.") \
F(ESTALE, "Reserved.") \
F(ETIME, " Stream ioctl() timeout.") \
F(ETIMEDOUT, "Connection timed out.") \
F(ETXTBSY, "Text file busy.") \
F(EXDEV, "Cross-device link.") 

#define ENUM_FUNC(F,G) \
  case (F): return #F ": " G;

namespace magnet {
  namespace detail {
    /*! \brief Helper function to convert standard error codes in
     *  errno.h and specified in IEEE Std 1003.1-2001, to string
     *  representations.
     *
     * This is a C++ reentrant and thread safe function of the
     * standard strerror function.
     *
     * \param errnum The error number to convert to a string representation
     */
    inline const char* strerror_enum(int errnum)
    {
      switch (errnum)
	{

	  ERRNO_ENUM_FACTORY(ENUM_FUNC)
#undef ERRNO_ENUM_FACTORY
#undef ENUM_FUNC
	default:
	  M_throw() << "Unhandled errno! " << errnum;
	}
    }
  }

  /*! \brief A reentrant, and C++ form of the strerror C function.
   * \sa detail::strerror_enum()
   */
  inline std::string strerror(int errnum)
  {
    return detail::strerror_enum(errnum);
  }
}

