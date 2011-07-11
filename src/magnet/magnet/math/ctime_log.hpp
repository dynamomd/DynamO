/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

namespace magnet {
  namespace math {
    namespace detail {
      //! \brief A worker for the \ref ctime_log class.
      template <int val, int base>
      struct ctime_log_worker
      {
	static const int value = ctime_log_worker<val / base, base>::value + 1;
      };
      
      template <int base>
      struct ctime_log_worker<0, base>
      {
	static const int value = 0;
      };
    }

    /*! \brief A compile time function to calculate the log of an
     * integer.
     *
     * This function actually returns \f${\textrm
     * floor}\left(\log_base(val)\right)\f$ as only integer
     * mathematics can be performed using templates.
     */
    template <int val, int base>
    struct ctime_log
    {
      static const int value = detail::ctime_log_worker<val / base, base>::value;
    };
    
    /*! \brief A specialization to produce an error for the log of zero.
     */
    template <int base> struct ctime_log<0, base> {};
  }
}
