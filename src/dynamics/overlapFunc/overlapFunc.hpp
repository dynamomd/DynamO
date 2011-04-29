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
/*! \file overlapFunc.hpp
 *  \brief Just a documentation header for the OverlapFunctions Namespace
 */

namespace dynamo
{
  /*! \brief Overlap functions between two specified shapes.
   *
   * Specifically these functions are used by CLocal event objects to
   * determine if they are within a given cell.
   *
   * As the cells are cubic there will probably be alot of cube vs *
   */
  namespace OverlapFunctions;
};
