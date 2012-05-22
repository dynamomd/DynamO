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
#include <magnet/math/vector.hpp>

namespace coil {
  namespace glprimatives { 
    //a is the head, b the tail
    inline void drawArrow(Vector b, Vector a)
    {
      Vector arrowAxis = a - b;
      Vector headpoint = b + arrowAxis * 0.75;
      Vector headaxis = (arrowAxis ^ Vector(1,0,0));
      double headaxisnorm = headaxis.nrm();

      if (headaxisnorm == 0)
	{
	  headaxis = ((a-b) ^ Vector(0,0,1));
	  headaxisnorm = headaxis.nrm();
	}
      
      headaxis *= 0.15 * arrowAxis.nrm() / headaxisnorm;
      
      GLfloat A[]={a(0), a(1), a(2)},
	B[]={b(0), b(1), b(2)},
	  C[]={headpoint(0) + headaxis(0), headpoint(1) + headaxis(1), headpoint(2) + headaxis(2)},
	    D[]={headpoint(0) - headaxis(0), headpoint(1) - headaxis(1), headpoint(2) - headaxis(2)};
	    
	    glBegin (GL_LINES);
	    glVertex3fv (A);
	    glVertex3fv (B);  
	    glVertex3fv (A);
	    glVertex3fv (C);
	    glVertex3fv (A);
	    glVertex3fv (D);
	    glEnd();
    }
  }
}
