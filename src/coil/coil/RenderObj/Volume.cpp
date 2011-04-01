/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "Volume.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>

namespace coil {
  RVolume::RVolume(std::string name):
    RQuads(name)
  {}

  void 
  RVolume::initOpenCL() 
  {
    {
      float vertices[] = {0,0,0, 1,0,0, 1,1,0, 0,1,0,
			  0,0,1, 0,1,1, 1,1,1, 1,0,1};
      
      std::vector<float> VertexPos(vertices, vertices + sizeof(vertices) / sizeof(float));
      setGLPositions(VertexPos);
    }
    
    {
      int elements[] = {0,1,2,3, 2,1,7,6, 6,7,4,5, 5,4,0,3, 5,3,2,6, 1,0,4,7}; 
      std::vector<int> ElementData(elements, elements + sizeof(elements) / sizeof(int));
      setGLElements(ElementData);
    }
  }  

  void 
  RVolume::glRender()
  {
    RQuads::glRender();

  }

}
