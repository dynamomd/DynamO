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
#pragma once
#include <magnet/GL/objects/cairo.hpp>
#include <sstream>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A quad textured with text and an optional background
          box.
       */
      class TextSurface: public CairoSurface
      {
      public:
	TextSurface(): _valid(false) {}
	
	/*! \brief Output operator for the text box.
	 */
	template<class T>
	inline TextSurface& operator<<(const T& value) 
	{
	  _os << value;
	  _valid = false;
	  return *this;
	}

	virtual void deinit() 
	{
	  CairoSurface::deinit();
	  clear();
	}

	virtual void resize(size_t width, size_t height)
	{
	  if ((width == _width) && (height == _height))
	    return;

	  std::string olddata = _os.str();
	  init(width, height, _alpha_testing);	  
	  _os << olddata;
	}


	inline void clear() 
	{ 
	  _os.str(""); 
	  _valid = false;
	}
	
	inline void 
	glRender(const GLMatrix& projection = GLMatrix::identity(),
		 const GLMatrix& modelview = GLMatrix::identity())
	{
	  if (_os.str().empty()) return;

	  if (!_valid) { redraw(); _valid = true; }

	  CairoSurface::glRender(projection, modelview);
	}

      protected:
	std::ostringstream _os;
	bool _valid;

	virtual void drawCommands()
	{
	  _cairoContext->scale(_width, _height);
	  _cairoContext->set_font_size(0.01);
	  _cairoContext->set_source_rgba(1.0, 0.1, 0.1, 1);
	  _cairoContext->move_to(0.5, 0.5);
	  std::string ostr = _os.str();
	  _cairoContext->show_text(ostr.c_str());
	}
      };
    }
  }
}
