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
#include <pangomm.h>

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

	virtual void init(size_t width, size_t height, size_t alpha_testing = 0)
	{
	  CairoSurface::init(width, height, alpha_testing);
	  _pango = Pango::Layout::create(_cairoContext);
	  Pango::FontDescription font("sans 12");
	  _pango->set_font_description(font);
	  
	  _pos[0] = 0.5 * _width;
	  _pos[1] = 0.5 * _height;
	}

	virtual void deinit()
	{
	  CairoSurface::deinit();
	  clear();
	  _pango.clear();
	}

	virtual void resize(size_t width, size_t height)
	{
	  if ((width == _width) && (height == _height))
	    return;

	  std::string olddata = _os.str();
	  init(width, height, _alpha_testing);
	  (*this) << olddata;
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

	void setPosition(double x, double y)
	{
	  if ((x == _pos[0]) && (y == _pos[1])) return;
	  _pos[0] = x;
	  _pos[1] = y;
	  _valid = false;
	}

      protected:
	std::ostringstream _os;
	bool _valid;
	Glib::RefPtr<Pango::Layout> _pango;
	double _pos[2];

	virtual void drawCommands()
	{
	  const double padding = 5;
	  _pango->set_text(_os.str());

	  //Fetch the box dimensions
	  _cairoContext->move_to(_pos[0], _pos[1]);
	  _pango->add_to_cairo_context(_cairoContext);
	  double topleft[2], bottomright[2];
	  _cairoContext->get_stroke_extents(topleft[0], topleft[1], bottomright[0], bottomright[1]);
	  topleft[0] += padding; topleft[1] += padding;
	  bottomright[0] += 3 * padding; bottomright[1] += 3 * padding;

	  _cairoContext->begin_new_path();	  
	  _cairoContext->set_source_rgba(0, 0, 0, 1.0);	 
	  _cairoContext->move_to(_pos[0] - 2 * padding, _pos[1]);
	  _cairoContext->line_to(_pos[0] + 2 * padding, _pos[1]);
	  _cairoContext->move_to(_pos[0], _pos[1] - 2 * padding);
	  _cairoContext->line_to(_pos[0], _pos[1] + 2 * padding);
	  _cairoContext->set_line_width(2.0);
	  _cairoContext->stroke();

	  //Make sure the box doesn't overlap the sides. The left hand
	  //side takes priority over the right
	  double dimensions[2] = {_width, _height};
	  for (size_t i(0); i < 2; ++i)
	    {
	      //right/bottom edge
	      double shift = std::min(0.0, dimensions[i] - bottomright[i]);
	      topleft[i] += shift;
	      bottomright[i] += shift;
	      //left/top edge
	      shift = std::max(-topleft[i], 0.0);
	      topleft[i] += shift;
	      bottomright[i] += shift;
	    }

	  //Background box
	  _cairoContext->begin_new_path();	  
	  _cairoContext->rectangle(topleft[0], topleft[1], 
				   bottomright[0] - topleft[0],
				   bottomright[1] - topleft[1]);

	  _cairoContext->set_source_rgba(0.5, 0.70588, 0.94118, 0.7);
	  _cairoContext->fill();

	  //Main text
	  _cairoContext->set_source_rgba(0, 0, 0, 1);
	  _cairoContext->move_to(topleft[0] + padding, topleft[1] + padding);
	  _pango->show_in_cairo_context(_cairoContext);
	}
      };
    }
  }
}
