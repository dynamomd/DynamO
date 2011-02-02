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

#pragma once
#include <gtkmm.h>
#include <magnet/color/colormaps.hpp>

namespace coil {
  namespace Gtk {
    class ColorMapSelector: public ::Gtk::ComboBox
    {
    public:
      ColorMapSelector()
      {
	//Icon data
	size_t width = 100;
	size_t height = 20;

	m_refTreeModel = ::Gtk::ListStore::create(m_Columns);
	set_model(m_refTreeModel);

	::Gtk::TreeModel::Row row = *(m_refTreeModel->append());
	row[m_Columns.m_col_id] = 0;
	row[m_Columns.m_col_name] = "HSV";

	{
	  Glib::RefPtr<Gdk::Pixbuf> tmpPixBuf(Gdk::Pixbuf::create(Gdk::COLORSPACE_RGB, true,8,width,height));
	  guint8* pixels = tmpPixBuf->get_pixels();
	  
	  for (size_t column(0); column < width; ++column)
	    {
	      magnet::color::HSVtoRGB((cl_uchar4&)(pixels[4*column]), (float)(column)/ width);
	      //magnet::color::SebastiantoRGB((cl_uchar4&)(pixels[4*column]), (float)(column)/ width);
	      //magnet::color::MarcustoRGB((cl_uchar4&)(pixels[4*column]), (float)(column)/ width);
	      for (size_t row(1); row < height; ++row)
		(int32_t&)(pixels[4*(column + width * row)]) = (int32_t&)(pixels[4*column]);
	    }

	  row[m_Columns.m_col_icon] = Glib::RefPtr< ::Gtk::Image>(new ::Gtk::Image(tmpPixBuf));
	}

	pack_start(m_Columns.m_col_name);
	pack_start(m_Columns.m_col_icon);
      }

    private:

      struct ColorMapColumns : public ::Gtk::TreeModel::ColumnRecord
      {
	ColorMapColumns() { add(m_col_id); add(m_col_name); add(m_col_icon); }
	::Gtk::TreeModelColumn<int> m_col_id;
	::Gtk::TreeModelColumn<Glib::ustring> m_col_name;
	::Gtk::TreeModelColumn<Glib::RefPtr< ::Gtk::Image> > m_col_icon;
      };

      ColorMapColumns m_Columns;
      Glib::RefPtr< ::Gtk::ListStore> m_refTreeModel;
    };
  }
}
