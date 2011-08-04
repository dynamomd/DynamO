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
#include <gtkmm.h>
#include <magnet/color/sebastian.hpp>
#include <magnet/color/marcus.hpp>
#include <magnet/color/HSV.hpp>

namespace magnet {
  namespace gtk {
    class ColorMapSelector: public ::Gtk::HBox
    {
    protected:
      ::Gtk::ComboBox _comboBox;
    public:

      typedef enum 
      {
	SEBASTIAN,
	HSV,
	MARCUS
      } Mode_t;

      void map(cl_uchar4& color, float val)
      {
	GLfloat floatcolor[4];
	map(floatcolor, val);
	for (size_t i(0); i < 4; ++i)
	  color.s[i] = 255 * floatcolor[i];
      }

      void map(float color[4], float val)
      {
	switch (_mode)
	  {
	  case SEBASTIAN:
	    magnet::color::SebastiantoRGB(color, val);
	    break;
	  case HSV:
	    magnet::color::HSVtoRGB(color, val);
	    break;
	  case MARCUS:
	    magnet::color::MarcustoRGB(color, val);
	    break;
	  default:
	    M_throw() << "Unknown color map mode";
	  }
      }

      ColorMapSelector()
      {
	m_refTreeModel = ::Gtk::ListStore::create(m_Columns);
	_comboBox.set_model(m_refTreeModel);
	
	buildEntry(SEBASTIAN, "Sebastian (Matlab-esque)");
	buildEntry(HSV, "HSV");
	buildEntry(MARCUS, "Gnuplot inspired (Good for grayscale)");
	
	_comboBox.pack_start(m_Columns.m_col_name, true);
	_comboBox.pack_start(m_Columns.m_col_icon, false);

	_comboBox.signal_changed().connect(sigc::mem_fun(*this, &ColorMapSelector::on_combobox_changed));
	_comboBox.set_active(0);
	_comboBox.show();
	pack_start(_comboBox, false, false, 5);
	show();
      }

      inline Mode_t getMode() const { return _mode; }

      ::sigc::signal<void> signal_changed() { return _signal_changed; }

    private:
      ::sigc::signal<void> _signal_changed;

      void buildEntry(Mode_t newMode, std::string name)
      {
	//Icon data
	size_t width = 100;
	size_t height = 20;
	::Gtk::TreeModel::Row row = *(m_refTreeModel->append());
	row[m_Columns.m_col_id] = _mode = newMode;
	row[m_Columns.m_col_name] = name;
	row[m_Columns.m_col_icon]
	  = Glib::RefPtr<Gdk::Pixbuf>(Gdk::Pixbuf::create(Gdk::COLORSPACE_RGB, true, 8, width, height));
	
	guint8* pixels = (static_cast<Glib::RefPtr<Gdk::Pixbuf> >(row[m_Columns.m_col_icon]))->get_pixels();
	
	for (size_t column(0); column < width; ++column)
	  {
	    map((cl_uchar4&)(pixels[4*column]), (float)(column)/ width);
	    for (size_t row(1); row < height; ++row)
	      (int32_t&)(pixels[4*(column + width * row)]) = (int32_t&)(pixels[4*column]);
	  }
      }

      void on_combobox_changed()
      {
	::Gtk::TreeModel::iterator iter = _comboBox.get_active();
	if (iter) _mode = ((*iter)[m_Columns.m_col_id]);
	_signal_changed.emit();
      }

      struct ColorMapColumns : public ::Gtk::TreeModel::ColumnRecord
      {
	ColorMapColumns() { add(m_col_id); add(m_col_name); add(m_col_icon); }
	::Gtk::TreeModelColumn<Mode_t> m_col_id;
	::Gtk::TreeModelColumn<Glib::ustring> m_col_name;
	::Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > m_col_icon;
      };

      ColorMapColumns m_Columns;
      Glib::RefPtr< ::Gtk::ListStore> m_refTreeModel;
      Mode_t _mode;
    };
  }
}
