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
#include <magnet/gtk/numericEntry.hpp>

namespace magnet {
  namespace gtk {
    class ColorMapSelector: public ::Gtk::HBox
    {
    public:
      enum Mode_t { SEBASTIAN, HSV, MARCUS };

      ColorMapSelector():
	_min(0), _max(1)
      {
	m_refTreeModel = ::Gtk::ListStore::create(m_Columns);
	_comboBox.set_model(m_refTreeModel);
	
	buildEntry(SEBASTIAN, "Heat");
	buildEntry(HSV, "HSV");
	buildEntry(MARCUS, "Grayscale safe");
	
	_comboBox.pack_start(m_Columns.m_col_name, true);
	_comboBox.pack_start(m_Columns.m_col_icon, false);
	_comboBox.signal_changed().connect(sigc::mem_fun(*this, &ColorMapSelector::on_combobox_changed));
	_comboBox.set_active(0);
	_comboBox.show();

	show();

	_minValue.set_text("0");
	_minValue.set_width_chars(5);
	_minValue.show();
	_minValue.signal_changed().connect(sigc::mem_fun(*this, &ColorMapSelector::on_limits_changed));

	_maxValue.set_text("1.0");
	_maxValue.set_width_chars(5);
	_maxValue.show();
	_maxValue.signal_changed().connect(sigc::mem_fun(*this, &ColorMapSelector::on_limits_changed));

	Gtk::Label* label = Gtk::manage(new Gtk::Label("Range")); label->show();
	pack_start(*label, false, false, 5);
	pack_start(_minValue, false, false, 5);
	label = Gtk::manage(new Gtk::Label(":")); label->show();
	pack_start(*label, false, false, 2);
	pack_start(_maxValue, false, false, 5);

	label = Gtk::manage(new Gtk::Label("Scale")); label->show(); label->set_alignment(0.95, 0.5);
	pack_start(*label, true, true, 5);
	pack_start(_comboBox, false, false, 5);
      }

      void map(cl_uchar4& color, float val)
      {
	GLfloat floatcolor[4];
	map(floatcolor, val);
	for (size_t i(0); i < 4; ++i)
	  color.s[i] = 255 * floatcolor[i];
      }

      void map(float color[4], float val)
      {
	val = (val - _min) / (_max - _min);
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

      inline Mode_t getMode() const { return _mode; }

      ::sigc::signal<void> signal_changed() { return _signal_changed; }

      inline void setRange(GLfloat min, GLfloat max)
      {
	_minValue.set_text(boost::lexical_cast<std::string>(min));
	_maxValue.set_text(boost::lexical_cast<std::string>(max));
      }
	
    protected:
      ::Gtk::ComboBox _comboBox;
      ::Gtk::Entry _minValue;
      ::Gtk::Entry _maxValue;
      ::sigc::signal<void> _signal_changed;
      GLfloat _min;
      GLfloat _max;

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

      void on_limits_changed()
      {
	magnet::gtk::forceNumericEntry(_minValue);
	try { _min = boost::lexical_cast<GLfloat>(_minValue.get_text()); } catch(...) {}
	magnet::gtk::forceNumericEntry(_maxValue);
	try { _max = boost::lexical_cast<GLfloat>(_maxValue.get_text()); } catch(...) {}
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
