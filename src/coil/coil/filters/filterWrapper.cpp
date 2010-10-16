/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <coil/clWindow.hpp>
#include "filterWrapper.hpp"

extern const char _binary_src_coil_coil_filters_SSAO_gladexml_start[];
extern const char _binary_src_coil_coil_filters_SSAO_gladexml_end[];

namespace coil 
{
  SSAOWrapper::SSAOWrapper()
  { 
    _filter.build(); 

    Glib::ustring glade_data
      (reinterpret_cast<const char *>(_binary_src_coil_coil_filters_SSAO_gladexml_start), 
       _binary_src_coil_coil_filters_SSAO_gladexml_end
       -_binary_src_coil_coil_filters_SSAO_gladexml_start);
    
    _refXml = Gtk::Builder::create_from_string(glade_data);
  }

  SSAOWrapper::~SSAOWrapper()
  { 
    std::cerr << "!!!!!!!!!!!!!DESTROYING THE FILTER!!!!!!!!!!!!";
    Gtk::Window* win;
    _refXml->get_widget("SSAOeditWindow", win);
    win->hide();
  }

  void SSAOWrapper::edit()
  {
    Gtk::Window* win;
    _refXml->get_widget("SSAOeditWindow", win);
    win->show();
  }
}
