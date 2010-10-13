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

#include <coil/filters/filter.hpp>


namespace coil 
{
//Format is F(enum,stringname,type)
//#define FilterFactory(F) \
//  F()

  void 
  filter::populateComboBox(Gtk::ComboBox* filterSelectBox)
  {
    struct filterSelectColumns : public Gtk::TreeModel::ColumnRecord
    {
      filterSelectColumns() { add(m_col_name); add(m_col_id); }
      Gtk::TreeModelColumn<int> m_col_id;
      Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    };
    
    filterSelectColumns vals;
    Glib::RefPtr<Gtk::ListStore> m_refTreeModel = Gtk::ListStore::create(vals);
    filterSelectBox->set_model(m_refTreeModel);
    
    {
      Gtk::TreeModel::Row row = *(m_refTreeModel->prepend());
      row[vals.m_col_id] = 1;
      row[vals.m_col_name] = "Testing the name";
    }
    
    filterSelectBox->pack_start(vals.m_col_name);
  }
}
