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
#include <coil/clWindow.hpp>
#include <coil/filters/DOF.hpp>
#include <coil/filters/SSAO.hpp>
#include <coil/filters/bilateralBlur.hpp>
#include <coil/filters/filter.hpp>
#include <coil/filters/filterWrapper.hpp>
#include <coil/filters/multiply.hpp>
#include <magnet/exception.hpp>

namespace coil {
void Filter::populateComboBox(Gtk::ComboBox *filterSelectBox) {
  Glib::RefPtr<Gtk::ListStore> m_refTreeModel =
      Gtk::ListStore::create(getSelectColumnsInstance());
  filterSelectBox->set_model(m_refTreeModel);

#define COMBO_BOX_GEN_FUNC(enumeration, stringname, type)                      \
  {                                                                            \
    Gtk::TreeModel::Row row = *(m_refTreeModel->prepend());                    \
    row[getSelectColumnsInstance().m_col_id] = enumeration;                    \
    row[getSelectColumnsInstance().m_col_name] = stringname;                   \
  }

  FILTER_FACTORY(COMBO_BOX_GEN_FUNC)

#undef COMBO_BOX_GEN_FUNC

  // filterSelectBox->set_active(0);
  filterSelectBox->pack_start(getSelectColumnsInstance().m_col_name);
}

Filter *Filter::createFromID(size_t type_id) {

#define CREATE_GEN_FUNC(enumeration, stringname, type)                         \
  case enumeration:                                                            \
    return new type;

  switch (type_id) {
    FILTER_FACTORY(CREATE_GEN_FUNC)
  default:
    M_throw() << "Bad id passed";
  }
#undef CREATE_GEN_FUNC
}

std::string Filter::getName(size_t type_id) {

#define NAME_GEN_FUNC(enumeration, stringname, type)                           \
  case enumeration:                                                            \
    return stringname;

  switch (type_id) {
    FILTER_FACTORY(NAME_GEN_FUNC)
  default:
    M_throw() << "Bad id passed";
  }
#undef NAME_GEN_FUNC
}
} // namespace coil
