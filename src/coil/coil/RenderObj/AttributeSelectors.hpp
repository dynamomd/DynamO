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
#include <coil/RenderObj/DataSet.hpp>

namespace coil {
class AttributeSelector : public Gtk::VBox {
public:
  AttributeSelector(bool enableDataFiltering = true);

  void buildEntries(std::string name, DataSet &ds, size_t minComponents,
                    size_t maxComponents, int typeMask, size_t components,
                    int defaultMask = 0);

  struct ModelColumns : Gtk::TreeModelColumnRecord {
    ModelColumns() {
      add(m_name);
      add(m_ptr);
    }

    Gtk::TreeModelColumn<Glib::ustring> m_name;
    Gtk::TreeModelColumn<std::shared_ptr<Attribute>> m_ptr;
  };

  magnet::GL::Buffer<GLfloat> &getBuffer();

  virtual void bindAttribute(size_t attrnum, size_t divisor = 1) {
    if (singleValueMode())
      setConstantAttribute(attrnum);
    else
      getBuffer().attachToAttribute(attrnum, divisor);
  }

  virtual std::vector<GLfloat> getValue(size_t id) {
    Gtk::TreeModel::iterator iter = _comboBox.get_active();
    if (!iter)
      return std::vector<GLfloat>();

    std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
    if (!ptr)
      return std::vector<GLfloat>();

    std::vector<GLfloat> retval(ptr->components());

    for (size_t i(0); i < ptr->components(); ++i)
      retval[i] = (*ptr)[id * ptr->components() + i];

    return retval;
  }

  virtual std::vector<GLfloat> getMin();
  virtual std::vector<GLfloat> getMax();

  ModelColumns _modelColumns;
  Gtk::ComboBox _comboBox;
  Gtk::ComboBoxText _componentSelect;
  Gtk::Label _label;
  Gtk::Label _singleValueLabel;
  Glib::RefPtr<Gtk::ListStore> _model;
  Gtk::Entry _scalarvalues[4];
  Gtk::HBox _selectorRow;

protected:
  Attribute *_lastAttribute;
  size_t _lastAttributeDataCount;
  int _lastComponentSelected;
  magnet::GL::Buffer<GLfloat> _filteredData;
  magnet::GL::Context::ContextPtr _context;
  size_t _components;
  bool _enableDataFiltering;

  inline bool singleValueMode() {
    Gtk::TreeModel::iterator iter = _comboBox.get_active();
    if (!iter)
      return true;
    std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
    return !ptr;
  }

  void generateFilteredData(std::vector<GLfloat> &scalardata,
                            const std::shared_ptr<Attribute> &ptr, int mode);

  void setConstantAttribute(size_t attr);

  virtual void updateGui();
};
} // namespace coil
