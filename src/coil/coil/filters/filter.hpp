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

#pragma once

#include <gtkmm.h>

#include <magnet/GL/laplacianFilter.hpp>
#include <magnet/GL/blur.hpp>
#include <magnet/GL/DoF.hpp>

namespace coil { template<class T> class magnetFilterWrapper; }

//Format is F(enumeration,stringname,type)
#define FILTER_FACTORY(F)						\
  F(0, "5x5 Laplacian", magnetFilterWrapper<magnet::GL::laplacianFilter5>) \
  F(1, "3x3 Laplacian A", magnetFilterWrapper<magnet::GL::laplacianFilter3A>) \
  F(2, "3x3 Laplacian B", magnetFilterWrapper<magnet::GL::laplacianFilter3B>) \
  F(3, "9x9 Laplacian of Gaussian", magnetFilterWrapper<magnet::GL::LoGFilter>)	\
  F(4, "5x5 Gaussian Blur", magnetFilterWrapper<magnet::GL::blurFilter>) \
  F(5, "5x5 Box Filter", magnetFilterWrapper<magnet::GL::boxFilter>) \
  F(6, "Depth of Field", magnetFilterWrapper<magnet::GL::DoF>) 


namespace coil 
{
  namespace detail 
  {

    ////////Generate a type trait for the filter enumerations (type id's)
    template<class T> struct filterEnum {};

#define ENUM_GEN_FUNC(enumeration,stringname,type) \
    template<> struct filterEnum<type> { static const size_t val = enumeration; };
    
    FILTER_FACTORY(ENUM_GEN_FUNC)

#undef ENUM_GEN_FUNC
  }

  class filter
  {
  public:
    //////////////Static members
    static void populateComboBox(Gtk::ComboBox*);
    
    static filter* createFromID(size_t type_id);

    static std::string getName(size_t type_id);

    struct filterSelectColumns : public Gtk::TreeModel::ColumnRecord
    {
      filterSelectColumns() { add(m_col_name); add(m_col_id); }
      Gtk::TreeModelColumn<int> m_col_id;
      Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    };

    inline static filterSelectColumns& getSelectColumnsInstance()
    { 
      static filterSelectColumns vals;
      return vals;
    }

    //////////////Virtual members
    virtual size_t type_id() = 0;
    virtual bool isEditable() = 0;
    virtual void edit() {}
    virtual void invoke(GLuint colorTextureUnit, GLuint depthTextureUnit, size_t width, size_t height) = 0;
    
  protected:
  };
}
