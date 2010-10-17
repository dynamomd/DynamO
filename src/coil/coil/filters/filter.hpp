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

namespace coil 
{ 
  template<class T, bool> class magnetFilterWrapper; 
  class SSAOWrapper; 
  class SSAOBlurWrapper; 
}

typedef coil::magnetFilterWrapper<magnet::GL::laplacianFilter5,false> lap5x5;
typedef coil::magnetFilterWrapper<magnet::GL::laplacianFilter3A,false> lap3x3A;
typedef coil::magnetFilterWrapper<magnet::GL::laplacianFilter3B,false> lap3x3B;
typedef coil::magnetFilterWrapper<magnet::GL::LoGFilter,false> lapgauss9x9;
typedef coil::magnetFilterWrapper<magnet::GL::blurFilter,false> gauss5x5;
typedef coil::magnetFilterWrapper<magnet::GL::boxFilter,false> box5x5;

//Format is F(enumeration,stringname,type)
#define FILTER_FACTORY(F)						\
  F(0, "5x5 Laplacian", lap5x5)						\
  F(1, "3x3 Laplacian A", lap3x3A)					\
  F(2, "3x3 Laplacian B", lap3x3B)					\
  F(3, "9x9 Laplacian of Gaussian", lapgauss9x9)			\
  F(4, "5x5 Gaussian Blur", gauss5x5)					\
  F(5, "5x5 Box Filter", box5x5)					\
  F(6, "SSAO: Shadow Gen (Pass 1)", SSAOWrapper)			\
  F(7, "SSAO: Blur (Pass 2)", SSAOBlurWrapper) 

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
    virtual inline ~filter() {}

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
    
    inline virtual bool needsNormalDepth() = 0;

    //////////////Virtual members
    virtual size_t type_id() = 0;
    virtual bool isEditable() = 0;
    virtual void edit() {}
    virtual void invoke(GLuint colorTextureUnit, size_t width, size_t height) = 0;
    
  protected:
  };
}
