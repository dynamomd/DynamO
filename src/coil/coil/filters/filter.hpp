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

#include <magnet/GL/camera.hpp>
#include <magnet/GL/shader/laplacianFilter.hpp>
#include <magnet/GL/shader/blur.hpp>
#include <magnet/GL/shader/hipass.hpp>
#include <gtkmm.h>

namespace coil 
{ 
  template<class T> class magnetFilterWrapper; 
  class SSAOWrapper; 
  class BilateralBlurWrapper;
  class MultiplyFilter;
  class DOFFilter;
  class FlushToOriginal;
}

typedef coil::magnetFilterWrapper<magnet::GL::shader::Laplacian5x5> lap5x5;
typedef coil::magnetFilterWrapper<magnet::GL::shader::Laplacian3x3A> lap3x3A;
typedef coil::magnetFilterWrapper<magnet::GL::shader::Laplacian3x3B> lap3x3B;
typedef coil::magnetFilterWrapper<magnet::GL::shader::LoG9x9> lapgauss9x9;
typedef coil::magnetFilterWrapper<magnet::GL::shader::Gaussian5x5Blur> gauss5x5;
typedef coil::magnetFilterWrapper<magnet::GL::shader::Box5x5Blur> box5x5;
typedef coil::magnetFilterWrapper<magnet::GL::shader::HiPass3x3> HiPass3x3;
typedef coil::magnetFilterWrapper<magnet::GL::shader::HiPass5x5> HiPass5x5;

//Format is F(enumeration,stringname,type)
#define FILTER_FACTORY(F)						\
  F(0, "3x3 Hi-Pass Filter", HiPass3x3)					\
  F(1, "5x5 Hi-Pass Filter", HiPass5x5)					\
  F(2, "5x5 Laplacian", lap5x5)						\
  F(3, "3x3 Laplacian A", lap3x3A)					\
  F(4, "3x3 Laplacian B", lap3x3B)					\
  F(5, "9x9 Laplacian of Gaussian", lapgauss9x9)			\
  F(6, "5x5 Gaussian Blur", gauss5x5)					\
  F(7, "5x5 Box Filter", box5x5)					\
  F(8, "5x5 Gaussian Bilateral Blur", BilateralBlurWrapper)		\
  F(9, "SSAO: Shadow (After, use a bilateral blur, then multiply)", SSAOWrapper)	\
  F(10, "Multiply with Original Image", MultiplyFilter)			\
  F(11, "DOF Filter (Make a fully blurred image first)", DOFFilter)      \
  F(12, "Flush Current (Do after a SSAO/DOF filter to allow further processing)", FlushToOriginal) 

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

  class Filter
  {
  public:
    inline Filter ():_active(true) {}
    virtual inline ~Filter() {}

    //////////////Static members
    static void populateComboBox(Gtk::ComboBox*);
    
    static Filter* createFromID(size_t type_id);

    static std::string getName(size_t type_id);

    struct FilterSelectColumns : public Gtk::TreeModel::ColumnRecord
    {
      FilterSelectColumns() { add(m_col_name); add(m_col_id); }
      Gtk::TreeModelColumn<int> m_col_id;
      Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    };

    inline static FilterSelectColumns& getSelectColumnsInstance()
    {
      static FilterSelectColumns vals;
      return vals;
    }
    
    //////////////Virtual members
    virtual size_t type_id() = 0;
    virtual void showControls(Gtk::ScrolledWindow*) {}
    virtual void invoke(GLint colorTextureUnit, size_t width, size_t height, const magnet::GL::Camera& vp) = 0;
    
    inline bool getActive() const { return _active; }
    inline void setActive(bool nv) { _active = nv; }

  protected:
    bool _active;
  };

  class FlushToOriginal: public Filter
  {
  public:
    inline virtual size_t type_id() { return detail::filterEnum<FlushToOriginal>::val; }
    inline virtual void invoke(GLint, size_t, size_t, const magnet::GL::Camera&) {}
  protected:
  };
}
