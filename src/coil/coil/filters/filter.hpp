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

//Format is F(enumeration,stringname,type)
#define FILTER_FACTORY(F) \
  F(0, "5x5 Laplacian", Laplacian5x5Filter)


namespace coil 
{
  ////////Generate prototypes for all of the filter classes
#define PROTOTYPE_GEN_FUNC(enumeration,stringname,type) \
  class type;

  FILTER_FACTORY(PROTOTYPE_GEN_FUNC)

#undef PROTOTYPE_GEN_FUNC

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
    static void populateComboBox(Gtk::ComboBox*);
    
    static filter* createFromID(size_t type_id);

    virtual size_t type_id() = 0;

    static std::string getName(size_t type_id);

  protected:
  };
}
