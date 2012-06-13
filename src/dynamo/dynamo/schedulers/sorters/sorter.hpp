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
#include <dynamo/schedulers/sorters/event.hpp>
#include <dynamo/base.hpp>
#include <dynamo/eventtypes.hpp>

namespace magnet { namespace xml { class Node; } } 
namespace xml { class XmlStream; } 

namespace dynamo {
  /*! \brief Future Event Lists (FEL) sort the Particle Event Lists
      (PEL) to determine the next event to occur.

      Classes Derived from this base class provide a mechanism to sort
      \ref Event s. These events are first pre-sorted using a Particle
      Event List before being sorted by these classes.
   */

  class FEL: public dynamo::SimBase_const
  {
  public:
    FEL(const dynamo::Simulation* const& SD, const char *aName);

    virtual ~FEL() {}
    virtual size_t size()                              const = 0;
    virtual bool   empty()                             const = 0;
    virtual void   resize(const size_t&)                     = 0;
    virtual void   clear()                                   = 0;
    virtual void   init()                                    = 0;
    //A slient version of init
    virtual void   rebuild()                                 = 0;
    virtual void   stream(const double&)                       = 0;
    virtual void   push(const Event&, const size_t&)       = 0;
    virtual void   update(const size_t&)                     = 0;
    virtual size_t next_ID()                           const = 0;
    //virtual PELHeap& next_Data()                               = 0;
    //virtual const PELHeap& next_Data()                   const = 0;
    //virtual const PELHeap& operator[](const size_t&)     const = 0;
    //virtual PELHeap& operator[](const size_t&)                 = 0;
    virtual double   next_dt()                           const = 0;
    virtual EEventType next_type() const                     = 0;
    virtual unsigned long next_collCounter2() const          = 0;
    virtual size_t next_p2() const                           = 0;

    virtual void   sort()                                    = 0;
    virtual void   rescaleTimes(const double&)                 = 0;
    virtual void   clearPEL(const size_t&)                   = 0;
    virtual void   popNextPELEvent(const size_t&)            = 0;
    virtual void   popNextEvent()                            = 0;
    virtual bool nextPELEmpty() const                        = 0;

    //! Fetch the next event in the list, 
    virtual Event   copyNextEvent() const               = 0;

    static shared_ptr<FEL>
    getClass(const magnet::xml::Node&, const dynamo::Simulation*);

    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const FEL&);

  private:
    virtual void outputXML(magnet::xml::XmlStream&) const = 0;
  
  };
}
