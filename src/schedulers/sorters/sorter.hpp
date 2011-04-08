/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include "../../base/is_base.hpp"

class CSSorter: public DYNAMO::SimBase_const
{
public:
  CSSorter(const DYNAMO::SimData* const& SD, const char *aName);

  virtual ~CSSorter() {}
  virtual size_t size()                              const = 0;
  virtual bool   empty()                             const = 0;
  virtual void   resize(const size_t&)                     = 0;
  virtual void   clear()                                   = 0;
  virtual void   init()                                    = 0;
  //A slient version of init
  virtual void   rebuild()                                 = 0;
  virtual void   stream(const double&)                       = 0;
  virtual void   push(const intPart&, const size_t&)       = 0;
  virtual void   update(const size_t&)                     = 0;
  virtual size_t next_ID()                           const = 0;
  //virtual pList& next_Data()                               = 0;
  //virtual const pList& next_Data()                   const = 0;
  //virtual const pList& operator[](const size_t&)     const = 0;
  //virtual pList& operator[](const size_t&)                 = 0;
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
  virtual intPart   copyNextEvent() const               = 0;
  virtual CSSorter* Clone()                          const = 0;

  static CSSorter* getClass(const magnet::xml::Node&, const DYNAMO::SimData*);

  friend xml::XmlStream& operator<<(xml::XmlStream&, const CSSorter&);

private:
  virtual void outputXML(xml::XmlStream&) const = 0;
  
};
