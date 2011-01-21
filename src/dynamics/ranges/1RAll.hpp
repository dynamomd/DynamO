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

#ifndef C1RAll_H
#define C1RAll_H

#include "1range.hpp"
#include "../../base/is_base.hpp"
#include "../../base/is_simdata.hpp"

class CRAll: public CRange, public DYNAMO::SimBase_const
{
public:
  CRAll(const DYNAMO::SimData* SimDat):
    SimBase_const(SimDat,"CRAll",IC_red){}

  CRAll(const XMLNode&, const DYNAMO::SimData*);

  virtual CRange* Clone() const { return new CRAll(*this); };

  virtual bool isInRange(const Particle&) const
  { return true; }

  //The data output classes
  virtual void operator<<(const XMLNode&);

  virtual unsigned long size() const { return Sim->particleList.size(); }

  virtual iterator begin() const { return CRange::iterator(0, this); }

  virtual iterator end() const { return CRange::iterator(Sim->particleList.size(), this); }

  virtual unsigned long operator[](unsigned long i) const  
  { return i; }

  virtual unsigned long at(unsigned long i) const 
  { 
    if (i >= Sim->particleList.size())
      M_throw() << "Bad array access value in range.at()";

    return i;
  }

protected:

  virtual const unsigned long& getIteratorID(const unsigned long &i) const { return i; }

  virtual void outputXML(xml::XmlStream&) const;
};

#endif
