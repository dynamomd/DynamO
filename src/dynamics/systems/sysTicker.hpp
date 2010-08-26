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

#ifndef CSTicker_HPP
#define CSTicker_HPP

#include "system.hpp"

class CSTicker: public System
{
public:
  CSTicker(DYNAMO::SimData*, Iflt, std::string);
  
  virtual System* Clone() const { return new CSTicker(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&) {}

  void setdt(Iflt);

  void increasedt(Iflt);

  void setTickerPeriod(const Iflt&);

  const Iflt& getPeriod() const { return period; }
protected:
  virtual void outputXML(xmlw::XmlStream&) const {}

  Iflt period;
};

#endif
