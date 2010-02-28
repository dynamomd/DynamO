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
#ifndef OPVTK_H
#define OPVTK_H

#include "ticker.hpp"

class OPVTK: public OPTicker
{
 public:
  OPVTK(const DYNAMO::SimData*, const XMLNode&);

  virtual OutputPlugin *Clone() const
  { return new OPVTK(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void output(xmlw::XmlStream&);

  void operator<<(const XMLNode&);
  
 protected:
  unsigned long frameCount;

  CVector<size_t> nBins;
  Vector  binWidth;
  Vector  invBinWidth;
  
  std::vector<Iflt> mVsquared;
  std::vector<unsigned long> SampleCounter;
  std::vector<Vector  > Momentum;

  size_t getCellID(Vector );

  unsigned long imageCounter;
};

#endif
