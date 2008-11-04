/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#ifndef COPVTK_H
#define COPVTK_H

#include "ticker.hpp"

class COPVTK: public COPTicker
{
 public:
  COPVTK(const DYNAMO::SimData*, const XMLNode&);

  virtual COutputPlugin *Clone() const
  { return new COPVTK(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void output(xmlw::XmlStream&);

  void operator<<(const XMLNode&);
  
 protected:
  unsigned long frameCount;

  CVector<size_t> nBins;
  CVector<> binWidth;
  CVector<> invBinWidth;
  
  std::vector<Iflt> mVsquared;
  std::vector<unsigned long> SampleCounter;
  std::vector<CVector<> > Momentum;

  size_t getCellID(CVector<>);

  unsigned long imageCounter;
};

#endif
