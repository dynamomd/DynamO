/*  dynamo:- Event driven molecular dynamics simulator
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

#include "chatteringCorrelator.hpp"
#include "../../dynamics/2particleEventData.hpp"
#include "../../dynamics/include.hpp"
#include <magnet/xmlwriter.hpp>

OPChatteringCorrelator::OPChatteringCorrelator(const dynamo::SimData* t1,
					       const magnet::xml::Node& XML):
  OP2PP(t1,"ChatteringCorrelator") {}

void
OPChatteringCorrelator::initialise()
{
  //Set the history size
  chatterTracker.resize(Sim->N, std::pair<double,double>(0,0));

  //Histogram in mean free times
  hist = C1DWeightHistogram(1);
}

void
OPChatteringCorrelator::A2ParticleChange(const PairEventData& PDat)
{
  size_t ID1 = PDat.particle1_.getParticle().getID();
  size_t ID2 = PDat.particle2_.getParticle().getID();


  if (chatterTracker[ID1].first == ID2)
  {
    chatterTracker[ID1].second++;
  }
  else
  {
    if(chatterTracker[ID1].second != 0)
    {
      hist.addVal(chatterTracker[ID1].second, chatterTracker[ID1].second);
    }

    chatterTracker[ID1].first = ID2;
    chatterTracker[ID1].second = 1;

  }

  if (chatterTracker[ID2].first == ID1)
  {
    chatterTracker[ID2].second++;
  }
  else
  {
    if(chatterTracker[ID2].second != 0)
    {
      hist.addVal(chatterTracker[ID2].second, chatterTracker[ID2].second);
    }

    chatterTracker[ID2].first = ID1;
    chatterTracker[ID2].second = 1;

  }
}

void
OPChatteringCorrelator::output(magnet::xml::XmlStream &XML)
{
  XML << magnet::xml::tag("ChatteringCorrelator");

  hist.outputHistogram(XML, 1.0);

  XML << magnet::xml::endtag("ChatteringCorrelator");
}
