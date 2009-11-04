/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "msdOrientationalCorrelator.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../../dynamics/systems/sysTicker.hpp"

OPMSDOrientationalCorrelator::OPMSDOrientationalCorrelator(const DYNAMO::SimData* tmp,
				   const XMLNode& XML):
  OPTicker(tmp,"MSDOrientationalCorrelator"),
  length(50),
  currCorrLength(0),
  ticksTaken(0),
  notReady(true)
{
  operator<<(XML);
}

void
OPMSDOrientationalCorrelator::operator<<(const XMLNode& XML)
{
  try
  {
    if (XML.isAttributeSet("Length"))
    {
      length = boost::lexical_cast<size_t> (XML.getAttribute("Length"));
    }
  }
  catch (boost::bad_lexical_cast &)
  {
    D_throw() << "Failed a lexical cast in OPMSDCorrelator";
  }
}

void
OPMSDOrientationalCorrelator::initialise()
{
  if (dynamic_cast<const LNOrientation*>(&Sim->dynamics.getLiouvillean()) == NULL)
  {
    D_throw() << "Plugin requires species to define an orientation";
  }

  I_cout() << "The length of the MSD orientational correlator is " << length;

  historicalData.resize(Sim->lN, boost::circular_buffer<RUpair>(length));

  stepped_data_parallel.resize(length, Iflt(0.0));
  stepped_data_perpendicular.resize(length, Iflt(0.0));

  currCorrLength = 1.0;

  const std::vector<LNOrientation::rotData>& initial_rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getCompleteRotData());

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
  {
    historicalData[part.getID()].push_front(RUpair(part.getPosition(), initial_rdat[part.getID()].orientation));
  }
}

void
OPMSDOrientationalCorrelator::ticker()
{
  const std::vector<LNOrientation::rotData>& current_rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getCompleteRotData());
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
  {
    historicalData[part.getID()].push_front(RUpair(part.getPosition(), current_rdat[part.getID()].orientation));
  }

  if (notReady)
    {
      if (++currCorrLength != length)
      {
	return;
      }

      notReady = false;
    }

  accPass();
}

void
OPMSDOrientationalCorrelator::accPass()
{
  ++ticksTaken;

  Iflt longitudinal_projection(0.0);
  Vector displacement_term(0,0,0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
  {
    for (size_t step(1); step < length; ++step)
    {
      displacement_term = historicalData[part.getID()][step].first - historicalData[part.getID()][0].first;
      longitudinal_projection = (displacement_term | historicalData[part.getID()][0].second);

      stepped_data_parallel[step] += std::pow(longitudinal_projection, 2);
      stepped_data_perpendicular[step] += (displacement_term - (longitudinal_projection * historicalData[part.getID()][0].second)).nrm2();

    }
  }
}

void
OPMSDOrientationalCorrelator::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("MSDOrientationalCorrelator");

  Iflt dt = dynamic_cast<const CSTicker&> (*Sim->dynamics.getSystem("SystemTicker")).getPeriod() / Sim->dynamics.units().unitTime();

  XML << xmlw::tag("Component")
      << xmlw::attr("Type") << "Parallel"
      << xmlw::chardata();

  for (size_t step(0); step < length; ++step)
  {
    XML << dt * step << "\t"
	<< stepped_data_parallel[step] / (static_cast<Iflt>(ticksTaken) * static_cast<Iflt>(Sim->lN) * Sim->dynamics.units().unitArea())
	<< "\n";
  }

  XML << xmlw::endtag("Component")
      << xmlw::tag("Component")
      << xmlw::attr("Type") << "Perpendicular"
      << xmlw::chardata();

  for (size_t step(0); step < length; ++step)
  {
    XML << dt * step << "\t"
	<< stepped_data_perpendicular[step] / (static_cast<Iflt>(ticksTaken) * static_cast<Iflt>(Sim->lN) * Sim->dynamics.units().unitArea())
	<< "\n";
  }

  XML << xmlw::endtag("Component")
      << xmlw::endtag("MSDOrientationalCorrelator");
}
