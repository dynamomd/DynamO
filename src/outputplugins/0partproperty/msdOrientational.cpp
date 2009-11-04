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

#include "msdOrientational.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include <vector>


OPMSDOrientational::OPMSDOrientational(const DYNAMO::SimData* tmp, const XMLNode&):
  OutputPlugin(tmp,"MSDOrientational")
{}

OPMSDOrientational::~OPMSDOrientational()
{}

void
OPMSDOrientational::initialise()
{
  if (dynamic_cast<const LNOrientation*>(&Sim->dynamics.getLiouvillean()) == NULL)
  {
    D_throw() << "Plugin requires species to define an orientation";
  }

  initialConfiguration.clear();
  initialConfiguration.resize(Sim->lN);

  const std::vector<LNOrientation::rotData>& rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getCompleteRotData());

  for (size_t ID = 0; ID < Sim->lN; ++ID)
  {
    initialConfiguration[ID] = RUpair(Sim->vParticleList[ID].getPosition(), rdat[ID].orientation);
  }
}

void
OPMSDOrientational::output(xmlw::XmlStream &XML)
{
  std::pair<Iflt,Iflt> MSDOrientational(calculate());

  XML << xmlw::tag("MSDOrientational")

      << xmlw::tag("Perpendicular")
      << xmlw::attr("val") << MSDOrientational.first
      << xmlw::attr("diffusionCoeff")
      << MSDOrientational.first * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xmlw::endtag("Perpendicular")

      << xmlw::tag("Parallel")
      << xmlw::attr("val") << MSDOrientational.second
      << xmlw::attr("diffusionCoeff")
      << MSDOrientational.second * Sim->dynamics.units().unitTime() / Sim->dSysTime
      << xmlw::endtag("Parallel")

      << xmlw::endtag("MSDOrientational");
}

std::pair<Iflt, Iflt>
OPMSDOrientational::calculate() const
{
  //Required to get the correct results
  Sim->dynamics.getLiouvillean().updateAllParticles();

  Iflt acc_perp(0.0), acc_parallel(0.0), longitudinal_projection(0.0);

  Vector displacement_term(0,0,0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
  {
    displacement_term = part.getPosition() - initialConfiguration[part.getID()].first;
    longitudinal_projection = (displacement_term | initialConfiguration[part.getID()].second);

    acc_perp += (displacement_term - (longitudinal_projection * initialConfiguration[part.getID()].second)).nrm2();
    acc_parallel += std::pow(longitudinal_projection, 2);
  }

  // In the N-dimensional case, the parallel component is 1-dimensonal and the perpendicular (N-1)
  acc_perp /= (initialConfiguration.size() * 2.0 * (NDIM - 1) * Sim->dynamics.units().unitArea());
  acc_parallel /= (initialConfiguration.size() * 2.0 * Sim->dynamics.units().unitArea());

  return std::pair<Iflt,Iflt> (acc_perp, acc_parallel);
}
