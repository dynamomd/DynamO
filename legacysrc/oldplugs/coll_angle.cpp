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

#include "coll_angle.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../dynamics/include.hpp"

#define binWidth 0.005

COPCollAngle::COPCollAngle(DYNAMO::SimData* tmp):
  COutputPlugin(tmp, "CollisionAngle")
{
  theta = C1DHistogram(binWidth);
  for (int i = 0; i < NDIM; i++)
    rhat[i] = C1DHistogram(binWidth);
}

COPCollAngle::~COPCollAngle()
{}

void 
COPCollAngle::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  theta.addVal(-preColl.rvdot/preColl.v12.length());

  CVector<> CVrhat = preColl.r12.unitVector();
  for (int i = 0; i < NDIM; i++)
    (rhat[i]).addVal(CVrhat[i]);
}

void
COPCollAngle::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("Collision_angle")
      << xmlw::tag("Theta")
      << xmlw::tag("Columns")
      << xmlw::attr("x") << "Theta"
      << xmlw::attr("y") << "f"
      << xmlw::endtag("Columns");
  
      theta.outputHistogram(XML, 1.0);

      XML << xmlw::endtag("Theta"); 

  for (int dim = 0; dim < NDIM; dim++)
    { 
      char dimnom = 'x' + dim;
      char name[2] = {dimnom, '\0'};
      XML << xmlw::tag(name)
	  << xmlw::tag("Columns")
	  << xmlw::attr("x") << "\\hat{r}_" <<  dimnom
	  << xmlw::attr("y") << "f"
	  << xmlw::endtag("Columns");

      rhat[dim].outputHistogram(XML, 1.0);

      XML << xmlw::endtag(name); 
    }

  XML << xmlw::endtag("Collision_angle");
}
