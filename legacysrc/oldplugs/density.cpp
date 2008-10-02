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

#include "density.hpp"

COPDensity::COPDensity(const std::vector<CParticle> &pList, const CDynamics * const dyn):
  COutputPlugin(pList,dyn),
  binwidth(1.0/bincount),
  samplecount(0)
{
  std::cout << "COPDensity: Loaded\n";

  for (int i = 0; i < bincount; i++)
    bin[i] = CVector<long>(0);
}

COPDensity::~COPDensity()
{
  std::cout << "COPDensity: Unloaded\n";
}

void 
COPDensity::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{}

void
COPDensity::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("density");

  for (int dim = 0; dim < NDIM; dim++)
    { 
      char *name;
      asprintf(&name, "dim%d", dim);
      XML << xmlw::tag(name)
	  << xmlw::tag("columns")
	  << xmlw::attr("x") << "r" 
	  << xmlw::attr("y") << "f" 
	  << xmlw::endtag("columns")
	  << xmlw::tag("data") << xmlw::chardata();
      
      for (int i = 0; i < bincount; i++)
	  XML << (((Iflt) i) + 0.5) * binwidth - 0.5 << " " << ((Iflt) (bin[i][dim]))/(binwidth*samplecount*((Iflt) particleList.size())) << "\n";
      
      XML << xmlw::endtag("data")
	  << xmlw::endtag(name); 
    }
}

void
COPDensity::periodicOutput()
{
  CVector<> pos;
  
  samplecount++;

  for (std::vector<CParticle>::const_iterator iPtr = particleList.begin(); iPtr != particleList.end(); iPtr++)
    {
      pos = iPtr->getPosition();

      dynamics->setPBC(pos);

      for (int dim = 0; dim < NDIM; dim++)
	{
	  int i = (int) ((pos[dim]+0.5)/binwidth);  
	  if (i >= bincount || i < 0)
	    throw CException() << "i out of range in density histogram, dim " 
			       << dim << ", val " << i << ", pos val " << pos[dim]*binwidth;
	  bin[i][dim]++;
	}
    }
}
