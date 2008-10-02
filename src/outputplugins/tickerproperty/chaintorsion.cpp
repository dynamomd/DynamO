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

#include "chaintorsion.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../dynamics/ranges/1range.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/topology/include.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

COPCTorsion::COPCTorsion(const DYNAMO::SimData* tmp):
  COPTicker(tmp,"Torsion")
{}

void 
COPCTorsion::initialise()
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(CTCdata(dynamic_cast<const CTChain*>(plugPtr.get_ptr()), 0.005, 0.005));
}

void 
COPCTorsion::changeSystem(COutputPlugin* plug)
{
  std::swap(Sim, static_cast<COPCTorsion*>(plug)->Sim);
  
#ifdef DYNAMO_DEBUG
  if (chains.size() != static_cast<COPCTorsion*>(plug)->chains.size())
    I_throw() << "CTorsion chain data size mismatch in replex exchange";
#endif

  std::list<CTCdata>::iterator iPtr1 = chains.begin(), 
    iPtr2 = static_cast<COPCTorsion*>(plug)->chains.begin();

  while(iPtr1 != chains.end())
    {

#ifdef DYNAMO_DEBUG
      if (iPtr1->chainPtr->getName() != iPtr2->chainPtr->getName())
	I_throw() << "Chain name mismatch when swapping chain plugins";
#endif

      std::swap(iPtr1->chainPtr, iPtr2->chainPtr);

      ++iPtr1;
      ++iPtr2;     
    }
}

void 
COPCTorsion::ticker()
{
  Sim->Dynamics.Liouvillean().updateAllParticles();

  BOOST_FOREACH(CTCdata& dat,chains)
    {
      Iflt sysGamma  = 0.0;
      long count = 0;
      BOOST_FOREACH(const smrtPlugPtr<CRange>& range,  dat.chainPtr->getMolecules())
	{
	  if (range->size() < 3)//Need three for curv and torsion
	    break;

#ifdef DYNAMO_DEBUG
	  if (NDIM != 3)
	    I_throw() << "Not implemented chain curvature in non 3 dimensional systems";
#endif
	  
	  CVector<> tmp;
	  std::vector<CVector<> > dr1;
	  std::vector<CVector<> > dr2;
	  std::vector<CVector<> > dr3;
	  std::vector<CVector<> > vec;

	  //Calc first and second derivatives
	  for (CRange::iterator it = range->begin() + 1; it != range->end() - 1; it++)
	    {
#ifdef DYNAMO_DEBUG
	      tmp = 0.5 * (Sim->vParticleList.at(*(it+1)).getPosition()
			   - Sim->vParticleList.at(*(it-1)).getPosition());
#else
	      tmp = 0.5 * (Sim->vParticleList[*(it+1)].getPosition()
			   - Sim->vParticleList[*(it-1)].getPosition());
#endif

	      dr1.push_back(tmp);
	      
#ifdef DYNAMO_DEBUG
	      tmp = Sim->vParticleList.at(*(it+1)).getPosition() 
		- (2.0 * Sim->vParticleList.at(*it).getPosition())
		+ Sim->vParticleList.at(*(it-1)).getPosition();
#else
	      tmp = Sim->vParticleList[*(it+1)].getPosition() 
		- (2.0 * Sim->vParticleList[*it].getPosition())
		+ Sim->vParticleList[*(it-1)].getPosition();
#endif
	      
	      dr2.push_back(tmp);
	      
	      vec.push_back(dr1.back().Cross(dr2.back()));
	    }
	  
	  //Create third derivative
	  for (std::vector<CVector<> >::iterator it2 = dr2.begin() + 1; it2 != dr2.end() - 1; it2++)
	    dr3.push_back(0.5 * (*(it2+1) - *(it2-1)));
	  
	  size_t derivsize =  dr3.size();

	  //Gamma Calc
	  Iflt gamma = 0.0;
	  Iflt torsion = 0, curvature = 0;
	  for (unsigned int i = 0; i < derivsize; i++)
	    {
#ifdef DYNAMO_DEBUG
	      torsion = ((vec.at(i+1)) % (dr3.at(i))) / (vec.at(i+1).square()); //Torsion
	      curvature = (vec.at(i+1).length()) / pow(dr1.at(i+1).length(), 3); //Curvature
#else
	      torsion = ((vec[i+1]) % (dr3[i])) / (vec[i+1].square()); //Torsion
	      curvature = (vec[i+1].length()) / pow(dr1[i+1].length(), 3); //Curvature
#endif
	      gamma += torsion / curvature;
	    }
	  gamma /= derivsize;
	  sysGamma += gamma;
	  ++count;
	  //Restrict the data collection to reasonable bounds
	  if (gamma < 10 && gamma > -10)
	    dat.gammaMol.addVal(gamma);
	}

      if (sysGamma < 10 && sysGamma > -10)
	dat.gammaSys.addVal(sysGamma/count);
    }
}

void 
COPCTorsion::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("ChainTorsion");
  
  BOOST_FOREACH(CTCdata& dat, chains)
    {
      XML << xmlw::tag(dat.chainPtr->getName().c_str())
	  << xmlw::tag("MolecularHistogram");
      
      dat.gammaMol.outputHistogram(XML, 1.0);
      
      XML << xmlw::endtag("MolecularHistogram")
	  << xmlw::tag("SystemHistogram");
      
      dat.gammaSys.outputHistogram(XML, 1.0);
      
      XML << xmlw::endtag("SystemHistogram")
	  << xmlw::endtag(dat.chainPtr->getName().c_str());
    }

  XML << xmlw::endtag("ChainTorsion");
}
