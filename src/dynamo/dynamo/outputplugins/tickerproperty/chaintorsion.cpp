/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/tickerproperty/chaintorsion.hpp>
#include <dynamo/include.hpp>
#include <dynamo/ranges/1range.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/BC/None.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
  OPCTorsion::OPCTorsion(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OPTicker(tmp,"Torsion")
  {}

  void 
  OPCTorsion::initialise()
  {
    BOOST_FOREACH(const shared_ptr<Topology>& plugPtr, Sim->topology)
      if (std::tr1::dynamic_pointer_cast<TChain>(plugPtr))
	chains.push_back(CTCdata(static_cast<const TChain*>(plugPtr.get()), 
				 0.005, 0.005, 0.01));

    if (!std::tr1::dynamic_pointer_cast<BCNone>(Sim->BCs))
      M_throw() << "Can only use this plugin with Null BC's"
		<< "\nPositions must be unwrapped";
  }

  void 
  OPCTorsion::changeSystem(OutputPlugin* plug)
  {
    std::swap(Sim, static_cast<OPCTorsion*>(plug)->Sim);
  
#ifdef DYNAMO_DEBUG
    if (chains.size() != static_cast<OPCTorsion*>(plug)->chains.size())
      M_throw() << "CTorsion chain data size mismatch in replex exchange";
#endif

    std::list<CTCdata>::iterator iPtr1 = chains.begin(), 
      iPtr2 = static_cast<OPCTorsion*>(plug)->chains.begin();

    while(iPtr1 != chains.end())
      {

#ifdef DYNAMO_DEBUG
	if (iPtr1->chainPtr->getName() != iPtr2->chainPtr->getName())
	  M_throw() << "Chain name mismatch when swapping chain plugins";
#endif

	std::swap(iPtr1->chainPtr, iPtr2->chainPtr);

	++iPtr1;
	++iPtr2;     
      }
  }

  void 
  OPCTorsion::ticker()
  {
    BOOST_FOREACH(CTCdata& dat,chains)
      {
	double sysGamma  = 0.0;
	long count = 0;
	BOOST_FOREACH(const shared_ptr<Range>& range,  dat.chainPtr->getMolecules())
	  {
	    if (range->size() < 3)//Need three for curv and torsion
	      break;

#ifdef DYNAMO_DEBUG
	    if (NDIM != 3)
	      M_throw() << "Not implemented chain curvature in non 3 dimensional systems";
#endif
	  
	    Vector tmp;
	    std::vector<Vector> dr1;
	    std::vector<Vector> dr2;
	    std::vector<Vector> dr3;
	    std::vector<Vector> vec;

	    //Calc first and second derivatives
	    for (Range::iterator it = range->begin() + 1; it != range->end() - 1; it++)
	      {
		tmp = 0.5 * (Sim->particles[*(it+1)].getPosition()
			     - Sim->particles[*(it-1)].getPosition());

		dr1.push_back(tmp);
	      
		tmp = Sim->particles[*(it+1)].getPosition() 
		  - (2.0 * Sim->particles[*it].getPosition())
		  + Sim->particles[*(it-1)].getPosition();
	      
		dr2.push_back(tmp);
	      
		vec.push_back(dr1.back() ^ dr2.back());
	      }
	  
	    //Create third derivative
	    for (std::vector<Vector  >::iterator it2 = dr2.begin() + 1; it2 != dr2.end() - 1; it2++)
	      dr3.push_back(0.5 * (*(it2+1) - *(it2-1)));
	  
	    size_t derivsize = dr3.size();

	    //Gamma Calc
	    double gamma = 0.0;
	    double fsum = 0.0;

	    for (unsigned int i = 0; i < derivsize; i++)
	      {
		double torsion = ((vec[i+1]) | (dr3[i])) / (vec[i+1].nrm2()); //Torsion
		double curvature = (vec[i+1].nrm()) / pow(dr1[i+1].nrm(), 3); //Curvature

		double instGamma = torsion / curvature;
		gamma += instGamma;

		double helixradius = 1.0/(curvature * (1.0+instGamma*instGamma));

		double minradius = HUGE_VAL;

		for (Range::iterator it1 = range->begin(); 
		     it1 != range->end(); it1++)
		  //Check this particle is not the same, or adjacent
		  if (*it1 != *(range->begin()+2+i)
		      && *it1 != *(range->begin()+1+i)
		      && *it1 != *(range->begin()+3+i))
		    for (Range::iterator it2 = range->begin() + 1; 
			 it2 != range->end() - 1; it2++)
		      //Check this particle is not the same, or adjacent to the studied particle
		      if (*it1 != *it2
			  && *it2 != *(range->begin()+2+i)
			  && *it2 != *(range->begin()+1+i)
			  && *it2 != *(range->begin()+3+i))
			{
			  //We have three points, calculate the lengths
			  //of the triangle sides
			  double a = (Sim->particles[*it1].getPosition() 
				      - Sim->particles[*it2].getPosition()).nrm(),
			    b = (Sim->particles[*(range->begin()+2+i)].getPosition() 
				 - Sim->particles[*it2].getPosition()).nrm(),
			    c = (Sim->particles[*it1].getPosition() 
				 - Sim->particles[*(range->begin()+2+i)].getPosition()).nrm();

			  //Now calc the area of the triangle
			  double s = (a + b + c) / 2.0;
			  double A = std::sqrt(s * (s - a) * (s - b) * (s - c));
			  double R = a * b * c / (4.0 * A);
			  if (R < minradius) minradius = R;			 
			}
		fsum += minradius / helixradius;
	      }

	    gamma /= derivsize;
	    sysGamma += gamma;
	    fsum /= derivsize;

	    ++count;
	    //Restrict the data collection to reasonable bounds
	    if (gamma < 10 && gamma > -10)
	      dat.gammaMol.addVal(gamma);

	    dat.f.addVal(fsum);
	  }

	if (sysGamma < 10 && sysGamma > -10)
	  dat.gammaSys.addVal(sysGamma/count);
      }
  }

  void 
  OPCTorsion::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("ChainTorsion");
  
    BOOST_FOREACH(CTCdata& dat, chains)
      {
	XML << magnet::xml::tag(dat.chainPtr->getName().c_str())
	    << magnet::xml::tag("MolecularHistogram");
      
	dat.gammaMol.outputHistogram(XML, 1.0);
      
	XML << magnet::xml::endtag("MolecularHistogram")
	    << magnet::xml::tag("SystemHistogram");
      
	dat.gammaSys.outputHistogram(XML, 1.0);
      
	XML << magnet::xml::endtag("SystemHistogram")
	    << magnet::xml::tag("FHistogram");
      
	dat.f.outputHistogram(XML, 1.0);
      
	XML << magnet::xml::endtag("FHistogram")
	    << magnet::xml::endtag(dat.chainPtr->getName().c_str());
      }

    XML << magnet::xml::endtag("ChainTorsion");
  }
}
