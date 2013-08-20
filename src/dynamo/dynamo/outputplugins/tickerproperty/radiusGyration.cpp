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

#include <dynamo/outputplugins/tickerproperty/radiusGyration.hpp>
#include <dynamo/include.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <vector>
#include <fstream>
#include <cmath>

namespace dynamo {
  OPRGyration::OPRGyration(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"GyrationRadius"),
    binwidth1(0.01),
    binwidth2(0.001),
    binwidth3(0.01)
  {
    operator<<(XML);
  }

  void 
  OPRGyration::operator<<(const magnet::xml::Node& XML)
  {
    if (XML.hasAttribute("binwidth1"))
      binwidth1 = XML.getAttribute("binwidth1").as<double>();

    if (XML.hasAttribute("binwidth2"))
      binwidth2 = XML.getAttribute("binwidth2").as<double>();

    if (XML.hasAttribute("binwidth3"))
      binwidth3 = XML.getAttribute("binwidth3").as<double>();
  }

  void 
  OPRGyration::initialise()
  {
    for (const shared_ptr<Topology>& plugPtr : Sim->topology)
      if (std::dynamic_pointer_cast<TChain>(plugPtr))
	chains.push_back(CTCdata(static_cast<const TChain*>(plugPtr.get()), 
				 binwidth1 * Sim->units.unitArea(), binwidth2, binwidth3));
  }

  void 
  OPRGyration::replicaExchange(OutputPlugin& plug)
  {
    std::swap(Sim, static_cast<OPRGyration&>(plug).Sim);

    std::list<CTCdata>::iterator iPtr1 = chains.begin(), iPtr2 = static_cast<OPRGyration&>(plug).chains.begin();

#ifdef DYNAMO_DEBUG
    if (chains.size() != static_cast<OPRGyration&>(plug).chains.size())
      M_throw() << "Size mismatch when exchanging!";
#endif

    while (iPtr1 != chains.end())
      {
#ifdef DYNAMO_DEBUG
	if (iPtr1->chainPtr->getName() != iPtr2->chainPtr->getName())
	  M_throw() << "Name mismatch while replexing!";
#endif
	std::swap(iPtr1->chainPtr, iPtr2->chainPtr);

	++iPtr1;
	++iPtr2;
      }
  }

  OPRGyration::molGyrationDat
  OPRGyration::getGyrationEigenSystem(const shared_ptr<IDRange>& range, const dynamo::Simulation* Sim)
  {
    //Determine the centre of mass. Watch for periodic images
    Vector  tmpVec;  
  
    molGyrationDat retVal;
    retVal.MassCentre = Vector (0,0,0);

    double totmass = Sim->species[Sim->particles[*(range->begin())]]->getMass(*(range->begin()));
    std::vector<Vector> relVecs;
    relVecs.reserve(range->size());
    relVecs.push_back(Vector(0,0,0));
  
    //Walk along the chain
    for (IDRange::iterator iPtr = range->begin()+1; iPtr != range->end(); iPtr++)
      {
	Vector currRelPos = Sim->particles[*iPtr].getPosition() 
	  - Sim->particles[*(iPtr - 1)].getPosition();

	Sim->BCs->applyBC(currRelPos);

	relVecs.push_back(currRelPos + relVecs.back());

	double mass = Sim->species[Sim->particles[*iPtr]]->getMass(*iPtr);

	retVal.MassCentre += relVecs.back() * mass;
	totmass += mass;
      }

    retVal.MassCentre /= totmass;

    //Now determine the inertia tensor
    Matrix inertiaTensor;
    inertiaTensor.zero();

    for (Vector & vec : relVecs)
      {
	vec -= retVal.MassCentre;
	inertiaTensor += Dyadic(vec, vec);
      }

    std::pair<std::array<Vector, 3>, std::array<double, 3> > result
      = inertiaTensor.symmetric_eigen_decomposition();

    for (size_t i = 0; i < NDIM; i++)
      {	
	retVal.EigenVal[i] = result.second[i] / range->size();

	//EigenVec Components
	for (size_t j = 0; j < NDIM; j++)
	  retVal.EigenVec[i][j] = result.first[i][j];
      }

    retVal.MassCentre += Sim->particles[*(range->begin())].getPosition();

    return retVal;
  }

  Vector 
  OPRGyration::NematicOrderParameter(const std::list<Vector  >& molAxis)
  {
    Matrix Q;
    Q.zero();

    for (const Vector & vec : molAxis)
      for (size_t i = 0; i < NDIM; i++)
	for (size_t j = i; j < NDIM; j++)
	  Q(i,j) += (3.0 * vec[i] * vec[j]) - (i==j ? 1 : 0);

    double Factor = 1.0 / (2.0 * molAxis.size());

    for (size_t i = 0; i < NDIM; i++)
      for (size_t j = i; j < NDIM; j++)
	Q(i,j) *= Factor;

    //Copy over the triangle matrix
    for (size_t i = 0; i < NDIM-1; i++)
      for (size_t j = i+1; j < NDIM; j++)
	Q(j,i) = Q(i,j);
    
    std::pair<std::array<Vector, 3>, std::array<double, 3> > result
      = Q.symmetric_eigen_decomposition();

    return Vector(result.second[0], result.second[1], result.second[2]);
  }

  void 
  OPRGyration::ticker()
  {
    for (CTCdata& dat : chains)
      {
	std::list<Vector  > molAxis;

	for (const shared_ptr<IDRange>& range : dat.chainPtr->getMolecules())
	  {
	    molGyrationDat vals = getGyrationEigenSystem(range, Sim);	  
	    //Take the largest eigenvector as the molecular axis
	    molAxis.push_back(vals.EigenVec[NDIM-1]);
	    //Now add the radius of gyration
	    for (size_t iDim = 0; iDim < NDIM; ++iDim)
	      dat.gyrationRadii[iDim].addVal(vals.EigenVal[iDim]);
	  }
      
	Vector  EigenVal = NematicOrderParameter(molAxis);
      
	for (size_t i = 0; i < NDIM; i++)
	  if (std::isnormal(EigenVal[i]))
	    dat.nematicOrder[i].addVal(EigenVal[i]);
      }
  }

  void 
  OPRGyration::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("ChainGyration");

    for (CTCdata& dat : chains)
      {
	XML << magnet::xml::tag("Chain") << magnet::xml::attr("Name") 
	    << dat.chainPtr->getName().c_str()
	    << magnet::xml::tag("GyrationRadii");
      
	for (size_t i = 0; i< NDIM; i++)
	  dat.gyrationRadii.at(i).outputHistogram(XML,1.0/Sim->units.unitArea());

	XML << magnet::xml::endtag("GyrationRadii")
	    << magnet::xml::tag("NematicOrderParameter");

	std::list<Vector  > molAxis;

	for (const shared_ptr<IDRange>& range : dat.chainPtr->getMolecules())
	  molAxis.push_back(getGyrationEigenSystem(range, Sim).EigenVec[NDIM-1]);

	Vector  EigenVal = NematicOrderParameter(molAxis);
            
	for (size_t i = 0; i < NDIM; i++)
	  if (std::isnormal(EigenVal[i]))
	    {
	      char lett[2] = {char('x' + i), '\0'};
	    
	      XML << magnet::xml::attr(lett)
		  << EigenVal[i];
	    }

	for (size_t i = 0; i<NDIM; i++)
	  dat.nematicOrder.at(i).outputHistogram(XML, 1.0);
      
	XML << magnet::xml::endtag("NematicOrderParameter")
	    << magnet::xml::endtag("Chain");
      }

    XML << magnet::xml::endtag("ChainGyration");
  
  }
}
