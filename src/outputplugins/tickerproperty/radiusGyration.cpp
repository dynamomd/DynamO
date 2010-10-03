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

#include "radiusGyration.hpp"
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../dynamics/ranges/1range.hpp"
#include <boost/foreach.hpp>
#include <vector>
#include "../../datatypes/vector.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/topology/include.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <fstream>
#include <cmath>

OPRGyration::OPRGyration(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"GyrationRadius"),
  binwidth1(0.01),
  binwidth2(0.001),
  binwidth3(0.01)
{
  operator<<(XML);
}

void 
OPRGyration::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth1"))
	binwidth1 = boost::lexical_cast<double>(XML.getAttribute("binwidth1"));

      if (XML.isAttributeSet("binwidth2"))
	binwidth2 = boost::lexical_cast<double>(XML.getAttribute("binwidth2"));

      if (XML.isAttributeSet("binwidth3"))
	binwidth3 = boost::lexical_cast<double>(XML.getAttribute("binwidth3"));
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPRGyration";
    }
  
}

void 
OPRGyration::initialise()
{
  BOOST_FOREACH(const magnet::ClonePtr<Topology>& plugPtr, Sim->dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(CTCdata(dynamic_cast<const CTChain*>
			       (plugPtr.get_ptr()), 
			       binwidth1 * Sim->dynamics.units().unitArea(), binwidth2, binwidth3));
}

void 
OPRGyration::changeSystem(OutputPlugin* plug)
{
  std::swap(Sim, static_cast<OPRGyration*>(plug)->Sim);

  std::list<CTCdata>::iterator iPtr1 = chains.begin(), iPtr2 = static_cast<OPRGyration*>(plug)->chains.begin();

#ifdef DYNAMO_DEBUG
  if (chains.size() != static_cast<OPRGyration*>(plug)->chains.size())
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
OPRGyration::getGyrationEigenSystem(const magnet::ClonePtr<CRange>& range, const DYNAMO::SimData* Sim)
{
  //Determine the centre of mass. Watch for periodic images
  Vector  tmpVec;  
  
  molGyrationDat retVal;
  retVal.MassCentre = Vector (0,0,0);

  double totmass = Sim->dynamics.getSpecies(Sim->particleList[*(range->begin())]).getMass();
  std::vector<Vector> relVecs;
  relVecs.reserve(range->size());
  relVecs.push_back(Vector(0,0,0));
  
  //Walk along the chain
  for (CRange::iterator iPtr = range->begin()+1; iPtr != range->end(); iPtr++)
    {
      Vector currRelPos = Sim->particleList[*iPtr].getPosition() 
	- Sim->particleList[*(iPtr - 1)].getPosition();

      Sim->dynamics.BCs().applyBC(currRelPos);

      relVecs.push_back(currRelPos + relVecs.back());

      retVal.MassCentre += relVecs.back() 
	* Sim->dynamics.getSpecies(Sim->particleList[*iPtr]).getMass();

      totmass += Sim->dynamics.getSpecies(Sim->particleList[*iPtr]).getMass();
    }

  retVal.MassCentre /= totmass;

  //Now determine the inertia tensor
  double data[NDIM * NDIM];
  for (size_t i = 0; i < NDIM; i++)
    for (size_t j = i; j < NDIM; j++)
      data[i + j * NDIM] = 0.0;
  
  BOOST_FOREACH(Vector & vec, relVecs)
    {
      vec -= retVal.MassCentre;

      for (size_t i = 0; i < NDIM; i++)
	for (size_t j = i; j < NDIM; j++)
	  data[i+j*NDIM] += vec[i] * vec[j];      
    }

  for (size_t i = 0; i < NDIM; ++i)
    for (size_t j = i+1; j < NDIM; ++j)
      data[j+i*NDIM] = data[i+j*NDIM];
  
  gsl_matrix_view m = gsl_matrix_view_array(data, NDIM, NDIM);
  gsl_vector *eval = gsl_vector_alloc (NDIM);
  gsl_matrix *evec = gsl_matrix_alloc (NDIM, NDIM);
  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (NDIM);
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free (w);
  //Sort Eigenvalues in ascending order
  gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_VAL_ASC);
  
  for (size_t i = 0; i < NDIM; i++)
    {
      retVal.EigenVal[i] = gsl_vector_get(eval, i) / range->size();

      gsl_vector_view evec_i 
	= gsl_matrix_column (evec, i);

      //EigenVec Components
      for (size_t j = 0; j < NDIM; j++)
	retVal.EigenVec[i][j] = gsl_vector_get(&evec_i.vector, j);
    }

  //Cleanup GSL
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  
  retVal.MassCentre += Sim->particleList[*(range->begin())].getPosition();

  return retVal;
}

Vector 
OPRGyration::NematicOrderParameter(const std::list<Vector  >& molAxis)
{
  double Q[NDIM * NDIM];
  
  for (size_t i = 0; i < NDIM; i++)
    for (size_t j = i; j < NDIM; j++)
      Q[i + j * NDIM] = 0.0;

  BOOST_FOREACH(const Vector & vec, molAxis)
    for (size_t i = 0; i < NDIM; i++)
      for (size_t j = i; j < NDIM; j++)
	Q[i + (j * NDIM)] += (3.0 * vec[i] * vec[j]) - (i==j ? 1 : 0);

  double Factor = 1.0 / (2.0 * molAxis.size());

  for (size_t i = 0; i < NDIM; i++)
    for (size_t j = i; j < NDIM; j++)
      Q[(j * NDIM) + i] *= Factor;

  //Copy over the triangle matrix
  for (size_t i = 0; i < NDIM-1; i++)
    for (size_t j = i+1; j < NDIM; j++)
      Q[(i * NDIM) + j] = Q[(j * NDIM) + i];
  
  gsl_matrix_view m = gsl_matrix_view_array(Q, NDIM, NDIM);
  gsl_vector *eval = gsl_vector_alloc (NDIM);
  gsl_matrix *evec = gsl_matrix_alloc (NDIM, NDIM);
  
  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (NDIM);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  
  gsl_eigen_symmv_free (w);
  
  gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_VAL_ASC);
  
  Vector  retval; 
  for (size_t iDim = 0; iDim < NDIM; iDim++)
    retval[iDim] = gsl_vector_get (eval, iDim);

  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  
  return retval;
}

double 
OPRGyration::CubaticOrderParameter(const std::list<Vector  >& molAxis)
{
  if (NDIM != 3)
    M_throw() << "Cubatic Order Parameter not implemented for non 3d sims!";

  //Cubatic Order Parameter
  double Q_cub[NDIM][NDIM][NDIM][NDIM];
  
  BOOST_FOREACH(const Vector & vec, molAxis)
    for (unsigned int iDim = 0; iDim < NDIM; iDim++)
      for (unsigned int jDim = 0; jDim < NDIM; jDim++)
	for (unsigned int kDim = 0; kDim < NDIM; kDim++)
	  for (unsigned int lDim = 0; lDim < NDIM; lDim++)
	    Q_cub[iDim][jDim][kDim][lDim] 
	      = (35.0/8.0) * vec[iDim] * vec[jDim] * vec[kDim] * vec[lDim]
	      - (5.0/8.0) * ( vec[iDim] * vec[jDim] * (kDim == lDim ? 1:0)
			      + vec[iDim] * vec[kDim] * (jDim == lDim ? 1 : 0)
			      + vec[iDim] * vec[lDim] * (jDim == kDim ? 1 : 0)
			      + vec[jDim] * vec[kDim] * (iDim == lDim ? 1 : 0)
			      + vec[jDim] * vec[lDim] * (iDim == kDim ? 1 : 0)
			      + vec[kDim] * vec[lDim] * (iDim == jDim ? 1 : 0))
	      + (1.0/8.0) * ( (iDim == jDim ? 1 : 0) * (kDim == lDim ? 1 : 0)
			      + (iDim == kDim ? 1 : 0) * (jDim == lDim ? 1 : 0)
			      + (iDim == lDim ? 1 : 0) * (jDim == kDim ? 1 : 0));
  
  double F[5*5];
  F[0 * 5 + 0] = Q_cub[0][0][0][0] - Q_cub[0][0][2][2];
  F[0 * 5 + 1] = Q_cub[0][0][0][1] + Q_cub[0][0][1][0];
  F[0 * 5 + 2] = Q_cub[0][0][0][2] + Q_cub[0][0][2][0];
  F[0 * 5 + 3] = Q_cub[0][0][1][1] - Q_cub[0][0][2][2];
  F[0 * 5 + 4] = Q_cub[0][0][1][2] + Q_cub[0][0][2][1];
  F[1 * 5 + 0] = Q_cub[0][1][0][0] - Q_cub[0][1][2][2];
  F[1 * 5 + 1] = Q_cub[0][1][0][1] + Q_cub[0][1][1][0];
  F[1 * 5 + 2] = Q_cub[0][1][0][2] + Q_cub[0][1][2][0];
  F[1 * 5 + 3] = Q_cub[0][1][1][1] - Q_cub[0][1][2][2];
  F[1 * 5 + 4] = Q_cub[0][1][1][2] + Q_cub[0][1][2][1];
  F[2 * 5 + 0] = Q_cub[0][2][0][0] - Q_cub[0][2][2][2];
  F[2 * 5 + 1] = Q_cub[0][2][0][1] + Q_cub[0][2][1][0];
  F[2 * 5 + 2] = Q_cub[0][2][0][2] + Q_cub[0][2][2][0];
  F[2 * 5 + 3] = Q_cub[0][2][1][1] - Q_cub[0][2][2][2];
  F[2 * 5 + 4] = Q_cub[0][2][1][2] + Q_cub[0][2][2][1];
  F[3 * 5 + 0] = Q_cub[1][1][0][0] - Q_cub[1][1][2][2];
  F[3 * 5 + 1] = Q_cub[1][1][0][1] + Q_cub[1][1][1][0];
  F[3 * 5 + 2] = Q_cub[1][1][0][2] + Q_cub[1][1][2][0];
  F[3 * 5 + 3] = Q_cub[1][1][1][1] - Q_cub[1][1][2][2];
  F[3 * 5 + 4] = Q_cub[1][1][1][2] + Q_cub[1][1][2][1];
  F[4 * 5 + 0] = Q_cub[1][2][0][0] - Q_cub[1][2][2][2];
  F[4 * 5 + 1] = Q_cub[1][2][0][1] + Q_cub[1][2][1][0];
  F[4 * 5 + 2] = Q_cub[1][2][0][2] + Q_cub[1][2][2][0];
  F[4 * 5 + 3] = Q_cub[1][2][1][1] - Q_cub[1][2][2][2];
  F[4 * 5 + 4] = Q_cub[1][2][1][2] + Q_cub[1][2][2][1];
  
  gsl_matrix_view m = gsl_matrix_view_array(F, 5, 5);
  gsl_vector *eval = gsl_vector_alloc (5);
  gsl_matrix *evec = gsl_matrix_alloc (5, 5);
  
  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (5);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);    
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_VAL_ASC);
  
  double retval = 8.0 * gsl_vector_get(eval, 4) / (7.0 * molAxis.size());
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  return retval;
}

void 
OPRGyration::ticker()
{
  BOOST_FOREACH(CTCdata& dat,chains)
    {
      std::list<Vector  > molAxis;

      BOOST_FOREACH(const magnet::ClonePtr<CRange>& range,  dat.chainPtr->getMolecules())
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
      
      dat.cubaticOrder.addVal(CubaticOrderParameter(molAxis));
    }
}

void 
OPRGyration::output(xml::XmlStream& XML)
{
  XML << xml::tag("ChainGyration");

  BOOST_FOREACH(CTCdata& dat, chains)
    {
      XML << xml::tag("Chain") << xml::attr("Name") 
	  << dat.chainPtr->getName().c_str()
	  << xml::tag("GyrationRadii");
      
      for (size_t i = 0; i< NDIM; i++)
	dat.gyrationRadii.at(i).outputHistogram(XML,1.0/Sim->dynamics.units().unitArea());

      XML << xml::endtag("GyrationRadii")
	  << xml::tag("NematicOrderParameter");

      std::list<Vector  > molAxis;

      BOOST_FOREACH(const magnet::ClonePtr<CRange>& range,  dat.chainPtr->getMolecules())
	molAxis.push_back(getGyrationEigenSystem(range, Sim).EigenVec[NDIM-1]);

      Vector  EigenVal = NematicOrderParameter(molAxis);
            
      for (size_t i = 0; i < NDIM; i++)
	if (std::isnormal(EigenVal[i]))
	  {
	    char lett[2] = {char('x' + i), '\0'};
	    
	    XML << xml::attr(lett)
		<< EigenVal[i];
	  }

      for (size_t i = 0; i<NDIM; i++)
	dat.nematicOrder.at(i).outputHistogram(XML, 1.0);
      
      XML << xml::endtag("NematicOrderParameter")
	  << xml::tag("CubaticOrderParameter")
	  << xml::attr("CurrentVal")
	  << CubaticOrderParameter(molAxis);

      dat.cubaticOrder.outputHistogram(XML, 1.0);

      XML << xml::endtag("CubaticOrderParameter")
	  << xml::endtag("Chain");
    }

  XML << xml::endtag("ChainGyration");
  
}
