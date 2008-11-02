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

#include "radiusGyration.hpp"
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <fstream>
#include <cmath>

COPRGyration::COPRGyration(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPTicker(tmp,"GyrationRadius"),
  binwidth1(0.01),
  binwidth2(0.001),
  binwidth3(0.01)
{
  operator<<(XML);
}

void 
COPRGyration::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("binwidth1"))
	binwidth1 = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth1"));

      if (XML.isAttributeSet("binwidth2"))
	binwidth2 = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth2"));

      if (XML.isAttributeSet("binwidth3"))
	binwidth3 = boost::lexical_cast<Iflt>(XML.getAttribute("binwidth3"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPRGyration";
    }
  
}

void 
COPRGyration::initialise()
{
  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      chains.push_back(CTCdata(dynamic_cast<const CTChain*>
			       (plugPtr.get_ptr()), 
			       binwidth1 * Sim->Dynamics.units().unitLength(), binwidth2, binwidth3));
}

void 
COPRGyration::changeSystem(COutputPlugin* plug)
{
  std::swap(Sim, static_cast<COPRGyration*>(plug)->Sim);

  std::list<CTCdata>::iterator iPtr1 = chains.begin(), iPtr2 = static_cast<COPRGyration*>(plug)->chains.begin();

#ifdef DYNAMO_DEBUG
  if (chains.size() != static_cast<COPRGyration*>(plug)->chains.size())
    D_throw() << "Size mismatch when exchanging!";
#endif

  while (iPtr1 != chains.end())
    {
#ifdef DYNAMO_DEBUG
      if (iPtr1->chainPtr->getName() != iPtr2->chainPtr->getName())
	D_throw() << "Name mismatch while replexing!";
#endif
      std::swap(iPtr1->chainPtr, iPtr2->chainPtr);

      ++iPtr1;
      ++iPtr2;
    }
}

COPRGyration::molGyrationDat
COPRGyration::getGyrationEigenSystem(const smrtPlugPtr<CRange>& range, const DYNAMO::SimData* Sim)
{
  //I'm assuming you Updated all particles positions first!

  //Determine the centre of mass. Watch for periodic images
  CVector<> tmpVec;  
  
  molGyrationDat retVal;
  retVal.MassCentre = CVector<>(0);  

  CVector<> currRelPos(0);
  Iflt totmass = Sim->Dynamics.getSpecies(Sim->vParticleList[*(range->begin())]).getMass();
  std::list<CVector<> > relVecs;
  relVecs.push_back(currRelPos);
  
  unsigned long ID;
  //Walk along the chain
  for (CRange::iterator iPtr = range->begin()+1; iPtr != range->end(); iPtr++)
    {
      ID = *iPtr;

      currRelPos = Sim->vParticleList[ID].getPosition() - Sim->vParticleList[*(iPtr - 1)].getPosition();
      Sim->Dynamics.BCs().setPBC(currRelPos);

      relVecs.push_back(currRelPos + relVecs.back());

      retVal.MassCentre += relVecs.back() * Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();
      totmass += Sim->Dynamics.getSpecies(Sim->vParticleList[ID]).getMass();

    }

  retVal.MassCentre /= totmass;

  //Now determine the inertia tensor
  double data[NDIM * NDIM];
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      data[i + j * NDIM] = 0.0;
  
  BOOST_FOREACH(CVector<>& vec, relVecs)
    {
      vec -= retVal.MassCentre;

      for (int i = 0; i < NDIM; i++)
	for (int j = i; j < NDIM; j++)
	  data[i+j*NDIM] += vec[i] * vec[j];      
    }

  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
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
  
  for (int i = 0; i < NDIM; i++)
    {
      retVal.EigenVal[i] = gsl_vector_get(eval, i);

      gsl_vector_view evec_i 
	= gsl_matrix_column (evec, i);

      //EigenVec Components
      for (int j = 0; j < NDIM; j++)
	retVal.EigenVec[i][j] = gsl_vector_get(&evec_i.vector, j);
    }

  //Cleanup GSL
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  
  retVal.MassCentre += Sim->vParticleList[*(range->begin())].getPosition();
  
  return retVal;
}

CVector<>
COPRGyration::NematicOrderParameter(const std::list<CVector<> >& molAxis)
{
  double Q[NDIM * NDIM];
  
  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      Q[i + j * NDIM] = 0.0;

  BOOST_FOREACH(const CVector<>& vec, molAxis)
    for (int i = 0; i < NDIM; i++)
      for (int j = i; j < NDIM; j++)
	Q[i + (j * NDIM)] += (3.0 * vec[i] * vec[j]) - (i==j ? 1 : 0);

  Iflt Factor = 1.0 / (2.0 * molAxis.size());

  for (int i = 0; i < NDIM; i++)
    for (int j = i; j < NDIM; j++)
      Q[(j * NDIM) + i] *= Factor;

  //Copy over the triangle matrix
  for (int i = 0; i < NDIM-1; i++)
    for (int j = i+1; j < NDIM; j++)
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
  
  CVector<> retval; 
  for (int iDim = 0; iDim < NDIM; iDim++)
    retval[iDim] = gsl_vector_get (eval, iDim);

  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  
  return retval;
}

Iflt 
COPRGyration::CubaticOrderParameter(const std::list<CVector<> >& molAxis)
{
  if (NDIM != 3)
    D_throw() << "Cubatic Order Parameter not implemented for non 3d sims!";

  //Cubatic Order Parameter
  Iflt Q_cub[NDIM][NDIM][NDIM][NDIM];
  
  BOOST_FOREACH(const CVector<>& vec, molAxis)
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
  
  Iflt retval = 8.0 * gsl_vector_get(eval, 4) / (7.0 * molAxis.size());
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  return retval;
}

void 
COPRGyration::ticker()
{
  BOOST_FOREACH(CTCdata& dat,chains)
    {
      std::list<CVector<> > molAxis;

      BOOST_FOREACH(const smrtPlugPtr<CRange>& range,  dat.chainPtr->getMolecules())
	{
	  //Update what you need
	  BOOST_FOREACH(const unsigned long& ID, *range)
	    Sim->Dynamics.Liouvillean().updateParticle(Sim->vParticleList[ID]);

	  molGyrationDat vals = getGyrationEigenSystem(range, Sim);	  
	  //Take the largest eigenvector as the molecular axis
	  molAxis.push_back(vals.EigenVec[NDIM-1]);
	  //Now add the radius of gyration
	  for (int iDim = 0; iDim < NDIM; ++iDim)
	    dat.gyrationRadii[iDim].addVal(vals.EigenVal[iDim]);
	}
      
      CVector<> EigenVal = NematicOrderParameter(molAxis);
      
      for (int i = 0; i < NDIM; i++)
	if (std::isnormal(EigenVal[i]))
	  dat.nematicOrder[i].addVal(EigenVal[i]);
      
      dat.cubaticOrder.addVal(CubaticOrderParameter(molAxis));
    }
}

void 
COPRGyration::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("ChainGyration");

  BOOST_FOREACH(CTCdata& dat, chains)
    {
      XML << xmlw::tag("Chain") << xmlw::attr("Name") 
	  << dat.chainPtr->getName().c_str()
	  << xmlw::tag("GyrationRadii");
      
      for (int i = 0; i<NDIM; i++)
	dat.gyrationRadii.at(i).outputHistogram(XML,1.0/Sim->Dynamics.units().unitArea());

      XML << xmlw::endtag("GyrationRadii")
	  << xmlw::tag("NematicOrderParameter");

      std::list<CVector<> > molAxis;

      BOOST_FOREACH(const smrtPlugPtr<CRange>& range,  dat.chainPtr->getMolecules())
	molAxis.push_back(getGyrationEigenSystem(range, Sim).EigenVec[NDIM-1]);

      CVector<> EigenVal = NematicOrderParameter(molAxis);
            
      for (int i = 0; i < NDIM; i++)
	if (std::isnormal(EigenVal[i]))
	  {
	    char lett[2] = {'x' + i, '\0'};
	    
	    XML << xmlw::attr(lett)
		<< EigenVal[i];
	  }

      for (int i = 0; i<NDIM; i++)
	dat.nematicOrder.at(i).outputHistogram(XML, 1.0);
      
      XML << xmlw::endtag("NematicOrderParameter")
	  << xmlw::tag("CubaticOrderParameter")
	  << xmlw::attr("CurrentVal")
	  << CubaticOrderParameter(molAxis);

      dat.cubaticOrder.outputHistogram(XML, 1.0);

      XML << xmlw::endtag("CubaticOrderParameter")
	  << xmlw::endtag("Chain");
    }

  XML << xmlw::endtag("ChainGyration");
  
}
