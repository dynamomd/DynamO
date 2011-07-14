#include "CubePlane.hpp"
#include <boost/foreach.hpp>


bool 
dynamo::OverlapFunctions::CubePlane(const Vector& CubeOrigin, 
				    const Vector& CubeDimensions,
				    const Vector& PlaneOrigin, 
				    const Vector& PlaneNormal,
				    const double tol)
{
  Vector  relpos(CubeOrigin - PlaneOrigin);
  
  size_t counter[3] = {0, 0, 0};
  
  while (counter[NDIM-1] < 2)
    {
      Vector  pointpos(relpos);
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	pointpos[iDim] += counter[iDim] * CubeDimensions[iDim];

      if ((pointpos | PlaneNormal) < tol) return true;
      
      ++counter[0];

      for (size_t iDim(0); iDim < NDIM-1; ++iDim)
	if (counter[iDim] > 1)
	  {
	    counter[iDim] = 0;
	    ++counter[iDim+1];
	  }
    }

  return false;
}
