#include "CubePlane.hpp"
#include "../../datatypes/vector.hpp"
#include <boost/foreach.hpp>


bool 
DYNAMO::OverlapFunctions::CubePlane(const Vector& CubeOrigin, 
				    const Vector& CubeDimensions,
				    const Vector& PlaneOrigin, 
				    const Vector& PlaneNormal)
{
  Vector  relpos(CubeOrigin - PlaneOrigin);
  
  //Get which side the cube origin is on 
  bool OriginSign(std::signbit(relpos | PlaneNormal));
  
  CVector<size_t> counter(0);
  
  while (counter[NDIM-1] < 2)
    {
      Vector  pointpos(relpos);
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	if (counter[iDim]) pointpos[iDim] += CubeDimensions[iDim];

      if (std::signbit(pointpos | PlaneNormal) != OriginSign)
	return true;
      
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
