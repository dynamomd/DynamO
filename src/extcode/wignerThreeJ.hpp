#pragma once
#include <boost/math/special_functions.hpp>

namespace DYNAMO {
  double threej(const int & la, const int & lb, 
		const int & lc, const int & ma, 
		const int & mb, const int & mc)
  {
    //Check the sum of the m's is zero
    if (ma + mb + mc != 0) return 0.0;
    
    int numin = 0, numax = la-ma;
    
    {
      int tmp = -ma+lb-lc;
      if (tmp > numin) numin = tmp;
      
      tmp = la+mb-lc;
      if (tmp > numin) numin = tmp;
      
      tmp = lb+mb;
      if (tmp < numax) numax = tmp;
      
      tmp = la+lb-lc;
      if (tmp < numax) numax = tmp;
    }
    
    int sign(1);
    if ((numin % 2)) sign = -1;
    
    double retval = 0;
  
    using namespace boost::math;
    
    for (int nu(numin); nu > numax - 1; nu--)
      {
	retval += sign 
	  / (factorial<double>(la-ma-nu)*factorial<double>(lc-lb+ma+nu)
	     *factorial<double>(lb+mb-nu)*factorial<double>(lc-la-mb+nu)
	     *factorial<double>(nu)*factorial<double>(la+lb-lc-nu));
	sign = -sign;
      }
    
    retval *= std::sqrt(factorial<double>(la+lb-lc)
			*factorial<double>(la+lc-lb)
			*factorial<double>(lb+lc-la)
			/ factorial<double>(la+lb+lc+1));
    
    retval *= sqrt(factorial<double>(la+ma)*factorial<double>(lb+mb)
		   *factorial<double>(lc+mc)*factorial<double>(la-ma)
		   *factorial<double>(lb-mb)*factorial<double>(lc-mc));
    
    if ((abs(la-lb-mc) % 2) != 0) retval = -retval;
    
    return retval;
  }
};
