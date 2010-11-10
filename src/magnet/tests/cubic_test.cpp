#include <iostream>
#include <magnet/math/quartic.hpp>
#include <vector>

int main()
{
  const size_t nroots = 9;
  double rootvals[nroots] = {-1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6 };

  for (size_t root1 = 0; root1 < nroots; ++root1)
    for (size_t root2 = 0; root2 < nroots; ++root2)
      for (size_t root3 = 0; root3 < nroots; ++root3)
	{
	  double a = -rootvals[root1] - rootvals[root2]  - rootvals[root3],
	    b = rootvals[root1] * rootvals[root2]
	    + rootvals[root1] * rootvals[root3]
	    + rootvals[root2] * rootvals[root3],
	    c = - rootvals[root1] * rootvals[root2] * rootvals[root3];
	  
	  std::vector<double> originals(3);
	  originals[0] = rootvals[root1];
	  originals[1] = rootvals[root2];
	  originals[2] = rootvals[root3];
	  
	  std::vector<double> roots(3);
	  
	  size_t rootcount = magnet::math::cubicSolve(a, b, c,
						      roots[0], roots[1], 
						      roots[2]);	    
	  sort(originals.begin(), originals.end());
	  sort(roots.begin(), roots.end());
	  
	  if (rootcount != 3)
	    {
	      rootcount = magnet::math::cubicSolve(a, b, c,
						      roots[0], roots[1], 
						      roots[2]);

	      std::cout << "\n\nRoot Count Failed"
			<< "\nActual             roots = " 
			<< originals[0] << ","
			<< originals[1] << ","
			<< originals[2]
			<< "\nAlgorithm found " << rootcount << ", roots = ";
	      for (size_t i = 0; i < rootcount; ++i)
		std::cout << roots[i] << ",";

	      std::cout.flush();

	      break;
	    }

	  for (size_t i = 0; i < rootcount; ++i)
	    if (std::abs((roots[i] - originals[i]) / originals[i]) > 0.000001)
	      {
		static size_t accuracyfails = 0;

		std::cout << "\n\nRoot accuracy Failure = " << ++accuracyfails
			  << "\nActual             roots = " 
			  << originals[0] << ","
			  << originals[1] << ","
			  << originals[2]
			  << "\nAlgorithm found " << rootcount << ", roots = ";
		for (size_t i = 0; i < rootcount; ++i)
		  std::cout << roots[i] << ",";

		std::cout.flush();
		break;
	      }
	}
  return 0;
}
