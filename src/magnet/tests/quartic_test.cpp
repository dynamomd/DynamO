#include <iostream>
#include <magnet/math/quartic.hpp>
#include <vector>

int main()
{
  const size_t nroots = 7;
  double rootvals[nroots] = {-1e-6, -100, -1, 0, 1, +100, 1e6 };

  for (size_t root1 = 0; root1 < nroots; ++root1)
    for (size_t root2 = 0; root2 < nroots; ++root2)
      for (size_t root3 = 0; root3 < nroots; ++root3)
	for (size_t root4 = 0; root4 < nroots; ++root4)
 	  {
 	    double a = rootvals[root1] + rootvals[root2] + rootvals[root3]
 	      + rootvals[root4],
 	      b = rootvals[root1] * rootvals[root2]
 	      + rootvals[root1] * rootvals[root3]
 	      + rootvals[root1] * rootvals[root4]
 	      + rootvals[root2] * rootvals[root3]
 	      + rootvals[root2] * rootvals[root4]
 	      + rootvals[root3] * rootvals[root4],
 	      c = rootvals[root1] * rootvals[root2] * rootvals[root3]
 	      + rootvals[root1] * rootvals[root2] * rootvals[root4]
 	      + rootvals[root1] * rootvals[root3] * rootvals[root4]
 	      + rootvals[root2] * rootvals[root3] * rootvals[root4],
 	      d = rootvals[root1] * rootvals[root2] * rootvals[root3]
 	      * rootvals[root4];
	    
 	    std::vector<double> originals(4);
 	    originals[0] = -rootvals[root1];
 	    originals[1] = -rootvals[root2];
 	    originals[2] = -rootvals[root3];
 	    originals[3] = -rootvals[root4];

 	    std::vector<double> roots(4);
	    
 	    size_t rootcount = magnet::math::quarticSolve(a, b, c, d,
 							  roots[0], roots[1], 
							  roots[2], roots[3]);
	    
 	    sort(originals.begin(), originals.end());
 	    sort(roots.begin(), roots.end());


	    std::cout << "\nActual             roots = " 
		      << originals[0] << ","
		      << originals[1] << ","
		      << originals[2] << ","
		      << originals[3];

	    //if (rootcount == 4) continue;

	    std::cout << "\nAlgorithm found " << rootcount << ", roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";
	    
	    rootcount = magnet::math::ferrariQuarticSolve(a, b, c, d,
							  roots[0], roots[1],
							  roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());

	    std::cout << "\nFerrari found " << rootcount << ",   roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    rootcount = magnet::math::yacfraidQuarticSolve(a, b, c, d,
							   roots[0], roots[1],
							   roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nYacfraid found " << rootcount << ",  roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    rootcount = magnet::math::descartesQuarticSolve(a, b, c, d,
							    roots[0], roots[1], 
							    roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nDescartes found " << rootcount << ", roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    rootcount = magnet::math::neumarkQuarticSolve(a, b, c, d,
							  roots[0], roots[1],
							  roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nNeumark found " << rootcount << ",   roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";
	    
	  }	  
  return 0;
}
