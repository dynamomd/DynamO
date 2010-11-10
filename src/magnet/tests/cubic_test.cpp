#include <iostream>
#include <magnet/math/quartic.hpp>
#include <vector>
#include "quartic.c"
#include <complex>

int main()
{

  //Test for triple roots
  const size_t nroots = 9;
  double rootvals[nroots] = {-1e6, -1e3, -100, -1, 0, 1, +100, 1e3, 1e6 };


  //Triple roots
//  for (size_t root1 = 0; root1 < nroots; ++root1)
//    for (size_t root2 = root1; root2 < nroots; ++root2)
//      for (size_t root3 = root2; root3 < nroots; ++root3)
//	{
//	  double a = -rootvals[root1] - rootvals[root2]  - rootvals[root3],
//	    b = rootvals[root1] * rootvals[root2]
//	    + rootvals[root1] * rootvals[root3]
//	    + rootvals[root2] * rootvals[root3],
//	    c = - rootvals[root1] * rootvals[root2] * rootvals[root3];
//	  
//	  std::vector<double> originals(3);
//	  originals[0] = rootvals[root1];
//	  originals[1] = rootvals[root2];
//	  originals[2] = rootvals[root3];
//	  
//	  std::vector<double> roots(3);
//	  
//	  size_t rootcount = magnet::math::cubicSolve(a, b, c,
//						      roots[0], roots[1], 
//						      roots[2]);	    
//	  sort(originals.begin(), originals.end());
//	  sort(roots.begin(), roots.end());
//	  
//	  if (rootcount != 3)
//	    {
//	      rootcount = magnet::math::cubicSolve(a, b, c,
//						   roots[0], roots[1], 
//						   roots[2]);
//
//	      std::cout << "\n\nRoot Count Failed"
//			<< "\nActual             roots = " 
//			<< originals[0] << ","
//			<< originals[1] << ","
//			<< originals[2]
//			<< "\nAlgorithm found " << rootcount << ", roots = ";
//	      for (size_t i = 0; i < rootcount; ++i)
//		std::cout << roots[i] << ",";
//
//	      double rts[4];
//	      int nroots = cubic(a,b,c,rts);
//
//	      std::cout << "\n\nOriginal found " << nroots
//			<< "\nOriginal           roots = ";
//	      for (size_t i = 0; i < nroots; ++i)
//		std::cout << roots[i] << ",";
//
//	      std::cout.flush();
//
//	      break;
//	    }
//
//	  for (size_t i = 0; i < rootcount; ++i)
//	    if (std::abs((roots[i] - originals[i]) / originals[i]) > 0.000001)
//	      {
//		static size_t accuracyfails = 0;
//
//		std::cout << "\n\nRoot accuracy Failure = " << ++accuracyfails
//			  << "\nActual             roots = " 
//			  << originals[0] << ","
//			  << originals[1] << ","
//			  << originals[2]
//			  << "\nAlgorithm found " << rootcount << ", roots = ";
//		for (size_t i = 0; i < rootcount; ++i)
//		  std::cout << roots[i] << ",";
//		
//		double rts[4];
//		int nroots = cubic(a,b,c,rts);
//		
//		std::cout << "\n\nOriginal found " << nroots
//			  << "\nOriginal           roots = ";
//		for (size_t i = 0; i < nroots; ++i)
//		  std::cout << roots[i] << ",";
//
//		std::cout.flush();
//		break;
//	      }
//	}

  //One real root and one pair of imaginary roots
  for (size_t root1 = 0; root1 < nroots; ++root1)
    for (size_t root2real = 0; root2real < nroots; ++root2real)
      for (size_t root2im = 0; root2im < nroots; ++root2im)
	{
	  if (rootvals[root2im] == 0) continue;

	  std::complex<double>
	    root1val(rootvals[root1],      0),	    
	    root2val(rootvals[root2real],  rootvals[root2im]),
	    root3val(rootvals[root2real], -rootvals[root2im]);
	  
	  double a = (-root1val - root2val  - root3val).real(),
	    b = (root1val * root2val
		 + root1val * root3val
		 + root2val * root3val).real(),
	    c = - (root1val * root2val * root3val).real();
	  
	  std::vector<double> roots(3);
	  
	  size_t rootcount = magnet::math::cubicSolve(a, b, c,
						      roots[0], roots[1], 
						      roots[2]);
	  if (rootcount != 1)
	    {
	      rootcount = magnet::math::cubicSolve(a, b, c,
						   roots[0], roots[1], 
						   roots[2]);

	      std::cout << "\n\nRoot Count Failed"
			<< "\nActual             roots = " 
			<< root1val.real() << ","
			<< root2val.real() << " + " << root2val.imag() << " i,"
			<< root3val.real() << " + " << root3val.imag() << " i,"
			<< "\nAlgorithm found " << rootcount << ", roots = ";
	      for (size_t i = 0; i < rootcount; ++i)
		std::cout << roots[i] << ",";

	      double rts[4];
	      int nroots = cubic(a,b,c,rts);

	      std::cout << "\n\nOriginal found " << nroots
			<< "\nOriginal           roots = ";
	      for (size_t i = 0; i < nroots; ++i)
		std::cout << roots[i] << ",";

	      std::cout.flush();

	      break;
	    }

	  for (size_t i = 0; i < rootcount; ++i)
	    if (std::abs((roots[0] - root1val.real()) / root1val.real()) > 0.000001)
	      {
		static size_t accuracyfails = 0;

		std::cout << "\n\nRoot accuracy Failure = " << ++accuracyfails
			  << "\nActual             roots = " 
			  << root2val.real() << " + " << root2val.imag() << " i,"
			  << root3val.real() << " + " << root3val.imag() << " i,"
			  << "\nAlgorithm found " << rootcount << ", roots = "
			  << root1val.real() << ",";
		
		double rts[4];
		int nroots = cubic(a,b,c,rts);
		
		std::cout << "\n\nOriginal found " << nroots
			  << "\nOriginal           roots = ";
		for (size_t i = 0; i < nroots; ++i)
		  std::cout << roots[i] << ",";

		std::cout.flush();
		break;
	      }
	}


  return 0;
}
