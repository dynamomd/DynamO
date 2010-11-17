#include <iostream>
#include <magnet/math/quartic.hpp>
#include <vector>
#include "quartic.c"
#include <complex>

int main()
{
  setcns();

  const size_t nroots = 9;
  double rootvals[nroots] = {-1e3, -100, -10, -1, 0, 1, 10, +100, 1e3};

  std::cout << "\n\n///////////////////////4 real roots//////////////////////////";

  size_t counter = 0, accuracy = 0, rootcountfail = 0;
  for (size_t root1 = 0; root1 < nroots; ++root1)
    for (size_t root2 = root1; root2 < nroots; ++root2)
      for (size_t root3 = root2; root3 < nroots; ++root3)
	for (size_t root4 = root3; root4 < nroots; ++root4)
 	  {
	    ++counter;

 	    double a = -rootvals[root1] - rootvals[root2] - rootvals[root3]
 	      - rootvals[root4],
 	      b = rootvals[root1] * rootvals[root2]
 	      + rootvals[root1] * rootvals[root3]
 	      + rootvals[root1] * rootvals[root4]
 	      + rootvals[root2] * rootvals[root3]
 	      + rootvals[root2] * rootvals[root4]
 	      + rootvals[root3] * rootvals[root4],
 	      c = -rootvals[root1] * rootvals[root2] * rootvals[root3]
 	      - rootvals[root1] * rootvals[root2] * rootvals[root4]
 	      - rootvals[root1] * rootvals[root3] * rootvals[root4]
 	      - rootvals[root2] * rootvals[root3] * rootvals[root4],
 	      d = rootvals[root1] * rootvals[root2] * rootvals[root3]
 	      * rootvals[root4];
	    
 	    std::vector<double> originals(4);
 	    originals[0] = rootvals[root1];
 	    originals[1] = rootvals[root2];
 	    originals[2] = rootvals[root3];
 	    originals[3] = rootvals[root4];

 	    std::vector<double> roots(4);
	    
 	    size_t rootcount = magnet::math::quarticSolve(a, b, c, d,
 							  roots[0], roots[1], 
							  roots[2], roots[3]);
	    
	    size_t unique_roots = 1;
	    if (root2!=root1) ++unique_roots;
	    if ((root3!=root1) && (root3!=root2)) ++unique_roots;
	    if ((root4!=root1) && (root4!=root2) && (root4!=root3)) ++unique_roots;

 	    sort(originals.begin(), originals.end());
 	    sort(roots.begin(), roots.end());

	    bool accuracyfail = false;

	    for (size_t i = 0; i < rootcount; ++i)
	      if (std::abs((roots[i] - originals[i]) / originals[i]) > 0.0002)
		{ ++accuracy; accuracyfail = true; }

	    if ((rootcount == 4) && !accuracyfail) continue;	    
	    if (rootcount != 4) ++rootcountfail;	   
		
//	    if (accuracyfail) 
//	      std::cout << "\n\nAccuracy fail";
//	    else
//	      std::cout << "\n\nRoot count fail";
//	      
//	    std::cout << "\nActual             roots = " 
//		      << originals[0] << ","
//		      << originals[1] << ","
//		      << originals[2] << ","
//		      << originals[3];
//	    
//
//	    std::cout << "\nAlgorithm found " << rootcount << ", roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    double rts[4];
//	    int origrootcount = quartic(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::ferrariQuarticSolve(a, b, c, d,
//	    						  roots[0], roots[1],
//	    						  roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    
//	    std::cout << "\nFerrari found " << rootcount << ",   roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = ferrari(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::yacfraidQuarticSolve(a, b, c, d,
//	    						   roots[0], roots[1],
//	    						   roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nYacfraid found " << rootcount << ",  roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = yacfraid(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//
//	    
//	    rootcount = magnet::math::descartesQuarticSolve(a, b, c, d,
//	    						    roots[0], roots[1], 
//	    						    roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nDescartes found " << rootcount << ", roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = descartes(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::neumarkQuarticSolve(a, b, c, d,
//	    						  roots[0], roots[1],
//	    						  roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nNeumark found " << rootcount << ",   roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = neumark(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
	  }	  
  std::cout << "\nTested " << counter << " root combinations";
  std::cout << "\nAccuracy Failed in " << accuracy << " root combinations";
  std::cout << "\nRoot count failed in " << rootcountfail << " root combinations";



  std::cout << "\n\n////////////////////2 imaginary and 2 real roots//////////////////////";

  counter = 0; accuracy = 0; rootcountfail = 0;
  for (size_t root1real = 0; root1real < nroots; ++root1real)
    for (size_t root3real   = 0; root3real < nroots; ++root3real)
      for (size_t root2real = 0; root2real < nroots; ++root2real)
	for (size_t root2im   = 0; root2im  < nroots; ++root2im)
 	  {
	    if (rootvals[root2im] == 0) continue; //skip real double roots

	    ++counter;
	    
	    std::complex<double>
	      root1val(rootvals[root1real], 0),	    
	      root2val(rootvals[root3real], 0),
	      root3val(rootvals[root2real],  rootvals[root2im]),
	      root4val(rootvals[root2real], -rootvals[root2im]);

 	    double a = (- root1val - root2val 
			- root3val - root4val).real(),
 	      b = (+ root1val * root2val
		   + root1val * root3val
		   + root1val * root4val
		   + root2val * root3val
		   + root2val * root4val
		   + root3val * root4val).real(),
 	      c = (- root1val * root2val * root3val
		   - root1val * root2val * root4val
		   - root1val * root3val * root4val
		   - root2val * root3val * root4val).real(),
 	      d = (root1val * root2val * root3val
		   * root4val).real();
	    
 	    std::vector<double> roots(4);
	    
 	    size_t rootcount = magnet::math::quarticSolve(a, b, c, d,
 							  roots[0], roots[1], 
							  roots[2], roots[3]);

	    std::vector<double> originals;
	    originals.push_back(root1val.real());
	    originals.push_back(root2val.real());
	    sort(originals.begin(), originals.end());
	    sort(roots.begin(), roots.begin()+2);

	    bool accuracyfail = false;

	    for (size_t i = 0; (i < rootcount) && (i < 2); ++i)
	      if (std::abs((roots[i] - originals[i]) / originals[i]) > 0.0002)
		{ ++accuracy; accuracyfail = true; }
	    if (rootcount != 2) ++rootcountfail;

	    if ((rootcount == 2) && !accuracyfail) continue;

 	    rootcount = magnet::math::quarticSolve(a, b, c, d,
						   roots[0], roots[1], 
						   roots[2], roots[3]);


	    std::cout << "\n\nActual             roots = " 
		      << root1val.real() <<  ","
		      << root2val.real() <<  ","
		      << root3val.real() << " + " << root3val.imag() <<  " i,"
		      << root4val.real() << " + " << root4val.imag() <<  " i,";
	    

	    std::cout << "\nAlgorithm found " << rootcount << ", roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    double rts[4];
	    int origrootcount = quartic(a,b,c,d,rts);
	    std::sort(rts, rts+4);
	    std::cout << "\n Original found " << origrootcount << ", roots = ";
	    for (size_t i = 0; i < origrootcount; ++i)
	      std::cout << rts[i] << ",";
	    
	    rootcount = magnet::math::ferrariQuarticSolve(a, b, c, d,
	    						  roots[0], roots[1],
	    						  roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    
	    std::cout << "\nFerrari found " << rootcount << ",   roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    origrootcount = ferrari(a,b,c,d,rts);
	    std::sort(rts, rts+4);

	    std::cout << "\n Original found " << origrootcount << ", roots = ";
	    for (size_t i = 0; i < origrootcount; ++i)
	      std::cout << rts[i] << ",";
	    
	    rootcount = magnet::math::yacfraidQuarticSolve(a, b, c, d,
	    						   roots[0], roots[1],
	    						   roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nYacfraid found " << rootcount << ",  roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    origrootcount = yacfraid(a,b,c,d,rts);
	    std::sort(rts, rts+4);

	    std::cout << "\n Original found " << origrootcount << ", roots = ";
	    for (size_t i = 0; i < origrootcount; ++i)
	      std::cout << rts[i] << ",";

	    
	    rootcount = magnet::math::descartesQuarticSolve(a, b, c, d,
	    						    roots[0], roots[1], 
	    						    roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nDescartes found " << rootcount << ", roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    origrootcount = descartes(a,b,c,d,rts);
	    std::sort(rts, rts+4);

	    std::cout << "\n Original found " << origrootcount << ", roots = ";
	    for (size_t i = 0; i < origrootcount; ++i)
	      std::cout << rts[i] << ",";
	    
	    rootcount = magnet::math::neumarkQuarticSolve(a, b, c, d,
	    						  roots[0], roots[1],
	    						  roots[2], roots[3]);
 	    sort(roots.begin(), roots.end());
	    std::cout << "\nNeumark found " << rootcount << ",   roots = ";
	    for (size_t i = 0; i < rootcount; ++i)
	      std::cout << roots[i] << ",";

	    origrootcount = neumark(a,b,c,d,rts);
	    std::sort(rts, rts+4);

	    std::cout << "\n Original found " << origrootcount << ", roots = ";
	    for (size_t i = 0; i < origrootcount; ++i)
	      std::cout << rts[i] << ",";
	    
	  }

  std::cout << "\nTested " << counter << " root combinations";
  std::cout << "\nAccuracy Failed in " << accuracy << " root combinations";
  std::cout << "\nRoot count failed in " << rootcountfail << " root combinations";


  std::cout << "\n\n///////////////////////4 imaginary roots//////////////////////////";

  counter = 0; rootcountfail = 0;
  for (size_t root1real = 0; root1real < nroots; ++root1real)
    for (size_t root1im   = 0; root1im < nroots; ++root1im)
      for (size_t root2real = 0; root2real < nroots; ++root2real)
	for (size_t root2im   = 0; root2im  < nroots; ++root2im)
 	  {
	    if ((rootvals[root1im] == 0) || rootvals[root2im] == 0) continue;

	    ++counter;
	    
	    std::complex<double>
	      root1val(rootvals[root1real],  rootvals[root1im]),	    
	      root2val(rootvals[root1real], -rootvals[root1im]),
	      root3val(rootvals[root2real],  rootvals[root2im]),
	      root4val(rootvals[root2real], -rootvals[root2im]);

 	    double a = (-root1val - root2val 
			- root3val - root4val).real(),
 	      b = (root1val * root2val
		   + root1val * root3val
		   + root1val * root4val
		   + root2val * root3val
		   + root2val * root4val
		   + root3val * root4val).real(),
 	      c = (- root1val * root2val * root3val
		   - root1val * root2val * root4val
		   - root1val * root3val * root4val
		   - root2val * root3val * root4val).real(),
 	      d = (root1val * root2val * root3val
		   * root4val).real();
	    
 	    std::vector<double> roots(4);
	    
 	    size_t rootcount = magnet::math::quarticSolve(a, b, c, d,
 							  roots[0], roots[1], 
							  roots[2], roots[3]);

	    if (rootcount == 0) 
	      continue;
	    else
	      ++rootcountfail;

// 	    rootcount = magnet::math::quarticSolve(a, b, c, d,
// 							  roots[0], roots[1], 
//							  roots[2], roots[3]);
//
// 	    sort(roots.begin(), roots.end());
//
//	    std::cout << "\n\nActual             roots = " 
//		      << root1val.real() << " + " << root1val.imag() <<  " i,"
//		      << root2val.real() << " + " << root2val.imag() <<  " i,"
//		      << root3val.real() << " + " << root3val.imag() <<  " i,"
//		      << root4val.real() << " + " << root4val.imag() <<  " i,";
//	    
//
//	    std::cout << "\nAlgorithm found " << rootcount << ", roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    double rts[4];
//	    int origrootcount = quartic(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::ferrariQuarticSolve(a, b, c, d,
//	    						  roots[0], roots[1],
//	    						  roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    
//	    std::cout << "\nFerrari found " << rootcount << ",   roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = ferrari(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::yacfraidQuarticSolve(a, b, c, d,
//	    						   roots[0], roots[1],
//	    						   roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nYacfraid found " << rootcount << ",  roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = yacfraid(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//
//	    
//	    rootcount = magnet::math::descartesQuarticSolve(a, b, c, d,
//	    						    roots[0], roots[1], 
//	    						    roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nDescartes found " << rootcount << ", roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = descartes(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
//	    rootcount = magnet::math::neumarkQuarticSolve(a, b, c, d,
//	    						  roots[0], roots[1],
//	    						  roots[2], roots[3]);
// 	    sort(roots.begin(), roots.end());
//	    std::cout << "\nNeumark found " << rootcount << ",   roots = ";
//	    for (size_t i = 0; i < rootcount; ++i)
//	      std::cout << roots[i] << ",";
//
//	    origrootcount = neumark(a,b,c,d,rts);
//	    std::sort(rts, rts+4);
//
//	    std::cout << "\n Original found " << origrootcount << ", roots = ";
//	    for (size_t i = 0; i < origrootcount; ++i)
//	      std::cout << rts[i] << ",";
//	    
	  }	  
  std::cout << "\nTested " << counter << " root combinations";
  std::cout << "\nRoot count failed in " << rootcountfail << " root combinations";

  return 0;
}
