#include <iostream>
#include <iomanip>

#define DYNAMO_double_precsision
#include "src/extcode/mathtemplates.hpp"

typedef size_t (*quarticSolverType)(const Iflt&, const Iflt&, const Iflt&, const Iflt&,
		   Iflt&, Iflt&, Iflt&, Iflt&);


Iflt hybroot1,hybroot2,hybroot3, hybroot4;

void printDetails(std::string name, quarticSolverType func, size_t totalRoots, 
		  const Iflt& a, const Iflt& b, const Iflt& c, const Iflt& d)
{

  boost::array<Iflt, 4> roots;
  size_t rootCount = (*func)(a, b, c, d, roots[0], roots[1], roots[2], roots[3]);
  
  
  if (rootCount != totalRoots)
    std::cout << name << " found a different number of roots " 
	      << rootCount << "v" << totalRoots << "\n"
	      << "\t a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
  else
    {
      Iflt error = quarticError(a,b,c,d,roots,rootCount);
      if (error > 1e-6)
	{
	  if ((rootCount >= 1) && (totalRoots >= 1))
	    std::cout << name << " hybroot1=" << hybroot1 << " root1=" << roots[0] 
		      <<  " dev=" << roots[0]/hybroot1-1
		      << " a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
	  
	  if ((rootCount >= 2) && (totalRoots >= 2))
	    std::cout << name << " hybroot2=" << hybroot2 << " root2=" << roots[1] 
		      <<  " dev=" << roots[1]/hybroot2-1
		      << " a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
	  
	  if ((rootCount >= 3) && (totalRoots >= 3))
	    std::cout << name << " hybroot3=" << hybroot3 << " root3=" << roots[2] 
		      <<  " dev=" << roots[2]/hybroot3-1
		      << " a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
	  
	  if ((rootCount >= 4) && (totalRoots >= 4))
	    std::cout << name << " hybroot4=" << hybroot4 << " root4=" << roots[3] 
		      <<  " dev=" << roots[3]/hybroot4-1
		      << " a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
	  std::cout << "\n";
	}
    }
}

int main()
{
  
  //  std::cout << "Testing cubic solutions\n";
//  Iflt p = 0,
//    q = -5000000,
//    r = 1;
//  Iflt root1,root2,root3;
//  size_t roots = cubicSolve(p, q, r, root1, root2, root3);
//  
//  std::cout << "root1 = " << root1 << "\n";
//  if (roots >1)
//    std::cout << "root2 = " << root2 << "\n"
//	      << "root3 = " << root3 << "\n";
//
//  
  std::cout << "Testing quartic solutions\n"
	    << std::setprecision(std::numeric_limits<Iflt>::digits10);

  boost::array<Iflt, 10> array = {1e8,1e4,1,1e-4,1e-8,1e8,-1e4,-1,-1e-4,-1e-8};
  
  for (boost::array<Iflt, 4>::const_iterator aPtr = array.begin();
       aPtr != array.end(); ++aPtr)
    for (boost::array<Iflt, 4>::const_iterator bPtr = array.begin();
	 bPtr != array.end(); ++bPtr)
      for (boost::array<Iflt, 4>::const_iterator cPtr = array.begin();
	   cPtr != array.end(); ++cPtr)
      for (boost::array<Iflt, 4>::const_iterator dPtr = array.begin();
	   dPtr != array.end(); ++dPtr)
	{
	  size_t roots = quarticSolve(*aPtr, *bPtr, *cPtr, *dPtr, hybroot1, hybroot2, hybroot3, hybroot4);
	  
	  printDetails("YacFraid ", yacfraidQuarticSolve, roots, *aPtr, *bPtr, *cPtr, *dPtr);
	  printDetails("Neumark  ", neumarkQuarticSolve, roots, *aPtr, *bPtr, *cPtr, *dPtr);
	  printDetails("Descartes", descartesQuarticSolve, roots, *aPtr, *bPtr, *cPtr, *dPtr);
	  printDetails("Ferrari  ", ferrariQuarticSolve, roots, *aPtr, *bPtr, *cPtr, *dPtr);
	  
	  std::cout << "\n";
	}
  return 0;
}
