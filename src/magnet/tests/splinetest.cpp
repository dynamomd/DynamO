/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <magnet/math/spline.hpp>
#include <fstream>
#include <magnet/color/transferFunction.hpp>
#include <iostream>


int main(int argc, char *argv[])
{
  using namespace magnet::math;
  Spline spline;
  //Add points to the spline in any order, they're sorted in ascending
  //x later. (If you want to spline a circle you'll need to change the
  //class)
  spline.addPoint(0,        0.0);
  spline.addPoint(40.0/255, 0.0);
  spline.addPoint(60.0/255, 0.2);
  spline.addPoint(63.0/255, 0.05);
  spline.addPoint(80.0/255, 0.0);
  spline.addPoint(82.0/255, 0.9);
  spline.addPoint(1.0, 1.0);

  //spline.addPoint(0, 0);
  //spline.addPoint(0.75, 0.5);
  //spline.addPoint(1, 1);  
  //spline.addPoint(0.2, 0.6);
  //spline.addPoint(0.5, 0.6);

  { //We can extract the original data points by treating the spline as
    //a read-only STL container.
    std::ofstream of("orig.dat");
    for (Spline::const_iterator iPtr = spline.begin();
	 iPtr != spline.end(); ++iPtr)
      of << iPtr->first << " " << iPtr->second << "\n";
  }
  
  { //A "natural spline" where the second derivatives are set to 0 at the boundaries.

    //Each boundary condition is set seperately

    //The following aren't needed as its the default setting. The 
    //zero values are the second derivatives at the spline points.
    spline.setLowBC(Spline::FIXED_2ND_DERIV_BC, 0);
    spline.setHighBC(Spline::FIXED_2ND_DERIV_BC, 0);

    //Note: We can calculate values outside the range spanned by the
    //points. The extrapolation is based on the boundary conditions
    //used.

    std::ofstream of("spline.natural.dat");
    for (double x(-0.2); x <= 1.2001; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

  { //A spline with a fixed first derivative at the boundaries
    
    //The zero values are the value of the gradient at that point
    spline.setLowBC(Spline::FIXED_1ST_DERIV_BC, 0);
    spline.setHighBC(Spline::FIXED_1ST_DERIV_BC, 0);

    std::ofstream of("spline.fixedy1.dat");
    for (double x(-0.2); x <= 1.2001; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

  { //A spline which turns into a parabola at the boundaries.
    
    //This BC does not need extra parameters, so just the condition is
    //enough.
    spline.setLowBC(Spline::PARABOLIC_RUNOUT_BC);
    spline.setHighBC(Spline::PARABOLIC_RUNOUT_BC);

    std::ofstream of("spline.parabolicrunout.dat");
    for (double x(-0.2); x <= 1.2001; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

  { //A spline which turns into a parabola at the boundaries.
    
    //This BC does not need extra parameters, so just the condition is
    //enough.
    spline.setLowBC(Spline::FIXED_1ST_DERIV_BC, 100);
    spline.setHighBC(Spline::PARABOLIC_RUNOUT_BC);

    std::ofstream of("spline.mixed.dat");
    for (double x(-0.2); x <= 1.2001; x += 0.005)
      of << x << " " << spline(x) << "\n";
  }

  {
    magnet::color::TransferFunction tf;
    tf.addKnot(0,        0.91, 0.7, 0.61, 0.0);
    tf.addKnot(40.0/255, 0.91, 0.7, 0.61, 0.0);
    tf.addKnot(60.0/255, 0.91, 0.7, 0.61, 0.2);
    tf.addKnot(63.0/255, 0.91, 0.7, 0.61, 0.05);
    tf.addKnot(80.0/255, 0.91, 0.7, 0.61, 0.0);
    tf.addKnot(82.0/255, 1.0,  1.0, 0.85, 0.9);
    tf.addKnot(1.0,      1.0,  1.0, 0.85, 1.0);
    
    //tf.getColorMap();
  }

}
