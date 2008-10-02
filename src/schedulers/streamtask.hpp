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

#ifndef CSStreamTask_H
#define CSStreamTask_H

template<class T = std::vector<CIntEvent> >
class CSStreamTask
{ 
public:
  typedef typename T::iterator T_it;

  CSStreamTask(const T_it nstart, const T_it nend, const Iflt ndt):
    start(nstart), end(nend), dt(ndt)
  {};
  
  void operator()() 
  {
    for (T_it iPtr = start; iPtr != end; iPtr++)
      iPtr->incrementTime(dt);
  };
  
  
private:   
  const T_it start;
  const T_it end;
  const Iflt dt;
};

#endif
