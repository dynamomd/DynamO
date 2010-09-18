/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSSGeneralCBT_H
#define CSSGeneralCBT_H
#include <vector>
#include <cmath>
#include "../../base/is_exception.hpp"

template <typename data>
class CSSGeneralCBT
{
private:
  std::vector<unsigned long> CBT;
  std::vector<unsigned long> Leaf;
  std::vector<data> Min;
  unsigned long NP, N;

public:  
  inline size_t size() const { return Min.size(); }
  inline bool empty() const { return Min.empty(); }

  void resize(size_t a)
  {
    clear();
    N = a;
    CBT.resize(2 * N);
    Leaf.resize(N + 1);
    Min.resize(N + 1);
  }

  void clear() 
  {
    CBT.clear();
    Leaf.clear();
    Min.clear();
    N = 0;
    NP = 0;
  }

  void init()
  {
    for (unsigned long i = 1; i <= N; i++)
      Insert(i);  
  }

  inline const data& operator[](size_t a) const 
  {
#ifdef DYNAMO_DEBUG 
    if (Min.empty())
      M_throw() << "Heap not yet sized";
#endif
    
    return Min[a+1]; 
  }
  
  inline data& operator[](size_t a) 
  {
#ifdef DYNAMO_DEBUG 
    if (Min.empty())
      M_throw() << "Heap not yet sized";
#endif

    return Min[a+1]; 
  }

  inline void update(size_t a) { UpdateCBT(a+1); }

  inline data& next_Data() { return Min[CBT[1]]; }

  inline unsigned long next_ID() { return CBT[1] - 1; }

private:
  inline void UpdateCBT(unsigned int i)
  {
    unsigned int f = Leaf[i]/2,l,r,w;
    
    for(; f > 0; f = f / 2) 
      {
	if (CBT[f] != i) break; /* jumps to the next "for" */
	l = CBT[f*2];
	r = CBT[f*2+1];
	if (Min[r] > Min[l] )
	  CBT[f] = l;
	else
	  CBT[f] = r;
      }
    //Walk up finding the winners till it doesn't change or you hit the top of the tree
    for( ; f>0; f=f/2) 
      {
	w = CBT[f]; /* old winner */
	l = CBT[f*2];
	r = CBT[f*2+1];
	if (Min[r] > Min[l])
	  CBT[f] = l;
	else
	  CBT[f] = r;
	if (CBT[f] == w ) return; /* end of the event time comparisons */
      }
  }

  inline void Insert(unsigned int i)
  {
    if (NP == 0) {CBT[1]=i; NP++; return;}
    int j = CBT[NP];
    CBT [NP*2] = j;
    CBT [NP*2+1] = i;
    Leaf[j] = NP*2;
    Leaf[i]= NP*2+1;
    ++NP;
    UpdateCBT (j);
  }
  
  inline void Delete(unsigned int i)
  {
    if (NP < 2) { CBT[1]=0; Leaf[0]=1; --NP; return; }
    int l = NP * 2 - 1;
    if (CBT[l-1] != i)
      {
	Leaf[CBT[l-1]] = l/2;
	CBT[l/2] = CBT[l-1];
	UpdateCBT(CBT[l-1]);
      }
    else 
      {
	Leaf[CBT[l]] = l/2;
	CBT[l/2] =CBT[l];
	UpdateCBT(CBT[l]);
	--NP;
	return;
      }

    if (CBT[l] != i) 
      {
	CBT[Leaf[i]] = CBT[l];
	Leaf[CBT[l]] = Leaf[i];
	UpdateCBT(CBT[l]);
      }
    
    --NP;
  }
};
#endif
