/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_colvar_IonDistanceBase_h
#define __PLUMED_colvar_IonDistanceBase_h
#include "Colvar.h"

namespace PLMD{

class NeighborList;

namespace colvar{

class IonDistanceBase : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  std::vector<PLMD::AtomNumber> list_a,list_b,list_c;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  bool invalidateList;
  bool firsttime;
  int  pn, qn, pq, qq, lambda;
  double d0, d1, r0, less_than, sum_exp;
  
public:
  explicit IonDistanceBase(const ActionOptions&);
  ~IonDistanceBase();
// active methods:
  virtual void calculate();
  virtual void prepare();
//  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const=0;
  //virtual double pairing(double arg,double&dfunc, double sw_cutoff, double sw_shift, int p, int q)=0;
  virtual double sw_func(double arg, double sw_cutoff, double sw_shift, int p, int q)=0;
  virtual double dsw_func(double arg, double sw_cutoff, double sw_shift, int p, int q)=0;
  static void registerKeywords( Keywords& keys );
};

}
}
#endif
