/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"

namespace PLMD {
namespace function {

class AltMin : public Function {
private:
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  explicit AltMin( const ActionOptions& ao );
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const ;
  void transformFinalValueAndDerivatives( const std::vector<double>& buf );
};

PLUMED_REGISTER_ACTION(AltMin,"ALT_MIN")

void AltMin::registerKeywords( Keywords& keys ) {
  Function::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

AltMin::AltMin( const ActionOptions& ao ) :
  Action(ao),
  Function(ao)
{
  for(unsigned i=0;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic() ) error("MIN is not a meaningful option for periodic variables");
  }
  parse("BETA",beta);
  log.printf("  value of beta is %f \n",beta);
  addValueWithDerivatives(); checkRead();
}

void AltMin::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  if( args.size()==1 ){
      plumed_dbg_assert( actionInChain() );
      double val=exp(-beta*args[0]); 
      // Numerator 
      addValue( 0, val, myvals ); addDerivative( 0, 0, -beta*val, myvals );
  } else {
      double s=0; std::vector<double> mind( args.size() );
      for(unsigned i=0;i<args.size();++i){
          double val = exp(-beta*args[i]); s += val; mind[i] = -beta*val;
      }
      double fval = - std::log( s ) / beta; addValue( 0, fval, myvals ); 
      if( !doNotCalculateDerivatives() ){
          double pref = -1.0 / (beta*s);
          for(unsigned i=0;i<args.size();++i){
              // Derivatives of min
              addDerivative( 0, i, pref*mind[i], myvals );
          }
      }
  }  
}

void AltMin::transformFinalValueAndDerivatives( const std::vector<double>& buf ) {
  if( !actionInChain() || getNumberOfArguments()>1 ) return;
  Value* val0 = getPntrToComponent(0); double val = val0->get(); 
  double fval = -std::log( val ) / beta; val0->set( fval );
  if( !doNotCalculateDerivatives() ){
      double pref = -1.0 /(beta*val);
      for(unsigned j=0;j<val0->getNumberOfDerivatives();++j){ 
          val0->setDerivative( j, pref*val0->getDerivative(j) );
      }
  }
}

}
}
