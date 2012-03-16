#include "ActionWithDistribution.h"
#include "DistributionFunctions.h"

namespace PLMD {

DistributionFunctionDocs::DistributionFunctionDocs(){
  std::string docs;
  min::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("MIN",docs) );
  sum::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("SUM",docs) );
  mean::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("AVERAGE",docs) );
  less_than::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("LESS_THAN",docs) );
  more_than::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("MORE_THAN",docs) );
  within::writeDocs( docs ); alldocs.insert( std::pair<std::string,std::string>("WITHIN",docs) );
}

void DistributionFunctionDocs::printDocs(const std::string& name){
  if( alldocs.count(name)>0 ) std::cout<<alldocs[name];
}

DistributionFunction::DistributionFunction( const std::vector<std::string>& parameters ):
fine(true)
{
}

DistributionFunction::~DistributionFunction(){
  for(unsigned i=0;i<accumulators.size();++i){ delete thesevalues[i]; delete accumulators[i]; }
}

void DistributionFunction::addAccumulator( const bool wderiv ){
  accumulators.push_back(new Value());
  thesevalues.push_back(new Value());
  hasDeriv.push_back(wderiv);
  plumed_assert( accumulators.size()==thesevalues.size() );
  plumed_assert( accumulators.size()==hasDeriv.size() );
}

void DistributionFunction::setNumberOfDerivatives( const unsigned nder ){
  for(unsigned i=0;i<accumulators.size();++i){
      if(hasDeriv[i]) accumulators[i]->resizeDerivatives(nder);
      else accumulators[i]->resizeDerivatives(0);
  }
}

void DistributionFunction::reset(){
  for(unsigned i=0;i<accumulators.size();++i){
      accumulators[i]->set(0); accumulators[i]->clearDerivatives();
  }
}

void DistributionFunction::copyDerivatives( const unsigned nn, Value* value_in  ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to copy derivatives to an accumulator with no derivatives");

  unsigned nder=value_in->getNumberOfDerivatives(); 
  if( nder!=thesevalues[nn]->getNumberOfDerivatives() ){ thesevalues[nn]->resizeDerivatives(nder); }
  thesevalues[nn]->clearDerivatives();
  for(unsigned i=0;i<value_in->getNumberOfDerivatives();++i) thesevalues[nn]->addDerivative( i, value_in->getDerivative(i) );
}

void DistributionFunction::extractDerivatives( const unsigned nn, Value* value_out  ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to copy derivatives from an accumulator with no derivatives");
  plumed_assert( accumulators[nn]->getNumberOfDerivatives()==value_out->getNumberOfDerivatives() );
  for(unsigned i=0;i<value_out->getNumberOfDerivatives();++i) value_out->addDerivative( i, accumulators[nn]->getDerivative(i) );
}

void DistributionFunction::setValue( const unsigned nn, const double f ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  thesevalues[nn]->set(f);
} 

void DistributionFunction::chainRule( const unsigned nn, const double f ){
  plumed_massert( nn<accumulators.size(), "not enough accumulators in distribution function");
  plumed_massert( hasDeriv[nn], "trying to do chain rule on an accumulator with no derivatives");
  thesevalues[nn]->chainRule(f);
} 

void DistributionFunction::mergeDerivatives( const unsigned kk, ActionWithDistribution& action ){
  for(unsigned i=0;i<accumulators.size();++i){
     accumulators[i]->add( thesevalues[i]->get() );
     if(hasDeriv[i]){ action.mergeDerivatives( kk, thesevalues[i], accumulators[i] ); }
  }
}

unsigned DistributionFunction::requiredBufferSpace() const {
  unsigned nbuf=0;
  for(unsigned i=0;i<accumulators.size();++i){
      nbuf+=1;
      if( hasDeriv[i] ) nbuf+=accumulators[i]->getNumberOfDerivatives();
  }
  return nbuf;
}

void DistributionFunction::copyDataToBuffers( unsigned& bufsize, std::vector<double>& buffers ) const {
  plumed_assert( ( bufsize+requiredBufferSpace() )<=buffers.size() );
  for(unsigned i=0;i<accumulators.size();++i){
      buffers[bufsize]=accumulators[i]->get(); bufsize++;
      if( hasDeriv[i] ){
          for(unsigned j=0;j<accumulators[i]->getNumberOfDerivatives();++j){
              buffers[bufsize]=accumulators[i]->getDerivative(j); bufsize++;
          }
      }
  }
}

void DistributionFunction::retrieveDataFromBuffers( unsigned& bufsize, const std::vector<double>& buffers ){
  plumed_assert( ( bufsize+requiredBufferSpace() )<=buffers.size() );
  for(unsigned i=0;i<accumulators.size();++i){
      accumulators[i]->set( buffers[bufsize] ); bufsize++;
      if( hasDeriv[i] ){
          accumulators[i]->clearDerivatives();
          for(unsigned j=0;j<accumulators[i]->getNumberOfDerivatives();++j){
              accumulators[i]->addDerivative( j, buffers[bufsize] ); bufsize++;
          }
      }
  }
}

}
