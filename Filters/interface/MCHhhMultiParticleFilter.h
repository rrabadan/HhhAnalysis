#ifndef MCHhhMultiParticleFilter_h
#define MCHhhMultiParticleFilter_h
// -*- C++ -*-
//
// Package:    MCHhhMultiParticleFilter
// Class:      MCHhhMultiParticleFilter
// 
/* 

 Description: Filter to select events with an arbitrary number of given particle(s).

 Implementation: derived from MCSingleParticleFilter
     
*/
//
// Original Author:  Paul Lujan
//         Created:  Wed Feb 29 04:22:16 CST 2012
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <HepMC/GenEvent.h>

namespace edm {
      class HepMCProduct;
}

//
// class declaration
//

class MCHhhMultiParticleFilter : public edm::EDFilter {
 public:
  explicit MCHhhMultiParticleFilter(const edm::ParameterSet&);
  ~MCHhhMultiParticleFilter();
  
 private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  void printChildren(const HepMC::GenParticle* genP);
  void printParents(const HepMC::GenParticle* genP);
  // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::HepMCProduct> srctoken_;
  bool debug_;
  int totalEvents_;                // counters
  int passedEvents_;
};
#endif
