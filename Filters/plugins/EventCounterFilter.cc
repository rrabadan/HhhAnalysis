// -*- C++ -*-
//
// Package:    HhhAnalysis/Filters
// Class:      EventCounterFilter
// 
/**\class EventCounterFilter EventCounterFilter.cc HhhAnalysis/Filters/plugins/EventCounterFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tao Huang
//         Created:  Wed, 29 Mar 2017 23:01:21 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
//
// class declaration
//

class EventCounterFilter : public edm::stream::EDFilter<> {
   public:
      explicit EventCounterFilter(const edm::ParameterSet&);
      ~EventCounterFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      int ievent ;
      //TString histname;
      TH1F *hevent;
      edm::Service< TFileService > fs;

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCounterFilter::EventCounterFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   ievent  = 0;

}


EventCounterFilter::~EventCounterFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EventCounterFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   ievent++;
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
EventCounterFilter::beginStream(edm::StreamID)
{
  hevent = fs->make<TH1F>("hevent_filter", "event counter before any cuts",10,0,10);
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
EventCounterFilter::endStream() {
  hevent->Fill(1, ievent);
}

// ------------ method called when starting to processes a run  ------------
/*
void
EventCounterFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
EventCounterFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
EventCounterFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
EventCounterFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventCounterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventCounterFilter);
