#include "HhhAnalysis/Filters/interface/MCHhhMultiParticleFilter.h"

MCHhhMultiParticleFilter::MCHhhMultiParticleFilter(const edm::ParameterSet& iConfig) :
  srctoken_(consumes<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("src", edm::InputTag(std::string("generator"),"unsmeared")))),
  debug_(iConfig.getUntrackedParameter<bool>("debug", false)),
  totalEvents_(0), passedEvents_(0)
{
  //here do whatever other initialization is needed
  
}

MCHhhMultiParticleFilter::~MCHhhMultiParticleFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to skim the data  ------------
bool MCHhhMultiParticleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(srctoken_, evt);
  
  totalEvents_++;
  
  bool htobb = false;
  bool htoWW = false;
  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end(); ++p ) {
      if (debug_){
	  if ((*p)->pdg_id() == 99926 or (*p)->pdg_id()==99927) printChildren(*p); 
	  if ((*p)->pdg_id() == 25) printChildren(*p); 
      }
      if ((*p)->pdg_id() == 25){
	  const HepMC::GenParticle *genP=*p;
	  bool getb1 = false; bool getb2 = false;
	  bool getw1 = false; bool getw2 = false;
	  if (genP->end_vertex()&& (*genP->end_vertex()->particles_begin(HepMC::children))){
	    for ( HepMC::GenVertex::particle_iterator child = genP->end_vertex()->particles_begin(HepMC::children);
		                                  child != genP->end_vertex()->particles_end(HepMC::children);  ++child ){
		if ((*child)->pdg_id()==5) getb1 = true;
		if ((*child)->pdg_id()==-5) getb2 = true;
		if ((*child)->pdg_id()==24) getw1 = true;
		if ((*child)->pdg_id()==-24) getw2 = true;
	    }

	  }
	  htobb = ((getb1 and getb2) or htobb);
	  htoWW = ((getw1 and getw2) or htoWW);
      }
  }
  
  /*
  if (nFound == numRequired_) {
    passedEvents_++;
    return true;
  } else {
    return false;
  }*/
  if (debug_){
      if (htobb and htoWW)
	std::cout <<"h2tohhtoWWbb "<< std::endl;
      else if (htobb and not(htoWW))
	  std::cout <<"only htobb" << std::endl;
      else if (not(htobb) and htoWW)
	  std::cout <<"only htoWW" << std::endl;
      else std::cout <<"no htobb nor htoWW "<< std::endl;
  }
  if (htobb and htoWW){
    passedEvents_++;
    return true;
  }else {
    return false;
  }
  
}

void 
MCHhhMultiParticleFilter::printParents(const HepMC::GenParticle* genP){


   std::cout << "Print the parents  of genP id " << genP->pdg_id() << std::endl;
    genP->print();
   if (genP->production_vertex()&& !genP->is_beam()){
   	for ( HepMC::GenVertex::particle_iterator mother = genP->production_vertex()->particles_begin(HepMC::parents);
		                                  mother != genP->production_vertex()->particles_end(HepMC::parents);  ++mother ){
                std::cout << "its mother "; (*mother)->print();
         }
   }

}



void 
MCHhhMultiParticleFilter::printChildren(const HepMC::GenParticle* genP){

   std::cout << "Print the children of genP id" << genP->pdg_id() << std::endl; 
   genP->print();
   if (genP->end_vertex()&& (*genP->end_vertex()->particles_begin(HepMC::children))){
   	for ( HepMC::GenVertex::particle_iterator child = genP->end_vertex()->particles_begin(HepMC::children);
		                                  child != genP->end_vertex()->particles_end(HepMC::children);  ++child ){
                std::cout << "its child "; (*child)->print();
         }

    }
}


// ------------ method called once each job just after ending the event loop  ------------
void MCHhhMultiParticleFilter::endJob() {
  edm::LogInfo("MCHhhMultiParticleFilter") << "=== Results of MCHhhMultiParticleFilter: passed "
                                        << passedEvents_ << "/" << totalEvents_ << " events" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCHhhMultiParticleFilter);

