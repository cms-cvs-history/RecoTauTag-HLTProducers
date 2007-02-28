// 
// // -*- :    TauDebug
// Class:      TauDebug
// 
/**\class TauDebug TauDebug.cc RecoTauTag/ConeIsolation/test/TauDebug.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Simone Gennai
//         Created:  Wed Apr 12 11:12:49 CEST 2006
// $Id: TauDebug.cc,v 1.2 2006/12/29 09:26:04 gennai Exp $
//
//


// user include files
// system include files
#include <memory>
#include <string>
#include <iostream>
#include <cmath>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/IsolatedTauTagInfo.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/Math/interface/Vector3D.h"

// Math
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include <DataFormats/VertexReco/interface/Vertex.h>

using namespace edm;
using namespace reco; 


class TauDebug : public edm::EDAnalyzer {
public:
  explicit TauDebug(const edm::ParameterSet&);
  ~TauDebug() {}

  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
private:
  int nEventsL1TauJets, nEventsL2TauJets,nEventsL25TauJets,nEventsL3TauJets;
  vector<float> nEventsUsed;
  vector<InputTag> jetTagSrc;
  vector<float> nEventsRiso;
  int nEventsL3TaggedJets,nEventsL25TaggedJets;
  InputTag l1ParticleMap;
  InputTag l2TauJetsSingle;
  InputTag l2TauJetsDouble;

};

TauDebug::TauDebug(const edm::ParameterSet& iConfig)
{
  nEventsL1TauJets = 0;
  nEventsL2TauJets = 0;
  nEventsL25TauJets = 0;
  nEventsL3TauJets = 0;
  nEventsL25TaggedJets=0;
  nEventsL3TaggedJets=0;

  jetTagSrc = iConfig.getParameter< vector<InputTag> >("L2TauJets");
  l1ParticleMap = iConfig.getParameter<InputTag>("L1ParticleMap");
  l2TauJetsSingle = iConfig.getParameter<InputTag>("L2TauJetsSingle");
  l2TauJetsDouble = iConfig.getParameter<InputTag>("L2TauJetsDouble");

  
}

void TauDebug::beginJob()
{


}

void TauDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
 using namespace l1extra;
  Handle<CaloJetCollection> caloJetL2SingleTauHandle;
  Handle<CaloJetCollection> caloJetL2DoubleTauHandle;

  int nL1TauJets =0;
  
  /*
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("pixelVertices",vertices);
  const reco::VertexCollection vertCollection = *(vertices.product());
  VertexCollection::const_iterator ci = vertCollection.begin();
  */
  
  Handle< L1ParticleMapCollection > mapColl ;
  iEvent.getByLabel( l1ParticleMap, mapColl );
  const L1ParticleMap& singleTauMap = ( *mapColl )[L1ParticleMap::kSingleTau ] ;
  const L1ParticleMap& doubleTauMap = ( *mapColl )[L1ParticleMap::kDoubleTau ] ;

  bool singleTauFired = singleTauMap.triggerDecision() ;
  bool doubleTauFired = doubleTauMap.triggerDecision() ;

  if(singleTauFired || doubleTauFired) {
    nEventsL1TauJets++;
    int nL2TauJets = 0;
    
   
    try{ iEvent.getByLabel(l2TauJetsSingle, caloJetL2SingleTauHandle);}catch(...) {;}
    try{ iEvent.getByLabel(l2TauJetsDouble, caloJetL2DoubleTauHandle);}catch(...) {;}
    
    if(caloJetL2SingleTauHandle.isValid()){
      nL2TauJets = nL2TauJets + caloJetL2SingleTauHandle->size();
    }
    if(caloJetL2DoubleTauHandle.isValid()){
      nL2TauJets = nL2TauJets + caloJetL2DoubleTauHandle->size();
    }
    
    if(nL2TauJets >1){
      nEventsL2TauJets++;
      int nTag=0;
      int it=0;
	    
      for( vector<InputTag>::const_iterator s = jetTagSrc.begin(); s != jetTagSrc.end(); ++ s ) {
	Handle<IsolatedTauTagInfoCollection> tauTagInfoHandle;
	try{ iEvent.getByLabel(*s, tauTagInfoHandle);}catch(...) {;}
	
	if(tauTagInfoHandle.isValid())
	  {
	    const IsolatedTauTagInfoCollection & tauTagInfo = *(tauTagInfoHandle.product());
	    //cout << "Found " << tauTagInfo.size() << " Tau candidates" << endl;
	    
	    IsolatedTauTagInfoCollection::const_iterator i = tauTagInfo.begin();
	    
	    for (; i != tauTagInfo.end(); ++i) {
	      //cout <<"Jet Number "<<it<<endl;
	      nTag++;
	      //To compute the efficiency as a function of the isolation cone (Riso)
	      
	      //Prints out some quantities
	      /*
	      const TrackRef leadTk= (i->leadingSignalTrack(0.1, 6.));
	      JetTracksAssociationRef     jetTracks = i->jetRef()->jtaRef();
	      
	      cout <<"PV z "<<ci->z()<<endl;
	      cout <<"Associated Tracks "<<jetTracks->val.size()<< endl;
	      edm::RefVector<reco::TrackCollection>::const_iterator it = jetTracks->val.begin();
	      for(;it!= jetTracks->val.end();it++)
		{
		  
		  cout <<"Track Pt " << (*it)->pt() <<" Chi2 "<<
		    (*it)->normalizedChi2() <<" TIP "<<
		    (*it)->d0() <<" recHits "<<
		    (*it)->recHitsSize() <<" PixelValid Hits "<<
		    (*it)->hitPattern().numberOfValidPixelHits()<<
		    " Zimp "<< (*it)->dz()<<endl; 
		  
		}
	      cout <<"Selected Tracks "<< i->selectedTracks().size()<< endl;
	      
	      if(!leadTk){
		cout <<" Discriminator = 0"<<endl;
		cout <<"No Leading Track "<<endl;
	      }else{
		cout <<"Leading Track pt "<<(*leadTk).pt()<<endl;
		math::XYZVector momentum = (*leadTk).momentum();
		cout <<"Number of SignalTracks = "<<(i->tracksInCone(momentum, 0.07,  1.)).size()<<endl;
		cout <<"New Parameters:signal cone = 0.07 and  isolation cone = 0.45 ";
		cout <<" Discriminator = " << i->discriminator(0.1,0.07,0.45,6.,1.) <<endl;
		if(i->discriminator(0.1,0.07,0.45,6.,1.)) nTag++;
	      }
	      */

	      if(nTag > 1) nEventsL25TauJets++;
	    }
	    
	    
	  }
      }
      
      Handle<JetTagCollection> tauTagHandle;
      try{ iEvent.getByLabel("isolatedL25DoubleTau", tauTagHandle);}catch(...) {;}
      if(tauTagHandle.isValid())
	{
	  const JetTagCollection & tauTagL25 = *(tauTagHandle.product());
	  if(tauTagL25.size() >1) nEventsL25TaggedJets++;
	}
      int nL3Tag=0;
      Handle<JetTagCollection> tauTagL3Handle;
      try{ iEvent.getByLabel("coneIsolationL3DoubleTau", tauTagL3Handle);}catch(...) {;}
      if(tauTagL3Handle.isValid())
	{
	  const JetTagCollection & tauTag = *(tauTagL3Handle.product());
	  if(tauTag.size() >1) nEventsL3TauJets++;
	  
	  JetTagCollection::const_iterator tau = tauTag.begin();
	  for(;tau != tauTag.end();tau++)
	    {
	      if(tau->discriminator() > 0) nL3Tag++;
	    }
	}
      if(nL3Tag > 1) nEventsL3TaggedJets++;



    }
  }
}
void TauDebug::endJob(){
  cout <<"L2TauJets Events "<<nEventsL2TauJets<<endl;
  cout <<"L25TauJets Events "<<nEventsL25TauJets<<endl;
  cout <<"L25TaggedJets Events "<<nEventsL25TaggedJets<<endl;
  cout <<"L3TauJets Events "<<nEventsL3TauJets<<endl;
  cout <<"L3TaggedJets Events "<<nEventsL3TaggedJets<<endl;

  double effL25 = 1.*nEventsL25TaggedJets/ nEventsL2TauJets;
  double errEffL25  = sqrt(effL25*(1-effL25)/nEventsL2TauJets);
  cout <<"L25 Efficiency "<<setprecision(2)<<effL25 <<" +- " <<setprecision(1)<<errEffL25<<endl;

  double effL3 = 1.*nEventsL3TaggedJets/ nEventsL25TaggedJets;
  double errEffL3  = sqrt(effL3*(1-effL3)/nEventsL25TaggedJets);
  cout <<"L3 Efficiency "<<setprecision(2)<<effL3 <<" +- " <<setprecision(1)<<errEffL3<<endl;

  double effHLT = 1.*nEventsL3TaggedJets/ nEventsL2TauJets;
  double errEffHLT  = sqrt(effHLT*(1-effHLT)/nEventsL2TauJets);
  cout <<"HLT Efficiency "<< setprecision(2) << effHLT <<" +- " <<setprecision(1)<<errEffHLT<<endl;

  
}
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TauDebug);
