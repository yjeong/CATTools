// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"

// Kinematic Reconstruction
#include "CATTools/CatAnalyzer/interface/TTbarFCNCFitter.h"

#include "TH1.h"
#include "TTree.h"

template<class T>
struct bigger_second : std::binary_function<T,T,bool>
{
   inline bool operator()(const T& lhs, const T& rhs)
   {
      return lhs.second > rhs.second;
   }
};
typedef std::pair<int,double> data_t;

using namespace edm;
using namespace reco;
using namespace cat;

class TopAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit TopAnalyzer(const edm::ParameterSet&);
  ~TopAnalyzer() {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  void clear();

  double transverseMass( const reco::Candidate::LorentzVector& lepton, const reco::Candidate::LorentzVector& met);

  edm::EDGetTokenT<cat::GenTopCollection>          genTopToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>    genToken_;
  edm::EDGetTokenT<cat::MuonCollection>            muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>        electronToken_;
  edm::EDGetTokenT<cat::JetCollection>             jetToken_;
  edm::EDGetTokenT<cat::METCollection>             metToken_;
  edm::EDGetTokenT<int>                            pvToken_;
  edm::EDGetTokenT<float>                          puWeight_;
  edm::EDGetTokenT<cat::GenWeights>                genWeightToken_;
  edm::EDGetTokenT<edm::TriggerResults>            triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

  // ---------- CSV weight ------------
  BTagWeightEvaluator csvWeight;

  // ----------member data ---------------------------

  TTree * tree;
  TH1F * tmp;

  int EVENT;
  int RUN;
  int LUMI;
  int events;

  float PUWeight;
  float GenWeight;
  float CSVWeight[19]; // 0 = central, 1-18 = systematics
  int NVertex; 

  double MET;
  double MET_Px;
  double MET_Py;

  const static int kMax = 100;

  int NMuon;

  float Muon_Pt[kMax];
  float Muon_Eta[kMax];
  float Muon_Phi[kMax];
  float Muon_E[kMax];
  float Muon_Iso03[kMax];
  float Muon_Iso04[kMax];
  float Muon_Charge[kMax];

  int NLooseMuon;

  float LooseMuon_Pt[kMax];
  float LooseMuon_Eta[kMax];
  float LooseMuon_Phi[kMax];
  float LooseMuon_E[kMax];
  float LooseMuon_Iso03[kMax];
  float LooseMuon_Iso04[kMax];
  float LooseMuon_Charge[kMax];

  int NElectron;

  float Electron_Pt[kMax];
  float Electron_Eta[kMax];
  float Electron_Phi[kMax];
  float Electron_E[kMax];
  float Electron_Iso03[kMax];
  float Electron_Iso04[kMax];
  float Electron_Charge[kMax];

  int NLooseElectron;

  float LooseElectron_Pt[kMax];
  float LooseElectron_Eta[kMax];
  float LooseElectron_Phi[kMax];
  float LooseElectron_E[kMax];
  float LooseElectron_Iso03[kMax];
  float LooseElectron_Iso04[kMax];
  float LooseElectron_Charge[kMax];


  int NJet;

  const static int Nptw = 6;
//  float NJet_ptw[Nptw];
//  int NJet_ptw50;
//  int NJet_ptw100;
//  int NJet_ptw150;
//  int NJet_ptw200;
//  int NJet_ptw250;
//  int NJet_ptw300;

  float Jet_Pt[kMax];
  float Jet_Eta[kMax];
  float Jet_Phi[kMax];
  float Jet_E[kMax];
  float Jet_partonFlavour[kMax];
  float Jet_hadronFlavour[kMax];
  float Jet_BTag[kMax];
  float Jet_bDiscriminator[kMax];
  float Jet_pfCombinedCvsLJetTags[kMax];
  float Jet_pfCombinedCvsBJetTags[kMax];
  float Jet_Pt_tot;
  float Jet_Pt_Sum[kMax];
  float Jet_Pt_Sum2[kMax];
  float Jet_mass_tot;
  float Jet_mass_sum[kMax];
  float Jet_HT;
  float Jet_dR_min[kMax];
  float Jet_dR_from_W[kMax];
  float Jet_dR_min_with_q[kMax];
  float Jet_dE_with_q[kMax];
  float Jet_dpt_with_q[kMax];
  float W_dR[kMax];
  float W_dE[kMax];
  float W_dpt[kMax];
  float Jet_mW[kMax];
  float Jet_mW_Flavour[kMax];
  int NW;
  float Jet_mt[kMax];
  int Nt;

  float Jet_JES_Up[kMax];
  float Jet_JES_Dw[kMax];

  int csvid[kMax];

  int NBJet;

  int DiLeptonic;
  int SemiLeptonic;
  int TTBJ;
  int TTBB;
  int TTCC;
  int TTJJ;

  int GenNJet20;
  int GenNBJet20;
  int GenNCJet20;
  int GenNAddJet20;
  int GenNAddBJet20;
  int GenNAddCJet20;

  float GenLepton1_Pt;
  float GenLepton1_Eta;
  float GenLepton2_Pt;
  float GenLepton2_Eta;

  float MT_MuonMET[kMax];
  float Phi_MuonMET[kMax];
  float MT_ElectronMET[kMax];
  float Phi_ElectronMET[kMax]; 

  float Kin_Hmass;
  float Kin_HdRbb;
  float Kin_Chi2;
  float Kin_TopMHc;
  float Kin_TopMWb;
  float Kin_Wmass;

  int IsMuonTrig;
  int IsElectronTrig; 
  int IsHadronTrig; 

  int nb;
  int nbbar;
  int nWp;
  int nWm;
  int nq;
  int nl;
  int id_from_W[kMax];
};

TopAnalyzer::TopAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  genTopToken_      = consumes<cat::GenTopCollection>       (iConfig.getParameter<edm::InputTag>("genTopLabel"));
  genToken_      = consumes<reco::GenParticleCollection>    (iConfig.getParameter<edm::InputTag>("genLabel"));
  muonToken_     = consumes<cat::MuonCollection>            (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_ = consumes<cat::ElectronCollection>        (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_      = consumes<cat::JetCollection>             (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_      = consumes<cat::METCollection>             (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_       = consumes<int>                            (iConfig.getParameter<edm::InputTag>("pvLabel"));
  puWeight_      = consumes<float>                          (iConfig.getParameter<edm::InputTag>("puWeight"));
  genWeightToken_  = consumes<cat::GenWeights>              (iConfig.getParameter<edm::InputTag>("genWeightLabel"));
  // Trigger  
  triggerBits_     = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerObjects_  = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjects"));


  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("events", "Tree for Top quark study");
  tmp = fs->make<TH1F>("EventSummary","EventSummary",2,0,2);

  // CSV re-shape 
  csvWeight.initCSVWeight(false, "csvv2");

}

void TopAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  const float mW0 = 80.385;
  const float mt0 = 172.44;
  tmp->Fill(0);

  clear();

  EVENT  = iEvent.id().event();
  RUN    = iEvent.id().run();
  LUMI   = iEvent.id().luminosityBlock();
  events++;

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByToken(genToken_, gens);

  int nb_gen_from_top=0, nbbar_gen_from_top=0, nWp_gen_from_top=0, nWm_gen_from_top=0;
  int nq_gen_from_W=0, nlep_gen_from_W=0;

  TLorentzVector vq[kMax];
  TLorentzVector vW[kMax];
  int nq_gen=0;
  int nW_gen=0;

  for (GenParticleCollection::const_iterator gen=gens->begin(); gen!=gens->end(); gen++) {
    int id = gen->pdgId();
    int ndau = gen->numberOfDaughters();
    //int status = gen->status();
    const Candidate * mom  = gen->mother();
    const Candidate *mom2=0;
    //const Candidate *mom3=0;
    //int momid=0, mom2id=0;//, mom3id=0;
    //if(mom!=0) {momid = mom->pdgId(); mom2 = (*mom).mother();} if(mom!=0 && mom2!=0) mom2id = mom2->pdgId();
    //if(mom2!=0) {mom2id = mom2->pdgId(); /*mom3 = (*mom2).mother();*/} //if(mom2!=0 && mom3!=0) mom3id = mom3->pdgId();
    ////if(mom3!=0) {mom3id = mom3->pdgId();}
    //if(events<=10) {cout<<"ev "<<events<<", "<<mom2id<<", "<<momid<<", "<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }
    const Candidate *mom3=0;
    int momid=0, mom2id=0, mom3id=0;
    if(mom!=0) {momid = mom->pdgId(); mom2 = (*mom).mother();} if(mom!=0 && mom2!=0) mom2id = mom2->pdgId();
    if(mom2!=0) {mom2id = mom2->pdgId(); mom3 = (*mom2).mother();} if(mom2!=0 && mom3!=0) mom3id = mom3->pdgId();
    if(mom3!=0) {mom3id = mom3->pdgId();}
    if(events<=10) {cout<<"ev "<<events<<", "<<mom3id<<", "<<mom2id<<", "<<momid<<", "<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }
//    if(EVENT==3357031) {cout<<"ev "<<events<<", "<<mom3id<<", "<<mom2id<<", "<<momid<<", "<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }

//    if((abs(id)==6 || abs(id)==24) && ndau>=2) {
    if((abs(id)==6 || abs(id)==24) && ndau==2) {
      for(int i=0;i<ndau;i++) {
        int daughter_id = gen->daughter(i)->pdgId();

        if(id==6&&daughter_id==5)    nb_gen_from_top++;
        if(id==-6&&daughter_id==-5)  nbbar_gen_from_top++;
        if(id==6&&daughter_id==24)   nWp_gen_from_top++;
        if(id==-6&&daughter_id==-24) nWm_gen_from_top++;

        if(abs(id)==6&&abs(daughter_id)==5) {
          vq[nq_gen].SetPtEtaPhiE(gen->daughter(i)->pt(),gen->daughter(i)->eta(),gen->daughter(i)->phi(),gen->daughter(i)->energy());
          nq_gen++;
        }

        if(abs(id)==6&&abs(daughter_id)==24) {
          vW[nq_gen].SetPtEtaPhiE(gen->daughter(i)->pt(),gen->daughter(i)->eta(),gen->daughter(i)->phi(),gen->daughter(i)->energy());
          nW_gen++;
        }

        if(abs(id)==24&& (abs(daughter_id)>=1&&abs(daughter_id)<=8)) {
          id_from_W[nq_gen_from_W] = daughter_id;
          nq_gen_from_W++;
          vq[nq_gen].SetPtEtaPhiE(gen->daughter(i)->pt(),gen->daughter(i)->eta(),gen->daughter(i)->phi(),gen->daughter(i)->energy());
          nq_gen++;
        }
        if(abs(id)==24&&abs(daughter_id)==11) nlep_gen_from_W++;
        if(abs(id)==24&&abs(daughter_id)==13) nlep_gen_from_W++;
        if(abs(id)==24&&abs(daughter_id)==15) nlep_gen_from_W++;
        
      }
    }
//    cout<<"ev "<<events<<", "<<mom3id<<", "<<mom2id<<", "<<momid<<", "<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; //cout<<endl;
//    cout<<"nb "<<nb_gen_from_top<<", nbbar "<<nbbar_gen_from_top<<", nWp "<<nWp_gen_from_top<<", nWm "<<nWm_gen_from_top<<endl;
  }
  //cout<<"ev "<<events<<", nb "<<nb_gen_from_top<<", nbbar "<<nbbar_gen_from_top<<", nW+ "<<nWp_gen_from_top<<", nW- "<<nWm_gen_from_top<<" from top // nq "<<nq_gen_from_W<<", nlep "<<nlep_gen_from_W<<" from W"<<endl;
  nb = nb_gen_from_top;
  nbbar = nbbar_gen_from_top;
  nWp = nWp_gen_from_top;
  nWm = nWm_gen_from_top;
  nq = nq_gen_from_W;
  nl = nlep_gen_from_W;

  if((nb+nbbar+nWp+nWm)!=4&&(nb+nbbar+nWp+nWm)!=8)
    cout<<"ev "<<events<<", nb "<<nb<<". nbbar "<<nbbar<<", nWp "<<nWp<<", nWm "<<nWm<<endl;

/*
  if((nb+nbbar+nWp+nWm)!=4&&(nb+nbbar+nWp+nWm)!=8)
  for (GenParticleCollection::const_iterator gen=gens->begin(); gen!=gens->end(); gen++) {
    int id = gen->pdgId();
    int ndau = gen->numberOfDaughters();
    const Candidate * mom  = gen->mother();
    const Candidate *mom2=0;
    const Candidate *mom3=0;
    int momid=0, mom2id=0, mom3id=0;
    if(mom!=0) {momid = mom->pdgId(); mom2 = (*mom).mother();} if(mom!=0 && mom2!=0) mom2id = mom2->pdgId();
    if(mom2!=0) {mom2id = mom2->pdgId(); mom3 = (*mom2).mother();} if(mom2!=0 && mom3!=0) mom3id = mom3->pdgId();
    if(mom3!=0) {mom3id = mom3->pdgId();}
    cout<<"ev "<<events<<", "<<mom3id<<", "<<mom2id<<", "<<momid<<", "<<id<<" : "; 
    for(int i=0;i<ndau;i++) 
      cout<<gen->daughter(i)->pdgId()<<", "; 
    cout<<"nb "<<nb<<". nbbar "<<nbbar<<", nWp "<<nWp<<", nWm "<<nWm<<endl;
    cout<<endl;
  }
*/

/*
  for (unsigned int i = 0; i < gens->size() ; i++) {
    //const cat::Muon & muon = muons->at(i);
    const reco::Candidate & gen = gens->at(i);
    int id = gen.pdgId();
    int ndau = gen.numberOfDaughters();
    cout<<"ev "<<events<<", pdgid "<<id<<", ndau "<<ndau<<endl;
//    const reco::Candidate * mom = gen->mother();
//    cout<<"ev "<<EVENT<<", momid "<<mom->pdgId()<<endl;
//    const reco::Candidate & mom = gen.mother();
//    cout<<"ev "<<EVENT<<", momid "<<mom.pdgId()<<endl;
    cout<<"ev "<<events<<", momid "<<gen.mother()->pdgId()<<endl;
//    cout<<"ev "<<EVENT<<", pdgid "<<gen.pdgId()<<", momid "<<mom->pdgId()<<endl;
  }
*/

  edm::Handle<cat::GenTopCollection> genTops;
  iEvent.getByToken(genTopToken_, genTops);

  if( genTops.isValid() ) {
    cat::GenTop catGenTop = genTops->at(0);
    DiLeptonic = catGenTop.diLeptonic(0);
    SemiLeptonic = catGenTop.semiLeptonic(0);

    GenNJet20 = catGenTop.NJets20();
    GenNBJet20 = catGenTop.NbJets20();
    GenNCJet20 = catGenTop.NcJets20();
    GenNAddJet20 = catGenTop.NaddJets20();
    GenNAddBJet20 = catGenTop.NaddbJets20();
    GenNAddCJet20 = catGenTop.NaddcJets20();

    if( GenNAddBJet20 == 1 ) TTBJ = 1;
    if( GenNAddBJet20 >= 2 ) TTBB = 1;
    if( GenNAddCJet20 >= 2 ) TTCC = 1;
    if( GenNAddJet20 >= 2 ) TTJJ = 1;

    if( catGenTop.lepton1().pt() > 0){
      GenLepton1_Pt = catGenTop.lepton1().pt();
      GenLepton1_Eta = catGenTop.lepton1().eta();
    }

    if( catGenTop.lepton2().pt() > 0){
      GenLepton2_Pt = catGenTop.lepton2().pt();
      GenLepton2_Eta = catGenTop.lepton2().eta();
    }

  }

  if(!iEvent.isRealData()) {
    edm::Handle<float> PileUpWeight;
    iEvent.getByToken(puWeight_, PileUpWeight);
    PUWeight = *PileUpWeight;

    edm::Handle<cat::GenWeights> genWeight;
    iEvent.getByToken(genWeightToken_, genWeight);
    GenWeight = genWeight->genWeight();
    tmp->Fill(1,GenWeight);
  }

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);
  //we don't use triggerObjects here: can be removed. 

  bool PassMuonTrigger = (trigHelper.triggerFired("HLT_IsoMu24_v") || trigHelper.triggerFired("HLT_IsoTkMu24_v"));
  if(PassMuonTrigger) IsMuonTrig = 1;
  bool PassElectronTrigger = (trigHelper.triggerFired("HLT_Ele32_eta2p1_WPTight_Gsf_v"));
  if(PassElectronTrigger) IsElectronTrig = 1;
  bool PassHadronTrigger = 0;
  //if(!iEvent.isRealData()) PassHadronTrigger = (trigHelper.triggerFired("HLT_PFHT450_SixJet40_PFBTagCSV") || trigHelper.triggerFired("HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV"));
  if(iEvent.isRealData()) PassHadronTrigger = (trigHelper.triggerFired("HLT_PFHT450_SixJet40_PFBTagCSV0p72_v") || trigHelper.triggerFired("HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v"));
  if(!iEvent.isRealData()) PassHadronTrigger = (trigHelper.triggerFired("HLT_PFHT450_SixJet40_BTagCSV_p056_v") || trigHelper.triggerFired("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v"));
  //if(iEvent.isRealData()) PassHadronTrigger = (trigHelper.triggerFired("HLT_PFHT450_SixJet40") || trigHelper.triggerFired("HLT_PFHT400_SixJet30"));
  if(PassHadronTrigger) IsHadronTrig = 1;


  if(events==1){
    cout<<"HLT paths : "<<endl;
    for(unsigned int i=0;i<triggerBits->size();i++){
      const std::string triggerName = triggerNames.triggerName(i);
      cout<<triggerName<<endl;
    }
  }


  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );

  NVertex = *pvHandle;

  Handle<cat::METCollection> METHandle;
  iEvent.getByToken(metToken_, METHandle);

  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  Handle<cat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  MET     = METHandle->begin()->pt();
  MET_Px  = METHandle->begin()->px();
  MET_Py  = METHandle->begin()->py();

  int nmuons = 0;
  int nloosemuons = 0;
  for (unsigned int i = 0; i < muons->size() ; i++) {
    const cat::Muon & muon = muons->at(i);

    bool passLooseMuon = muon.pt() > 15 && fabs(muon.eta()) < 2.4 && muon.isLooseMuon();  
    //bool passID = muon.pt() > 30 && fabs(muon.eta()) < 2.1 && muon.isTightMuon(); 

    //if( !passLooseMuon ) continue;
    if( passLooseMuon ) {

      LooseMuon_Pt[nloosemuons] = muon.pt();
      LooseMuon_Eta[nloosemuons] = muon.eta();
      LooseMuon_Phi[nloosemuons] = muon.phi();
      LooseMuon_E[nloosemuons] = muon.energy();
      LooseMuon_Iso03[nloosemuons] = muon.relIso(0.3);
      LooseMuon_Iso04[nloosemuons] = muon.relIso(0.4);
      LooseMuon_Charge[nloosemuons] = muon.charge();
      nloosemuons++;
    }

    bool passTightMuon = muon.pt() > 30 && fabs(muon.eta()) < 2.1 && muon.isTightMuon() && muon.relIso(0.4) < 0.15;

    if( passTightMuon ) {

      Muon_Pt[nmuons] = muon.pt(); 
      Muon_Eta[nmuons] = muon.eta(); 
      Muon_Phi[nmuons] = muon.phi(); 
      Muon_E[nmuons] = muon.energy();
      Muon_Iso03[nmuons] = muon.relIso(0.3);
      Muon_Iso04[nmuons] = muon.relIso(0.4);
      Muon_Charge[nmuons] = muon.charge();

      MT_MuonMET[nmuons] = transverseMass( muon.p4(), METHandle->begin()->p4() );
      Phi_MuonMET[nmuons] = fabs(deltaPhi( muon.phi(), METHandle->begin()->p4().phi()));
      nmuons++;
    }
  }

  NMuon = nmuons;
  NLooseMuon = nloosemuons;
  int nelectrons = 0;
  int nlooseelectrons = 0;
  for (unsigned int i = 0; i < electrons->size() ; i++) {
    const cat::Electron & electron = electrons->at(i);

    bool passLooseElectron = electron.pt() > 15 && fabs(electron.eta()) < 2.4 && electron.electronID("cutBasedElectronID-Summer16-80X-V1-veto") > 0;
    //bool passID = electron.pt() > 30 && fabs(electron.eta()) < 2.1 && electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") > 0;

    //if ( !passLooseElectron ) continue;
    if ( passLooseElectron ) {

      LooseElectron_Pt[nlooseelectrons] = electron.pt();
      LooseElectron_Eta[nlooseelectrons] = electron.eta();
      LooseElectron_Phi[nlooseelectrons] = electron.phi();
      LooseElectron_E[nlooseelectrons] = electron.energy();
      LooseElectron_Iso03[nlooseelectrons] = electron.relIso(0.3);
      LooseElectron_Iso04[nlooseelectrons] = electron.relIso(0.4);
      LooseElectron_Charge[nlooseelectrons] = electron.charge();
      nlooseelectrons++;
    }

    bool passTightElectron = electron.pt() > 30 && fabs(electron.eta()) < 2.1 && electron.electronID("cutBasedElectronID-Summer16-80X-V1-tight") > 0;
    //bool passIso = electron.relIso() < 0.12;

    //if ( !passTightElectron ) continue;
    if ( passTightElectron ) {

      Electron_Pt[nelectrons] = electron.pt();
      Electron_Eta[nelectrons] = electron.eta();
      Electron_Phi[nelectrons] = electron.phi();
      Electron_E[nelectrons] = electron.energy();
      Electron_Iso03[nelectrons] = electron.relIso(0.3);
      Electron_Iso04[nelectrons] = electron.relIso(0.4);
      Electron_Charge[nelectrons] = electron.charge();

      MT_ElectronMET[nelectrons] = transverseMass( electron.p4(), METHandle->begin()->p4() );
      Phi_ElectronMET[nelectrons] = fabs(deltaPhi( electron.phi(), METHandle->begin()->p4().phi()));
      nelectrons++;
    }

  }

  NElectron = nelectrons;
  NLooseElectron = nlooseelectrons;

  // CSV re-shape 
  // Initialize SF_btag
  float Jet_SF_CSV[19];
  for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

  //for CSV ordering
  std::map<int,double> mapJetBDiscriminator;

  int nJets = 0;
  float nJets_ptw[Nptw] = {0,};
//  float nJets_ptw50 = 0;
//  float nJets_ptw100 = 0;
//  float nJets_ptw150 = 0;
//  float nJets_ptw200 = 0;
//  float nJets_ptw250 = 0;
//  float nJets_ptw300 = 0;
  int nbJets = 0;

  std::vector<cat::Jet> selectedJets;
  TLorentzVector VJet_sum[12], VJet_tot;

  for (unsigned int i = 0; i < jets->size() ; i++) {

    const cat::Jet & jet = jets->at(i);

    bool pass = std::abs(jet.eta()) < 2.4 && jet.pt() > 30 && jet.LooseId() ;
    //if (!pass ) continue; 

    if(pass) {
      double dr = 999.9;
      TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());
      VJet_tot += vjet;
      if(nJets<12) {
        VJet_sum[nJets] = VJet_tot;
        //Jet_mass_sum[nJets] = VJet[nJets].M();
        //Jet_mass_sum[nJets] = VJet_tot.M();
        Jet_mass_sum[nJets] = VJet_sum[nJets].M();
      }

      double dR_min = 10000;
      for(int j = 0 ; j < nq_gen ; j++){
        double dR_temp = vq[j].DeltaR(vjet);
        if(dR_min>dR_temp) {
          dR_min=dR_temp;
          Jet_dE_with_q[nJets] = (vjet.E() - vq[j].E())/vjet.E();
          Jet_dpt_with_q[nJets] = (vjet.Pt() - vq[j].Pt())/vjet.Pt();
        }
      }
      Jet_dR_min_with_q[nJets] = dR_min;

      for(int j = 0 ; j < NMuon ; j++){ 
        TLorentzVector vlep;
        vlep.SetPtEtaPhiE(Muon_Pt[j], Muon_Eta[j], Muon_Phi[j], Muon_E[j]);
        dr = vjet.DeltaR(vlep);
        if( dr < 0.4 ) break;
      }
      //if( dr < 0.4 ) continue;

      for(int j = 0 ; j < NElectron ; j++){
        TLorentzVector vlep;
        vlep.SetPtEtaPhiE(Electron_Pt[j], Electron_Eta[j], Electron_Phi[j], Electron_E[j]);
        dr = vjet.DeltaR(vlep);
        if( dr < 0.4 ) break; 
      }
      //if( dr < 0.4) continue;

      Jet_Pt[nJets] = jet.pt();
      Jet_Eta[nJets] = jet.eta();
      Jet_Phi[nJets] = jet.phi();
      Jet_E[nJets] = jet.energy();

      Jet_partonFlavour[nJets] = jet.partonFlavour();
      Jet_hadronFlavour[nJets] = jet.hadronFlavour();

      Jet_JES_Up[nJets] = jet.shiftedEnUp();
      Jet_JES_Dw[nJets] = jet.shiftedEnDown();

      Jet_Pt_tot += jet.pt();
      Jet_HT += jet.et();

      for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] *= csvWeight.getSF(jet, iu);

      double bDiscriminator = jet.bDiscriminator(BTAG_CSVv2);
      double pfCombinedCvsLJetTags = jet.bDiscriminator("pfCombinedCvsLJetTags");
      double pfCombinedCvsBJetTags = jet.bDiscriminator("pfCombinedCvsBJetTags");

      Jet_bDiscriminator[nJets] = bDiscriminator;
      Jet_pfCombinedCvsLJetTags[nJets] = pfCombinedCvsLJetTags;
      Jet_pfCombinedCvsBJetTags[nJets] = pfCombinedCvsBJetTags;

      if( bDiscriminator > WP_BTAG_CSVv2M) {
        nbJets++;
        Jet_BTag[nJets] = 1;
      }else{
        Jet_BTag[nJets] = 0;
      }

      mapJetBDiscriminator[nJets] = bDiscriminator;

      selectedJets.push_back( jet ); 

      nJets++;
      //nJets_ptw50 += jet.pt()/50;
      for(int np=0;np<Nptw;np++) nJets_ptw[np] += jet.pt()/(50.*(np+1.));
    }
  }

  NJet = nJets;
  NBJet = nbJets;

//  for(int np=0;np<Nptw;np++) NJet_ptw[np] = nJets_ptw[np];

  if(NJet>=1) {
    Jet_Pt_Sum[0] = Jet_Pt[0];
    if(NJet>=12) Jet_Pt_Sum2[0] = Jet_Pt[11];
    else Jet_Pt_Sum2[0] = 0;

    for(int i=1;i<12;i++) {
      if(i<NJet) Jet_Pt_Sum[i] = Jet_Pt_Sum[i-1] + Jet_Pt[i];
      else Jet_Pt_Sum[i] = Jet_Pt_Sum[NJet-1];

      if(NJet>=12) Jet_Pt_Sum2[i] = Jet_Pt_Sum2[i-1] + Jet_Pt[11-i];
      if(NJet<12) {
        if(i<(12-NJet)) Jet_Pt_Sum2[i] = 0;
        else Jet_Pt_Sum2[i] = Jet_Pt_Sum2[i-1] + Jet_Pt[11-i];
      }
    }

//    cout<<"NJet "<<NJet<<endl;
//    for(int i=0;i<12;i++) cout<<Jet_Pt_Sum[i]<<", "; cout<<endl;
//    for(int i=0;i<12;i++) cout<<Jet_Pt_Sum2[i]<<", "; cout<<endl;
  }

  Jet_mass_tot = VJet_tot.M();
  if(NJet<12) {
    for(int nj=NJet;nj<12;nj++) {
      Jet_mass_sum[nj] = Jet_mass_sum[NJet-1];
    }
  }

  //const int NJET = NJet;
//  const static int NJET = NJet;
//  int num_mW_final[NJET]={0,};
  int num_mW_final[kMax]={0,};
//  int num_temp[NJet]={0,};
//  int num_temp[kMax]={0,};

  if(NJet>=2) {
    //for (unsigned int i = 0; i < selectedJets.size() ; i++) {
    for (unsigned int i = 0; i < selectedJets.size() ; i++) if(Jet_BTag[i]==0) {
      double dR_min = 10000;
      double mW_final = 10000;
      double dmW = 10000;
      double mW_final_Flavour = 10000;
      double dmW_Flavour = 10000;
      //int num_mW_final = 0;
      TLorentzVector j1(selectedJets[i].px(), selectedJets[i].py(), selectedJets[i].pz(), selectedJets[i].energy());
      for (unsigned int j = 0; j < selectedJets.size() ; j++) {
        //if(i!=j) {
        if(i!=j && Jet_BTag[j]==0) {
          TLorentzVector j2(selectedJets[j].px(), selectedJets[j].py(), selectedJets[j].pz(), selectedJets[j].energy());
          double dR = j1.DeltaR(j2);
          if(dR_min>dR) dR_min = dR;

          double mW = (j1+j2).M();
          double dmW_temp = fabs(mW0 - mW);
          if(dmW>dmW_temp) {
            dmW = dmW_temp;
            mW_final = mW;
            //num_mW_final = j;
            num_mW_final[i] = j;
            //num_temp[i] = j;
          }
          if(dmW_Flavour>dmW_temp && Jet_partonFlavour[i]==Jet_partonFlavour[j]) {
            dmW_Flavour = dmW_temp;
            mW_final_Flavour = mW;
          }
        }
      }
      Jet_dR_min[i] = dR_min;
      Jet_mW[i] = mW_final;
      Jet_mW_Flavour[i] = mW_final_Flavour;
      if(fabs(Jet_mW[i] - mW0)<20) NW++;

      //int k = num_mW_final;
      //k = num_temp[i];
      int k = num_mW_final[i];
      TLorentzVector j3(selectedJets[k].px(), selectedJets[k].py(), selectedJets[k].pz(), selectedJets[k].energy());
      Jet_dR_from_W[i] = j1.DeltaR(j3);

      double dR_W_min = 10000;
      double dE_W = 0;
      double dpt_W = 0;
      for(int nw=0;nw<nW_gen;nw++) {
        double dR_W = (j1+j3).DeltaR(vW[nw]);
        if(dR_W_min > dR_W) {
          dR_W_min = dR_W;
          dE_W = ((j1+j3).E() - vW[nw].E())/(j1+j3).E();
          dpt_W = ((j1+j3).Pt() - vW[nw].Pt())/(j1+j3).Pt();
        }
      }
      W_dR[i] = dR_W_min;
      W_dE[i] = dE_W;
      W_dpt[i] = dpt_W;
//    cout<<dR_min<<", ";
//      cout<<mW_final<<", ";
    }
    NW = NW/2;
  }
//  cout<<endl;


  //int num_mt_final[kMax] = {0,};

  if(NJet>=3) {
    for (unsigned int i = 0; i < selectedJets.size() ; i++) {
      if(Jet_BTag[i]) {
        double dmt = 10000;
        double mt_final = 10000;
        TLorentzVector j1(selectedJets[i].px(), selectedJets[i].py(), selectedJets[i].pz(), selectedJets[i].energy());
        for (unsigned int j = 0; j < selectedJets.size() ; j++) {
          //if(i!=j) {
          if(i!=j && Jet_BTag[j]==0) {
            TLorentzVector j2(selectedJets[j].px(), selectedJets[j].py(), selectedJets[j].pz(), selectedJets[j].energy());
            int k = num_mW_final[j];
            TLorentzVector j3(selectedJets[k].px(), selectedJets[k].py(), selectedJets[k].pz(), selectedJets[k].energy());
            double mt = (j1+j2+j3).M();
            double dmt_temp = fabs(mt0 - mt);
            if(dmt>dmt_temp) {
              dmt = dmt_temp;
              mt_final = mt;
              //num_mt_final[i] = j;
            }
          }
        }
        Jet_mt[i] = mt_final;
        //if(Jet_mt[i]>140&&Jet_mt[i]<210) Nt++;
        if(Jet_mt[i]>130&&Jet_mt[i]<250) Nt++;
      }
    }
  }

/*
  float Jet_Pt_temp[NJet]={0,};
  for(int i=0;i<NJet;i++) {
//    cout<<i<<" "<<Jet_Pt[i]<<", ";
    Jet_Pt_temp[i] = Jet_Pt[i];
  }
//  cout<<endl;

  sort(Jet_Pt_temp,Jet_Pt_temp+NJet);
  int Jet_Pt_index[NJet]={0,};

  for(int i=0;i<NJet;i++) {
//    cout<<i<<" "<<Jet_Pt_temp[i]<<", ";
    for(int j=0;j<NJet;j++) {
      if(Jet_Pt_temp[i]==Jet_Pt[j]) Jet_Pt_index[i]=j;
    }
  }
//  cout<<endl;
//  for(int i=0;i<NJet;i++) cout<<Jet_Pt_index[i]<<", "; cout<<endl;
*/

/*
  for (unsigned int iu=0; iu<19; iu++) CSVWeight[iu] = Jet_SF_CSV[iu];

  //csv order
  std::vector< std::pair<int,double> > vecJetBDisc(mapJetBDiscriminator.begin(), mapJetBDiscriminator.end());
  std::sort(vecJetBDisc.begin(), vecJetBDisc.end(), bigger_second<data_t>());
  int ncsvid = 0;
  for( std::vector< std::pair<int,double> >::iterator it = vecJetBDisc.begin() ; it != vecJetBDisc.end(); ++it){
    csvid[ncsvid] = (*it).first;
    ncsvid++;
  }

  //---------------------------------------------------------------------------
  // Kinematic Reconstruction
  //---------------------------------------------------------------------------
  TLorentzVector Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit;
  Kinnu.SetPtEtaPhiE(0,0,0,0);
  Kinblrefit.SetPtEtaPhiE(0,0,0,0);
  Kinbjrefit.SetPtEtaPhiE(0,0,0,0);
  Kinj1refit.SetPtEtaPhiE(0,0,0,0);
  Kinj2refit.SetPtEtaPhiE(0,0,0,0);

  std::vector<int> KinBestIndices;
  KinBestIndices.push_back(-999);
  KinBestIndices.push_back(-999);
  KinBestIndices.push_back(-999);
  KinBestIndices.push_back(-999);
  float bestchi2 = 0;

  if(NJet > 3){

    TLorentzVector leptonp4;
    if( NMuon == 1){
      leptonp4.SetPtEtaPhiE(Muon_Pt[0],Muon_Eta[0],Muon_Phi[0],Muon_E[0]);
    }else if (NElectron == 1){
      leptonp4.SetPtEtaPhiE(Electron_Pt[0],Electron_Eta[0],Electron_Phi[0],Electron_E[0]);
    }

    const cat::MET & catmet = METHandle->at(0);
    bool usebtaginfo = true; 
    fcnc::FindHadronicTop(leptonp4, selectedJets, catmet, usebtaginfo, csvid,  KinBestIndices, bestchi2, Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit);

    if( bestchi2 < 1.0e6 ){

      TLorentzVector Higgs = Kinj1refit + Kinj2refit;
      TLorentzVector TopHc = Higgs + Kinbjrefit ;
      TLorentzVector W = leptonp4 + Kinnu;
      TLorentzVector TopWb = W + Kinblrefit;

      Kin_Hmass = Higgs.M(); 
      Kin_HdRbb = Kinj1refit.DeltaR( Kinj2refit );
      Kin_Chi2 = bestchi2; 
      Kin_TopMHc =TopHc.M();
      Kin_TopMWb =TopWb.M();
      Kin_Wmass = W.M();

    }
  
  }
*/
  
  
  

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

//   if ( (NMuon + NElectron) == 1 ) {
     tree->Fill();
//   }
}


// ------------ method called once each job just before starting event loop  ------------
void TopAnalyzer::beginJob()
{
  events=0;

  tree->Branch("EVENT",&EVENT,"EVENT/i");
  tree->Branch("RUN",&RUN,"RUN/i");
  tree->Branch("LUMI",&LUMI,"LUMI/i");
  tree->Branch("PUWeight",&PUWeight,"PUWeight/F");
  tree->Branch("GenWeight",&GenWeight,"GenWeight/F");
  tree->Branch("CSVWeight",CSVWeight,"CSVWeight_[19]/F");
  tree->Branch("NVertex",&NVertex,"NVertex/i");

  tree->Branch("MET",&MET,"MET/d");
  tree->Branch("MET_Px",&MET_Px,"MET_Px/d");
  tree->Branch("MET_Py",&MET_Py,"MET_Py/d");

  tree->Branch("NMuon",&NMuon,"NMuon/I");
  tree->Branch("Muon_Pt",Muon_Pt,"Muon_Pt[NMuon]/F");
  tree->Branch("Muon_Eta",Muon_Eta,"Muon_Eta[NMuon]/F");
  tree->Branch("Muon_Phi",Muon_Phi,"Muon_Phi[NMuon]/F");
  tree->Branch("Muon_E",Muon_E,"Muon_E[NMuon]/F");
  tree->Branch("Muon_Iso03",Muon_Iso03,"Muon_Iso03[NMuon]/F");
  tree->Branch("Muon_Iso04",Muon_Iso04,"Muon_Iso04[NMuon]/F");
  tree->Branch("Muon_Charge",Muon_Charge,"Muon_Charge[NMuon]/F");

  tree->Branch("NLooseMuon",&NLooseMuon,"NLooseMuon/I");
  tree->Branch("LooseMuon_Pt",LooseMuon_Pt,"LooseMuon_Pt[NLooseMuon]/F");
  tree->Branch("LooseMuon_Eta",LooseMuon_Eta,"LooseMuon_Eta[NLooseMuon]/F");
  tree->Branch("LooseMuon_Phi",LooseMuon_Phi,"LooseMuon_Phi[NLooseMuon]/F");
  tree->Branch("LooseMuon_E",LooseMuon_E,"LooseMuon_E[NLooseMuon]/F");
  tree->Branch("LooseMuon_Iso03",LooseMuon_Iso03,"LooseMuon_Iso03[NLooseMuon]/F");
  tree->Branch("LooseMuon_Iso04",LooseMuon_Iso04,"LooseMuon_Iso04[NLooseMuon]/F");
  tree->Branch("LooseMuon_Charge",LooseMuon_Charge,"LooseMuon_Charge[NLooseMuon]/F");

  tree->Branch("NElectron",&NElectron,"NElectron/I");
  tree->Branch("Electron_Pt",Electron_Pt,"Electron_Pt[NElectron]/F");
  tree->Branch("Electron_Eta",Electron_Eta,"Electron_Eta[NElectron]/F");
  tree->Branch("Electron_Phi",Electron_Phi,"Electron_Phi[NElectron]/F");
  tree->Branch("Electron_E",Electron_E,"Electron_E[NElectron]/F");
  tree->Branch("Electron_Iso03",Electron_Iso03,"Electron_Iso03[NElectron]/F");
  tree->Branch("Electron_Iso04",Electron_Iso04,"Electron_Iso04[NElectron]/F");
  tree->Branch("Electron_Charge",Electron_Charge,"Electron_Charge[NElectron]/F");

  tree->Branch("NLooseElectron",&NLooseElectron,"NLooseElectron/I");
  tree->Branch("LooseElectron_Pt",LooseElectron_Pt,"LooseElectron_Pt[NLooseElectron]/F");
  tree->Branch("LooseElectron_Eta",LooseElectron_Eta,"LooseElectron_Eta[NLooseElectron]/F");
  tree->Branch("LooseElectron_Phi",LooseElectron_Phi,"LooseElectron_Phi[NLooseElectron]/F");
  tree->Branch("LooseElectron_E",LooseElectron_E,"LooseElectron_E[NLooseElectron]/F");
  tree->Branch("LooseElectron_Iso03",LooseElectron_Iso03,"LooseElectron_Iso03[NLooseElectron]/F");
  tree->Branch("LooseElectron_Iso04",LooseElectron_Iso04,"LooseElectron_Iso04[NLooseElectron]/F");
  tree->Branch("LooseElectron_Charge",LooseElectron_Charge,"LooseElectron_Charge[NLooseElectron]/F");

  tree->Branch("NJet",&NJet,"NJet/i");
//  tree->Branch("NJet_ptw",NJet_ptw,"NJet_ptw[6]/F");
//  tree->Branch("NJet_ptw50",&NJet_ptw50,"NJet_ptw50/F");
//  tree->Branch("NJet_ptw100",&NJet_ptw100,"NJet_ptw100/F");
//  tree->Branch("NJet_ptw150",&NJet_ptw150,"NJet_ptw150/F");
//  tree->Branch("NJet_ptw200",&NJet_ptw200,"NJet_ptw200/F");
//  tree->Branch("NJet_ptw250",&NJet_ptw250,"NJet_ptw250/F");
//  tree->Branch("NJet_ptw300",&NJet_ptw300,"NJet_ptw300/F");
  tree->Branch("Jet_Pt",Jet_Pt,"Jet_Pt[NJet]/F");
  tree->Branch("Jet_Pt_tot",&Jet_Pt_tot,"Jet_Pt_tot/F");
//  tree->Branch("Jet_Pt_Sum",Jet_Pt_Sum,"Jet_Pt_Sum[NJet]/F");
//  tree->Branch("Jet_Pt_Sum2",Jet_Pt_Sum2,"Jet_Pt_Sum2[NJet]/F");
//  tree->Branch("Jet_mass_sum",Jet_mass_sum,"Jet_mass_sum[NJet]/F");
  tree->Branch("Jet_Pt_Sum",Jet_Pt_Sum,"Jet_Pt_Sum[12]/F");
  tree->Branch("Jet_Pt_Sum2",Jet_Pt_Sum2,"Jet_Pt_Sum2[12]/F");
  tree->Branch("Jet_mass_sum",Jet_mass_sum,"Jet_mass_sum[12]/F");
  tree->Branch("Jet_mass_tot",&Jet_mass_tot,"Jet_mass_tot/F");
  tree->Branch("Jet_dR_min",Jet_dR_min,"Jet_dR_min[NJet]/F");
  tree->Branch("Jet_dR_from_W",Jet_dR_from_W,"Jet_dR_from_W[NJet]/F");
  tree->Branch("Jet_dR_min_with_q",Jet_dR_min_with_q,"Jet_dR_min_with_q[NJet]/F");
  tree->Branch("Jet_dE_with_q",Jet_dE_with_q,"Jet_dE_with_q[NJet]/F");
  tree->Branch("Jet_dpt_with_q",Jet_dpt_with_q,"Jet_dpt_with_q[NJet]/F");
  tree->Branch("W_dR",W_dR,"W_dR[NJet]/F");
  tree->Branch("W_dE",W_dE,"W_dE[NJet]/F");
  tree->Branch("W_dpt",W_dpt,"W_dpt[NJet]/F");
  tree->Branch("Jet_mW",Jet_mW,"Jet_mW[NJet]/F");
  tree->Branch("Jet_mt",Jet_mt,"Jet_mt[NJet]/F");
  tree->Branch("Jet_mW_Flavour",Jet_mW_Flavour,"Jet_mW_Flavour[NJet]/F");
  tree->Branch("NW",&NW,"NW/I");
  tree->Branch("Jet_mt",Jet_mt,"Jet_mt[NJet]/F");
  tree->Branch("Nt",&Nt,"Nt/I");
  tree->Branch("Jet_HT",&Jet_HT,"Jet_HT/F");
  tree->Branch("Jet_Eta",Jet_Eta,"Jet_Eta[NJet]/F");
  tree->Branch("Jet_Phi",Jet_Phi,"Jet_Phi[NJet]/F");
  tree->Branch("Jet_E",Jet_E,"Jet_E[NJet]/F");
  tree->Branch("Jet_partonFlavour",Jet_partonFlavour,"Jet_partonFlavour[NJet]/F");
  tree->Branch("Jet_hadronFlavour",Jet_hadronFlavour,"Jet_hadronFlavour[NJet]/F");
  tree->Branch("Jet_BTag",Jet_BTag,"Jet_BTag[NJet]/F");
  tree->Branch("Jet_bDiscriminator",Jet_bDiscriminator,"Jet_bDiscriminator[NJet]/F"); 
  tree->Branch("Jet_pfCombinedCvsLJetTags",Jet_pfCombinedCvsLJetTags,"Jet_pfCombinedCvsLJetTags[NJet]/F"); 
  tree->Branch("Jet_pfCombinedCvsBJetTags",Jet_pfCombinedCvsBJetTags,"Jet_pfCombinedCvsBJetTags[NJet]/F"); 

  tree->Branch("Jet_JES_Up",Jet_JES_Up,"Jet_JES_Up[NJet]/F");
  tree->Branch("Jet_JES_Dw",Jet_JES_Dw,"Jet_JES_Dw[NJet]/F");

  tree->Branch("csvid",csvid,"csvid[NJet]/i");

  tree->Branch("NBJet",&NBJet,"NBJet/i");
  tree->Branch("DiLeptonic",&DiLeptonic,"DiLeptonic/i");
  tree->Branch("SemiLeptonic",&SemiLeptonic,"SemiLeptonic/i");

  tree->Branch("TTBJ",&TTBJ,"TTBJ/i");
  tree->Branch("TTBB",&TTBB,"TTBB/i");
  tree->Branch("TTCC",&TTCC,"TTCC/i");
  tree->Branch("TTJJ",&TTJJ,"TTJJ/i");

  tree->Branch("GenNJet20",&GenNJet20, "GenNJet20/i");
  tree->Branch("GenNBJet20",&GenNBJet20, "GenNBJet20/i");
  tree->Branch("GenNCJet20",&GenNCJet20, "GenNCJet20/i");
  tree->Branch("GenNAddJet20",&GenNAddJet20, "GenNAddJet20/i");
  tree->Branch("GenNAddBJet20",&GenNAddBJet20, "GenNAddBJet20/i");
  tree->Branch("GenNAddCJet20",&GenNAddCJet20, "GenNAddCJet20/i");   

  tree->Branch("GenLepton1_Pt",&GenLepton1_Pt, "GenLepton1_Pt/f");
  tree->Branch("GenLepton1_Eta",&GenLepton1_Eta, "GenLepton1_Eta/f");
  tree->Branch("GenLepton2_Pt",&GenLepton2_Pt, "GenLepton2_Pt/f");
  tree->Branch("GenLepton2_Eta",&GenLepton2_Eta, "GenLepton2_Eta/f");

  tree->Branch("MT_MuonMET",MT_MuonMET,"MT_MuonMET[NMuon]/F"); 
  tree->Branch("Phi_MuonMET",Phi_MuonMET,"Phi_MuonMET[NMuon]/F"); 
  tree->Branch("MT_ElectronMET",MT_ElectronMET,"MT_ElectronMET[NElectron]/F"); 
  tree->Branch("Phi_ElectronMET",Phi_ElectronMET,"Phi_ElectronMET[NElectron]/F"); 

  tree->Branch("Kin_Hmass",&Kin_Hmass,"Kin_Hmass/F");
  tree->Branch("Kin_HdRbb",&Kin_HdRbb,"Kin_HdRbb/F");
  tree->Branch("Kin_Chi2",&Kin_Chi2,"Kin_Chi2/F"); 
  tree->Branch("Kin_TopMHc",&Kin_TopMHc,"Kin_TopMHc/F"); 
  tree->Branch("Kin_TopMWb",&Kin_TopMWb,"Kin_TopMWb/F"); 
  tree->Branch("Kin_Wmass",&Kin_Wmass,"Kin_Wmass/F"); 

  tree->Branch("IsMuonTrig",&IsMuonTrig,"IsMuonTrig/i"); 
  tree->Branch("IsElectronTrig",&IsElectronTrig,"IsElectronTrig/i"); 
  tree->Branch("IsHadronTrig",&IsHadronTrig,"IsHadronTrig/i"); 

  tree->Branch("nb",&nb,"nb/i");
  tree->Branch("nbbar",&nbbar,"nbbar/i");
  tree->Branch("nWp",&nWp,"nWp/i");
  tree->Branch("nWm",&nWm,"nWm/i");
  tree->Branch("nq",&nq,"nq/i");
  tree->Branch("nl",&nl,"nl/i");
  //tree->Branch("id_from_W",id_from_W,"id_from_W[nq]/i");
}

void TopAnalyzer::clear()
{
  PUWeight = 1.0;
  GenWeight = 1.0;
  NVertex = -1;
  NMuon = -1;
  NLooseMuon = -1;
  NElectron = -1;
  NLooseElectron = -1;
  NJet = -1;
  NBJet = -1;
  DiLeptonic = -1;
  SemiLeptonic = -1;
  TTBJ = -1; 
  TTBB = -1; 
  TTCC = -1; 
  TTJJ = -1; 

  GenNJet20 = -1; 
  GenNBJet20 = -1;
  GenNCJet20 = -1;
  GenNAddJet20 = -1;
  GenNAddBJet20 = -1;
  GenNAddCJet20 = -1;

  GenLepton1_Pt = -9.0;
  GenLepton1_Eta = -9.0;
  GenLepton2_Pt = -9.0;
  GenLepton2_Eta = -9.0;

  Kin_Hmass = -1.0;
  Kin_HdRbb = -1.0;
  Kin_Chi2 = -1.0;
  Kin_TopMHc = -1.0;
  Kin_TopMWb = -1.0;
  Kin_Wmass = -1.0;

  IsMuonTrig = 0;
  IsElectronTrig = 0;
  IsHadronTrig = 0;

  nb = 0;
  nbbar = 0;
  nWp = 0;
  nWm = 0;
  nq = 0;
  nl = 0;

  Jet_Pt_tot = 0;
  NW = 0;
  Nt = 0;
  Jet_HT = 0;
}

double TopAnalyzer::transverseMass( const reco::Candidate::LorentzVector& lepton,
                                    const reco::Candidate::LorentzVector& met)
{
  reco::Candidate::LorentzVector leptonT(lepton.Px(),lepton.Py(),0.,lepton.E()*sin(lepton.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+met;
  return std::sqrt(sumT.M2());
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TopAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TopAnalyzer);
