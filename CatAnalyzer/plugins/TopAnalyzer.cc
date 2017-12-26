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

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Kinematic Reconstruction
#include "CATTools/CatAnalyzer/interface/LepJetsFitterFCNH.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"

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
  TH1D *hskim0, *hskim1;
  TH1D *htrig0, *htrig1;

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

  int NLooseLep2;

  int NJet;
  float NJetW; // pt weighted jet multiplicity
  float NJetW2;

  float T; // Thrust of jets
  float T2;
  float T3;
  float T4;
  float T5;
  float T6;

  float DR_jc; // mean dR between jet & center of jets
  float DR_jc_sigma; // standard variance of DR_jc

//  const static int Nptw = 6;
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
  float Jet_Et[kMax];
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
  float Jet_H;
  float HTb;
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
  bool dTau = false;
  int nTau;
  int id_from_W[kMax];
  int Ngp[kMax];
  int Ngp_noc[kMax];
  int Ngp_other[kMax];
  float dRmax_gp[kMax];
  int nqjet;
  int nqjet2;
  float dR_qqjet[kMax];

  int nq_all;
  float mindRqq[kMax];
  int nq_in_a_jet[kMax];

  //const static int NWset = 105;
  int ij[3][105][2];
  //int W4[3][105][4][2];
  int W4[3][200][4][2];
  //int ij[1000][1000][1000];
  //int W4[1000][1000][1000][1000];

  int nw;
  int nw2;
  float MW[4];
  float MWm;
  int nt;
  int nt2;
  float Mt[4];
  float Mtm;
  float M8j;
  float M12j;
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
  hskim0 = fs->make<TH1D>("hskim0","skim info",9,0,9);
  hskim1 = fs->make<TH1D>("hskim1","skim info",9,0,9);
  htrig0 = fs->make<TH1D>("htrig0","trig info",9,0,9);
  htrig1 = fs->make<TH1D>("htrig1","trig info",9,0,9);

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
  if(events%1000==0) cout<<"events "<<events<<endl;

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByToken(genToken_, gens);

  int nb_gen_from_top=0, nbbar_gen_from_top=0, nWp_gen_from_top=0, nWm_gen_from_top=0;
  int nq_gen_from_W=0, nlep_gen_from_W=0, ntau_gen_from_W=0;

  TLorentzVector vq[kMax];
  TLorentzVector vqW[kMax];
  TLorentzVector vW[kMax];
  //TLorentzVector gp[kMax];
  int nq_gen=0;
  int nqW_gen=0;
  int nW_gen=0;
  //int ngp=0;
  int ngp_noc=0;
  const static int NQ=11, Nq=2;
  int NGP[NQ][Nq];
  for(int i=0;i<NQ;i++) for(int j=0;j<Nq;j++) NGP[i][j] = 0;
  int Ngp_temp[kMax]={0,};
  TLorentzVector vGP[kMax][kMax];

  //TLorentzVector vFP[12];
  //for(int i=0;i<12;i++) vFP[i].SetPtEtaPhiE(0,0,0,0);

  //int NGP[NQ][Nq] = {0,};
  //int NGP[11][2]={0,};
  //vector<int> momID;
  //vector<vector<int>> tops(0, vector<int>(2,0));

  for (GenParticleCollection::const_iterator gen=gens->begin(); gen!=gens->end(); gen++) {
    int id = gen->pdgId();
    int ndau = gen->numberOfDaughters();
    int status = gen->status();

    const static int Nmom = 15;
    const Candidate * MOM[Nmom] = {0,};
    int MOMID[Nmom] = {0,};
    MOM[0] = (*gen).mother();
    for(int i=0;i<Nmom-1;i++) {
      if(MOM[i]!=0) {MOMID[i] = MOM[i]->pdgId(); MOM[i+1] = (*MOM[i]).mother();} 
      if(MOM[i]!=0 && MOM[i+1]!=0) MOMID[i+1] = MOM[i+1]->pdgId();
    }

    //if(events<=10) {cout<<"ev "<<events<<", "; for(int i=Nmom-1;i>=0;i--) cout<<MOMID[i]<<", "; cout<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }

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

        if(abs(id)==24&& (abs(daughter_id)>=1&&abs(daughter_id)<=6)) {
          id_from_W[nq_gen_from_W] = daughter_id;
          nq_gen_from_W++;
          vq[nq_gen].SetPtEtaPhiE(gen->daughter(i)->pt(),gen->daughter(i)->eta(),gen->daughter(i)->phi(),gen->daughter(i)->energy());
          nq_gen++;
          vqW[nqW_gen].SetPtEtaPhiE(gen->daughter(i)->pt(),gen->daughter(i)->eta(),gen->daughter(i)->phi(),gen->daughter(i)->energy());
          nqW_gen++;
        }
        if(abs(id)==24&&abs(daughter_id)==11) nlep_gen_from_W++;
        if(abs(id)==24&&abs(daughter_id)==13) nlep_gen_from_W++;
        if(abs(id)==24&&abs(daughter_id)==15) nlep_gen_from_W++;
      }
    }
    if(abs(MOMID[0])==24&&abs(id)==15) ntau_gen_from_W++;
    for(int i=0; i<ndau; i++){
     int daughter_id = gen->daughter(i)->pdgId();
      if(abs(MOMID[0])==24&&abs(id)==15){
	if(abs(daughter_id)==11) dTau = true;
	if(abs(daughter_id)==13) dTau = true;
	//if(abs(daughter_id)==15) dTau = true;
         //cout<<gen->daughter(i)->pdgId()<<endl;
         //cout<<"daughter of tau: "<<dTau<<endl;
      }
     //if(dTau == 1)cout<<"number of tau: "<<ntau_gen_from_W<<endl;
     //if(nlep_gen_from_W == 2)cout<<"number of tau: "<<nq_gen_from_W<<endl;
    }
    bool momid_pass = false;
    bool momid_pass5 = false;
    bool momid_pass0 = false;
    bool momid_pass6 = false;
    //bool momid_pass24 = false;
    //for(int i=0;i<Nmom;i++) if(abs(MOMID[i])<=5) momid_pass = true;
    //for(int i=0;i<Nmom;i++) if(abs(MOMID[i])<=5&&MOMID[i]!=0) momid_pass = true;
    //for(int i=0;i<Nmom;i++) if(abs(MOMID[i])<=5 && MOMID[i]!=0 && abs(MOMID[i])==6) momid_pass = true;
    for(int i=0;i<Nmom;i++) {
      if(abs(MOMID[i])>=1&&abs(MOMID[i])<=5) momid_pass5 = true;
      if(abs(MOMID[i])!=0) momid_pass0 = true;
      if(abs(MOMID[i])==6) momid_pass6 = true;
    }
    if(momid_pass5 && momid_pass0 && momid_pass6) momid_pass = true;
    //if(ndau==0 && status==1 && gen->pt()>10 && gen->charge()!=0 && (abs(momid)<=5||abs(mom2id)<=5||abs(mom3id)<=5||abs(mom4id)<=5||abs(mom5id)<=5) && abs(momid)>0&&abs(mom2id)>0&&abs(mom3id)>0&&abs(mom4id)>0&&abs(mom5id)>0)
    //if(ndau==0 && status==1 && (abs(momid)<=5||abs(mom2id)<=5||abs(mom3id)<=5||abs(mom4id)<=5||abs(mom5id)<=5) && abs(momid)>0&&abs(mom2id)>0&&abs(mom3id)>0&&abs(mom4id)>0&&abs(mom5id)>0)
    if(ndau==0 && status==1 && momid_pass)
    {
      //if(gen->pt()>10 && gen->charge()!=0) {
      //  gp[ngp].SetPtEtaPhiE(gen->px(), gen->py(), gen->pz(), gen->energy());
      //  ngp++;
      //}
      //if(gen->charge()!=0) ngp++;
      ngp_noc++;
      int MOMid=0;
      int MOMID_first_index=0;
      //if(momid!=0)  MOMid = momid;
      //if(mom2id!=0) MOMid = mom2id;
      //if(mom3id!=0) MOMid = mom3id;
      //if(mom4id!=0) MOMid = mom4id;
      //if(mom5id!=0) MOMid = mom5id;
      for(int i=0;i<Nmom;i++) {
        if(abs(MOMID[i])<=5 && MOMID[i]!=0) {
          //if(MOMid==0) MOMID_first_index=i;
          MOMID_first_index=i;
          MOMid = MOMID[i];
        }
      }
//      for(int i=0;i<Nmom;i++) MOMID_first_index = i;
      NGP[MOMid+5][0] = MOMid;
      NGP[MOMid+5][1]++;
      //cout<<"MOMid "<<MOMid<<endl;
      //cout<<"ev "<<events<<", MOMID_first_index "<<MOMID_first_index<<", MOMid "<<MOMid<<endl;
      //if(MOMid!=0 && MOMID_first_index!=0) cout<<"ev "<<events<<", MOMID_first_index "<<MOMID_first_index<<", MOMid "<<MOMid<<endl;
      //for(int i=0;i<nq_gen;i++) {
      //  double dpt_q = -1;
      //  for(int j=Nmom-1;j>=0;j--) if(abs(MOMID[j]<=5)&&MOMID[j]!=0) {
      //    //TLorentzVector Q(MOM[j]->px(), MOM[j]->py(), MOM[j]->pz(), MOM[j]->energy());
      //    dpt_q = vq_gen[i].pt()-MOM[j]->pt();
      //    cout<<dpt_q<<" ";
      //  }
      //  cout<<endl;
      //  if(dpt_q==0) break;
      //}
      //for(int i=Nmom-1;i>=0;i--) if(abs(MOMID[i]<=5)&&MOMID[i]!=0) {
        double dpt_q = -1;
        double dR_q = -1;
        int I = MOMID_first_index;
        //TLorentzVector Q
        TLorentzVector Q(gen->px(), gen->py(), gen->pz(), gen->energy());
        for(int j=0;j<nq_gen;j++) {
          //if(MOM[I]->pt()>0) Q.SetPtEtaPhiE(gen->px(), gen->py(), gen->pz(), gen->energy());
          //if(MOM[I]->pt()>0) dR_q = vq[j].DeltaR(Q);
          if(gen->pt()>0) dR_q = vq[j].DeltaR(Q);
          if(MOMID_first_index<Nmom) dpt_q = vq[j].Pt()-MOM[I]->pt();
          //cout<<dpt_q<<" "<<dR_q<<", ";
          //if(dpt_q==0 || dR_q==0) cout<<"momid "<<MOM[I]->pdgId()<<endl;;
          //if(dpt_q==0 || dR_q==0) Ngp_noc[j]++; 
          if(dpt_q==0 || dR_q==0) {
            vGP[j][Ngp_temp[j]] = Q;
            Ngp_temp[j]++;
            //vFP[j] += Q;
            break;
          }
        }
        //if(dpt_q!=0) cout<<endl;
        //if(!(dpt_q==0 || dR_q==0)) cout<<endl;
        //break;
      //}
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
  nTau = ntau_gen_from_W;
  nq_all = nq_gen;

  //if(ngp_noc>0) cout<<"ev "<<events<<" : ";
  //for(int i=0;i<NQ;i++) {
  //  //if(abs(NGP[i][0])>0) cout<<NGP[i][0]<<" "<<NGP[i][1]<<", ";
  //  if(NGP[i][1]>0) cout<<NGP[i][0]<<" "<<NGP[i][1]<<", ";  
  //}
  //if(ngp_noc>0) cout<<endl;

  double eta_gp[nq_all];
  double phi_gp[nq_all];
  TLorentzVector qjet_gp[nq_all];
  TLorentzVector qjet[nq_all];
  int Ngp_ot[nq_all];
  for(int i=0;i<nq_gen;i++) {
  //  if(i<NQ && NGP[i][1]>0) Ngp_noc[i] = NGP[i][1];
  //  else Ngp_noc[i] = 0;
    Ngp_noc[i] = Ngp_temp[i];
    //if(vFP[i].Pt()>0) cout<<vq[i].Pt()<<" "<<vFP[i].Pt()<<" "<<vq[i].M()<<" "<<vFP[i].M()<<" / ";
    eta_gp[i] = 0;
    phi_gp[i] = 0;
    Ngp_ot[i] = 0;
    //if(Ngp_noc[i]==0) vGP[i][0].SetPtEtaPhiE(0,0,0,0);
    //qjet_gp[i].SetPtEtaPhiE(0,0,0,0);
    qjet[i].SetPtEtaPhiE(0,0,0,0);

    double dR_max=0;
    for(int j=0;j<Ngp_temp[i];j++) {
      double dR_temp=-1;
      qjet[i] += vGP[i][j];
      for(int k=0;k<Ngp_temp[i];k++) if(j!=k) {
        if(vGP[i][j].Pt()>0 && vGP[i][k].Pt()>0) dR_temp = vGP[i][j].DeltaR(vGP[i][k]);
        if(dR_max < dR_temp) {
          dR_max = dR_temp;
          eta_gp[i] = (vGP[i][j].Eta() + vGP[i][k].Eta())/2.;
          phi_gp[i] = (vGP[i][j].Phi() + vGP[i][k].Phi())/2.;
        }
      }
    }
    if(Ngp_temp[i]==1) eta_gp[i] = vGP[i][0].Eta();
    if(Ngp_temp[i]==1) phi_gp[i] = vGP[i][0].Phi();
    dRmax_gp[i] = dR_max;
    qjet_gp[i].SetPtEtaPhiE(1,eta_gp[i],phi_gp[i],1);

    if(Ngp_temp[i]>=1) nqjet++;
    if(Ngp_temp[i]>=1 && qjet[i].Pt()>0) {
      dR_qqjet[nqjet2] = vq[i].DeltaR(qjet[i]);
      nqjet2++;
    }

    //if(qjet[i].Pt()>0) cout<<"ev "<<events<<", "<<vq[i].Pt()-qjet[i].Pt()<<", "<<vq[i].DeltaR(qjet[i])<<endl;
    //if(nl==0 && dR_max>0) cout<<dR_max<<" ";
    for(int j=0;j<nq_gen;j++) if(i!=j) {
      for(int k=0;k<Ngp_temp[j];k++) {
        double dR_temp = 10000;
        if(vGP[j][k].Pt()>0) dR_temp = qjet_gp[i].DeltaR(vGP[j][k]);
        if((dRmax_gp[i]>0.5 && dR_temp<dRmax_gp[i]) || (dRmax_gp[i]<0.5 && dR_temp<0.5)) Ngp_ot[i]++;
      }
    }
    Ngp_other[i] = Ngp_ot[i];
  }

  if((nb+nbbar+nWp+nWm)!=4&&(nb+nbbar+nWp+nWm)!=8)
    cout<<"ev "<<events<<", nb "<<nb<<". nbbar "<<nbbar<<", nWp "<<nWp<<", nWm "<<nWm<<endl;

  
  //for(int i=0;i<nqW_gen/2;i+=2) cout<<vqW[i].DeltaR(vqW[i+1])<<" "; cout<<endl;
  for(int i=0;i<nq_gen;i++) {
    double mindR=10000;
    for(int j=0;j<nq_gen;j++) if(i!=j) {
      if(mindR>vq[i].DeltaR(vq[j])) mindR=vq[i].DeltaR(vq[j]);
    }
    mindRqq[i] = mindR;
//    cout<<mindR<<" ";

//    for (GenParticleCollection::const_iterator gen2=gens->begin(); gen2!=gens->end(); gen2++) {
//      int id = gen2->pdgId();
//      //int ndau = gen->numberOfDaughters();
//      if(abs(id)>8 && abs(id)!=2212 && gen2->pt()>0) {
//        TLorentzVector GP(gen2->px(), gen2->py(), gen2->pz(), gen2->energy());
//        if(GP.DeltaR(vq[i])==0) cout<<"ev "<<events<<", id "<<id<<endl;
 //     }
//    }
  }



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


//  if(events==1){
//    cout<<"HLT paths : "<<endl;
//    for(unsigned int i=0;i<triggerBits->size();i++){
//      const std::string triggerName = triggerNames.triggerName(i);
//      cout<<triggerName<<endl;
//    }
//  }


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

  TLorentzVector vllep[kMax];

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
      vllep[nloosemuons].SetPtEtaPhiE(muon.pt(), muon.eta(), muon.phi(), muon.energy());
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
      vllep[nloosemuons + nlooseelectrons].SetPtEtaPhiE(electron.pt(), electron.eta(), electron.phi(), electron.energy());
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

/*
  std::vector<reco::Jet> GJets;
  TLorentzVector vgj[kMax];
  int ngjet=0;

  for (unsigned int i = 0; i < gjets->size() ; i++) {
    const reco::GenJet & gjet = gjets->at(i);
    //cout<<"emEnergy "<<gjet.emEnergy()<<endl;
    //cout<<"hadEnergy "<<gjet.hadEnergy()<<endl;
    //cout<<"invisibleEnergy "<<gjet.invisibleEnergy()<<endl;
    //cout<<"energy "<<gjet.energy()<<endl;
    //cout<<"emEnergy+hadEnergy "<<gjet.emEnergy()+gjet.hadEnergy()<<endl;
    //cout<<"emEnergy+hadEnergy+invisibleEnergy "<<gjet.emEnergy()+gjet.hadEnergy()+gjet.invisibleEnergy()<<endl;
    //vlep.SetPtEtaPhiE(Muon_Pt[j], Muon_Eta[j], Muon_Phi[j], Muon_E[j]);
    if(gjet.pt()>30 && fabs(gjet.eta())<2.4) {
      GJets.push_back( gjet ); 
      vgj[ngjet].SetPtEtaPhiE(gjet.pt(),gjet.eta(),gjet.phi(),gjet.emEnergy()+gjet.hadEnergy());
      ngjet++;
    }
    TLorentzVector j1(gjet.px(), gjet.py(), gjet.pz(), gjet.energy());
    //double mindM = 1000;
    //double mW_best = 0;
    //for (unsigned int j = 0; j < gjets->size() ; j++) if(i!=j) {
    //  const reco::GenJet & gjet2 = gjets->at(j);
    //  TLorentzVector j2(gjet2.px(), gjet2.py(), gjet2.pz(), gjet2.energy());
    //  if(mindM > fabs(mW0 - (j1+j2).M())) {
    //    mindM = fabs(mW0 - (j1+j2).M());
    //    mW_best = (j1+j2).M();
    //  }
    //}
cout<<"check  99"<<endl;
    std::vector <const GenParticle*> mcparts = gjet.getGenConstituents();
cout<<"check 100"<<endl;
    //for (unsigned j = 0; j < gjet.getGenConstituents().size(); j++) {
    for (unsigned j = 0; j < mcparts.size(); j++) {
cout<<"check 101"<<endl;
      //reco::GenParticle const* ptl = gjet.getGenConstituent(j);
      const GenParticle* mcpart = mcparts[i];
cout<<"check 102"<<endl;
      //reco::Candidate const* ptl = gjet.getGenConstituent(j);
      //const Candidate * mom  = gen->mother();
      //cout<<ptl->pdgId()<<" ";
      //cout<<ptl->pt()<<" ";
      cout<<mcpart->pt()<<" ";
cout<<"check 103"<<endl;
    }
cout<<"check 104"<<endl;
    cout<<endl;
    //cout<<mW_best<<" ";
    //double mindR = 1000;
    //double energy_ratio = 0;
    //for(int i=0;i<nqW_gen;i++) {
    //  if(mindR > j1.DeltaR(vqW[i])) {
    //    mindR = j1.DeltaR(vqW[i]);
    //    energy_ratio = fabs(j1.E() - vqW[i].E())/j1.E();
    //  }
    //}
    //cout<<mindR<<"("<<energy_ratio<<")"<<" ";
  }
  cout<<endl;
*/

  //for(int i=0;i<ngjet;i++) {
  //  double mindR = 1000;
  //  for(int j=0;j<ngjet;j++) if(i!=j) {
  //    if(mindR > vgj[i].DeltaR(vgj[j])) mindR = vgj[i].DeltaR(vgj[j]);
  //  }
  //  cout<<mindR<<" ";
  //}
  //cout<<endl;

  // CSV re-shape 
  // Initialize SF_btag
  float Jet_SF_CSV[19];
  for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSV[iu] = 1.0;

  //for CSV ordering
  std::map<int,double> mapJetBDiscriminator;

  int nJets = 0;
//  float nJets_ptw[Nptw] = {0,};
//  float nJets_ptw50 = 0;
//  float nJets_ptw100 = 0;
//  float nJets_ptw150 = 0;
//  float nJets_ptw200 = 0;
//  float nJets_ptw250 = 0;
//  float nJets_ptw300 = 0;
  int nbJets = 0;
//  int nWJets = 0;

  std::vector<cat::Jet> selectedJets;
  std::vector<cat::Jet> WJets;
  std::vector<cat::Jet> BJets;
  std::vector<cat::Jet> TJets;
  TLorentzVector VJet_sum[12], VJet_tot;
  TLorentzVector vjet[kMax];

  for (unsigned int i = 0; i < jets->size() ; i++) {

    const cat::Jet & jet = jets->at(i);

    //bool pass = std::abs(jet.eta()) < 2.4 && jet.pt() > 30 && jet.LooseId() ;
    bool pass = std::fabs(jet.eta()) < 2.4 && jet.pt() > 30 && jet.LooseId() ;
    //if (!pass ) continue; 

    if(pass) {
      double dr = 999.9;
      //TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());
      vjet[nJets].SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      VJet_tot += vjet[nJets];
      if(nJets<12) {
        VJet_sum[nJets] = VJet_tot;
        //Jet_mass_sum[nJets] = VJet[nJets].M();
        //Jet_mass_sum[nJets] = VJet_tot.M();
        Jet_mass_sum[nJets] = VJet_sum[nJets].M();
      }

      double dR_min = 10000;
      for(int j = 0 ; j < nq_gen ; j++){
        double dR_temp = vq[j].DeltaR(vjet[nJets]);
        if(dR_min>dR_temp) {
          dR_min=dR_temp;
          Jet_dE_with_q[nJets] = (vjet[nJets].E() - vq[j].E())/vjet[nJets].E();
          Jet_dpt_with_q[nJets] = (vjet[nJets].Pt() - vq[j].Pt())/vjet[nJets].Pt();
        }
      }
      Jet_dR_min_with_q[nJets] = dR_min;

      for(int j = 0 ; j < NMuon ; j++){ 
        TLorentzVector vlep;
        vlep.SetPtEtaPhiE(Muon_Pt[j], Muon_Eta[j], Muon_Phi[j], Muon_E[j]);
        dr = vjet[nJets].DeltaR(vlep);
        if( dr < 0.4 ) break;
      }
      //if( dr < 0.4 ) continue;

      for(int j = 0 ; j < NElectron ; j++){
        TLorentzVector vlep;
        vlep.SetPtEtaPhiE(Electron_Pt[j], Electron_Eta[j], Electron_Phi[j], Electron_E[j]);
        dr = vjet[nJets].DeltaR(vlep);
        if( dr < 0.4 ) break; 
      }
      //if( dr < 0.4) continue;

      Jet_Pt[nJets] = jet.pt();
      Jet_Eta[nJets] = jet.eta();
      Jet_Phi[nJets] = jet.phi();
      Jet_E[nJets] = jet.energy();
      Jet_Et[nJets] = jet.et();

      Jet_partonFlavour[nJets] = jet.partonFlavour();
      Jet_hadronFlavour[nJets] = jet.hadronFlavour();

      Jet_JES_Up[nJets] = jet.shiftedEnUp();
      Jet_JES_Dw[nJets] = jet.shiftedEnDown();

      Jet_Pt_tot += jet.pt();
      Jet_HT += jet.et();
      Jet_H += jet.energy();

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
        HTb += jet.et();
      }else{
        Jet_BTag[nJets] = 0;
      }

      mapJetBDiscriminator[nJets] = bDiscriminator;

      selectedJets.push_back( jet ); 

      if(Jet_BTag[nJets]==0&&WJets.size()<8) {
        WJets.push_back( jet ); 
        TJets.push_back( jet ); 
      }
      if(Jet_BTag[nJets]==1&&BJets.size()<4) {
        BJets.push_back( jet ); 
        TJets.push_back( jet ); 
      }

      nJets++;
      //nJets_ptw50 += jet.pt()/50;
      //for(int np=0;np<Nptw;np++) nJets_ptw[np] += jet.pt()/(50.*(np+1.));
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
//  const int NWJet = NJet - NBJet;
//  const int NW_set = 0;

//  TLorentzVector vW[1000];
//  int NvW=0;

  double NjW=0, NjW2=0;
  double ptmin=30, ptmax=125;
  double pt1=ptmin, pt2=0;
  for(int i=NJet-1;i>0;i--) {
    //double pt1 = selectedJets[i].pt();
    //double pt2 = selectedJets[i-1].pt();
    pt2 = selectedJets[i].pt();
    NjW2 += (i+1)*(pt2*pt2 - pt1*pt1)/2.;
    if(pt2<ptmax) {
      //NjW += (i+1)*(pt2*pt2 - pt1*pt1)/2.;
      NjW = NjW2;
    }
    pt1 = selectedJets[i].pt();
  }
  double factor = (ptmax*ptmax - ptmin*ptmin)/2.;
  NjW2 /= factor;
  NjW /= factor;
  NJetW2 = NjW2;
  NJetW = NjW;

  if(NJet>=4) {
    TLorentzVector trijet1 = vjet[0]+vjet[1]+vjet[2];
    double Thrust=0, Thrust2=0;
    double sumpt=0, sumpt2=0;
    for(int i=0;i<NJet;i++) {
      if(i>=3) Thrust += vjet[i].Pt() * fabs(cos( trijet1.DeltaPhi(vjet[i]) ));
      if(i>=3) sumpt += vjet[i].Pt();
      Thrust2 += vjet[i].Pt() * fabs(cos( trijet1.DeltaPhi(vjet[i]) ));
      sumpt2 += vjet[i].Pt();
    }
    if(sumpt>0) Thrust /= sumpt;
    if(sumpt2>0) Thrust2 /= sumpt2;
    T = Thrust;
    T2 = Thrust2;
    //cout<<"ev "<<events<<", T "<<Thrust<<endl;
  }

  double nlep_loose_wo_jet=0;
  for(int i=0;i<(NLooseMuon+NLooseElectron);i++) {
    bool matching = false;
    for(int j=0;j<NJet;j++) if(vllep[i].DeltaR(vjet[j])<0.4) matching = true;
    if(matching==false) nlep_loose_wo_jet++;
  }
  NLooseLep2 = nlep_loose_wo_jet;

  double moment_eta=0, moment_phi=0;
  double total_pt=0;
  for(int i=0;i<NJet;i++) {
    moment_eta = vjet[i].Pt() * vjet[i].Eta();
    moment_phi = vjet[i].Pt() * vjet[i].Phi();
    total_pt += vjet[i].Pt();
  }
  if(total_pt>0) moment_eta /= total_pt;
  if(total_pt>0) moment_phi /= total_pt;
  TLorentzVector COJ; // centor of jets
  COJ.SetPtEtaPhiE(1, moment_eta, moment_phi, 1);

  double Thrust6=0;
  for(int i=0;i<NJet;i++) {
    Thrust6 += vjet[i].Pt() * fabs(cos( COJ.DeltaPhi(vjet[i]) ));
  }
  if(total_pt>0) Thrust6 /= total_pt;
  T6 = Thrust6;

  int NTJet = TJets.size();
  double moment_eta5=0, moment_phi5=0;
  double total_pt5=0;
  for(int i=0;i<NTJet;i++) {
    moment_eta = TJets[i].pt() * TJets[i].eta();
    moment_phi = TJets[i].pt() * TJets[i].phi();
    total_pt5 += TJets[i].pt();
  }
  if(total_pt5>0) moment_eta5 /= total_pt5;
  if(total_pt5>0) moment_phi5 /= total_pt5;
  TLorentzVector COJ5; // centor of jets
  COJ5.SetPtEtaPhiE(1, moment_eta5, moment_phi5, 1);

  double Thrust3=-1;
  double Thrust4=-1;
  double Thrust5=-1;

  double dR_jc=0;
  double dR_jc_sigma=0;

  for(int i=0;i<NTJet;i++) {
    double Thrust=0;
    double sumpt=0;
    TLorentzVector j1(TJets[i].px(), TJets[i].py(), TJets[i].pz(), TJets[i].energy());
    for(int j=0;j<NTJet;j++) if(i!=j) {
      TLorentzVector j2(TJets[j].px(), TJets[j].py(), TJets[j].pz(), TJets[j].energy());
      Thrust += j2.Pt() * fabs(cos( j1.DeltaPhi(j2) ));
      sumpt += j2.Pt();
    }
    if(sumpt>0) Thrust /= sumpt;
    if(i==0) Thrust3 = Thrust;
    if(Thrust4 < Thrust) Thrust4 = Thrust;

    Thrust5 += j1.Pt() * fabs(cos( COJ5.DeltaPhi(j1) ));
    if(total_pt5>0) Thrust5 /= total_pt5;

    dR_jc += COJ5.DeltaR(j1);
    dR_jc_sigma += pow(COJ5.DeltaR(j1),2);
  }
  T3 = Thrust3;
  T4 = Thrust4;
  T5 = Thrust5;

  if(NTJet>0) dR_jc /= NTJet;
  if(NTJet>0) dR_jc_sigma /= NTJet;
  dR_jc_sigma = dR_jc_sigma - dR_jc*dR_jc;
  dR_jc_sigma = sqrt(dR_jc_sigma);

  DR_jc = dR_jc;
  DR_jc_sigma = dR_jc_sigma;

  //cout<<"NJet "<<NJet<<", NGJet "<<ngjet<<endl;
  //for (unsigned int i = 0; i < selectedJets.size() ; i++) {
    //double mindR = 10000;
    //double energy_ratio = 0;
    //TLorentzVector j1(selectedJets[i].px(), selectedJets[i].py(), selectedJets[i].pz(), selectedJets[i].energy());
    //for (int j=0;j<ngjet;j++) {
    //  if(mindR > j1.DeltaR(vgj[j])) {
    //    mindR = j1.DeltaR(vgj[j]);
    //    //energy_ratio = fabs(selectedJets[i].energy() - vgj[j].E())/selectedJets[i].energy();
    //  }
    //}
    //cout<<mindR<<"("<<energy_ratio<<")"<<" ";
  //}
  //cout<<endl;


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
      TLorentzVector j1(selectedJets[i].px(), selectedJets[i].py(), selectedJets[i].pz(), selectedJets[i].energy());
      if(Jet_BTag[i]) {
        double dmt = 10000;
        double mt_final = 10000;
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
      //double mindR=10000;
      int nq_temp=0;
      for(int j=0;j<nq_gen;j++) {
        //if(mindR > j1.DeltaR(vq[j]) mindR = j1.DeltaR(vq[j]);
        if(j1.DeltaR(vq[j]) < 0.5) nq_temp++;
      }
      nq_in_a_jet[i] = nq_temp;

//cout<<"check 100, ev "<<events<<", nl "<<nl<<", ngp "<<ngp<<endl;
//      int ngp_temp=0;
//cout<<"check 101, i "<<i<<endl;
//      for(int j=0;j<ngp;j++) {
//cout<<"check 102, j "<<j<<", gp pt "<<gp[j].Pt()<<endl;
//        if(j1.DeltaR(gp[j]) < 0.5) ngp_temp++;
//cout<<"check 103, ngp_temp "<<ngp_temp<<endl;
//      }
//cout<<"check 104"<<endl;
      //cout<<ngp_temp<<" ";
    }
//    cout<<endl;
  }

  int NWJet = WJets.size();
//cout<<"ev "<<events<<", NWJet "<<NWJet<<endl;
//  int NlW = NWJet*(NWJet-1)/2;
//  TLorentzVector lW[NlW];
//  for(int i=0;i<NlW;i++) if(NWJet>=6&&NWJet<=8){
//    TLorentzVector j1(WJets[ij[NWJet-6][i][0]].px(),WJets[ij[NWJet-6][i][0]].py(),WJets[ij[NWJet-6][i][0]].pz(),WJets[ij[NWJet-6][i][0]].energy());
//    TLorentzVector j2(WJets[ij[NWJet-6][i][1]].px(),WJets[ij[NWJet-6][i][1]].py(),WJets[ij[NWJet-6][i][1]].pz(),WJets[ij[NWJet-6][i][1]].energy());
//    lW[i] = j1+j2;
//  }

  int NWset=105;
  if(NWJet==6) NWset=15;
  if(NWJet==6||NWJet==7) nw=3;
  if(NWJet==8) nw=4;
  double MW_mean[105]={0,};
  typedef std::pair<double, int> MvsI;
  vector<MvsI> iW;
  TLorentzVector lW[NWset][4];
  TLorentzVector vw[nw];
  if(NWJet>=6&&NWJet<=8 && IsHadronTrig) {
    for(int i=0;i<NWset;i++) {
      //for(int j=0;j<4;j++) if(j<3 || (j==3&&NWJet>=8)){
      for(int j=0;j<nw;j++) {
        int n1 = W4[NWJet-6][i][j][0] - 1;
        int n2 = W4[NWJet-6][i][j][1] - 1;
        if(n1>=8) cout<<"n1 "<<n1<<endl;
        if(n2>=8) cout<<"n2 "<<n2<<endl;
        TLorentzVector j1(WJets[n1].px(),WJets[n1].py(),WJets[n1].pz(),WJets[n1].energy());
        TLorentzVector j2(WJets[n2].px(),WJets[n2].py(),WJets[n2].pz(),WJets[n2].energy());
        lW[i][j] = j1+j2;
        MW_mean[i] += lW[i][j].M();
      }
      if(NWJet==8) MW_mean[i] /= 4;
      if(NWJet==6||NWJet==7) MW_mean[i] /= 3;
      iW.push_back(make_pair(fabs(mW0-MW_mean[i]),i));
//      cout<<MW_mean[i]<<" ";
    }
//    cout<<"//"<<endl;
    sort(iW.begin(), iW.end(), less<MvsI>());
    //for(int i=0;i<NWset/5;i++) cout<<MW_mean[iW[i].second]<<" "; cout<<endl;
    //for(int i=0;i<NWset;i++) cout<<MW_mean[iW[i].second]<<" "; cout<<endl;

    //sort(MW_mean, MW_mean+NWset, greater<double>());
    //sort(MW_mean, MW_mean+NWset);
    //for(int i=0;i<NWset;i++) cout<<MW_mean[i]<<" "; cout<<endl;

    int best_lW = iW[0].second;
    MWm = MW_mean[best_lW];
    TLorentzVector m8j;
    //for(int i=0;i<4;i++) if(i<3 || (i==3&&NWJet>=8)) MW[i] = lW[iW[0].second][i].M();
    for(int i=0;i<nw;i++) 
    {
      vw[i] = lW[best_lW][i];
      MW[i] = vw[i].M();
      if(MW[i]>70&&MW[i]<95) nw2++;
      m8j += vw[i];
    }
    M8j = m8j.M()/nw;

    int nbj = BJets.size();
    //TLorentzVector vtc[nbj][nw];
    //int ntops = nbj*nw;
    //int ntop = TMath::Min(nbj,nw);
    //TLorentzVector tops[ntops];
    //TLorentzVector top[ntop];
    //for(int i=0;i<nbj;i++) {
    //  for(int j=0;j<nw;j++) { 
    //    vtc[i][j] = BJets[i]+lW[iW[0].second][i];
    //    vtc[i][j] = BJets[i]+vw[j];
    //  }
    //}

    if(nbj>=2) {
      const int ntops = nbj*nw;
      const int nt_min = TMath::Min(nbj,nw);
      const int nt_max = TMath::Max(nbj,nw);

      vector<vector<int>> tops(0, vector<int>(2,0));
      TLorentzVector vb[nbj];

      int ntN=1;
      for(int i=0;i<nt_min;i++) ntN *= nt_max - i;
      //cout<<"ntN "<<ntN<<endl;

      int Nt=0;
      for(int i=0;i<nbj;i++)
      {
        vb[i].SetPtEtaPhiE(BJets[i].pt(), BJets[i].eta(), BJets[i].phi(), BJets[i].energy());
        for(int j=0;j<nw;j++)
        {
          tops.push_back(vector<int>({i,j}));
          //cout<<tops[Nt][0]<<tops[Nt][1]<<" ";
          Nt++;
        }
      }
      //cout<<endl;

      int Ntop=0;
      vector<vector<int>> top(0, vector<int>(nt_min*2,0));
      for(int i=0;i<ntops;i++)
      {
        for(int j=i+1;j<ntops;j++)
        if(tops[j][0]!=tops[i][0]&&tops[j][1]!=tops[i][1])
        {
          if(nt_min==2)
          {
             top.push_back(vector<int>({tops[i][0],tops[i][1],tops[j][0],tops[j][1]}));
             Ntop++;
          }
          for(int k=j+1;k<ntops;k++)
          if(tops[k][0]!=tops[i][0]&&tops[k][1]!=tops[i][1]
          && tops[k][0]!=tops[j][0]&&tops[k][1]!=tops[j][1])
          {
            if(nt_min==3)
            {
              top.push_back(vector<int>({tops[i][0],tops[i][1],tops[j][0],tops[j][1],tops[k][0],tops[k][1]}));
              Ntop++;
              //cout<<tops[i][0]<<tops[i][1]<<" "<<tops[j][0]<<tops[j][1]<<" "<<tops[k][0]<<tops[k][1]<<endl;
            }
            for(int l=k+1;l<ntops;l++)
            if(tops[l][0]!=tops[i][0]&&tops[l][1]!=tops[i][1]
            && tops[l][0]!=tops[j][0]&&tops[l][1]!=tops[j][1]
            && tops[l][0]!=tops[k][0]&&tops[l][1]!=tops[k][1])
            {
              if(nt_min==4)
              {
                top.push_back(vector<int>({tops[i][0],tops[i][1],tops[j][0],tops[j][1],tops[k][0],tops[k][1],tops[l][0],tops[l][1]}));
                Ntop++;
                //cout<<tops[i][0]<<tops[i][1]<<" "<<tops[j][0]<<tops[j][1]<<" "<<tops[k][0]<<tops[k][1]<<" "<<tops[l][0]<<tops[l][1]<<endl;
              }
            }
          }
        }
      }
      //cout<<"Ntop "<<Ntop<<endl;
            
      //double Mtop[Ntop]={0,};
      double Mt_mean[Ntop]={0,};
      vector<MvsI> iT;
      //TLorentzVector lt[NWset][4];
      TLorentzVector vt[nt_min];
      int itop[Ntop][nt_min][2];

      int j3=0;
      for(vector<vector<int>>::iterator i_it = top.begin(); i_it != top.end(); i_it++)
      {       
        int j2=0;
        int i=0, j=0;
        for(vector<int>::iterator j_it = i_it->begin(); j_it != i_it->end(); j_it++)
        {
          if(j2%2==0) i = *j_it;
          if(j2%2==1) j = *j_it;
          //cout<<*j_it;
          //if(j2%2==1) cout<<" ";
          if(j2%2==1)
          {
            Mt_mean[j3] += (vb[i]+vw[j]).M();
            itop[j3][j2/2][0] = i;
            itop[j3][j2/2][1] = j;
            //cout<<i<<j<<" ";
          }
          j2++;
        }
        //Mtop[j3] = (vb[i]+vw[j]).M();
        Mt_mean[j3] /= nt_min;
        //cout<<Mt_mean[j3]<<" ";
        iT.push_back(make_pair(fabs(mt0-Mt_mean[j3]),j3));
        j3++;
      }               
      //cout<<endl;

      //sort(Mt_mean,Mt_mean+Ntop);
      //for(int i=0;i<Ntop;i++) cout<<Mt_mean[i]<<" "; cout<<endl;
      //sort(Mtop,Mtop+Ntop);
      //for(int i=0;i<Ntop;i++) cout<<Mtop[i]<<" "; cout<<endl;

      sort(iT.begin(), iT.end(), less<MvsI>());
      //for(int i=0;i<Ntop;i++) cout<<Mt_mean[iT[i].second]<<" "; cout<<"//"<<endl;

      nt = nt_min;
      Mtm = Mt_mean[iT[0].second];

      TLorentzVector m12j;
      for(int i=0;i<nt_min;i++) {
        int n1 = itop[iT[0].second][i][0];
        int n2 = itop[iT[0].second][i][1];
        vt[i] = vb[n1] + vw[n2];
        //Mt[i] = vt[i].M();
        if(vt[i].M()>150&&vt[i].M()<200) nt2++;
        m12j += vt[i];
      }
      M12j = m12j.M()/nt_min;
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
//     tree->Fill();
//   }
     //tree->Fill();
     if(nb==2&&nbbar==2&&nWp==2&&nWm==2) hskim0->Fill(nl);
     else if(nb==1&&nbbar==1&&nWp==1&&nWm==1) hskim0->Fill(nl+5);
     else hskim0->Fill(8);
     if(NJet>=6&&NBJet>=2) {
       if(nb==2&&nbbar==2&&nWp==2&&nWm==2) hskim1->Fill(nl);
       else if(nb==1&&nbbar==1&&nWp==1&&nWm==1) hskim1->Fill(nl+5);
       else hskim1->Fill(8);
       if(nb==2&&nbbar==2&&nWp==2&&nWm==2) htrig0->Fill(nl);
       else if(nb==1&&nbbar==1&&nWp==1&&nWm==1) htrig0->Fill(nl+5);
       else htrig0->Fill(8);
     }
     if(NJet>=6&&NBJet>=2 && IsHadronTrig) {
       if(nb==2&&nbbar==2&&nWp==2&&nWm==2) htrig1->Fill(nl);
       else if(nb==1&&nbbar==1&&nWp==1&&nWm==1) htrig1->Fill(nl+5);
       else htrig1->Fill(8);
       tree->Fill();
     }
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
//  tree->Branch("Muon_Pt",Muon_Pt,"Muon_Pt[NMuon]/F");
//  tree->Branch("Muon_Eta",Muon_Eta,"Muon_Eta[NMuon]/F");
//  tree->Branch("Muon_Phi",Muon_Phi,"Muon_Phi[NMuon]/F");
//  tree->Branch("Muon_E",Muon_E,"Muon_E[NMuon]/F");
//  tree->Branch("Muon_Iso03",Muon_Iso03,"Muon_Iso03[NMuon]/F");
//  tree->Branch("Muon_Iso04",Muon_Iso04,"Muon_Iso04[NMuon]/F");
//  tree->Branch("Muon_Charge",Muon_Charge,"Muon_Charge[NMuon]/F");

  tree->Branch("NLooseMuon",&NLooseMuon,"NLooseMuon/I");
//  tree->Branch("LooseMuon_Pt",LooseMuon_Pt,"LooseMuon_Pt[NLooseMuon]/F");
//  tree->Branch("LooseMuon_Eta",LooseMuon_Eta,"LooseMuon_Eta[NLooseMuon]/F");
//  tree->Branch("LooseMuon_Phi",LooseMuon_Phi,"LooseMuon_Phi[NLooseMuon]/F");
//  tree->Branch("LooseMuon_E",LooseMuon_E,"LooseMuon_E[NLooseMuon]/F");
//  tree->Branch("LooseMuon_Iso03",LooseMuon_Iso03,"LooseMuon_Iso03[NLooseMuon]/F");
//  tree->Branch("LooseMuon_Iso04",LooseMuon_Iso04,"LooseMuon_Iso04[NLooseMuon]/F");
//  tree->Branch("LooseMuon_Charge",LooseMuon_Charge,"LooseMuon_Charge[NLooseMuon]/F");

  tree->Branch("NElectron",&NElectron,"NElectron/I");
//  tree->Branch("Electron_Pt",Electron_Pt,"Electron_Pt[NElectron]/F");
//  tree->Branch("Electron_Eta",Electron_Eta,"Electron_Eta[NElectron]/F");
//  tree->Branch("Electron_Phi",Electron_Phi,"Electron_Phi[NElectron]/F");
//  tree->Branch("Electron_E",Electron_E,"Electron_E[NElectron]/F");
//  tree->Branch("Electron_Iso03",Electron_Iso03,"Electron_Iso03[NElectron]/F");
//  tree->Branch("Electron_Iso04",Electron_Iso04,"Electron_Iso04[NElectron]/F");
//  tree->Branch("Electron_Charge",Electron_Charge,"Electron_Charge[NElectron]/F");

  tree->Branch("NLooseElectron",&NLooseElectron,"NLooseElectron/I");
//  tree->Branch("LooseElectron_Pt",LooseElectron_Pt,"LooseElectron_Pt[NLooseElectron]/F");
//  tree->Branch("LooseElectron_Eta",LooseElectron_Eta,"LooseElectron_Eta[NLooseElectron]/F");
//  tree->Branch("LooseElectron_Phi",LooseElectron_Phi,"LooseElectron_Phi[NLooseElectron]/F");
//  tree->Branch("LooseElectron_E",LooseElectron_E,"LooseElectron_E[NLooseElectron]/F");
//  tree->Branch("LooseElectron_Iso03",LooseElectron_Iso03,"LooseElectron_Iso03[NLooseElectron]/F");
//  tree->Branch("LooseElectron_Iso04",LooseElectron_Iso04,"LooseElectron_Iso04[NLooseElectron]/F");
//  tree->Branch("LooseElectron_Charge",LooseElectron_Charge,"LooseElectron_Charge[NLooseElectron]/F");

  tree->Branch("NLooseLep2",&NLooseLep2,"NLooseLep2/I");

  tree->Branch("NJet",&NJet,"NJet/i");
  tree->Branch("NJetW",&NJetW,"NJetW/F");
  tree->Branch("NJetW2",&NJetW2,"NJetW2/F");
//  tree->Branch("T",&T,"T/F");
//  tree->Branch("T2",&T2,"T2/F");
//  tree->Branch("T3",&T3,"T3/F");
//  tree->Branch("T4",&T4,"T4/F");
  tree->Branch("T5",&T5,"T5/F");
//  tree->Branch("T6",&T6,"T6/F");

//  tree->Branch("DR_jc",&DR_jc,"DR_jc/F");
//  tree->Branch("DR_jc_sigma",&DR_jc_sigma,"DR_jc_sigma/F");

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

//  tree->Branch("Jet_Pt_Sum",Jet_Pt_Sum,"Jet_Pt_Sum[12]/F");
//  tree->Branch("Jet_Pt_Sum2",Jet_Pt_Sum2,"Jet_Pt_Sum2[12]/F");
//  tree->Branch("Jet_mass_sum",Jet_mass_sum,"Jet_mass_sum[12]/F");
  tree->Branch("Jet_mass_tot",&Jet_mass_tot,"Jet_mass_tot/F");
//  tree->Branch("Jet_dR_min",Jet_dR_min,"Jet_dR_min[NJet]/F");
//  tree->Branch("Jet_dR_from_W",Jet_dR_from_W,"Jet_dR_from_W[NJet]/F");
//  tree->Branch("Jet_dR_min_with_q",Jet_dR_min_with_q,"Jet_dR_min_with_q[NJet]/F");
//  tree->Branch("Jet_dE_with_q",Jet_dE_with_q,"Jet_dE_with_q[NJet]/F");
//  tree->Branch("Jet_dpt_with_q",Jet_dpt_with_q,"Jet_dpt_with_q[NJet]/F");
//  tree->Branch("W_dR",W_dR,"W_dR[NJet]/F");
//  tree->Branch("W_dE",W_dE,"W_dE[NJet]/F");
//  tree->Branch("W_dpt",W_dpt,"W_dpt[NJet]/F");
  tree->Branch("Jet_mW",Jet_mW,"Jet_mW[NJet]/F");
  tree->Branch("Jet_mt",Jet_mt,"Jet_mt[NJet]/F");
  tree->Branch("Jet_mW_Flavour",Jet_mW_Flavour,"Jet_mW_Flavour[NJet]/F");
  tree->Branch("NW",&NW,"NW/I");
  tree->Branch("Nt",&Nt,"Nt/I");
  tree->Branch("Jet_HT",&Jet_HT,"Jet_HT/F");
  tree->Branch("Jet_H",&Jet_H,"Jet_H/F");
//  tree->Branch("Jet_Eta",Jet_Eta,"Jet_Eta[NJet]/F");
//  tree->Branch("Jet_Phi",Jet_Phi,"Jet_Phi[NJet]/F");
  tree->Branch("Jet_E",Jet_E,"Jet_E[NJet]/F");
  tree->Branch("Jet_Et",Jet_Et,"Jet_Et[NJet]/F");
  tree->Branch("Jet_partonFlavour",Jet_partonFlavour,"Jet_partonFlavour[NJet]/F");
  tree->Branch("Jet_hadronFlavour",Jet_hadronFlavour,"Jet_hadronFlavour[NJet]/F");
  tree->Branch("Jet_BTag",Jet_BTag,"Jet_BTag[NJet]/F");
  tree->Branch("Jet_bDiscriminator",Jet_bDiscriminator,"Jet_bDiscriminator[NJet]/F"); 
//  tree->Branch("Jet_pfCombinedCvsLJetTags",Jet_pfCombinedCvsLJetTags,"Jet_pfCombinedCvsLJetTags[NJet]/F"); 
//  tree->Branch("Jet_pfCombinedCvsBJetTags",Jet_pfCombinedCvsBJetTags,"Jet_pfCombinedCvsBJetTags[NJet]/F"); 

//  tree->Branch("Jet_JES_Up",Jet_JES_Up,"Jet_JES_Up[NJet]/F");
//  tree->Branch("Jet_JES_Dw",Jet_JES_Dw,"Jet_JES_Dw[NJet]/F");

//  tree->Branch("csvid",csvid,"csvid[NJet]/i");

  tree->Branch("NBJet",&NBJet,"NBJet/i");
  tree->Branch("HTb",&HTb,"HTb/F");
//  tree->Branch("DiLeptonic",&DiLeptonic,"DiLeptonic/i");
//  tree->Branch("SemiLeptonic",&SemiLeptonic,"SemiLeptonic/i");

//  tree->Branch("TTBJ",&TTBJ,"TTBJ/i");
//  tree->Branch("TTBB",&TTBB,"TTBB/i");
//  tree->Branch("TTCC",&TTCC,"TTCC/i");
//  tree->Branch("TTJJ",&TTJJ,"TTJJ/i");

//  tree->Branch("GenNJet20",&GenNJet20, "GenNJet20/i");
//  tree->Branch("GenNBJet20",&GenNBJet20, "GenNBJet20/i");
//  tree->Branch("GenNCJet20",&GenNCJet20, "GenNCJet20/i");
//  tree->Branch("GenNAddJet20",&GenNAddJet20, "GenNAddJet20/i");
//  tree->Branch("GenNAddBJet20",&GenNAddBJet20, "GenNAddBJet20/i");
//  tree->Branch("GenNAddCJet20",&GenNAddCJet20, "GenNAddCJet20/i");   

//  tree->Branch("GenLepton1_Pt",&GenLepton1_Pt, "GenLepton1_Pt/f");
//  tree->Branch("GenLepton1_Eta",&GenLepton1_Eta, "GenLepton1_Eta/f");
//  tree->Branch("GenLepton2_Pt",&GenLepton2_Pt, "GenLepton2_Pt/f");
//  tree->Branch("GenLepton2_Eta",&GenLepton2_Eta, "GenLepton2_Eta/f");

//  tree->Branch("MT_MuonMET",MT_MuonMET,"MT_MuonMET[NMuon]/F"); 
//  tree->Branch("Phi_MuonMET",Phi_MuonMET,"Phi_MuonMET[NMuon]/F"); 
//  tree->Branch("MT_ElectronMET",MT_ElectronMET,"MT_ElectronMET[NElectron]/F"); 
//  tree->Branch("Phi_ElectronMET",Phi_ElectronMET,"Phi_ElectronMET[NElectron]/F"); 

//  tree->Branch("Kin_Hmass",&Kin_Hmass,"Kin_Hmass/F");
//  tree->Branch("Kin_HdRbb",&Kin_HdRbb,"Kin_HdRbb/F");
//  tree->Branch("Kin_Chi2",&Kin_Chi2,"Kin_Chi2/F"); 
//  tree->Branch("Kin_TopMHc",&Kin_TopMHc,"Kin_TopMHc/F"); 
//  tree->Branch("Kin_TopMWb",&Kin_TopMWb,"Kin_TopMWb/F"); 
//  tree->Branch("Kin_Wmass",&Kin_Wmass,"Kin_Wmass/F"); 

  tree->Branch("IsMuonTrig",&IsMuonTrig,"IsMuonTrig/i"); 
  tree->Branch("IsElectronTrig",&IsElectronTrig,"IsElectronTrig/i"); 
  tree->Branch("IsHadronTrig",&IsHadronTrig,"IsHadronTrig/i"); 

  tree->Branch("nb",&nb,"nb/i");
  tree->Branch("nbbar",&nbbar,"nbbar/i");
  tree->Branch("nWp",&nWp,"nWp/i");
  tree->Branch("nWm",&nWm,"nWm/i");
  tree->Branch("nq",&nq,"nq/i");
  tree->Branch("nl",&nl,"nl/i");
  tree->Branch("nTau",&nTau,"nTau/i");
  tree->Branch("dTau",&dTau,"dTau/O");
  //tree->Branch("id_from_W",id_from_W,"id_from_W[nq]/i");
  tree->Branch("nq_all",&nq_all,"nq_all/i");
  tree->Branch("mindRqq",mindRqq,"mindRqq[nq_all]/F");
  tree->Branch("nq_in_a_jet",nq_in_a_jet,"nq_in_a_jet[NJet]/i");
  tree->Branch("Ngp_noc",Ngp_noc,"Ngp_noc[nq_all]/i");
  tree->Branch("Ngp_other",Ngp_other,"Ngp_other[nq_all]/i");
  tree->Branch("dRmax_gp",dRmax_gp,"dRmax_gp[nq_all]/F");
  tree->Branch("nqjet",&nqjet,"nqjet/i");
  tree->Branch("nqjet2",&nqjet2,"nqjet2/i");
  tree->Branch("dR_qqjet",dR_qqjet,"dR_qqjet[nqjet2]/F");

  tree->Branch("nw",&nw,"nw/i");
  tree->Branch("nw2",&nw2,"nw2/i");
  //tree->Branch("MW",MW,"MW[nw]/F");
  tree->Branch("MWm",&MWm,"MWm/F");
  tree->Branch("nt",&nt,"nt/i");
  tree->Branch("nt2",&nt2,"nt2/i");
  //tree->Branch("Mt,Mt,"Mt[nw]/F");
  tree->Branch("Mtm",&Mtm,"Mtm/F");
  tree->Branch("M8j",&M8j,"M8j/F");
  tree->Branch("M12j",&M12j,"M12j/F");

  //int ij[3][105][2]={0,};
  int NW[3]={0,};
  for(int i=0;i<3;i++) {
    int njet=i+6;
    for(int j=1;j<=njet;j++) {
      for(int k=j+1;k<=njet;k++) {
        ij[i][NW[i]][0] = j;
        ij[i][NW[i]][1] = k;
        cout<<ij[i][NW[i]][0]<<ij[i][NW[i]][1]<<" ";
	NW[i]++;
      }
    }
    cout<<endl;
    cout<<"NJet "<<njet<<", NW cand"<<njet*(njet-1)/2<<", NW "<<NW[i]<<endl;
  }

  int NW4[3]={0,};
  for(int i=0;i<3;i++) {
    for(int j=0;j<NW[i];j++) {
      for(int k=j+1;k<NW[i];k++) {
        if(ij[i][k][0]!=ij[i][j][0]&&ij[i][k][0]!=ij[i][j][1]&&ij[i][k][1]!=ij[i][j][0]&&ij[i][k][1]!=ij[i][j][1])
        for(int l=k+1;l<NW[i];l++) {
          if(ij[i][l][0]!=ij[i][j][0]&&ij[i][l][0]!=ij[i][j][1]&&ij[i][l][1]!=ij[i][j][0]&&ij[i][l][1]!=ij[i][j][1]
          && ij[i][l][0]!=ij[i][k][0]&&ij[i][l][0]!=ij[i][k][1]&&ij[i][l][1]!=ij[i][k][0]&&ij[i][l][1]!=ij[i][k][1]) {
            if(NW4[i]>105) {cout<<"i"<<i<<", j"<<j<<", k"<<k<<", l"<<l<<", NW4 "<<NW4[i]<<endl; break;}
            W4[i][NW4[i]][0][0] = ij[i][j][0];
            W4[i][NW4[i]][0][1] = ij[i][j][1];
            W4[i][NW4[i]][1][0] = ij[i][k][0];
            W4[i][NW4[i]][1][1] = ij[i][k][1];
            W4[i][NW4[i]][2][0] = ij[i][l][0];
            W4[i][NW4[i]][2][1] = ij[i][l][1];
            if(i<=1) NW4[i]++;
//            if(i<=1) cout<<"i"<<i<<", "<<ij[i][j][0]<<ij[i][j][1]<<" "<<ij[i][k][0]<<ij[i][k][1]<<" "<<ij[i][l][0]<<ij[i][l][1]<<", NW4 "<<NW4[i]<<endl;
            for(int n=l+1;n<NW[i];n++) {
              if(ij[i][n][0]!=ij[i][j][0]&&ij[i][n][0]!=ij[i][j][1]&&ij[i][n][1]!=ij[i][j][0]&&ij[i][n][1]!=ij[i][j][1]
              && ij[i][n][0]!=ij[i][k][0]&&ij[i][n][0]!=ij[i][k][1]&&ij[i][n][1]!=ij[i][k][0]&&ij[i][n][1]!=ij[i][k][1]
              && ij[i][n][0]!=ij[i][l][0]&&ij[i][n][0]!=ij[i][l][1]&&ij[i][n][1]!=ij[i][l][0]&&ij[i][n][1]!=ij[i][l][1]) {
                if(NW4[i]>105) {cout<<"i"<<i<<", j"<<j<<", k"<<k<<", l"<<l<<", n"<<n<<", NW4 "<<NW4[i]<<endl; break;}
                W4[i][NW4[i]][3][0] = ij[i][n][0];
                W4[i][NW4[i]][3][1] = ij[i][n][1];
                NW4[i]++;
//                if(i==2) cout<<"i"<<i<<", "<<ij[i][j][0]<<ij[i][j][1]<<" "<<ij[i][k][0]<<ij[i][k][1]<<" "<<ij[i][l][0]<<ij[i][l][1]<<" "<<ij[i][n][0]<<ij[i][n][1]<<", NW4 "<<NW4[i]<<endl;
              }
            }
          }
        }
      }
    }
//    cout<<"NW4["<<i<<"] "<<NW4[i]<<endl;
  }


}

void TopAnalyzer::clear()
{
  PUWeight = 1.0;
  GenWeight = 1.0;
  NVertex = -1;
  NMuon = -1;
  NLooseMuon = -1;
  NElectron = -1;
  NLooseLep2 = -1;
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
  nTau = 0;

  nqjet = 0;
  nqjet2 = 0;

  Jet_Pt_tot = 0;
  NW = 0;
  Nt = 0;
  Jet_HT = 0;
  Jet_H = 0;
  HTb = 0;

  nw = 0;
  nw2 = 0;
  nt = 0;
  nt2 = 0;

  T = -1;
  T2 = -1;
  T3 = -1;
  T4 = -1;
  T5 = -1;
  T6 = -1;
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
