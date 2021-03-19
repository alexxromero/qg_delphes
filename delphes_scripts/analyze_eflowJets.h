//#include "process_caloJets.h"
#include <stdio.h>
#include <iostream>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Recluster.hh"

#else
class ExRootTreeReader;
class ExRootResult;
#endif

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

class analyze_eflowJets  {
  public:
    analyze_eflowJets(string inputFile, string outputFile);
    ~analyze_eflowJets();
    void build_ttree(TTree *ttree);
    void analyze_delphesevents(string inputFile, string outputFile);
  private:
    TFile *File, *outRootFile;
    TTree *ttree;
    TTree *outTTree;
    TH2F *EmImage, *HadImage, *TrackImage;
    TH2F *EmImage_multi, *HadImage_multi, *TrackImage_multi;
    Jet *jet;
    GenParticle *particle;
    TObject *object;
    Track *track;
    Tower *tower;

    const double pi=acos(-1.0);

    // ----- unprocessed tower features ----- //
    Float_t out_tau1_05;
    Float_t out_tau1_10;
    Float_t out_tau1_20;
    Float_t out_jetPT;
    Float_t out_jetEta;
    Float_t out_jetPhi;
    Float_t out_jetMass;
    Int_t out_nparticles;
    Int_t out_ntowers;
    Int_t out_nEmtowers;
    Int_t out_nHadtowers;
    Int_t out_ntracks;
    Float_t out_towers_pt[200];
    Float_t out_towers_eta[200];
    Float_t out_towers_phi[200];
    Int_t out_towers_isEcal[200];  // flag -- is 1 if tower is from the Ecal
    Float_t out_tracks_pt[200];
    Float_t out_tracks_eta[200];
    Float_t out_tracks_phi[200];

    // ----- binned tower features ----- //
    Float_t out_tau1_05_im;
    Float_t out_tau1_10_im;
    Float_t out_tau1_20_im;
    Float_t out_jetPT_im;
    Float_t out_jetEta_im;
    Float_t out_jetPhi_im;
    Float_t out_jetMass_im;
    Int_t out_nparticles_im;
    Int_t out_ntowers_im;
    Int_t out_nEmtowers_im;
    Int_t out_nHadtowers_im;
    Int_t out_ntracks_im;
    Int_t out_multiple_Emtowers_perpix;
    Int_t out_multiple_Hadtowers_perpix;
    Float_t out_towers_pt_im[200];
    Float_t out_towers_eta_im[200];
    Float_t out_towers_phi_im[200];
    Int_t out_towers_isEcal_im[200];  // flag -- is 1 if tower is from the Ecal
    Float_t out_tracks_pt_im[200];
    Float_t out_tracks_eta_im[200];
    Float_t out_tracks_phi_im[200];
    Int_t out_tracks_multi_im[200];
};

// class IsEcal: public PseudoJet::UseInfoBase {
//   IsEcal(int id) : _ecal_id(id);
//   int ecal_id() const {return _ecal_id;}
//   int _ecal_id;
// };
