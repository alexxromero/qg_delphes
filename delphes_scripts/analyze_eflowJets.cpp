#include "analyze_eflowJets.h"

// JetDefinition is needed if calculating observables that require
// reclustering. In that case, use the same JetDefinition defined in
// the Delphes card.
fastjet::JetDefinition *jet_def = new fastjet::JetDefinition(
    fastjet::antikt_algorithm, 0.4);

float deltaPhi(float this_phi, float that_phi) {
  float dphi = this_phi - that_phi;
  if (dphi < -pi)
    dphi += 2*pi;
  if (dphi > pi)
    dphi -= 2*pi;
  return dphi;
}

float Nsubjettiness_Nbeta(PseudoJet Jet, int N, float beta) {
    Nsubjettiness nSub(N, KT_Axes(), UnnormalizedMeasure(beta));
    return nSub(Jet);
}

void analyze_eflowJets::build_ttree(TTree *ttree)  {
  // jets before any preprocessing -- used mainly for control
  ttree->Branch("jet_PT", &out_jetPT, "jet_PT/F");
  ttree->Branch("jet_Eta", &out_jetEta, "jet_Eta/F");
  ttree->Branch("jet_Phi", &out_jetPhi, "jet_Phi/F");
  ttree->Branch("jet_Mass", &out_jetMass, "jet_Mass/F");
  ttree->Branch("tau1_05", &out_tau1_05, "tau1_05/F");
  ttree->Branch("tau1_10", &out_tau1_10, "tau1_10/F");
  ttree->Branch("tau1_20", &out_tau1_20, "tau1_20/F");
  ttree->Branch("nparticles", &out_nparticles, "nparticles/I");
  ttree->Branch("ntowers", &out_ntowers, "ntowers/I");
  ttree->Branch("nEmtowers", &out_nEmtowers, "nEmtowers/I");
  ttree->Branch("nHadtowers", &out_nHadtowers, "nHadtowers/I");
  ttree->Branch("ntracks", &out_ntracks, "ntracks/I");
  ttree->Branch("towers_PT", &out_towers_pt, "towers_PT[200]/F");
  ttree->Branch("towers_Eta", &out_towers_eta, "towers_Eta[200]/F");
  ttree->Branch("towers_Phi", &out_towers_phi, "towers_Phi[200]/F");
  ttree->Branch("towers_IsEcal", &out_towers_isEcal, "towers_IsEcal[200]/I");
  ttree->Branch("tracks_PT", &out_tracks_pt, "tracks_PT[200]/F");
  ttree->Branch("tracks_Eta", &out_tracks_eta, "tracks_Eta[200]/F");
  ttree->Branch("tracks_Phi", &out_tracks_phi, "tracks_Phi[200]/F");
  // jets after centering and binning
  ttree->Branch("jet_PT_image", &out_jetPT_im, "jet_PT_image/F");
  ttree->Branch("jet_Eta_image", &out_jetEta_im, "jet_Eta_image/F");
  ttree->Branch("jet_Phi_image", &out_jetPhi_im, "jet_Phi_image/F");
  ttree->Branch("jet_Mass_image", &out_jetMass_im, "jet_Mass_image/F");
  ttree->Branch("tau1_05_image", &out_tau1_05_im, "tau1_05_image/F");
  ttree->Branch("tau1_10_image", &out_tau1_10_im, "tau1_10_image/F");
  ttree->Branch("tau1_20_image", &out_tau1_20_im, "tau1_20_image/F");
  ttree->Branch("nparticles_image", &out_nparticles_im, "nparticles_image/I");
  ttree->Branch("ntowers_image", &out_ntowers_im, "ntowers_image/I");
  ttree->Branch("nEmtowers_image", &out_nEmtowers_im, "nEmtowers_image/I");
  ttree->Branch("nHadtowers_image", &out_nHadtowers_im, "nHadtowers_image/I");
  ttree->Branch("ntracks_image", &out_ntracks_im, "ntracks_image/I");
  ttree->Branch("multiple_Emtowers_perpix", &out_multiple_Emtowers_perpix, "multiple_Emtowers_perpix/I");
  ttree->Branch("multiple_Hadtowers_perpix", &out_multiple_Hadtowers_perpix, "multiple_Hadtowers_perpix/I");
  ttree->Branch("towers_PT_image", &out_towers_pt_im, "towers_PT_image[200]/F");
  ttree->Branch("towers_Eta_image", &out_towers_eta_im, "towers_Eta_image[200]/F");
  ttree->Branch("towers_Phi_image", &out_towers_phi_im, "towers_Phi_image[200]/F");
  ttree->Branch("towers_IsEcal_image", &out_towers_isEcal_im, "towers_IsEcal_image[200]/I");
  ttree->Branch("tracks_PT_image", &out_tracks_pt_im, "tracks_PT_image[200]/F");
  ttree->Branch("tracks_Eta_image", &out_tracks_eta_im, "tracks_Eta_image[200]/F");
  ttree->Branch("tracks_Phi_image", &out_tracks_phi_im, "tracks_Phi_image[200]/F");
  ttree->Branch("tracks_Multi_image", &out_tracks_multi_im, "tracks_Multi_image[200]/I");
}

void analyze_eflowJets::analyze_delphesevents(string inputFile,
                                              string outputFile)  {
  gSystem->Load("libDelphes");
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile.c_str());
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchMissingHT = treeReader->UseBranch("ScalarHT");

  Float_t deta_pix = 0.025;  // eta tower granularity - see Delphes card
  Float_t dphi_pix = pi/126.0;  // phi tower granularity - see Delphes card
  Float_t eta_npix = 33;  // odd No. of pixels is preferred to have a
  Float_t phi_npix = 33;  // central pixel
  Float_t eta_jetR = eta_npix * deta_pix / 2;  // 0.4125
  Float_t phi_jetR = phi_npix * dphi_pix / 2;  // 0.4125

  EmImage = new TH2F("ECAL", "ECAL Image",
                     eta_npix, -eta_jetR, +eta_jetR,
                     phi_npix, -phi_jetR, +phi_jetR);

  EmImage_multi = new TH2F("ECAL Multiplicity", "ECAL Image Multiplicity",
                           eta_npix, -eta_jetR, +eta_jetR,
                           phi_npix, -phi_jetR, +phi_jetR);

  HadImage = new TH2F("EHAD", "EHAD Image",
                      eta_npix, -eta_jetR, +eta_jetR,
                      phi_npix, -phi_jetR, +phi_jetR);

  HadImage_multi = new TH2F("EHAD Multiplicity", "EHAD Image Multiplicity",
                            eta_npix, -eta_jetR, +eta_jetR,
                            phi_npix, -phi_jetR, +phi_jetR);

  // using the same granularity as the calorimeters for the track image,
  // paired with a multiplicity image
  TrackImage = new TH2F("Tracks", "Rec Tracks Image",
                        eta_npix, -eta_jetR, +eta_jetR,
                        phi_npix, -phi_jetR, +phi_jetR);

  TrackImage_multi = new TH2F("Tracks Multiplicity", "Track Image Multiplicity",
                              eta_npix, -eta_jetR, +eta_jetR,
                              phi_npix, -phi_jetR, +phi_jetR);

  outRootFile = new TFile(outputFile.c_str(), "RECREATE");
  ttree = new TTree("EFlowJets", "Delphes EFlow jets");
  build_ttree(ttree);

  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events **" << endl;

  Float_t ptmin = 500;  // GeV
  Float_t ptmax = 550;  // GeV
  Int_t njets_accepted = 0;
  for (Int_t entry = 0; entry < allEntries; entry++) {
    treeReader->ReadEntry(entry);
    if (branchJet->GetEntries()>0) {
      jet = (Jet*) branchJet->At(0);  // Keep only the leading jet
      if ((ptmin <= jet->PT) && (jet->PT <= ptmax) && (abs(jet->Eta) <= 1.5)) {
        vector<PseudoJet> tower_collection;
        vector<PseudoJet> track_collection;
        vector<PseudoJet> constit_collection;  // towers + tracks
        Int_t nparticles = 0;  // Sanity check. should be zero
        Int_t ntracks = 0;
        Int_t nEmtowers = 0;
        Int_t nHadtowers = 0;
        // Initialize all histograms to zero pT (or multiplicity)
        for (Int_t i = 1; i <= eta_npix; i++)  {
          for (Int_t j = 1; j <= phi_npix; j++) {
            EmImage->SetBinContent(i, j, 0);
            EmImage_multi->SetBinContent(i, j, 0);
            HadImage->SetBinContent(i, j, 0);
            HadImage_multi->SetBinContent(i, j, 0);
            TrackImage->SetBinContent(i, j, 0);
            TrackImage_multi->SetBinContent(i, j, 0);
          }
        }
        for (Int_t j = 0; j < jet->Constituents.GetEntriesFast(); j++) {
          object = jet->Constituents.At(j);
          if(object->IsA() == GenParticle::Class())
            nparticles++;
          else if(object->IsA() == Tower::Class()) {
            tower = (Tower*) object;
            TLorentzVector momentum = tower->P4();
            Float_t dEta = tower->Eta - jet->Eta;
            Float_t dPhi = deltaPhi(tower->Phi, jet->Phi);
            PseudoJet pj = PseudoJet(momentum.Px(), momentum.Py(),
                                     momentum.Pz(), momentum.E());
            PseudoJet pj_centered;
            pj_centered.reset_PtYPhiM(momentum.Pt(), dEta,
                                      dPhi, momentum.M());
            if (tower->Eem > 0)  {
              EmImage->Fill(pj_centered.pseudorapidity(),
                            pj_centered.phi_std(),
                            pj_centered.perp());
              EmImage_multi->Fill(pj_centered.pseudorapidity(),
                                  pj_centered.phi_std(),
                                  1);
              pj.set_user_index(1);
              tower_collection.push_back(pj);
              constit_collection.push_back(pj);
              nEmtowers++;
            }
            else if (tower->Ehad > 0)  {
              HadImage->Fill(pj_centered.pseudorapidity(),
                             pj_centered.phi_std(),
                             pj_centered.perp());
              HadImage_multi->Fill(pj_centered.pseudorapidity(),
                                   pj_centered.phi_std(),
                                   1);
              pj.set_user_index(0);
              tower_collection.push_back(pj);
              constit_collection.push_back(pj);
              nHadtowers++;
            }
          }
          else if (object->IsA() == Track::Class())  {
            track = (Track*) object;
            TLorentzVector momentum = track->P4();
            Float_t dEta = track->Eta - jet->Eta;
            Float_t dPhi = deltaPhi(track->Phi, jet->Phi);
            PseudoJet pj = PseudoJet(momentum.Px(), momentum.Py(),
                                     momentum.Pz(), momentum.E());
            PseudoJet pj_centered;
            pj_centered.reset_PtYPhiM(momentum.Pt(), dEta,
                                      dPhi, momentum.M());
            TrackImage->Fill(pj_centered.pseudorapidity(),
                             pj_centered.phi_std(),
                             pj_centered.perp());
            TrackImage_multi->Fill(pj_centered.pseudorapidity(),
                                   pj_centered.phi_std(),
                                   1);
            pj.set_user_index(-1);  // dummy tag
            track_collection.push_back(pj);
            constit_collection.push_back(pj);
            ntracks++;
          }
          else continue;
        }

        // Recluster the jet constituents. If interested in observables of
        // the towers only, use tower_collection as the input to the
        // ClusterSequence.
        ClusterSequence ReclusterConsts(constit_collection, *jet_def);
        vector<PseudoJet> InclusiveJets = sorted_by_pt(ReclusterConsts.inclusive_jets(ptmin));
        if (InclusiveJets.size() == 0) continue;
        PseudoJet leading_jet = InclusiveJets[0];
        Float_t tau1_05 = Nsubjettiness_Nbeta(leading_jet, 1, 0.5);
        Float_t tau1_1 = Nsubjettiness_Nbeta(leading_jet, 1, 1);
        Float_t tau1_2 = Nsubjettiness_Nbeta(leading_jet, 1, 2);
        out_tau1_05 = tau1_05;
        out_tau1_10 = tau1_1;
        out_tau1_20 = tau1_2;
        out_jetPT = leading_jet.pt();
        out_jetEta = leading_jet.pseudorapidity();
        out_jetPhi = leading_jet.phi_std();
        out_jetMass = leading_jet.m();
        out_nparticles = nparticles;
        out_ntowers = nEmtowers + nHadtowers;
        out_nEmtowers = nEmtowers;
        out_nHadtowers = nHadtowers;
        out_ntracks = ntracks;

        if (jet->Constituents.GetEntriesFast() != leading_jet.constituents().size())  // if clustering constituents
          cout << "WARNING: mismatch in jet nconstituents.  ";

        // unfold histograms to get the index of the non-zero pixels
        vector<PseudoJet> tower_collection_image;
        vector<PseudoJet> track_collection_image;
        vector<PseudoJet> constit_collection_image;
        vector<Int_t> track_multi_image;
        Int_t ntracks_image = 0;
        Int_t nEmtowers_image = 0;
        Int_t nHadtowers_image = 0;
        Int_t multiple_Emtowers_perpix = 0;  // sanity check -- should be 0
        Int_t multiple_Hadtowers_perpix = 0;  // sanity check -- should be 0
        for (Int_t i = 1; i <= eta_npix; i++) {
          for (Int_t j = 1; j <= phi_npix; j++) {
            Float_t EmBinPT = EmImage->GetBinContent(i, j);
            if (EmBinPT > 0) {
              PseudoJet pEm;
              pEm.reset_momentum_PtYPhiM(EmImage->GetBinContent(i, j),
                                         EmImage->GetXaxis()->GetBinCenter(i),
                                         EmImage->GetYaxis()->GetBinCenter(j),
                                         0.0);
              pEm.set_user_index(1);
              tower_collection_image.push_back(pEm);
              constit_collection_image.push_back(pEm);
              nEmtowers_image++;
              if (EmImage_multi->GetBinContent(i, j) > 1)
                multiple_Emtowers_perpix++;
            }
            Float_t HadBinPT = HadImage->GetBinContent(i, j);
            if (HadBinPT > 0) {
              PseudoJet pHad;
              pHad.reset_momentum_PtYPhiM(HadImage->GetBinContent(i, j),
                                          HadImage->GetXaxis()->GetBinCenter(i),
                                          HadImage->GetYaxis()->GetBinCenter(j),
                                          0.0);
              pHad.set_user_index(0);
              tower_collection_image.push_back(pHad);
              constit_collection_image.push_back(pHad);
              nHadtowers_image++;
              if (HadImage_multi->GetBinContent(i, j) > 1)
                multiple_Hadtowers_perpix++;
            }
            Float_t TrackBinPT = TrackImage->GetBinContent(i, j);
            if (TrackBinPT > 0) {
              PseudoJet pTrack;
              pTrack.reset_momentum_PtYPhiM(TrackImage->GetBinContent(i, j),
                                            TrackImage->GetXaxis()->GetBinCenter(i),
                                            TrackImage->GetYaxis()->GetBinCenter(j),
                                            0.0);
              pTrack.set_user_index(-1);
              track_collection_image.push_back(pTrack);
              constit_collection_image.push_back(pTrack);
              ntracks_image++;
              track_multi_image.push_back(TrackImage_multi->GetBinContent(i, j));
            }
          }
        }

        // Jet features of after centering and binning -- //
        ClusterSequence ReclusterConsts_im(constit_collection_image, *jet_def);
        vector<PseudoJet> InclusiveJets_im = sorted_by_pt(ReclusterConsts_im.inclusive_jets(ptmin));
        if (InclusiveJets_im.size() == 0) continue;
        PseudoJet leading_jet_im = InclusiveJets_im[0];
        Float_t tau1_05_im = Nsubjettiness_Nbeta(leading_jet_im, 1, 0.5);
        Float_t tau1_1_im = Nsubjettiness_Nbeta(leading_jet_im, 1, 1);
        Float_t tau1_2_im = Nsubjettiness_Nbeta(leading_jet_im, 1, 2);
        out_tau1_05_im = tau1_05_im;
        out_tau1_10_im = tau1_1_im;
        out_tau1_20_im = tau1_2_im;
        out_jetPT_im = leading_jet_im.pt();
        out_jetEta_im = leading_jet_im.pseudorapidity();
        out_jetPhi_im = leading_jet_im.phi_std();
        out_jetMass_im = leading_jet_im.m();
        out_nparticles_im = 0;
        out_ntowers_im = nEmtowers_image + nHadtowers_image;
        out_nEmtowers_im = nEmtowers_image;
        out_nHadtowers_im = nHadtowers_image;
        out_ntracks_im = ntracks_image;
        out_multiple_Emtowers_perpix = multiple_Emtowers_perpix;
        out_multiple_Hadtowers_perpix = multiple_Hadtowers_perpix;

        // save the 3-M of the towers -- before and after preprocessing --
        // in padded branch arrays
        for (Int_t i=0; i<200; i++) {
          if (i < tower_collection.size()) {
            out_towers_pt[i] = tower_collection.at(i).pt();
            out_towers_eta[i] = tower_collection.at(i).pseudorapidity();
            out_towers_phi[i] = tower_collection.at(i).phi_std();
            out_towers_isEcal[i] = tower_collection.at(i).user_index();
          }
          else {
            out_towers_pt[i] = 0;
            out_towers_eta[i] = 0;
            out_towers_phi[i] = 0;
            out_towers_isEcal[i] = -999;
          }
          if (i < track_collection.size()) {
            out_tracks_pt[i] = track_collection.at(i).pt();
            out_tracks_eta[i] = track_collection.at(i).pseudorapidity();
            out_tracks_phi[i] = track_collection.at(i).phi_std();
          }
          else {
            out_tracks_pt[i] = 0;
            out_tracks_eta[i] = 0;
            out_tracks_phi[i] = 0;
          }
          if (i < tower_collection_image.size()) {
            out_towers_pt_im[i] = tower_collection_image.at(i).pt();
            out_towers_eta_im[i] = tower_collection_image.at(i).pseudorapidity();
            out_towers_phi_im[i] = tower_collection_image.at(i).phi_std();
            out_towers_isEcal_im[i] = tower_collection_image.at(i).user_index();
          }
          else {
            out_towers_pt_im[i] = 0;
            out_towers_eta_im[i] = 0;
            out_towers_phi_im[i] = 0;
            out_towers_isEcal_im[i] = -999;
          }
          if (i < track_collection_image.size()) {
            out_tracks_pt_im[i] = track_collection_image.at(i).pt();
            out_tracks_eta_im[i] = track_collection_image.at(i).pseudorapidity();
            out_tracks_phi_im[i] = track_collection_image.at(i).phi_std();
            out_tracks_multi_im[i] = track_multi_image.at(i);
          }
          else {
            out_tracks_pt_im[i] = 0;
            out_tracks_eta_im[i] = 0;
            out_tracks_phi_im[i] = 0;
            out_tracks_multi_im[i] = 0;
          }
        }

        ttree->Fill();
        njets_accepted++;
      }  // if leading jet within cuts
    }  // if event not empty
  }  // event loop
  ttree->Write();
  EmImage->Reset();
  HadImage->Reset();
  cout << "total of " << njets_accepted << " jets generated "  << endl;
}

analyze_eflowJets::analyze_eflowJets(string inputFile, string outputFile) {
  cout << "Analyzing " << inputFile << " ..." << endl;
  analyze_delphesevents(inputFile, outputFile);
}

analyze_eflowJets::~analyze_eflowJets() {
  cout << "...End of analysis :)" << endl;
}
