#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include "HelixExtrapolator.h"
#include "HelixPropagator.h"
// Add default constructor to fix std::map value type initialization
// struct HitPoint {
//     double theta = 0.0;
//     double phi = 0.0;
//     double z = 0.0;

//     HitPoint() = default;
//     HitPoint(double t, double p, double zz) : theta(t), phi(p), z(zz) {}

//     static std::vector<double> extractTheta(const std::vector<HitPoint>& hits) {
//         std::vector<double> out(hits.size());
//         for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].theta;
//         return out;
//     }

//     static std::vector<double> extractPhi(const std::vector<HitPoint>& hits) {
//         std::vector<double> out(hits.size());
//         for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].phi;
//         return out;
//     }

//     static std::vector<double> extractZ(const std::vector<HitPoint>& hits) {
//         std::vector<double> out(hits.size());
//         for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].z;
//         return out;
//     }
// };

void HelixFitter(char *INFILE, char *TNAME, char *OFILE, char *ONAME) {
  TFile *f = TFile::Open(INFILE);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file " << INFILE << std::endl;
    return;
  }

  TTree *t = dynamic_cast<TTree*>(f->Get(TNAME));
  if (!t) {
    std::cerr << "Cannot find TTree " << TNAME << " in file." << std::endl;
    return;
  }

  //Calorimeter and tracker values:
  double ENERES_stoch = 0.04; //4%/sqrt(E)
  double sigma_ene = 0.0; //to be calculated
  double z_cal_thick = 1.0e-3, z_track_thick = 0.3e-3, xy_track_size = 100e-6;
  double sigma_tht = 10e-6, sigma_phi = 100e-5, sigma_z_track = z_track_thick/std::sqrt(12.0);
  double sigma_z_cal = z_cal_thick/std::sqrt(12.0), sigma_xy_track = xy_track_size/std::sqrt(12.0);
  double phi_cal = 0.0; // to be calculated
  double z_cal_dist = 2.5, z_track_dist_first = 0.22, z_track_dist_last = 2.221;
  int n_tracker_layers = 5;
  double layer_spacing = (z_track_dist_last - z_track_dist_first) / (n_tracker_layers - 1);
  double z_tolerance = 0.10 * layer_spacing;  // 10% tolerance in cm

  //Reconstruction values
  double mean_ele_Lhood = -12000.0, std_ele_Lhood = 2000.0, e_score = 0.0;
  
  std::vector<double> *tht = nullptr, *phi = nullptr, *layz = nullptr;
  std::vector<double> *track_tht = nullptr, *track_phi = nullptr, *track_z = nullptr;
  std::vector<double> *ene = nullptr, *track_ene = nullptr;
  std::vector<double> *v_beamTH = nullptr;
  std::vector<double> *v_beamPH = nullptr;
  std::vector<double> *v_beamEne = nullptr;

  t->SetBranchAddress("TrackTht", &track_tht);
  t->SetBranchAddress("beamEnergy", &v_beamEne);
  t->SetBranchAddress("beamTH", &v_beamTH);
  t->SetBranchAddress("beamPH", &v_beamPH);
  t->SetBranchAddress("TrackPhi", &track_phi);
  t->SetBranchAddress("TrackZ", &track_z);
  t->SetBranchAddress("TrackEne", &track_ene);
  t->SetBranchAddress("LayerX", &tht);
  t->SetBranchAddress("LayerY", &phi);
  t->SetBranchAddress("LayerZ", &layz);
  t->SetBranchAddress("LayerEdep", &ene);

  // Output file and tree
  TFile *fout = new TFile(OFILE, "RECREATE");
  TTree *tout = new TTree(ONAME, "Track hits with extrapolated points");

  std::vector<double> extr_theta, extr_phi, extr_z;
  Int_t event = 0;
  double R,z0,phi0,dR,dz0,dphi0,p,dp,beamTH,beamPH,beamEne;
  double MIP,dMIP,Slope,dSlope,Curv,dCurv;
  double eLhood_helix, eLhood_both;
  //tout->Branch("ExTheta", &extr_theta);
  //tout->Branch("ExPhi", &extr_phi);
  tout->Branch("beamEne", &beamEne);
  tout->Branch("beamTH", &beamTH);
  tout->Branch("beamPH", &beamPH);
  //tout->Branch("ExZ", &extr_z);
  tout->Branch("Event", &event);
  //tout->Branch("R", &R);tout->Branch("dR", &dR);
  //tout->Branch("z0", &z0);tout->Branch("dz0", &dz0);
  //tout->Branch("phi0", &phi0);tout->Branch("dphi0", &dphi0);
  //tout->Branch("p", &p);tout->Branch("dp", &dp);
  //tout->Branch("MIP",&MIP);tout->Branch("dMIP",&dMIP);
  //tout->Branch("Slope",&Slope);tout->Branch("dSlope",&dSlope);
  //tout->Branch("Curv",&Curv);tout->Branch("dCurv",&dCurv);
  tout->Branch("eLhood_helix",&eLhood_helix);
  tout->Branch("eLhood_both",&eLhood_both);
  tout->Branch("e_score",&e_score);
  
  Long64_t nev = t->GetEntries();
  for (Long64_t i = 0; i < nev; ++i) {
    t->GetEntry(i);

    // Import event beam values
    beamEne = v_beamEne->at(0);
    beamTH = v_beamTH->at(0);
    beamPH = v_beamPH->at(0);

    // Calculate energy resolution and phi at calorimeter (2.5 m)
    sigma_ene = ENERES_stoch /std::sqrt(beamEne);
    phi_cal = beamPH + (0.3 * 3.5)/(beamEne * std::sin(beamTH)) * 2.5/std::tan(beamTH); //assuming calorimeter is measured at its start of 2.5 meters

    // Convert tracker hits to HitPoint list
    std::vector<HitPoint> hits;
    for (size_t j = 0; j < tht->size(); ++j) {
      double zval = (*track_z)[j];
      double thtval = (*track_tht)[j];
      double phival = (*track_phi)[j];
      double eneval = (*track_ene)[j];
      if (zval <= 10.0 || zval >= 240.0 || abs( thtval - beamTH) > 0.05 || abs( phival - beamPH ) > 0.05) continue;
      hits.emplace_back(thtval, phival, zval/100.0, eneval); //add to hitpoint object, convert to meters. 
    }

    //Whole and/or photon compensator code
    //Here we check if there are tracker layers with no hits
    //If there are, we add a zero energy hit

    // Track which layers are already hit
    std::vector<bool> has_hit(n_tracker_layers, false);

    // Mark hits by z-layer proximity
    for (const auto& hit : hits) {
      for (int k = 0; k < n_tracker_layers; ++k) {
        double expected_z_cm = z_track_dist_first + k * layer_spacing;
        if (std::abs(hit.z * 100.0 - expected_z_cm) < z_tolerance) {
	  has_hit[k] = true;
	  break;
        }
      }
    }

    // Insert zero energy hits where missing
    for (int k = 0; k < n_tracker_layers; ++k) {
      if (!has_hit[k]) {
        double missing_z = (z_track_dist_first + k * layer_spacing) / 100.0;  // m
        hits.emplace_back(beamTH, beamPH, missing_z, 0.0);  // zero-energy dummy hit
      }
    }

    
    // Apply extrapolation
    //MultipleScatteringRotationCorrection(hits);
    //MultipleScatteringSymmetricCorrection(hits);
    //ConstantZExtrapolation(hits);
    //DifferentZExtrapolation(hits);

    //Get LL for only helix
    eLhood_helix = ComputeTrackerHitPositionLogLikelihood(hits,beamTH,beamPH,beamEne,
							  z_cal_dist,sigma_tht,sigma_phi,
							  sigma_ene,sigma_z_cal,sigma_z_track,
							  xy_track_size,xy_track_size);

    //Get LL for helix and hit energy combined
    eLhood_both = ComputeTrackerPosAndEneLogLikelihood(hits,beamTH,beamPH,beamEne,
							  z_cal_dist,sigma_tht,sigma_phi,
							  sigma_ene,sigma_z_cal,sigma_z_track,
							  xy_track_size,xy_track_size);

    //Convert likelihood to electron score under assumption that electrons have gaussian mean and spread defined above
    e_score = 1.0 / (1.0 + std::exp((eLhood_helix - mean_ele_Lhood) / std_ele_Lhood));
    event = i;
    
    tout->Fill();
  }

  fout->cd();
  tout->Write();
  fout->Close();
  std::cout << "Wrote electron helix log-likelihood to " << OFILE << std::endl;
}


