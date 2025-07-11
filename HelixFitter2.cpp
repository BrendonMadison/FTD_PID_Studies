#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include "HelixExtrapolator.h"
#include "HelixPropagator.h"
#include "ShowerFitterHelper.h"
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
  double mean_gam_Lhood = -32.5, std_gam_Lhood = 20.0, g_score = 0.0;
  
  std::vector<double> *tht = nullptr, *phi = nullptr, *layz = nullptr;
  std::vector<double> *track_tht = nullptr, *track_phi = nullptr, *track_z = nullptr;
  std::vector<double> *ene = nullptr, *track_ene = nullptr;
  std::vector<double> *v_beamTH = nullptr;
  std::vector<double> *v_beamPH = nullptr;
  std::vector<double> *v_beamEne = nullptr;
  std::vector<double> *v_pdgID = nullptr;

  t->SetBranchAddress("TrackTht", &track_tht);
  t->SetBranchAddress("beamEnergy", &v_beamEne);
  t->SetBranchAddress("beamTH", &v_beamTH);
  t->SetBranchAddress("beamPH", &v_beamPH);
  t->SetBranchAddress("pdgID", &v_pdgID);
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
  double R,z0,phi0,dR,dz0,dphi0,p,dp,beamTH,beamPH,beamEne,pdgID;
  double eLhood_helix, eLhood_both;
  double xinter, xinter_err, slope, slope_err, chi2;
  double gamLhood_conv, gamLhood_ene;
  //tout->Branch("ExTheta", &extr_theta);
  //tout->Branch("ExPhi", &extr_phi);
  tout->Branch("beamEne", &beamEne);
  tout->Branch("beamTH", &beamTH);
  tout->Branch("beamPH", &beamPH);
  tout->Branch("pdgID", &pdgID);
  //tout->Branch("ExZ", &extr_z);
  tout->Branch("Event", &event);
  //tout->Branch("R", &R);tout->Branch("dR", &dR);
  //tout->Branch("z0", &z0);tout->Branch("dz0", &dz0);
  //tout->Branch("phi0", &phi0);tout->Branch("dphi0", &dphi0);
  //tout->Branch("p", &p);tout->Branch("dp", &dp);
  //tout->Branch("MIP",&MIP);tout->Branch("dMIP",&dMIP);
  //tout->Branch("Slope",&Slope);tout->Branch("dSlope",&dSlope);
  //tout->Branch("Curv",&Curv);tout->Branch("dCurv",&dCurv);
  tout->Branch("xinter",&xinter);tout->Branch("xinter_err",&xinter_err);
  tout->Branch("slope",&slope);tout->Branch("slope_err",&slope_err);
  tout->Branch("eLhood_helix",&eLhood_helix);
  tout->Branch("eLhood_both",&eLhood_both);
  tout->Branch("e_score",&e_score);
  tout->Branch("gamLhood_conv",&gamLhood_conv);
  tout->Branch("gamLhood_ene",&gamLhood_ene);
  tout->Branch("g_score",&g_score);
  
  // Generate path lengths and then
  // load in the sum of 5 electron MIP PDF
  //
  //
  // *****************************************************************************
  // *****************************************************************************
  // ********THIS ASSUMES THE ENTIRE DATASET HAS ONE THETA VALUE , ***************
  // ********OTHERWISE YOU NEED TO REMAKE THE SUMMED MIP PDF EACH TIME!!! ********
  // *****************************************************************************
  // *****************************************************************************
  //
  //
  t->GetEntry(0);
  beamTH = v_beamTH->at(0); // load one event in so we get beamTH
  std::vector<double> path_lengths = { 300e-6 / std::cos(beamTH) , 300e-6 / std::cos(beamTH) , 300e-6 / std::cos(beamTH) , 300e-6 / std::cos(beamTH) , 300e-6 / std::cos(beamTH) }; 
  TH1D* old = (TH1D*)gDirectory->Get("landau_sum_pdf");
  if (old) delete old;
  TH1D* MIP_PDF = GenerateLandauSumPDF(path_lengths);
  
  Long64_t nev = t->GetEntries();
  for (Long64_t i = 0; i < nev; ++i) {
    t->GetEntry(i);

    // Import event beam values
    beamEne = v_beamEne->at(0);
    beamTH = v_beamTH->at(0);
    beamPH = v_beamPH->at(0);
    pdgID = v_pdgID->at(0);
    
    // Calculate energy resolution and phi at calorimeter (2.5 m)
    sigma_ene = ENERES_stoch /std::sqrt(beamEne);
    phi_cal = beamPH + (0.3 * 3.5)/(beamEne * std::sin(beamTH)) * 2.5/std::tan(beamTH); //assuming calorimeter is measured at its start of 2.5 meters
    
    // Convert tracker hits to HitPoint list
    std::vector<HitPoint> hits;
    for (size_t j = 0; j < track_tht->size(); ++j) {
      double zval = (*track_z)[j];
      double thtval = (*track_tht)[j];
      double phival = (*track_phi)[j];
      double eneval = (*track_ene)[j];
      //if (zval <= 10.0 || zval >= 245.0 || abs( thtval - beamTH) > 0.001 || abs( phival - beamPH ) > 0.001) continue;
      if (zval <= 10.0 || zval >= 245.0 || abs( thtval - beamTH) > 0.001) continue;
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

    // Fill new hits vector to include layer values with tracker values
    std::vector<HitPoint> hits_lay;
    for (const auto& hit : hits) {
      hits_lay.emplace_back(hit.GetTheta(),hit.GetPhi(),hit.GetZ(),hit.GetEne());
    }
    // Layer values
    for (size_t j = 0; j < track_tht->size(); ++j) {
      double zval = (*track_z)[j];
      double thtval = (*track_tht)[j];
      double phival = (*track_phi)[j];
      double eneval = (*track_ene)[j];
      if (zval <= 245.0 || zval > 253.8) continue; //reject tracker layers and later calorimeter layers
      //we only want the first 2 X0 of the calorimeter
      //if (abs( thtval - beamTH) > 0.001 || abs( phival - beamPH ) > 0.001)
      if (abs( thtval - beamTH) > 0.001)
        {
          // this layer fails the angle cut, so put zero energy
          hits_lay.emplace_back(thtval,phival,zval/100.0,0.0);
        }
      else
        {
          hits_lay.emplace_back(thtval, phival, zval/100.0, eneval); //add to hitpoint object, convert to meters.
        }
    }

    // Compute photon values!
    int n_shower_vals = 3;
    // Compute the Shower start summer (X0 of start, hit energies, integrated radiation lengths, path lengths)
    ShowerStartSummary summary = FindShowerStartFromHits(hits_lay, n_shower_vals,beamTH);

    std::cout << "*" << std::endl;
    for(Int_t i = 0; i < summary.int_rad_vals.size(); i++)
      {
	std::cout << summary.energies[i] << " , " << summary.int_rad_vals[i] << " , " << summary.path_length[i] << "\n";
      }
    // Compute the log-likelihood of photon convertion
    // First have to check if the convertion occurs later than 2 X0, at which point we give this a value of 0.0, or 100% chance a photon;
    if (summary.int_rad_vals.size() < 1) gamLhood_conv = 0.0;
    else gamLhood_conv = ComputePhotonConversionLogLikelihood(summary.int_rad_vals[0], summary.int_rad_vals_err[0]);

    // Generate the MIP histogram for the potentially different path length shower start
    TH1D* s_old = (TH1D*)gDirectory->Get("shower_sum_pdf");
    if (s_old) delete s_old;
    TH1D* SHOWER_MIP_PDF = GenerateShowerLandauSumPDF(summary);

    // Now we can compute the individual hit and sum of hit energy likelihood , under assumption of 2 MIPs
    gamLhood_ene = ComputeLayered2MIPEnergyLogLikelihood(summary, SHOWER_MIP_PDF, beamEne);

    //Convert likelihood to gamma score under assumption that electrons have gaussian mean and spread defined above
    g_score = 1.0 - 1.0 / (1.0 + std::exp((gamLhood_conv + gamLhood_ene - mean_gam_Lhood) / std_gam_Lhood));
    
    auto profile = ComputeEnergyProfileFromHits(hits_lay); // Compute energy profile
    ShowerFitResult result = FitPhotonShowerStartWithTMinuit(profile); // Do initial shower fit
    xinter = result.x_intercept;xinter_err = result.x_intercept_err;
    slope = result.slope;slope_err = result.slope_err;
    
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
    // eLhood_both = ComputeTrackerPosAndEneLogLikelihood(hits,beamTH,beamPH,beamEne,
    // 							  z_cal_dist,sigma_tht,sigma_phi,
    // 							  sigma_ene,sigma_z_cal,sigma_z_track,
    // 							  xy_track_size,xy_track_size);

    //Get LL for helix, hit energy and hit energy sum combined
    eLhood_both = ComputeTrackerPosAndEneAndEneSumLogLikelihood(hits,MIP_PDF,
								beamTH,beamPH,beamEne,
								z_cal_dist,sigma_tht,sigma_phi,
								sigma_ene,sigma_z_cal,sigma_z_track,
								xy_track_size,xy_track_size);
    //Convert likelihood to electron score under assumption that electrons have gaussian mean and spread defined above
    e_score = 1.0 - 1.0 / (1.0 + std::exp((eLhood_helix - mean_ele_Lhood) / std_ele_Lhood));
    event = i;
    
    tout->Fill();
  }

  fout->cd();
  tout->Write();
  fout->Close();
  std::cout << "Wrote electron helix log-likelihood to " << OFILE << std::endl;
}


