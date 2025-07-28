#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include "StructContainer.h"
#include "ElectronFunctions.h"
#include "PhotonFunctions.h"
#include <fstream>
#include "TString.h"  // for Form()

double CleanValue(double x, double default_val = 0.0) {
    if (std::isnan(x) || std::isinf(x))
        return default_val;
    return x;
}

void HelixFitter(char *INFILE, char *TNAME, char *OFILE, char *ONAME, int i_PDGID) {
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

  // Open file for writing
  std::ofstream bdt_file(Form("%s.data",OFILE));

  if (!bdt_file.is_open()) {
    std::cerr << "Error: Could not open file for writing." << std::endl;
    return;
  }
  
  //Define the input tracker geometry
  //0.0005595970900951316 X0 from G10
  //0.0032017075773745993 X0 from Si
  //Standard ILD FTD values:
  TrackerGeo In_Trk(9.370,300.0,3.5,0.22,2.212,0.00376,5); //for 5 layers, 300 um Si, 100 um PCB
  In_Trk.SetMIPValues(80.0,15.0,45.6); //MPV (keV), Width (keV), nominal energy (GeV)
  In_Trk.SetMSValues(11e-6 * 45.6); // Multiple scattering (radians) * energy that the value is at (GeV)
  In_Trk.SetTimeRes(500.0e-12); //Set timing resolution for tracker
  //Define the input calorimeter geometry
  //Assuming 1mm PCB, 1mm Si, 0.521mm W, 0.5 mm Air
  CalGeo In_Cal(9.370,0.3504,1000.0,521.0,3.5,2.5,3.021e-3,1.0/6.0,240);
  In_Cal.SetMIPValues(80.0*10.0/3.0,15.0*std::sqrt(10.0/3.0),45.6);
  In_Cal.SetMSValues(11e-6 * 45.6);

  //50 layers 30 micron MAPS tracker
  //0.00032017075 X0 from Si , no G10
  //TrackerGeo In_Trk(9.370,30.0,3.5,0.22,2.212,0.000376,50); //for 5 layers, 300 um Si, 100 um PCB
  //In_Trk.SetMIPValues(80.0,15.0,45.6); //MPV (keV), Width (keV), nominal energy (GeV)
  //In_Trk.SetMSValues(11e-6 * 45.6 / 10.0); // Multiple scattering (radians) * energy that the value is at (GeV)
  //In_Trk.SetTimeRes(500.0e-12); //Set timing resolution for tracker
  //Define the input calorimeter geometry
  //Assuming 1mm PCB, 1mm Si, 0.521mm W, 0.5 mm Air
  //CalGeo In_Cal(9.370,0.3504,1000.0,521.0,3.5,2.5,3.021e-3,1.0/6.0,240);
  //In_Cal.SetMIPValues(80.0*10.0/3.0,15.0*std::sqrt(10.0/3.0),45.6);
  //In_Cal.SetMSValues(11e-6 * 45.6);
  
  //Combine into one detector geometry struct
  DetGeo In_Det(In_Trk,In_Cal);
  
  //Calorimeter and tracker values:
  // double ENERES_stoch = 0.04; //4%/sqrt(E)
  // double sigma_ene = 0.0; //to be calculated
  // double z_cal_thick = 1.0e-3, z_track_thick = 0.3e-3, xy_track_size = 100e-6;
  // double sigma_tht = 10e-6, sigma_phi = 100e-5, sigma_z_track = z_track_thick/std::sqrt(12.0);
  // double sigma_z_cal = z_cal_thick/std::sqrt(12.0), sigma_xy_track = xy_track_size/std::sqrt(12.0);
  // double phi_cal = 0.0; // to be calculated
  // double z_cal_dist = 2.5, z_track_dist_first = 0.22, z_track_dist_last = 2.221;
  // int n_tracker_layers = 5;
  // double layer_spacing = (z_track_dist_last - z_track_dist_first) / (n_tracker_layers - 1);
  // double z_tolerance = 0.10 * layer_spacing;  // 10% tolerance in cm

  //Reconstruction values
  double mean_ele_Lhood = -12000.0, std_ele_Lhood = 2000.0, e_score = 0.0;
  double mean_gam_Lhood = -32.5, std_gam_Lhood = 20.0, g_score = 0.0;
  
  std::vector<double> *tht = nullptr, *phi = nullptr,
    *layz = nullptr, *ene = nullptr;
  std::vector<double> *track_tht = nullptr, *track_phi = nullptr,
    *track_z = nullptr, *track_time = nullptr, *track_ene = nullptr;
  std::vector<double> *v_beamTH = nullptr, *v_beamPH = nullptr,
    *v_beamEne = nullptr, *v_pdgID = nullptr;

  t->SetBranchAddress("TrackTht", &track_tht);
  t->SetBranchAddress("TrackTime", &track_time);
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
  std::vector<double> v_deltaR, v_std_deltaR;
  std::vector<double> v_ene,v_pz,v_pt;
  double qene_sum, qpz_sum, qpt_sum;
  Int_t event = 0;
  double R,z0,phi0,dR,dz0,dphi0,p,dp,beamTH,beamPH,beamEne,pdgID;
  double eLhood_helix, eLhood_ene, eLhood_cnt, eLhood_tof;
  double xinter, xinter_err, slope, slope_err, chi2;
  double gamLhood_conv, gamLhood_ene;
  double track_enesum_tot , track_enesum_late;
  double enesum_mip_track_LL , track_latest_hit;
  double cal_enesum_tot , enesum_mip_cal_LL;
  double shower_ini, std_shower_ini;
  double deltaR, std_deltaR;
  double cal_LR_Asym, cal_UD_Asym, UD_LL, LR_LL;
  //tout->Branch("ExTheta", &extr_theta);
  //tout->Branch("ExPhi", &extr_phi);
  //tout->Branch("ExZ", &extr_z);
  tout->Branch("beamEne", &beamEne);tout->Branch("beamTH", &beamTH);
  tout->Branch("beamPH", &beamPH);tout->Branch("pdgID", &pdgID);
  tout->Branch("Event", &event);
  //tout->Branch("R", &R);tout->Branch("dR", &dR);
  //tout->Branch("z0", &z0);tout->Branch("dz0", &dz0);
  //tout->Branch("phi0", &phi0);tout->Branch("dphi0", &dphi0);
  //tout->Branch("p", &p);tout->Branch("dp", &dp);
  //tout->Branch("MIP",&MIP);tout->Branch("dMIP",&dMIP);
  //tout->Branch("Slope",&Slope);tout->Branch("dSlope",&dSlope);
  //tout->Branch("Curv",&Curv);tout->Branch("dCurv",&dCurv);
  tout->Branch("shower_ini",&shower_ini);tout->Branch("std_shower_ini",&std_shower_ini);
  tout->Branch("deltaR",&deltaR);tout->Branch("std_deltaR",&std_deltaR);
  tout->Branch("v_deltaR",&v_deltaR);tout->Branch("v_std_deltaR",&v_std_deltaR);
  tout->Branch("v_ene",&v_ene);tout->Branch("v_pz",&v_pz);
  tout->Branch("v_pt",&v_pt);
  tout->Branch("qene",&qene_sum);tout->Branch("qpz",&qpz_sum);tout->Branch("qpt",&qpt_sum);
  tout->Branch("cal_LR_Asym",&cal_LR_Asym);tout->Branch("cal_UD_Asym",&cal_UD_Asym);
  tout->Branch("LR_LL",&LR_LL);tout->Branch("UD_LL",&UD_LL);
  //tout->Branch("xinter",&xinter);tout->Branch("xinter_err",&xinter_err);
  //tout->Branch("slope",&slope);tout->Branch("slope_err",&slope_err);
  tout->Branch("eLhood_helix",&eLhood_helix);tout->Branch("eLhood_ene",&eLhood_ene);
  tout->Branch("eLhood_cnt",&eLhood_cnt);tout->Branch("eLhood_tof",&eLhood_tof);
  tout->Branch("e_score",&e_score);
  tout->Branch("gamLhood_conv",&gamLhood_conv);tout->Branch("gamLhood_ene",&gamLhood_ene);
  tout->Branch("g_score",&g_score);
  tout->Branch("track_enesum_tot",&track_enesum_tot);
  tout->Branch("track_enesum_late",&track_enesum_late);  
  tout->Branch("enesum_mip_track_LL",&enesum_mip_track_LL);
  tout->Branch("cal_enesum_tot",&cal_enesum_tot);
  tout->Branch("enesum_mip_cal_LL",&enesum_mip_cal_LL);
  tout->Branch("track_latest_hit",&track_latest_hit); //in ns
  
  t->GetEntry(0);
  beamTH = v_beamTH->at(0); // load one event in so we get beamTH
  beamEne = v_beamEne->at(0); // load one event in so we get beamTH

  //Set the detector resolution values
  //
  //
  // *****************************************************************************
  // *****************************************************************************
  // **********************DETECTOR RESOLUTION SETTINGS***************************
  // *****************************************************************************
  // *****************************************************************************
  //
  //
  In_Det.SetResolutions(beamEne * 0.04, 10.0e-6, 100.0e-6,
			In_Trk.Tracker_Thickness_um * 1e-6 / std::sqrt(12.0),
			1.0e-3,1.0e-3 / std::sqrt(12.0));
  //
  //
  
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
  //std::vector<double> path_lengths;
  //for(Int_t i = 0; i < In_Trk.Tracker_Num_Layers; i++) path_lengths.push_back(In_Trk.Tracker_Thickness_um/std::cos(beamTH));
  //TH1D* old = (TH1D*)gDirectory->Get("landau_sum_pdf");
  //if (old) delete old;
  //TH1D* MIP_PDF = GenerateLandauSumPDF(path_lengths);

  //
  // Load in the likelihood model for an electron
  // so that we can compute an electron score "e_score"
  //
  LikelihoodModel model = LoadLikelihoodModelSimple("Electron30uModel");
  
  Long64_t nev = t->GetEntries();
  for (Long64_t i = 0; i < nev; ++i) {
    t->GetEntry(i);

    // Import event beam values
    beamEne = v_beamEne->at(0);
    beamTH = v_beamTH->at(0);
    beamPH = v_beamPH->at(0);
    pdgID = v_pdgID->at(0);
    // Reset values
    track_enesum_late = 0.0;track_enesum_tot = 0.0;
    enesum_mip_track_LL = 0.0;track_latest_hit = 0.0;
    cal_enesum_tot = 0.0;enesum_mip_cal_LL = 0.0;
    qene_sum = 0.0;qpz_sum = 0.0;qpt_sum = 0.0;

    //reset deltaR vector values
    v_deltaR.clear();
    v_std_deltaR.clear();
    v_ene.clear();v_pt.clear();v_pz.clear();
    
    // Calculate energy resolution and phi at calorimeter (2.5 m)
    //sigma_ene = ENERES_stoch /std::sqrt(beamEne);
    //phi_cal = beamPH + (0.3 * 3.5)/(beamEne * std::sin(beamTH)) * 2.5/std::tan(beamTH); //assuming calorimeter is measured at its start of 2.5 meters
    
    // Convert tracker hits to HitPoint list
    DataLoader Data(track_tht,track_phi,track_z,track_time,track_ene,beamTH,beamPH,In_Det);

    // Compute the log-likelihood of photon convertion
    // First have to check if the convertion occurs later than 2 X0, at which point we give this a value of 0.0, or 100% chance a photon;
    //if (Data.Shower_Summary.int_rad_vals.size() < 1) gamLhood_conv = 0.0;
    //else gamLhood_conv = ComputePhotonConversionLogLikelihood(Data.Shower_Summary.int_rad_vals[0], Data.Shower_Summary.int_rad_vals_err[0]);

    // Generate the MIP histogram for the potentially different path length shower start
    //TH1D* s_old = (TH1D*)gDirectory->Get("shower_sum_pdf");
    //if (s_old) delete s_old;
    //TH1D* SHOWER_MIP_PDF = GenerateShowerLandauSumPDF(Data.Shower_Summary);

    // Now we can compute the individual hit and sum of hit energy likelihood , under assumption of 2 MIPs
    //gamLhood_ene = ComputeLayered2MIPEnergyLogLikelihood(Data.Shower_Summary, SHOWER_MIP_PDF, beamEne);

    //Convert likelihood to gamma score under assumption that electrons have gaussian mean and spread defined above
    //g_score = 1.0 - 1.0 / (1.0 + std::exp((gamLhood_conv + gamLhood_ene - mean_gam_Lhood) / std_gam_Lhood));

    //std::cout << i << std::endl;

    //Get LL for only helix , electron hypothesis
    eLhood_helix = ComputeTrackerHitPositionLogLikelihood(Data.Track_Hits,beamTH,
							  beamPH,beamEne,In_Det);

    //Get LL for electron 1 MIP hit track hypothesis
    eLhood_ene = ComputeEleMIPEnergyLogLikelihood(Data.Track_Hits,beamTH,
						  beamEne,In_Det);

    //Get LL for electron track count hit being equal to number of tracker layers
    //with penalties for missing hits
    //penalties for energies being too small or too large (10% of MIP and 1000% of MIP)
    //penalties for multiple hits within the acceptance
    eLhood_cnt = ComputeEleCountLogLikelihood(Data.Track_Hits,beamTH,
					      beamPH,beamEne,In_Det);

    //Get LL for speed of light consistency -- hypothesis being that the tracker hits must be consistent with speed of light
    eLhood_tof = ComputeSpeedOfLightLogLikelihood(Data.Track_Hits,beamTH,In_Det);

    //Get shower initial X0 given a MIP fraction threshold
    shower_ini = ComputeShowerStartX0(Data.All_Hits, In_Det, beamTH, 10.0);
    std_shower_ini = In_Det.Input_Cal.Cal_X0_Per_Layer/std::sqrt(12.0);

    //Get calorimeter and tracker delta-R
    auto [tmpR,tmpSR] = ComputeTrackClusterDeltaR(Data.All_Hits, beamTH, In_Det, true);
    tmpR = CleanValue(tmpR);
    tmpSR = CleanValue(tmpSR);
    deltaR = tmpR;
    std_deltaR = tmpSR;
    v_deltaR.push_back(tmpR);
    v_std_deltaR.push_back(tmpSR);

    //Get the quadrant deltaR and energy/momentum values
    for(Int_t ir = 1; ir < 5; ir++)
      {
	auto [qtmpR,qtmpSR] = ComputeQuadrantDeltaR(Data.All_Hits, beamTH, beamPH, ir, In_Det, true);
	qtmpR = CleanValue(qtmpR);
	qtmpSR = CleanValue(qtmpSR);
	v_deltaR.push_back(qtmpR);
	v_std_deltaR.push_back(qtmpSR);

	auto [qtmpE,qtmpPz,qtmpPt] = ComputeQuadrantEneStats(Data.All_Hits, beamTH, beamPH, ir, In_Det, true);
	qtmpE = CleanValue(qtmpE);
	qtmpPz = CleanValue(qtmpPz);
	qtmpPt = CleanValue(qtmpPt);
	v_ene.push_back(qtmpE);
	v_pz.push_back(qtmpPz);
	v_pt.push_back(qtmpPt);
	qene_sum += qtmpE;qpz_sum += qtmpPz;qpt_sum += qtmpPt;
      }

    //Compute the calorimeter shower asymmetry values
    //We average over two permutations with no diagonal sums
    //Left-right asymmetry
    cal_LR_Asym = 0.5*( (v_deltaR[4] - v_deltaR[1]) + (v_deltaR[3] - v_deltaR[2]) );
    //Up-down asymmetry
    cal_UD_Asym = 0.5*( (v_deltaR[1] - v_deltaR[2]) + (v_deltaR[4] - v_deltaR[3]) );

    //Compute the log-likelihood under the Molier (MS) shower spread theory
    LR_LL = ComputeMCSAngularLogLikelihood(cal_LR_Asym,beamEne,2.0);
    UD_LL = ComputeMCSAngularLogLikelihood(cal_UD_Asym,beamEne,2.0);
      
    //Print data out to the BDT input file
    bdt_file << Form("%.4e , %.4e , %.4e , %.4e , %.4e , %.4e , %.4e , %.4e , %.4e , %.4e , %i.\n", eLhood_helix, eLhood_ene, eLhood_cnt, eLhood_tof, shower_ini, deltaR, cal_LR_Asym, cal_UD_Asym, LR_LL , UD_LL , i_PDGID);
    
    //Convert likelihood to electron score under assumption that electrons
    //have a multidimensional mean, variance, covar, defined by
    //the previous four LL values.
    //You must first read in the electron model values from a csv!
    //Generally you will run this code first and generate the csv
    //(using idealized electron values)
    //then run again with the csv to get e_score values.
    std::vector<double> LLvec = {eLhood_ene, eLhood_helix, eLhood_cnt, eLhood_tof};
    e_score = ComputeMultivariateGaussianProb(LLvec, model);
    //e_score = ComputeIndependentGaussianProb(LLvec, model);
    if(e_score != 0) e_score = 1.0 - e_score;

    event = i;
    
    tout->Fill();
  }

  bdt_file.close();
  
  fout->cd();
  tout->Write();
  fout->Close();
  std::cout << "Wrote electron helix log-likelihood to " << OFILE << std::endl;
}

int main(int argc, char* argv[]) {
    // Example usage:
    // ./HelixFitter input.root treeName output.root outputTreeName 11

    if (argc < 6) {
      std::cerr << "Found " << argc-1 << " arguments! Need 5!\n";
      for(Int_t i = 1; i < argc; i++)
	{
	  std::cerr << argv[i] << "\n";
	}
        std::cerr << "Usage: " << argv[0] << " INFILE TNAME OFILE ONAME\n";
        return 1;
    }

    char* INFILE = argv[1];
    char* TNAME  = argv[2];
    char* OFILE  = argv[3];
    char* ONAME  = argv[4];
    int i_PDGID  = atoi(argv[5]);

    // Call your original logic here
    HelixFitter(INFILE, TNAME, OFILE, ONAME, i_PDGID);

    return 0;
}
