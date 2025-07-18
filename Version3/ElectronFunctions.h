#pragma once
#include "StructContainer.h"  // for HitPoint
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <tuple>
#include "TMath.h"
#include "TRandom3.h"
//#include <ROOT/RMath.hxx>  // ROOT::Math::landau_pdf

std::tuple<double, double, double, double> ExtrapolateXYAtZFromCalorimeter(
    double theta, double phi, double energy, DetGeo Detector, double z_tracker)
{
	//Import detector values
	double z_cal = Detector.Input_Cal.Cal_First_Pos_m;
	double sigma_theta = Detector.Sigma_Tht;
	double layers_passed = Detector.Input_Tracker.GetLayersPassed(z_tracker,theta);
	double sigma_z_tracker = Detector.Input_Tracker.Tracker_Thickness_um * 1e-6/std::sqrt(12.0);
	double sigma_z_cal = Detector.Input_Cal.Active_Thickness_um * 1e-6/std::sqrt(12.0);
	double sigma_energy = Detector.Sigma_Ene;
	double sigma_phi = Detector.Sigma_Phi;
	
	//Import physics values
	double B_TESLA = Detector.Input_Tracker.Mag_Field_tesla;
	double MULTIPLE_SCATTERING_THETA_RAD_PER_LAYER = Detector.Input_Tracker.Tracker_MS_Per_Layer * Detector.Input_Tracker.Tracker_Nominal_Ene;
	double TRACKER_LAYER_SPACING_M = Detector.Input_Tracker.Tracker_Layer_Spacing_m;
	double nominal_MPV_keV = Detector.Input_Tracker.Tracker_Nominal_MPV_keV;
	double nominal_sigma_keV = Detector.Input_Tracker.Tracker_Nominal_Wid_keV;
	double silicon_thickness_um = Detector.Input_Tracker.Tracker_Thickness_um;
	
	// All lengths should be in meters!
    // Step 1: Compute pT, R, omega
    double pt = energy * std::sin(theta);
    double R = pt / (0.3 * B_TESLA);
    double omega = 1.0 / R;

    // Step 2: Extrapolate to phi0 at origin using z_cal
    double dz = z_cal;
    double s_cal = dz / std::tan(theta);
    double delta_phi_cal = omega * s_cal;
    double phi0 = phi - delta_phi_cal;

    // Step 3: Extrapolate to tracker z
    double s_tracker = z_tracker / std::tan(theta);
    double delta_phi = omega * s_tracker;
    double phi_total = phi0 + delta_phi;

    double x = z_tracker * std::tan(theta) * std::cos(phi_total);
    double y = z_tracker * std::tan(theta) * std::sin(phi_total);

    // Step 4: Propagate uncertainties (assuming no higher order or covariances)

	// for the dx/dphi_total we need to expand phi_total = phi0 + delta_phi -> 
	// -> delta_phi = z_tracker / tan(theta) * omega = z_tracker / tan(theta) * 0.3 * B_TESLA / (energy * sin(theta))
	// -> phi0 = phi - delta_phi_cal = phi - z_cal / tan(theta) * 0.3 * B_TESLA / (energy * sin(theta))
	// so phi_total = phi + (z_tracker - z_cal) * 0.3 * B_TESLA / (energy * sin(theta) * tan(theta))
	// phi_total = phi + (z_tracker - z_cal) * DTERM
	// We simplify things:
	double ATERM = z_tracker * std::tan(theta);
	double BTERM = z_tracker-z_cal;
	double CTERM = 0.3 * B_TESLA / (energy);
	double DTERM = CTERM / (std::sin(theta) * std::tan(theta));
	double FTERM = DTERM * BTERM * energy;
	
	// dx/dz_tracker is equal to:
	double dx_dz_tracker = std::tan(theta) * std::sin(phi + DTERM * (BTERM) ) + 
						   DTERM * z_tracker * std::cos(phi + DTERM * (BTERM) );
	// dy/dz_tracker is equal to:
	double dy_dz_tracker = std::tan(theta) * std::cos(phi + DTERM * (BTERM) ) - 
						   DTERM * z_tracker * std::sin(phi + DTERM * (BTERM) );

    // dx/dz_cal is equal to:
	double dx_dz_cal = DTERM * ATERM * std::sin(phi + DTERM * (BTERM));
    // dx/dz_cal is equal to:
	double dy_dz_cal = -1.0 * DTERM * ATERM * std::sin(phi + DTERM * (BTERM));

    // dx/dtheta is equal to:
	double dx_dtheta = z_tracker * sec(theta) * (sec(theta) * std::cos(CTERM * BTERM * cot(theta) * csc(theta) + phi) +
						     CTERM*BTERM*( sqr(cot(theta)) + sqr(csc(theta)))*std::sin(CTERM * BTERM*cot(theta)*csc(theta) + phi));
    // dy/dtheta is equal to:
	double dy_dtheta = z_tracker * sec(theta) * (sec(theta) * std::sin(CTERM * BTERM * cot(theta) * csc(theta) + phi) -
						     CTERM*BTERM*( sqr(cot(theta)) + sqr(csc(theta)))*std::cos(CTERM * BTERM*cot(theta)*csc(theta) + phi));
	
	// dx/dphi is equal to:
	double dx_dphi = -1.0 * ATERM * std::sin(phi + DTERM * (BTERM));
	// dy/dphi is equal to:
	double dy_dphi = ATERM * std::cos(phi + DTERM * (BTERM));
	
	// dx/dE is equal to:
	double dx_dE = ATERM * FTERM * sin(phi + FTERM/energy) / sqr(energy);
	// dy/dE is equal to:
	double dy_dE = -1.0 * ATERM * FTERM * cos(phi + FTERM/energy) / sqr(energy);
	
	// Step 5: Combined variances with their associated sigma
	double var_x_noMS = sqr(dx_dE) * sqr(sigma_energy) + sqr(dx_dphi) * sqr(sigma_phi) + 
	                    sqr(dx_dtheta) * sqr(sigma_theta) + sqr(dx_dz_cal) * sqr(sigma_z_cal) + 
                        sqr(dx_dz_tracker) * sqr(sigma_z_tracker);
	double var_y_noMS = sqr(dy_dE) * sqr(sigma_energy) + sqr(dy_dphi) * sqr(sigma_phi) + 
	                    sqr(dy_dtheta) * sqr(sigma_theta) + sqr(dy_dz_cal) * sqr(sigma_z_cal) + 
                        sqr(dy_dz_tracker) * sqr(sigma_z_tracker);

	// Step 6: Add multiple scattering contribution for number of tracker layers passed
	// we assume that, for each layer passed, we add sqrt(layers_passed * (amount of theta deflection by MS)^2)
	double sigma_ms = MULTIPLE_SCATTERING_THETA_RAD_PER_LAYER / energy;
	// then, if the layers are spaced out an amount of layer_spacing (in m)
	// we expect an amount of change in x and y from MS
	// we approximate this to only depend on leading order in theta!
	// The following are propagated errors for one layer
	double dx_dms = TRACKER_LAYER_SPACING_M * std::cos(phi_total) * sqr(sec(theta));
	double dy_dms = TRACKER_LAYER_SPACING_M * std::sin(phi_total) * sqr(sec(theta));
	// Then we include for number of layers passed. So when layers_passed = 0 then its the no MS case
	double var_x_MS = var_x_noMS + sqr(dx_dms) * sqr(sigma_ms) * layers_passed;
	double var_y_MS = var_y_noMS + sqr(dy_dms) * sqr(sigma_ms) * layers_passed;

    return std::make_tuple(x, y, std::sqrt(var_x_MS), std::sqrt(var_y_MS));
}

double ComputeTrackerHitPositionLogLikelihood(
    const std::vector<HitPoint>& hits,
    double theta, double phi, double energy, DetGeo Detector)
{
	//Import detector values
	double z_cal = Detector.Input_Cal.Cal_First_Pos_m;
	double sigma_theta = Detector.Sigma_Tht;
	double sigma_z_tracker = Detector.Input_Tracker.Tracker_Thickness_um * 1e-6/std::sqrt(12.0);
	double sigma_z_cal = Detector.Input_Cal.Active_Thickness_um * 1e-6/std::sqrt(12.0);
	double sigma_energy = Detector.Sigma_Ene;
	double sigma_phi = Detector.Sigma_Phi;
	double cell_size_x = Detector.XY_size;
	double cell_size_y = Detector.XY_size;
	
	
    double logL = 0.0;
    double sigma_meas_x = cell_size_x / std::sqrt(12.0);
    double sigma_meas_y = cell_size_y / std::sqrt(12.0);

    // Step 1: Extract unique z-layer values and sort them
    std::set<double> unique_z_set;
    for (const auto& hit : hits)
        unique_z_set.insert(hit.z);

    std::vector<double> unique_z(unique_z_set.begin(), unique_z_set.end());
    std::sort(unique_z.begin(), unique_z.end());  // ensure ordered

    // Step 2: Loop over hits
    for (const auto& hit : hits) {
        double z = hit.z;
        double x_meas = hit.x();
        double y_meas = hit.y();
	double e_meas = hit.GetEne();
	if(e_meas < 1.0) //To do: change this for an actual threshold value, instead of using 1 keV
	  {
	    logL += std::log(1e-4);
	    continue;
	  }//Assign a 10^-4 chance (0.01%) that an electron hit is too low energy to be above threshold...

        // Step 3: Extrapolate predicted (x, y) and uncertainties
        auto [x_pred, y_pred, sigma_pred_x, sigma_pred_y] =
            ExtrapolateXYAtZFromCalorimeter(
                theta, phi, energy, Detector, z);

        // Step 4: Combine uncertainties
        double sigma_x2 = sigma_pred_x * sigma_pred_x + sigma_meas_x * sigma_meas_x;
        double sigma_y2 = sigma_pred_y * sigma_pred_y + sigma_meas_y * sigma_meas_y;

        // Step 5: Residuals and likelihood
        double dx = x_meas - x_pred;
        double dy = y_meas - y_pred;
	
        logL += -0.5 * (dx * dx / (sigma_x2 + 1e-12) + dy * dy / (sigma_y2 + 1e-12));
    }

    return logL;
}

double ComputeMIPEnergyLogLikelihood(
    double hit_energy_keV,
    double theta,
    double beam_energy_GeV,
	DetGeo Detector)
{
    // Step 1: Path length correction
    double path_length_um = Detector.Input_Tracker.Tracker_Thickness_um / std::cos(theta);  // assuming flat layers

    // Step 2: Scale MPV and width
    double mpv = Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta);
    double width = Detector.Input_Tracker.Tracker_Nominal_Wid_keV * std::sqrt(1.0 / std::cos(theta));  // approximate

    // Step 3: Evaluate Landau log-likelihood
    double pdf_val = TMath::Landau(hit_energy_keV, mpv, width, true);  // normalized = true
    //double pdf_val = ROOT::Math::landau_pdf(hit_energy_keV, width, mpv);
    double logL = std::log(pdf_val + 1e-12);  // prevent log(0)
    
    return logL;
}

double ComputeEleMIPEnergyLogLikelihood(
	const std::vector<HitPoint>& hits,
    double theta,
    double beam_energy_GeV,
	DetGeo Detector)
{
	double ene_LL = 0.0;
    for (const auto& hit : hits) {
		double e_meas = hit.GetEne();
	if(e_meas < 1.0) //To do: change this for an actual threshold value, instead of using 1 keV
	  {
	    ene_LL += std::log(1e-4);
	    continue;
	  }//Assign a 10^-4 chance (0.01%) that an electron hit is too low energy to be above threshold...

	// Compute the corresponding log-likelihood for a valid hit to be a 1 MIP electron signal
	ene_LL += ComputeMIPEnergyLogLikelihood(e_meas,theta,beam_energy_GeV,Detector);
	}
    
    return ene_LL;
}

double ComputeEleCountLogLikelihood(
    const std::vector<HitPoint>& hits,
    double theta, double phi, double energy, DetGeo Detector)
{
    const double z_tol = Detector.Input_Tracker.Tracker_Layer_Spacing_m * 0.10;  // 10% tolerance
    const int num_layers = Detector.Input_Tracker.Tracker_Num_Layers;
    const double first_z = Detector.Input_Tracker.Tracker_First_Pos_m;
    const double spacing = Detector.Input_Tracker.Tracker_Layer_Spacing_m;

    // Step 1: Group hits by expected layer index
    std::vector<std::vector<const HitPoint*>> layer_hits(num_layers);

    for (const auto& hit : hits) {
        double z_cm = hit.z;  // convert to cm
        for (int i = 0; i < num_layers; ++i) {
            double expected_z_cm = first_z + i * spacing;
	    //std::cout << z_cm << " , " << expected_z_cm << std::endl;
            if (std::abs(z_cm - expected_z_cm) < z_tol) {
                layer_hits[i].push_back(&hit);
                break;
            }
        }
    }

    // Step 2: Loop over each layer and compute log-likelihood
    double logL = 0.0;
    for (int i = 0; i < num_layers; ++i) {
        const auto& hits_in_layer = layer_hits[i];

        if (hits_in_layer.empty()) {
            // No hit → assign penalty
            logL += std::log(1e-16);  // 0.01% probability
            continue;
        }

        if (hits_in_layer.size() == 1) {
	  //std::cout << "One hit in layer!" << std::endl;
            double e_keV = hits_in_layer[0]->GetEne();
            if (e_keV < 0.1 * Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta)
		|| e_keV > 10.0 * Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta))
	      { // energy threshold of 10% of MIP energy
		// and upper threshold of 10 MIP energy
                logL += std::log(1e-16);
	      }
	    else {
                logL += ComputeMIPEnergyLogLikelihood(e_keV, theta, energy, Detector);
            }
        } else {
	  //std::cout << "Multiple hits!" << std::endl;
            // More than one hit: choose best one, penalize others
            double bestLL = -1e9;
            for (const auto* h : hits_in_layer) {
                double e_keV = h->GetEne();
                if (e_keV < Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta)
		    || e_keV > 10.0 * Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta)) continue;
                double ll = ComputeMIPEnergyLogLikelihood(e_keV, theta, energy, Detector);
                bestLL = std::max(bestLL, ll);
	    }

            if (bestLL == -1e9) {
                logL += std::log(1e-16);  // none were valid hits
            } else {
                logL += bestLL + std::log(8.0);  // penalty for multiple hits
            }
        }
    }
    return logL;
}

double ComputeSpeedOfLightLogLikelihood(
    const std::vector<HitPoint>& hits,
    double theta, DetGeo Detector)
{
    const double c_cm_per_ns = 29.9792458;
    const double mu = 1.0 / (c_cm_per_ns * std::cos(theta));  // expected t/z

    TRandom3 *rt = new TRandom3();
    TRandom3 *rz = new TRandom3();

    double sigma_time_ns = Detector.Input_Tracker.Tracker_Time_Res_ns;
    double sigma_z_cm = Detector.Input_Tracker.Tracker_Thickness_um * 1e-4/std::sqrt(12.0);
    double sigma_theta_rad = Detector.Sigma_Tht;
    
    // Theta error propagation (affects uncertainty on mu)
    double sigma_mu2 = 0.0;
    if (sigma_theta_rad > 0.0) {
        double dmu_dtheta = std::sin(theta) / (c_cm_per_ns * std::cos(theta) * std::cos(theta));
        sigma_mu2 = dmu_dtheta * dmu_dtheta * sigma_theta_rad * sigma_theta_rad;
    }

    double logL = 0.0;

    // Step 1: Loop over hits
    for (const auto& hit : hits) {
      // Step 2: Import hit values and blur them
      double z = hit.z + rz->Gaus(0.0,sigma_z_cm);
      double time = hit.GetTime() + rt->Gaus(0.0,sigma_time_ns);
      double e_meas = hit.GetEne();
      // Step 3: Reject empty or too-low energy hits
      if(e_meas < 0.1 * Detector.Input_Tracker.Tracker_Nominal_MPV_keV / std::cos(theta))
	{
	  logL += std::log(1e-4);
	  continue;
	}//Assign a 10^-4 chance (0.01%) that an electron hit is too low energy to be above threshold...	
	
      double ratio = time / z; //The inverse speed of light like ratio

      // Error propagation for ratio
      double var = (sigma_time_ns * sigma_time_ns) / (z * z)
	+ (time * time * sigma_z_cm * sigma_z_cm) / (z * z * z * z)
	+ sigma_mu2;

      if (var < 1e-12) var = 1e-12; // To counter erroneously small values
	
      double dratio = ratio - mu;
      logL += -0.5 * (dratio * dratio / var + std::log(var));
    }

    return logL;
}
