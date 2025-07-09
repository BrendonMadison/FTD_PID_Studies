#pragma once
#include "HelixExtrapolator.h"  // for HitPoint
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <tuple>
#include "TMath.h"
//#include <ROOT/RMath.hxx>  // ROOT::Math::landau_pdf

constexpr double silicon_thickness_um = 300.0;
constexpr double nominal_MPV_keV = 75.0;  // for 300 µm normal incidence
constexpr double nominal_sigma_keV = 15.0; // estimated FWHM/2.35

constexpr double B_TESLA = 3.5;
constexpr double MULTIPLE_SCATTERING_THETA_RAD_PER_LAYER = 4e-6 * 125.0;  // rad * GeV , divide by energy to get radians
constexpr double TRACKER_LAYER_SPACING_M = (2.2212 - 0.22) / 5.0;

double csc(double x)
{
	return 1.0 / std::sin(x);
}

double sec(double x)
{
	return 1.0 / std::cos(x);
}

double sqr(double x)
{
	return x*x;
}

double cot(double x)
{
	return 1.0 / std::tan(x);
}

std::tuple<double, double, double, double> ExtrapolateXYAtZFromCalorimeter(
    double theta, double phi, double energy, double z_cal,
    double sigma_theta, double sigma_phi, double sigma_energy, double sigma_z_cal,
    double z_tracker, double sigma_z_tracker,
    int layers_passed)
{
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
    double theta, double phi, double energy, double z_cal,
    double sigma_theta, double sigma_phi, double sigma_energy, double sigma_z_cal,
    double sigma_z_tracker,
    double cell_size_x, double cell_size_y)
{
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

        // Determine how many layers are passed before this z
        int layers_passed = 0;
        for (double z_layer : unique_z) {
            if (z_layer < z) layers_passed++;
            else break;
        }

        // Step 3: Extrapolate predicted (x, y) and uncertainties
        auto [x_pred, y_pred, sigma_pred_x, sigma_pred_y] =
            ExtrapolateXYAtZFromCalorimeter(
                theta, phi, energy, z_cal,
                sigma_theta, sigma_phi, sigma_energy, sigma_z_cal,
                z, sigma_z_tracker,
                layers_passed);

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
    double beam_energy_GeV)
{
    // Step 1: Path length correction
    double path_length_um = silicon_thickness_um / std::cos(theta);  // assuming flat layers

    // Step 2: Scale MPV and width
    double mpv = nominal_MPV_keV * (path_length_um / silicon_thickness_um);
    double width = nominal_sigma_keV * std::sqrt(path_length_um / silicon_thickness_um);  // approximate

    // Step 3: Evaluate Landau log-likelihood
    double pdf_val = TMath::Landau(hit_energy_keV, mpv, width, true);  // normalized = true
    //double pdf_val = ROOT::Math::landau_pdf(hit_energy_keV, width, mpv);
    double logL = std::log(pdf_val + 1e-12);  // prevent log(0)
    
    return logL;
}

double ComputeTrackerPosAndEneLogLikelihood(
    const std::vector<HitPoint>& hits,
    double theta, double phi, double energy, double z_cal,
    double sigma_theta, double sigma_phi, double sigma_energy, double sigma_z_cal,
    double sigma_z_tracker,
    double cell_size_x, double cell_size_y)
{
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

	// Compute the corresponding log-likelihood for a valid hit to be a 1 MIP electron signal
	double ene_LL = ComputeMIPEnergyLogLikelihood(e_meas,theta,energy);
	
        // Determine how many layers are passed before this z
        int layers_passed = 0;
        for (double z_layer : unique_z) {
            if (z_layer < z) layers_passed++;
            else break;
        }

        // Step 3: Extrapolate predicted (x, y) and uncertainties
        auto [x_pred, y_pred, sigma_pred_x, sigma_pred_y] =
            ExtrapolateXYAtZFromCalorimeter(
                theta, phi, energy, z_cal,
                sigma_theta, sigma_phi, sigma_energy, sigma_z_cal,
                z, sigma_z_tracker,
                layers_passed);

        // Step 4: Combine uncertainties
        double sigma_x2 = sigma_pred_x * sigma_pred_x + sigma_meas_x * sigma_meas_x;
        double sigma_y2 = sigma_pred_y * sigma_pred_y + sigma_meas_y * sigma_meas_y;

        // Step 5: Residuals and likelihood
        double dx = x_meas - x_pred;
        double dy = y_meas - y_pred;

        logL += -0.5 * (dx * dx / (sigma_x2 + 1e-12) + dy * dy / (sigma_y2 + 1e-12));
	// Add the energy LL
	logL += ene_LL;
    }

    return logL;
}
