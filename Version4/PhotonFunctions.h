#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "StructContainer.h"

double ComputeShowerStartX0(const std::vector<HitPoint>& calo_hits, DetGeo Detector, double theta, double MIP_fraction) {
    // Bin energy by layer (Z → integrated X0)
    std::map<double, double> energy_by_x0;

    for (const auto& hit : calo_hits) {
      // screen out tracker hits and empty energy hits
      if (hit.ene <= 0.0 || hit.z < Detector.Input_Tracker.Tracker_Last_Pos_m + 0.1) continue;
        double x0 = Detector.GetIntegratedRadLen(hit.z, theta);  // in X0
        energy_by_x0[x0] += hit.ene;
    }

    if (energy_by_x0.empty()) return -1.0; // return -1 if no threshold passed

    // Sliding sum: find first X0 where cumulative energy > threshold
    //double total_energy = 0.0;
    //for (const auto& [x0, e] : energy_by_x0) total_energy += e;

    double cumulative = 0.0;
    for (const auto& [x0, e] : energy_by_x0) {
        cumulative += e;
	// if passes MIP fraction threshold
        if (cumulative > MIP_fraction * Detector.Input_Cal.Cal_Nominal_MPV_keV) return x0;
    }

    return -1.0;  // If there is an erroroneous state it will reach here and return -1
}

std::tuple<double, double, double, double> ComputeTrackerAngularStats(
    const std::vector<HitPoint>& All_hits,
    DetGeo Detector)
{
//auto [mean_tht, std_tht, mean_phi, std_phi] =
//     ComputeTrackerAngularStats(calo_hits, detector);

// std::cout << "Mean theta: " << mean_tht
//           << ", Std theta: " << std_tht << std::endl;
// std::cout << "Mean phi: " << mean_phi
//           << ", Std phi: " << std_phi << std::endl;
    std::vector<double> theta_vals;
    std::vector<double> phi_vals;

    const double zmin = Detector.Input_Tracker.Tracker_First_Pos_m-0.1;
    const double zmax = Detector.Input_Tracker.Tracker_Last_Pos_m+0.1;

    for (const auto& hit : All_hits) {
        if (hit.z >= zmin && hit.z <= zmax && hit.ene > 0.0) {
            theta_vals.push_back(hit.theta);
            phi_vals.push_back(hit.phi);
        }
    }

    auto mean_std = [](const std::vector<double>& values) -> std::pair<double, double> {
        if (values.empty()) return {0.0, 0.0};

        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double sq_sum = std::accumulate(values.begin(), values.end(), 0.0,
            [mean](double acc, double val) { return acc + (val - mean) * (val - mean); });

        double stdev = std::sqrt(sq_sum / values.size());
        return {mean, stdev};
    };

    auto [mean_tht, std_tht] = mean_std(theta_vals);
    auto [mean_phi, std_phi] = mean_std(phi_vals);

    return std::make_tuple(mean_tht, std_tht, mean_phi, std_phi);
}

std::tuple<double, double, double, double> ComputeCalAngularStats(
    const std::vector<HitPoint>& All_hits,
    DetGeo Detector)
{
//auto [mean_tht, std_tht, mean_phi, std_phi] =
//     ComputeCalAngularStats(calo_hits, detector);

// std::cout << "Mean theta: " << mean_tht
//           << ", Std theta: " << std_tht << std::endl;
// std::cout << "Mean phi: " << mean_phi
//           << ", Std phi: " << std_phi << std::endl;
    std::vector<double> theta_vals;
    std::vector<double> phi_vals;

    const double zmin = Detector.Input_Cal.Cal_First_Pos_m-0.1;
    const double zmax = zmin + Detector.Input_Cal.GetX0Pos(2.0)+0.02; // Position of 2 X0

    for (const auto& hit : All_hits) {
        if (hit.z >= zmin && hit.z <= zmax && hit.ene > 0.0) {
            theta_vals.push_back(hit.theta);
            phi_vals.push_back(hit.phi);
        }
    }

    auto mean_std = [](const std::vector<double>& values) -> std::pair<double, double> {
        if (values.empty()) return {0.0, 0.0};

        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double sq_sum = std::accumulate(values.begin(), values.end(), 0.0,
            [mean](double acc, double val) { return acc + (val - mean) * (val - mean); });

        double stdev = std::sqrt(sq_sum / values.size());
        return {mean, stdev};
    };

    auto [mean_tht, std_tht] = mean_std(theta_vals);
    auto [mean_phi, std_phi] = mean_std(phi_vals);

    return std::make_tuple(mean_tht, std_tht, mean_phi, std_phi);
}

double DeltaR_Calculation(double deltaEta, double deltaPhi, double theta_rad, bool isDeltaR)
{
  if (isDeltaR == true)
    {
      //Compute deltaR in the typical way
      return std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
    }
  else
    {
      //Compute angular separation instead
      // Avoid division by zero
      if (theta_rad < 1e-6) return 0.0;

      double sin_theta = std::sin(theta_rad);
      double cos_theta = std::cos(theta_rad);
    
      // dθ/dη = (2 / (1 + cosθ)) * (1 / sinθ)
      double dtheta_deta = (2.0 / (1.0 + cos_theta)) / sin_theta;

      double deltaTheta = dtheta_deta * deltaEta;
      double deltaOmega = std::sqrt(deltaTheta * deltaTheta + (sin_theta * deltaPhi) * (sin_theta * deltaPhi));
    
      return deltaOmega; // in radians
    }
}


std::tuple<double, double> ComputeTrackClusterDeltaR(
    const std::vector<HitPoint>& All_hits,
    double theta,
    DetGeo Detector,
    bool AngleOption)
{
//auto [mean_dR, std_dR] =
//     ComputeTrackClusterDeltaR(all_hits, detector);
// std::cout << "Mean Delta-R: " << mean_dR
//           << ", Std Delta-R: " << std_dR << std::endl;
    // Step 0: Compute angular means and standard deviations
    auto [trk_mean_tht, trk_std_tht, trk_mean_phi, trk_std_phi] =
        ComputeTrackerAngularStats(All_hits, Detector);

    auto [cal_mean_tht, cal_std_tht, cal_mean_phi, cal_std_phi] =
        ComputeCalAngularStats(All_hits, Detector);

    // Step 1: Compute Δη and Δφ between cluster and track means
    double eta_trk = -std::log(std::tan(trk_mean_tht / 2.0));
    double eta_cal = -std::log(std::tan(cal_mean_tht / 2.0));
    double deta = eta_trk - eta_cal;

    double dphi = std::fabs(cal_mean_phi - trk_mean_phi);
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    
    //double deltaR = std::sqrt(deta * deta + dphi * dphi);

    double deltaR = DeltaR_Calculation(deta,dphi,theta,AngleOption);

    // Step 2: Propagate uncertainties to std(ΔR) via quadrature
    // Assume errors on Δη and Δφ are independent
    double eta_trk_std = 0.5 * trk_std_tht / std::sin(trk_mean_tht);
    double eta_cal_std = 0.5 * cal_std_tht / std::sin(cal_mean_tht);
    double sigma_deta = std::sqrt(eta_trk_std * eta_trk_std + eta_cal_std * eta_cal_std);
    double sigma_dphi = std::sqrt(trk_std_phi * trk_std_phi + cal_std_phi * cal_std_phi);

    double dr_std = 0.0;
    if (deltaR > 1e-6) {
        dr_std = std::sqrt((deta * deta * sigma_deta * sigma_deta +
                            dphi * dphi * sigma_dphi * sigma_dphi) / (deltaR * deltaR));
    }

    return std::make_tuple(deltaR, dr_std);
}

std::tuple<double, double, double, double> ComputeCalAngularQuadrantStats(
    const std::vector<HitPoint>& All_hits,
    double theta,
    double phi,
    int quadrant,
    DetGeo Detector)
{
  //quadrant is defined in phi
  //such that 1 is top right (phi=0 to phi = pi/2)
  //and then we rotate clockwise (pi/2 to pi , pi to 3pi/2 , 3pi/2 to 2pi)
  //
  //We are computing the mean theta, phi for an angular quadrant (in phi)
  //that is outside the core of the shower (around theta)
  //such that we are sensitive to how the shower is smearing and skewing
  //instead of the core of the shower.
  //
  if (quadrant < 1 || quadrant > 4)
    {
      std::cerr << "Invalid value of quadrant entered! Must be 1 to 4." << std::endl;
      return std::make_tuple(9999,9999,9999,9999);
    }
//auto [mean_tht, std_tht, mean_phi, std_phi] =
//     ComputeCalAngularQuadrantDeltaR(calo_hits, beamTH, 1, detector);

// std::cout << "Mean theta: " << mean_tht
//           << ", Std theta: " << std_tht << std::endl;
// std::cout << "Mean phi: " << mean_phi
//           << ", Std phi: " << std_phi << std::endl;
    std::vector<double> theta_vals;
    std::vector<double> phi_vals;

    double PIVAL = 3.14159265;

    const double zmin = Detector.Input_Cal.Cal_First_Pos_m-0.1;
    const double zmax = zmin + Detector.Input_Cal.GetX0Pos(2.0)+0.02; // Position of 2 X0

    //our new origin points in x,y,z
    double z_prime = Detector.Input_Cal.Cal_First_Pos_m-0.1;
    double x_prime = z_prime * std::tan(theta) * std::cos(phi);
    double y_prime = z_prime * std::tan(theta) * std::sin(phi);

    double phimin = (quadrant-1.0)*PIVAL/2.0;
    double phimax = quadrant*PIVAL/2.0;
    //double theta_window = 0.001745; //window of 1/10 of a degree
    //double theta_window = 0.001745/10.0; //window of 1/100 of a degree
    double theta_window = 0.0;
    
    for (const auto& hit : All_hits) {
      if (hit.z >= zmin && hit.z <= zmax && hit.ene > 0.0)
	{
	  //check that hit is valid
	  //then compute the translated coordinates
	  x_prime = (hit.z - z_prime) * std::tan(theta) * std::cos(phi);
	  y_prime = (hit.z - z_prime) * std::tan(theta) * std::sin(phi);
	  double xval = x_prime - hit.z * std::tan(hit.theta) * std::cos(hit.phi);
	  double yval = y_prime - hit.z * std::tan(hit.theta) * std::sin(hit.phi);
	  double zval = z_prime - hit.z;

	  //double phival = yval/abs(yval) * std::acos(xval / std::sqrt(xval*xval + yval*yval));
	  double phival = std::atan2(yval, xval);  // returns in (-π, π]
	  if (phival < 0) phival += 2*PIVAL;       // wrap to [0, 2π)
	  double thtval = std::acos(zval / std::sqrt(xval*xval + yval*yval + zval*zval));
	  if(phival >= phimin && phival < phimax &&
	     thtval >= theta_window)
	    {//accept the hit if its in the correct quadrant and outside the theta window
	      theta_vals.push_back(hit.theta);
	      phi_vals.push_back(hit.phi);	      
	    }
	}
    }

    auto mean_std = [](const std::vector<double>& values) -> std::pair<double, double> {
        if (values.empty()) return {0.0, 0.0};

        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double sq_sum = std::accumulate(values.begin(), values.end(), 0.0,
            [mean](double acc, double val) { return acc + (val - mean) * (val - mean); });

        double stdev = std::sqrt(sq_sum / values.size());
        return {mean, stdev};
    };

    auto [mean_tht, std_tht] = mean_std(theta_vals);
    auto [mean_phi, std_phi] = mean_std(phi_vals);

    return std::make_tuple(mean_tht, std_tht, mean_phi, std_phi);
}

std::tuple<double, double, double, double, double, double> ComputeCalEnergyQuadrantStats(
    const std::vector<HitPoint>& All_hits,
    double theta,
    double phi,
    int quadrant,
    DetGeo Detector)
{
  //quadrant is defined in phi
  //such that 1 is top right (phi=0 to phi = pi/2)
  //and then we rotate clockwise (pi/2 to pi , pi to 3pi/2 , 3pi/2 to 2pi)
  //
  //We are computing the mean theta, phi for an angular quadrant (in phi)
  //that is outside the core of the shower (around theta)
  //such that we are sensitive to how the shower is smearing and skewing
  //instead of the core of the shower.
  //
  if (quadrant < 1 || quadrant > 4)
    {
      std::cerr << "Invalid value of quadrant entered! Must be 1 to 4." << std::endl;
      return std::make_tuple(9999,9999,9999,9999,9999,9999);
    }
//auto [mean_tht, std_tht, mean_phi, std_phi] =
//     ComputeCalAngularQuadrantDeltaR(calo_hits, beamTH, 1, detector);

// std::cout << "Mean theta: " << mean_tht
//           << ", Std theta: " << std_tht << std::endl;
// std::cout << "Mean phi: " << mean_phi
//           << ", Std phi: " << std_phi << std::endl;
    std::vector<double> ene_vals;
    std::vector<double> pt_vals;
    std::vector<double> pz_vals;

    double PIVAL = 3.14159265;

    const double zmin = Detector.Input_Cal.Cal_First_Pos_m-0.1;
    const double zmax = zmin + Detector.Input_Cal.GetX0Pos(2.0)+0.02; // Position of 2 X0

    //our new origin points in x,y,z
    double z_prime = Detector.Input_Cal.Cal_First_Pos_m-0.1;
    double x_prime = z_prime * std::tan(theta) * std::cos(phi);
    double y_prime = z_prime * std::tan(theta) * std::sin(phi);

    double phimin = (quadrant-1.0)*PIVAL/2.0;
    double phimax = quadrant*PIVAL/2.0;
    //double theta_window = 0.001745; //window of 1/10 of a degree
    //double theta_window = 0.001745/10.0; //window of 1/100 of a degree
    double theta_window = 0.0;
    
    for (const auto& hit : All_hits) {
      if (hit.z >= zmin && hit.z <= zmax && hit.ene > 0.0)
	{
	  //check that hit is valid
	  //then compute the translated coordinates
	  x_prime = (hit.z - z_prime) * std::tan(theta) * std::cos(phi);
	  y_prime = (hit.z - z_prime) * std::tan(theta) * std::sin(phi);
	  double xval = x_prime - hit.z * std::tan(hit.theta) * std::cos(hit.phi);
	  double yval = y_prime - hit.z * std::tan(hit.theta) * std::sin(hit.phi);
	  double zval = z_prime - hit.z;

	  //double phival = yval/abs(yval) * std::acos(xval / std::sqrt(xval*xval + yval*yval));
	  double phival = std::atan2(yval, xval);  // returns in (-π, π]
	  if (phival < 0) phival += 2*PIVAL;       // wrap to [0, 2π)
	  double thtval = std::acos(zval / std::sqrt(xval*xval + yval*yval + zval*zval));
	  if(phival >= phimin && phival < phimax &&
	     thtval >= theta_window)
	    {//accept the hit if its in the correct quadrant and outside the theta window
	      ene_vals.push_back(hit.ene);
	      pt_vals.push_back(hit.ene * std::sin(thtval));
	      pz_vals.push_back(hit.ene * std::cos(thtval));
	    }
	}
    }

    auto mean_std = [](const std::vector<double>& values) -> std::pair<double, double> {
        if (values.empty()) return {0.0, 0.0};

        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double sq_sum = std::accumulate(values.begin(), values.end(), 0.0,
            [mean](double acc, double val) { return acc + (val - mean) * (val - mean); });

        double stdev = std::sqrt(sq_sum / values.size());
        return {mean, stdev};
    };

    auto [mean_ene, std_ene] = mean_std(ene_vals);
    auto [mean_pt, std_pt] = mean_std(pt_vals);
    auto [mean_pz, std_pz] = mean_std(pz_vals);
    
    return std::make_tuple(mean_ene, std_ene, mean_pz, std_pz, mean_pt, std_pt);
}

std::tuple<double, double> ComputeQuadrantDeltaR(
    const std::vector<HitPoint>& All_hits,
    double theta,
    double phi,
    int quadrant,
    DetGeo Detector,
    bool AngleOption)
{
  if (quadrant < 1 || quadrant > 4)
    {
      std::cerr << "Invalid value of quadrant entered! Must be 1 to 4." << std::endl;
      return std::make_tuple(9999,9999);
    }
  //auto [mean_dR, std_dR] =
//     ComputeTrackClusterDeltaR(all_hits, detector);
// std::cout << "Mean Delta-R: " << mean_dR
//           << ", Std Delta-R: " << std_dR << std::endl;
    // Step 0: Compute angular means and standard deviations
    auto [trk_mean_tht, trk_std_tht, trk_mean_phi, trk_std_phi] =
        ComputeTrackerAngularStats(All_hits, Detector);

    auto [tru_cal_mean_tht, tru_cal_std_tht, tru_cal_mean_phi, tru_cal_std_phi] =
        ComputeCalAngularStats(All_hits, Detector);

    //This is under the assumption that you just measure this well in the calorimeter
    //auto [cal_mean_tht, cal_std_tht, cal_mean_phi, cal_std_phi] =
    //ComputeCalAngularQuadrantStats(All_hits, theta, phi, quadrant, Detector);

    //Under the assumption you do a unweighted average to get theta,phi in cal and tracker
    auto [cal_mean_tht, cal_std_tht, cal_mean_phi, cal_std_phi] =
      ComputeCalAngularQuadrantStats(All_hits, tru_cal_mean_tht, tru_cal_mean_phi, quadrant, Detector);

    // Step 1: Compute Δη and Δφ between cluster and track means
    double eta_trk = -std::log(std::tan(trk_mean_tht / 2.0));
    double eta_cal = -std::log(std::tan(cal_mean_tht / 2.0));
    double deta = eta_trk - eta_cal;

    double dphi = std::fabs(cal_mean_phi - trk_mean_phi);
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;

    //double deltaR = std::sqrt(deta * deta + dphi * dphi);
    //Under the assumption of perfect theta measurement
    //double deltaR = DeltaR_Calculation(deta,dphi,theta,AngleOption);
    //Under the unweighted average
    double deltaR = DeltaR_Calculation(deta,dphi,tru_cal_mean_tht,AngleOption);
    
    // Step 2: Propagate uncertainties to std(ΔR) via quadrature
    // Assume errors on Δη and Δφ are independent
    double eta_trk_std = 0.5 * trk_std_tht / std::sin(trk_mean_tht);
    double eta_cal_std = 0.5 * cal_std_tht / std::sin(cal_mean_tht);
    double sigma_deta = std::sqrt(eta_trk_std * eta_trk_std + eta_cal_std * eta_cal_std);
    double sigma_dphi = std::sqrt(trk_std_phi * trk_std_phi + cal_std_phi * cal_std_phi);

    double dr_std = 0.0;
    if (deltaR > 1e-6) {
        dr_std = std::sqrt((deta * deta * sigma_deta * sigma_deta +
                            dphi * dphi * sigma_dphi * sigma_dphi) / (deltaR * deltaR));
    }

    //If the deltaR value is 0 or uncertainty is zero, then set
    //it to an errorstate value
    if (deltaR < 1e-9 || dr_std < 1e-9)
      {
	deltaR = -1000.0;
	dr_std = -1000.0;
      }

    return std::make_tuple(deltaR, dr_std);
}

std::tuple<double, double, double> ComputeQuadrantEneStats(
    const std::vector<HitPoint>& All_hits,
    double theta,
    double phi,
    int quadrant,
    DetGeo Detector,
    bool AngleOption)
{
  if (quadrant < 1 || quadrant > 4)
    {
      std::cerr << "Invalid value of quadrant entered! Must be 1 to 4." << std::endl;
      return std::make_tuple(9999,9999,9999);
    }
  //auto [mean_dR, std_dR] =
//     ComputeTrackClusterDeltaR(all_hits, detector);
// std::cout << "Mean Delta-R: " << mean_dR
//           << ", Std Delta-R: " << std_dR << std::endl;
    // Step 0: Compute angular means and standard deviations
    auto [trk_mean_tht, trk_std_tht, trk_mean_phi, trk_std_phi] =
        ComputeTrackerAngularStats(All_hits, Detector);

    auto [tru_cal_mean_tht, tru_cal_std_tht, tru_cal_mean_phi, tru_cal_std_phi] =
        ComputeCalAngularStats(All_hits, Detector);

    //This is under the assumption that you just measure this well in the calorimeter
    //auto [cal_mean_tht, cal_std_tht, cal_mean_phi, cal_std_phi] =
    //ComputeCalAngularQuadrantStats(All_hits, theta, phi, quadrant, Detector);

    //Under the assumption you do a unweighted average to get theta,phi in cal and tracker
    auto [mean_ene, std_ene, mean_pz, std_pz, mean_pt, std_pt] =
      ComputeCalEnergyQuadrantStats(All_hits, tru_cal_mean_tht, tru_cal_mean_phi, quadrant, Detector);

    return std::make_tuple(mean_ene, mean_pz, mean_pt);
}


// Return the log-likelihood (natural log) of an angular separation
// under the Gaussian MCS model
double ComputeMCSAngularLogLikelihood(double angle_rad, double energy_GeV, double X0)
{
    if (energy_GeV <= 0.0 || X0 <= 0.0)
      std::cerr << "Energy and radiation length must be positive." << std::endl;

    // Compute sigma_theta from Molière theory
    double sigma_theta = (13.6e-3 / energy_GeV) * std::sqrt(X0) * (1.0 + 0.038 * std::log(X0));

    // Gaussian log-likelihood
    // ln P = -0.5 * ((x - μ)^2 / σ^2 + ln(2πσ^2)), μ = 0
    double logL = -0.5 * ( (angle_rad * angle_rad) / (sigma_theta * sigma_theta) +
                           std::log(2.0 * M_PI * sigma_theta * sigma_theta) );

    return logL;
}
