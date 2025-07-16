#pragma once
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <map>
#include <numeric>
#include <unordered_map>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
//constexpr double mrad_tol = 0.005;  // 5 mrad
//constexpr double silicon_thickness_um = 300.0;
//constexpr double nominal_MPV_keV = 75.0;  // for 300 µm normal incidence
//constexpr double nominal_sigma_keV = 15.0; // estimated FWHM/2.35

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

double compute_pt(double p, double theta) {
    return p * std::sin(theta);
}

double compute_radius(double pt, double B = 3.5) {
    return pt / (0.3 * B);
}

struct TrackerGeo {
//Struct for defining tracker geometry
	double Tracker_Rad_Len_cm;
	double Tracker_Thickness_um;
	double Mag_Field_tesla;
	double Tracker_Layer_Spacing_m;
	double Tracker_First_Pos_m;
	double Tracker_Last_Pos_m;
	double Tracker_X0_Per_Layer;
	double Tracker_Nominal_MPV_keV;
	double Tracker_Nominal_Wid_keV;
	double Tracker_Nominal_Ene;//Sets the nominal energy that the MPV and Width are at
	double Tracker_MS_Per_Layer;
        double Tracker_Time_Res_ns;
	int Tracker_Num_Layers;
	
	TrackerGeo(double trl, double ttu, double mft, double tfpm, double tlpm, double txp, int tnl) : Tracker_Rad_Len_cm(trl), Tracker_Thickness_um(ttu), Mag_Field_tesla(mft), Tracker_First_Pos_m(tfpm), Tracker_Last_Pos_m(tlpm), Tracker_X0_Per_Layer(txp), Tracker_Num_Layers(tnl) {
		if ( mft <= 0.0 ) std::cerr << "Warning, tracker initialized with negative or zero magnetic field!" << std::endl;
		Tracker_Layer_Spacing_m = (Tracker_Last_Pos_m - Tracker_First_Pos_m)/(Tracker_Num_Layers-1);
	} //7 variable Constructor

        void SetTimeRes(double tres)
        {
	  Tracker_Time_Res_ns = tres;
        }
  
	void SetMIPValues(double mpv, double wid, double ene)
	{
		//For setting the 1 MIP mpv and wid and nominal energy
		Tracker_Nominal_MPV_keV = mpv;
		Tracker_Nominal_Wid_keV = wid;
		Tracker_Nominal_Ene = ene;
	}
	
	void SetMSValues(double mspl)
	{
		Tracker_MS_Per_Layer = mspl;
	}
	
	double GetLayersPassed(double zpos, double theta)
	{
		//Returns the effective number of layers passed from
		//the z position and the theta angle
		if (zpos < Tracker_First_Pos_m) return 0.0;
		if (zpos > Tracker_First_Pos_m && zpos < Tracker_Last_Pos_m+0.1)
		  {
			int n_track_layers = int(std::floor((zpos - Tracker_First_Pos_m)/Tracker_Layer_Spacing_m) + 1);
			return n_track_layers / std::cos(theta);
		  }
		//If it is in the calorimeter
		return Tracker_Num_Layers / std::cos(theta);
	}
};

struct CalGeo {
//Struct for defining calorimeter geometry
	double Active_Rad_Len_cm;
	double Passive_Rad_Len_cm;
	double Active_Thickness_um;
	double Passive_Thickness_um;
	double Mag_Field_tesla;
	double Cal_AirGap_m;
	double Cal_Layer_Spacing_m;
	double Cal_First_Pos_m;
	double Cal_Tot_Thickness_m;//including spacing
	double Cal_X0_Per_Layer;
	double Cal_Nominal_MPV_keV;
	double Cal_Nominal_Wid_keV;
	double Cal_Nominal_Ene;//Sets the nominal energy that the MPV and Width are at
	double Cal_MS_Per_Layer;
	int Cal_Num_Layers;
	
	CalGeo(double arl, double prl, double atu, double ptu, double mft, double cfp, double ctt, double cxp, int cnl) : Active_Rad_Len_cm(arl), Passive_Rad_Len_cm(prl), Active_Thickness_um(atu), Passive_Thickness_um(ptu), Mag_Field_tesla(mft), Cal_First_Pos_m(cfp), Cal_Tot_Thickness_m(ctt), Cal_X0_Per_Layer(cxp), Cal_Num_Layers(cnl) {
		if ( mft <= 0.0 ) std::cerr << "Warning, tracker initialized with negative or zero magnetic field!" << std::endl;
		//We assume a PCB of 1mm thickness
		Cal_AirGap_m = Cal_Tot_Thickness_m - 1.0e-3 - Passive_Thickness_um*1e-6 - Active_Thickness_um*1e-6;
	    Cal_Layer_Spacing_m = Cal_Tot_Thickness_m;
	} //9 variable Constructor
	
	void SetMIPValues(double mpv, double wid, double ene)
	{
		//For setting the 1 MIP mpv and wid and nominal energy
		Cal_Nominal_MPV_keV = mpv;
		Cal_Nominal_Wid_keV = wid;
		Cal_Nominal_Ene = ene;
	}
	
	void SetMSValues(double mspl)
	{
		Cal_MS_Per_Layer = mspl;
	}
	
};

struct DetGeo {
    //Struct for defining entire detector geometry	
	TrackerGeo Input_Tracker;
	CalGeo Input_Cal;
	double Sigma_Ene;
	double Sigma_Tht;
	double Sigma_Phi;
	double Sigma_Z;
	double XY_size;
	double Sigma_XY;
	
	DetGeo( TrackerGeo it, CalGeo ic ) : Input_Tracker(it), Input_Cal(ic) {} //double struct constructor

	void SetResolutions(double se, double st, double sp, double sz, double xy, double sxy)
	{
		//Transcribe the values over
		Sigma_Ene = se;
		Sigma_Tht = st;
		Sigma_Phi = sp;
		Sigma_Z = sz;
		XY_size = xy;
		Sigma_XY = sxy;
	}

	// Function for computing integrated radiation length
	// given input z position and theta (polar angle)
	double GetIntegratedRadLen(double zpos, double theta)
	{
		//zpos must be in meters and theta must be in radians
		if (zpos < Input_Tracker.Tracker_First_Pos_m) return 0.0;
		if (zpos > Input_Tracker.Tracker_First_Pos_m && zpos < Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			int n_track_layers = int(std::floor((zpos - Input_Tracker.Tracker_First_Pos_m)/Input_Tracker.Tracker_Layer_Spacing_m) + 1);
			return n_track_layers * Input_Tracker.Tracker_X0_Per_Layer / std::cos(theta);
		  }
		if (zpos > Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			int n_cal_layers = int(std::floor((zpos - Input_Cal.Cal_First_Pos_m)/Input_Cal.Cal_Tot_Thickness_m)+1);
			return Input_Tracker.Tracker_Num_Layers * Input_Tracker.Tracker_X0_Per_Layer / std::cos(theta) + n_cal_layers * Input_Cal.Cal_X0_Per_Layer / std::cos(theta);
		  }
		  //Error case, return 0.0
		return 0.0;
	}
	
	// Function for getting current radiation length
	// given input z position and theta (polar angle)
	double GetConvertionRadLen(double zpos, double theta)
	{
		//zpos must be in meters and theta must be in radians
		if (zpos < Input_Tracker.Tracker_First_Pos_m) return 0.0;
		if (zpos > Input_Tracker.Tracker_First_Pos_m && zpos < Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			  //in the tracker
			return Input_Tracker.Tracker_X0_Per_Layer / std::cos(theta);
		  }
		if (zpos > Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			  //in the calorimeter
			return Input_Cal.Cal_X0_Per_Layer / std::cos(theta);
		  }
		  //Error case, return 0.0
		return 0.0;
	}
	
    // Function for getting the path length through
	// the current active layer
	// given input z position and theta (polar angle)
	double GetLayerPathLength(double zpos, double theta)
	{
		//zpos must be in meters and theta must be in radians
		if (zpos < Input_Tracker.Tracker_First_Pos_m) return 0.0;
		if (zpos > Input_Tracker.Tracker_First_Pos_m && zpos < Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			  //in the tracker
			return Input_Tracker.Tracker_Thickness_um / std::cos(theta);
		  }
		if (zpos > Input_Tracker.Tracker_Last_Pos_m+0.1)
		  {
			  //in the calorimeter
			return Input_Cal.Active_Thickness_um / std::cos(theta);
		  }
		  //Error case, return 0.0
		return 0.0;
	}
};

struct HitPoint {
    double theta;
    double phi;
    double z;
    double time;
    double ene = 0.0;

  HitPoint(double t, double p, double zz) : theta(t), phi(p), z(zz) {time = 0.0; ene = 0.0;} //3 var constructor
  HitPoint(double t, double p, double zz, double tt) : theta(t), phi(p), z(zz), time(tt) {ene = 0.0;} //4 var constructor
  HitPoint(double t, double p, double zz, double tt, double ee) : theta(t), phi(p), z(zz), time(tt), ene(ee) {} //5 var constructor
  
  double x() const { return std::sin(theta) * std::cos(phi); }
  double y() const { return std::sin(theta) * std::sin(phi); }

  double GetTheta() const { return theta; }
  double GetPhi() const { return phi; }
  double GetZ() const { return z; }
  double GetTime() const { return time; }
  double GetTZRat() const { return time/z; }
  double GetEne() const { return ene; }

  static std::vector<double> extractTheta(const std::vector<HitPoint>& hits) {
        std::vector<double> out(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].theta;
        return out;
    }

    static std::vector<double> extractPhi(const std::vector<HitPoint>& hits) {
        std::vector<double> out(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].phi;
        return out;
    }

    static std::vector<double> extractZ(const std::vector<HitPoint>& hits) {
        std::vector<double> out(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].z;
        return out;
    }

  static std::vector<double> extractEne(const std::vector<HitPoint>& hits) {
    std::vector<double> out(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) out[i] = hits[i].ene;
    return out;
    }
};

struct HitPointErrors {
    double theta, phi, z;
    double sigma_theta, sigma_phi, sigma_z;

    HitPointErrors(double t, double p, double zz, double st, double sp, double sz)
        : theta(t), phi(p), z(zz), sigma_theta(st), sigma_phi(sp), sigma_z(sz) {}

    double x() const { return std::sin(theta) * std::cos(phi); }
    double y() const { return std::sin(theta) * std::sin(phi); }
};

struct ShowerStartSummary {
    std::vector<double> energies;     // size = nlayers
    std::vector<double> energy_errs;  // same size
    std::vector<double> int_rad_vals;
    std::vector<double> int_rad_vals_err;
    std::vector<double> path_length; //path length traveled in this layer in m
};

ShowerStartSummary FindShowerStartFromHits(const std::vector<HitPoint>& hits, DetGeo Detector, int nlayers, double theta) {
    std::map<double, std::vector<double>> zLayerMap;

    for (const auto& h : hits) {
      if (h.ene > 0.0)
	    {
	      double z_cm = h.z * 100.0;
	      zLayerMap[z_cm].push_back(h.ene);
	    }
    }

    ShowerStartSummary result{};

    result.energies.reserve(nlayers);
    result.energy_errs.reserve(nlayers);
    result.int_rad_vals.reserve(nlayers);
    result.int_rad_vals_err.reserve(nlayers);
    result.path_length.reserve(nlayers);

    if (zLayerMap.empty() || nlayers <= 0)
        return result;

    auto it = zLayerMap.begin();  // start from first layer

    int count = 0;
    for (; it != zLayerMap.end() && count < nlayers; ++it, ++count) {
        const auto& enes = it->second;
        double e_sum = std::accumulate(enes.begin(), enes.end(), 0.0);
        double e_err = 0.008 * e_sum;

        result.energies.push_back(e_sum);
        result.energy_errs.push_back(e_err);

        double zval = it->first/100.0;
	// Store total integrated X0 and error at shower start
        result.int_rad_vals.push_back(Detector.GetIntegratedRadLen(zval, theta));
        result.int_rad_vals_err.push_back(Detector.GetConvertionRadLen(zval, theta) / std::sqrt(12.0));
	result.path_length.push_back(Detector.GetLayerPathLength(zval,theta));
    }

    return result;
}

struct DataLoader{
	std::vector<HitPoint> Track_Hits;
	std::vector<HitPoint> All_Hits;
	ShowerStartSummary Shower_Summary;
	
	DataLoader(std::vector<double>* track_tht,
             std::vector<double>* track_phi,
             std::vector<double>* track_z,
	     std::vector<double>* track_time,
	     std::vector<double>* track_ene,
             double beamTH,
             double beamPH,
             DetGeo Detector,
             int n_shower_vals = 3)
  {

	  double Cal_2X0 = Detector.Input_Cal.Cal_Tot_Thickness_m * 2.0/Detector.Input_Cal.Cal_X0_Per_Layer + Detector.Input_Cal.Cal_First_Pos_m;
	  double z_tolerance = Detector.Input_Tracker.Tracker_Layer_Spacing_m*0.1; //10% tracker layer spacing
	  
    // Step 1: Convert tracker hits to HitPoint list
    for (size_t j = 0; j < track_tht->size(); ++j) {
      double zval = (*track_z)[j]/100.0;
      double thtval = (*track_tht)[j];
      double phival = (*track_phi)[j];
      double eneval = (*track_ene)[j];
      double tval = (*track_time)[j];
      
      if (zval <= Detector.Input_Tracker.Tracker_First_Pos_m-0.1 || zval >= Detector.Input_Tracker.Tracker_Last_Pos_m+0.1 || std::abs(thtval - beamTH) > 0.001) continue;

      Track_Hits.emplace_back(thtval, phival, zval, tval, eneval);  // meters
    }

    // Step 2: Tracker layer completeness compensation
    std::vector<bool> has_hit(Detector.Input_Tracker.Tracker_Num_Layers, false);
    for (const auto& hit : Track_Hits) {
      for (int k = 0; k < Detector.Input_Tracker.Tracker_Num_Layers; ++k) {
        double expected_z_cm = Detector.Input_Tracker.Tracker_First_Pos_m + k * Detector.Input_Tracker.Tracker_Layer_Spacing_m;
	//std::cout << " COM " << hit.z << " , " << expected_z_cm << std::endl;
        if (std::abs(hit.z - expected_z_cm) < z_tolerance) {
          has_hit[k] = true;
          break;
        }
      }
    }

    for (int k = 0; k < Detector.Input_Tracker.Tracker_Num_Layers; ++k) {
      if (!has_hit[k]) {
        double missing_z = (Detector.Input_Tracker.Tracker_First_Pos_m + k * Detector.Input_Tracker.Tracker_Layer_Spacing_m);
        Track_Hits.emplace_back(beamTH, beamPH, missing_z, 9999.9, 0.0); //Give it large time and zero energy
	//	std::cout << Detector.Input_Tracker.Tracker_First_Pos_m << " , " << Detector.Input_Tracker.Tracker_Layer_Spacing_m << " Insert empty hit at " << missing_z << std::endl;
      }
    }

    // Troubleshooting -- print out the values of Track_Hits
    //for (size_t j = 0; j < Track_Hits.size(); j++)
    //  {
    //	std::cout << Track_Hits[j].theta << " , " << Track_Hits[j].z << " , " << Track_Hits[j].GetEne() << " , " << Track_Hits[j].GetTime() << std::endl;
    //  }
    
    // Step 3: Fill All_Hits with tracker hits (now complete) + calorimeter front layers
    for (const auto& hit : Track_Hits) {
      All_Hits.emplace_back(hit.GetTheta(), hit.GetPhi(), hit.GetZ(), hit.GetEne());
    }

    for (size_t j = 0; j < track_tht->size(); ++j) {
      double zval = (*track_z)[j]/100.0;
      double thtval = (*track_tht)[j];
      double phival = (*track_phi)[j];
      double eneval = (*track_ene)[j];
      double tval = (*track_time)[j];

      if (zval <= Detector.Input_Tracker.Tracker_Last_Pos_m+0.1 || zval > Cal_2X0 ) continue; // only first 2 X0 of calorimeter

      if (std::abs(thtval - beamTH) > 0.001) {
        All_Hits.emplace_back(thtval, phival, zval, 9999.9, 0.0);  // angular mismatch → zero energy
      } else {
        All_Hits.emplace_back(thtval, phival, zval, tval, eneval);
      }
    }

    // Step 4: Compute shower start summary
    Shower_Summary = FindShowerStartFromHits(All_Hits, Detector, n_shower_vals, beamTH);
  }
	
};

struct LikelihoodModel {
    std::vector<std::string> labels;         // Variable names
    std::vector<double> mean;                // Mean vector
    std::vector<std::vector<double>> cov;    // Covariance matrix
    std::vector<std::vector<double>> inv_cov; // Inverse covariance matrix
    double det = 0.0;                        // Determinant of covariance
};

LikelihoodModel LoadLikelihoodModelSimple(const std::string& base_name) {
    LikelihoodModel model;
    std::ifstream f_labels(base_name + "_labels.csv");
    std::ifstream f_means(base_name + "_means.csv");
    std::ifstream f_vars(base_name + "_var.csv");
    std::ifstream f_covar(base_name + "_covar.csv");

    std::string line;
    while (std::getline(f_labels, line)) model.labels.push_back(line);

    double val;
    while (f_means >> val) model.mean.push_back(val);
    while (f_vars >> val) model.det += val;  // also store total variance if needed

    const int N = model.labels.size();
    model.cov.resize(N, std::vector<double>(N));
    model.inv_cov.resize(N, std::vector<double>(N));

    for (int i = 0; i < N; ++i) {
        std::getline(f_covar, line);
        std::stringstream ss(line);
        std::string item;
        int j = 0;
        while (std::getline(ss, item, ',') && j < N) {
            model.cov[i][j++] = std::stod(item);
        }
    }

    // Dummy inverse and det for now
    model.det = 1.0;
    for (int i = 0; i < N; ++i) model.inv_cov[i][i] = 1.0 / model.cov[i][i];

    return model;
}
 //to be used with the loglikelihood csv writer program!

double ComputeMultivariateGaussianProb(
    const std::vector<double>& x,
    const LikelihoodModel& model)
{
    const int N = x.size();
    if (model.mean.size() != N || model.inv_cov.size() != N) {
        std::cerr << "Dimension mismatch in input vector or model" << std::endl;
        return 0.0;
    }

    std::vector<double> diff(N, 0.0);
    for (int i = 0; i < N; ++i)
        diff[i] = x[i] - model.mean[i];

    // Compute Mahalanobis distance: d² = (x - μ)^T Σ⁻¹ (x - μ)
    double mahal = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            mahal += diff[i] * model.inv_cov[i][j] * diff[j];
        }
    }

    double norm = std::pow(2 * M_PI, -0.5 * N) / std::sqrt(model.det + 1e-12);
    double prob = norm * std::exp(-0.5 * mahal);
    return prob;
} //assumes ordering from the csv file
// this ordering is generally ene, helix, cnt, tof

double ComputeIndependentGaussianProb(
    const std::vector<double>& x,
    const LikelihoodModel& model)
{
    const int N = x.size();
    if (model.mean.size() != N || model.cov.size() != N) {
        std::cerr << "Size mismatch in model" << std::endl;
        return 0.0;
    }

    double logP = 0.0;

    for (int i = 0; i < N; ++i) {
        double var = (i < model.cov[i].size()) ? model.cov[i][i] : 0.0;
        if (var < 1e-12) var = 1e-12;  // avoid division by zero

        double diff = x[i] - model.mean[i];
        logP += -0.5 * (diff * diff / var + std::log(2.0 * M_PI * var));
    }

    return std::exp(logP);
}
