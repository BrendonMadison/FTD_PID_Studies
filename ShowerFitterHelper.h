#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <tuple>
#include <TMinuit.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <map>
#include <limits>
#include <algorithm>

double RadiationLengthAtZ(double z) {
    if (z < 245.0) return 0.0032;        // first 5 layers of thin silicon
    else return 1.0 / 6.0;               // rest are thick (1/6 X0 per layer)
}

struct EnergyProfilePoint {
    double X0;     // integrated radiation length
    double energy; // mean energy per hit
    double error;  // standard error
};

struct ShowerFitResult {
    double slope;
    double slope_err;
    double x_intercept;
    double x_intercept_err;
    double chi2;
    int ndf;
};

namespace {
    std::vector<double> g_X, g_Y, g_Yerr;

    void LinearChi2Function(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t flag) {
        double a = par[0];
        double b = par[1];
        double chi2 = 0.0;
        for (size_t i = 0; i < g_X.size(); ++i) {
            double yfit = a * g_X[i] + b;
            double resid = (g_Y[i] - yfit) / g_Yerr[i];
            chi2 += resid * resid;
        }
        fval = chi2;
    }
}

ShowerFitResult FitPhotonShowerStartWithTMinuit(const std::vector<EnergyProfilePoint>& profile) {
    g_X.clear();
    g_Y.clear();
    g_Yerr.clear();

    for (const auto& p : profile) {
        if (p.energy > 0 && p.error > 0) {
            g_X.push_back(p.X0);
            g_Y.push_back(std::log(p.energy));
            g_Yerr.push_back(p.error / (p.energy+1e-12));
        }
    }

    if (g_X.size() < 2) return {0, 0, 0, 0, -1, 0};

    TMinuit minuit(2);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(LinearChi2Function);

    double arglist[10];
    int ierflg = 0;

    minuit.DefineParameter(0, "slope", -1.0, 0.01, -10000.0, 10000.0);
    minuit.DefineParameter(1, "intercept", 0.0, 0.01, -10000.0, 10000.0);
    
    minuit.Migrad();

    double a, aerr, b, berr;
    minuit.GetParameter(0, a, aerr);
    minuit.GetParameter(1, b, berr);

    double x_int = -b / a;
    double x_int_err = std::sqrt((berr * berr) / (a * a) + (b * b * aerr * aerr) / (a * a * a * a));

    double chi2, edm, errdef;
    int nvpar, nparx;
    minuit.mnstat(chi2, edm, errdef, nvpar, nparx, ierflg);

    return {a, aerr, x_int, x_int_err, chi2, static_cast<int>(g_X.size()) - 2};
}

std::vector<EnergyProfilePoint> ComputeEnergyProfileFromHits(const std::vector<HitPoint>& hits) {
    std::map<double, std::vector<double>> zToEnergies;

    // Group hits by z position
    for (const auto& hit : hits) {
        zToEnergies[hit.z].push_back(hit.ene);  // z is already in meters
    }

    std::vector<EnergyProfilePoint> profile;
    double integratedX0 = 0.0;

    for (const auto& [z, energies] : zToEnergies) {
        double dX0 = RadiationLengthAtZ(z * 100.0);  // Convert z back to cm
        integratedX0 += dX0;

        if (energies.empty()) continue;

        double sum = 0.0;
        for (double e : energies) sum += e;
        double mean = sum / energies.size();

        double variance = 0.0;
        for (double e : energies) variance += (e - mean) * (e - mean);

        double stderr = 0.0;
        if (energies.size() == 1) {
            stderr = 0.1 * mean;  // Assign 10% error if only one entry
        } else {
            stderr = std::sqrt(variance / (energies.size() - 1)) / std::sqrt(energies.size());
        }

        profile.push_back({integratedX0, mean, stderr});
    }

    return profile;
}

struct ShowerStartSummary {
    std::vector<double> energies;     // size = nlayers
    std::vector<double> energy_errs;  // same size
    std::vector<double> int_rad_vals;
    std::vector<double> int_rad_vals_err;
    std::vector<double> path_length; //path length traveled in this layer in m
};

double GetIntegratedRadLen(double zpos,double theta)
  {
    // Function for calculation amount of integrated interaction length given z position
    // We also correct for the fact that angular deviation in theta results in more radiation length
    double LASTPOS = 2.212;
    double FIRSPOS = 0.22;
    double NLAY = 5.0;
    double TRACKSPACING = (LASTPOS - FIRSPOS)/(NLAY-1);
    double TRACK_X0 = 300.0e-6 / 9.37e-2 + 100.0e-6 / 17.87e-2; //for 300um Si and 100um PCB
    double CAL_X0 = 1.0/6.0;
    double CALFIRS = 2.501;
    double CALSPACING = 0.1+0.05+0.0521+0.1; //0.5 mm between layers
    if (zpos < FIRSPOS) return 0.0;
    if (zpos > FIRSPOS && zpos < LASTPOS+0.1)
      {
	int n_track_layers = int(std::floor((zpos-FIRSPOS)/TRACKSPACING) + 1);
	return n_track_layers * 1.0 * TRACK_X0 / std::cos(theta);
      }
    if (zpos > LASTPOS+0.1)
      {
	int n_cal_layers = int(std::floor((zpos-CALFIRS)/CALSPACING)+1);
	return (NLAY * TRACK_X0 + n_cal_layers * CAL_X0) / std::cos(theta);
      }
  }

double GetConvertionRadLen(double zpos,double theta)
  {
    // Function for calculating the radiation length for the layer convertion happens in
    double LASTPOS = 2.212;
    double FIRSPOS = 0.22;
    double NLAY = 5.0;
    double TRACKSPACING = (LASTPOS - FIRSPOS)/(NLAY-1);
    double TRACK_X0 = 300.0e-6 / 9.37e-2 + 100.0e-6 / 17.87e-2; //for 300um Si and 100um PCB
    double CAL_X0 = 1.0/6.0;
    double CALFIRS = 2.501;
    double CALSPACING = 0.1+0.05+0.0521+0.1; //0.5 mm between layers
    if (zpos < FIRSPOS) return 0.0;
    if (zpos > FIRSPOS && zpos < LASTPOS+0.1)
      {
	return TRACK_X0/std::cos(theta);
      }
    if (zpos > LASTPOS)
      {
	return CAL_X0/std::cos(theta);
      }
  }

double GetLayerPathLength(double zpos,double theta)
  {
    // Function for calculating the radiation length for the layer convertion happens in
    double LASTPOS = 2.212;
    double FIRSPOS = 0.22;
    double NLAY = 5.0;
    double TRACKSPACING = (LASTPOS - FIRSPOS)/(NLAY-1);
    double TRACK_X0 = 300.0e-6 / 9.37e-2 + 100.0e-6 / 17.87e-2; //for 300um Si and 100um PCB
    double CAL_X0 = 1.0/6.0;
    double CALFIRS = 2.501;
    double CALSPACING = 0.1+0.05+0.0521+0.1; //0.5 mm between layers
    double CAL_THICK = 1000.0e-6;
    double TRA_THICK = 300.0e-6;
    if (zpos < FIRSPOS) return 0.0;
    if (zpos > FIRSPOS && zpos < LASTPOS+0.1)
      {
	return TRA_THICK/std::cos(theta);
      }
    if (zpos > LASTPOS)
      {
	return CAL_THICK/std::cos(theta);
      }
  }

ShowerStartSummary FindShowerStartFromHits(const std::vector<HitPoint>& hits, int nlayers, double theta) {
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

	std::cout << e_sum << std::endl;
        result.energies.push_back(e_sum);
        result.energy_errs.push_back(e_err);

        double zval = it->first/100.0;
	// Store total integrated X0 and error at shower start
        result.int_rad_vals.push_back(GetIntegratedRadLen(zval, theta));
        result.int_rad_vals_err.push_back(GetConvertionRadLen(zval, theta) / std::sqrt(12.0));
	result.path_length.push_back(GetLayerPathLength(zval,theta));
    }

    return result;
}


double ComputePhotonConversionLogLikelihood(double x, double sigma_x) {
  //where x is the integrated radiation length at the convertion point
    const double lambda = 7.0 / 9.0;  // decay constant of exponential
    const double norm = lambda / 2.0;

    double exp_term = std::exp((lambda * lambda * sigma_x * sigma_x) / 2.0 - lambda * x);
    double arg = (lambda * sigma_x * sigma_x - x) / (std::sqrt(2.0) * sigma_x);
    double erfc_term = TMath::Erfc(arg);

    double likelihood = norm * exp_term * erfc_term;

    // To avoid log(0) or negative values due to numerical underflow:
    if (likelihood <= 1e-300) return -700.0;  // log(1e-300)
    return std::log(likelihood);
}

constexpr double s_silicon_thickness_um = 300.0;
constexpr double s_nominal_MPV_keV = 75.0;  // for 300 µm normal incidence
constexpr double s_nominal_sigma_keV = 15.0; // estimated FWHM/2.35

double ComputeLayered2MIPEnergyLogLikelihood(const ShowerStartSummary& summary, const TH1D* pdf_histogram, double beam_energy_GeV) {
  double log_likelihood = 0.0;
  double ene_sum = 0.0;
  int mip_factor = 2;
  
  // Error catch for if we don't have hits or a histogram
  if (!pdf_histogram || summary.energies.size() < 1)
    return std::log(1e-4);  // fallback for missing hit

  int nlayers = summary.energies.size();
  // Compute the log-likelihoods for the individual hits
  for (int i = 0; i < nlayers; ++i) {
    ene_sum += summary.energies[i];

    double expected = s_nominal_MPV_keV * mip_factor * summary.path_length[i]/(s_silicon_thickness_um*1e-6);
    double sigma = s_nominal_sigma_keV * std::sqrt(summary.path_length[i]/(s_silicon_thickness_um*1e-6));

    double pdf_val = TMath::Landau(summary.energies[i], expected, sigma, true);  // normalized = true
    log_likelihood += std::log(pdf_val + 1e-12);  // prevent log(0)
  }
	
  // Compute the log-likelihood for the summed energy
  double pdf = pdf_histogram->Interpolate(ene_sum);
  log_likelihood += std::log(pdf + 1e-12);

  return log_likelihood;
}

TH1D* GenerateShowerLandauSumPDF(
    const ShowerStartSummary& summary,  // one per layer
    int n_samples = 1e6,
    double mpv_300 = 75.0,     // keV
    double sigma_300 = 15.0,   // keV
    int n_bins = 1000,
    double hist_min = 0.0,
    double hist_max = 1000.0)
{
    TRandom3 rng(42);
    TH1D* hist = new TH1D("shower_sum_pdf", "Sum of Landau deposits", n_bins, hist_min, hist_max);

    int n_layers = summary.path_length.size();
    for (int i = 0; i < n_samples; ++i) {
        double total = 0.0;
        for (int j = 0; j < n_layers; ++j) {
            double scale = summary.path_length[j] / 300.0e-6;
            double mpv = mpv_300 * scale;
            double sigma = sigma_300 * std::sqrt(scale);  // approximate width scaling
            total += rng.Landau(mpv, sigma);
        }
        hist->Fill(total);
    }

    hist->Scale(1.0 / hist->Integral("width"));  // Normalize to a PDF
    return hist;
}
