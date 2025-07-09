// HelixExtrapolator.h
#pragma once
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <map>
#include <numeric>
#include <unordered_map>
#include <random>

constexpr double mrad_tol = 0.005;  // 5 mrad

struct HitPoint {
    double theta;
    double phi;
    double z;
    double ene = 0.0;

  HitPoint(double t, double p, double zz) : theta(t), phi(p), z(zz) {} //3 var constructor
  HitPoint(double t, double p, double zz, double ee) : theta(t), phi(p), z(zz), ene(ee) {} //4 var constructor
  
  double x() const { return std::sin(theta) * std::cos(phi); }
  double y() const { return std::sin(theta) * std::sin(phi); }

  double GetTheta() const { return theta; }
  double GetPhi() const { return phi; }
  double GetZ() const { return z; }
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


// Generate all combinations of k elements from input vector
void generateCombinations(const std::vector<HitPoint>& points, size_t k,
                          std::vector<std::vector<HitPoint>>& combinations) {
    std::vector<bool> bitmask(k, true);
    bitmask.resize(points.size(), false);
    do {
        std::vector<HitPoint> comb;
        for (size_t i = 0; i < points.size(); ++i)
            if (bitmask[i]) comb.push_back(points[i]);
        combinations.push_back(comb);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

// ConstantZExtrapolation implementation
void ConstantZExtrapolation(std::vector<HitPoint>& hits) {
    std::map<double, std::vector<HitPoint>> zGroups;
    for (const auto& h : hits)
        zGroups[h.z].push_back(h);

    std::vector<HitPoint> newHits;
    for (const auto& [z, group] : zGroups) {
        if (group.size() < 2) continue;
        for (size_t k = 2; k <= group.size(); ++k) {
            std::vector<std::vector<HitPoint>> combos;
            generateCombinations(group, k, combos);
            for (const auto& combo : combos) {
                double tAvg = 0, pAvg = 0;
                for (const auto& pt : combo) {
                    tAvg += pt.theta;
                    pAvg += pt.phi;
                }
                tAvg /= combo.size();
                pAvg /= combo.size();
                newHits.emplace_back(tAvg, pAvg, z);
            }
        }
    }
    hits.insert(hits.end(), newHits.begin(), newHits.end());
}

// Return distance in (theta, phi) space
inline double angularDist(const HitPoint& a, const HitPoint& b) {
    return std::sqrt((a.theta - b.theta) * (a.theta - b.theta) +
                     (a.phi - b.phi) * (a.phi - b.phi));
}

// Approximate helical interpolation from sorted hits
std::vector<HitPoint> interpolateHelix(const std::vector<HitPoint>& sortedHits) {
    std::vector<HitPoint> result;
    for (size_t i = 0; i + 1 < sortedHits.size(); ++i) {
        const auto& p1 = sortedHits[i];
        const auto& p2 = sortedHits[i + 1];

        double midZ = 0.5 * (p1.z + p2.z);
        double dz = p2.z - p1.z;
        double dTheta = p2.theta - p1.theta;
        double dPhi = p2.phi - p1.phi;

        result.emplace_back(p1.theta + dTheta * 0.5, p1.phi + dPhi * 0.5, midZ);
    }
    return result;
}

// DifferentZExtrapolation implementation
void DifferentZExtrapolation(std::vector<HitPoint>& hits) {
    std::vector<HitPoint> newHits;
    for (size_t i = 0; i < hits.size(); ++i) {
        std::vector<HitPoint> candidates;
        for (size_t j = 0; j < hits.size(); ++j) {
            if (i == j) continue;
            if (std::abs(hits[i].z - hits[j].z) < 1e-3) continue;
            if (angularDist(hits[i], hits[j]) < mrad_tol)
                candidates.push_back(hits[j]);
        }
        // Include the central hit
        candidates.push_back(hits[i]);
        std::sort(candidates.begin(), candidates.end(), [](const HitPoint& a, const HitPoint& b) {
            return a.z < b.z;
        });

        // Remove duplicates in z
        std::map<double, HitPoint> uniqueZ;
        for (const auto& pt : candidates)
	  {
	  auto it = uniqueZ.find(pt.z);
	  if (it == uniqueZ.end()) uniqueZ.emplace(pt.z, pt);
	  }
	//uniqueZ[pt.z] = pt;

        std::vector<HitPoint> uniqueHits;
        for (const auto& [_, pt] : uniqueZ)
            uniqueHits.push_back(pt);

        if (uniqueHits.size() < 5) continue;

        // Take 5-z-point subsequences and interpolate
        for (size_t start = 0; start + 5 <= uniqueHits.size(); ++start) {
            std::vector<HitPoint> slice(uniqueHits.begin() + start, uniqueHits.begin() + start + 5);
            auto midpoints = interpolateHelix(slice);
            newHits.insert(newHits.end(), midpoints.begin(), midpoints.end());
        }
    }
    hits.insert(hits.end(), newHits.begin(), newHits.end());
}

// Extractors for final use
void extractVectors(const std::vector<HitPoint>& hits,
                    std::vector<double>& theta,
                    std::vector<double>& phi,
                    std::vector<double>& z) {
    theta = HitPoint::extractTheta(hits);
    phi = HitPoint::extractPhi(hits);
    z = HitPoint::extractZ(hits);
}

void MultipleScatteringRotationCorrection(std::vector<HitPoint>& hits) {
    if (hits.size() < 6) return; // Not enough layers to iterate

    // Group hits by unique z layer values
    std::map<double, std::vector<HitPoint>> layerMap;
    for (const auto& pt : hits)
        layerMap[pt.z].push_back(pt);

    std::vector<double> sortedZ;
    for (const auto& entry : layerMap)
        sortedZ.push_back(entry.first);

    std::mt19937 rng(42); // deterministic seed for reproducibility
    std::uniform_real_distribution<double> dist(-0.015, 0.015); // ±150 microns

    for (size_t i = 0; i + 3 < sortedZ.size(); ++i) {
        double z0 = sortedZ[i];
        double z1 = sortedZ[i + 1];
        double z2 = sortedZ[i + 2];
        double z3 = sortedZ[i + 3];

        // Prepare linear fit for theta(z) and phi(z)
        std::vector<double> z_vals, theta_vals, phi_vals;
        for (double z : {z0, z1, z2}) {
            for (const auto& pt : layerMap[z]) {
                z_vals.push_back(z);
                theta_vals.push_back(pt.theta);
                phi_vals.push_back(pt.phi);
            }
        }
        if (z_vals.size() < 3) continue;

        // Least-squares linear fit: y = a*z + b
        auto linear_fit = [](const std::vector<double>& z, const std::vector<double>& y) {
            double S = 0, Sz = 0, Sy = 0, Szz = 0, Szy = 0;
            size_t n = z.size();
            for (size_t i = 0; i < n; ++i) {
                S += 1;
                Sz += z[i];
                Sy += y[i];
                Szz += z[i] * z[i];
                Szy += z[i] * y[i];
            }
            double denom = S * Szz - Sz * Sz;
            if (std::abs(denom) < 1e-10) return std::pair<double, double>{0, 0};
            double a = (S * Szy - Sz * Sy) / denom;
            double b = (Szz * Sy - Sz * Szy) / denom;
            return std::make_pair(a, b);
        };

        auto [a_tht, b_tht] = linear_fit(z_vals, theta_vals);
        auto [a_phi, b_phi] = linear_fit(z_vals, phi_vals);

        double delta_z = dist(rng);
        double zrand = z3 + delta_z;
        double tht_target = a_tht * zrand + b_tht;
        double phi_target = a_phi * zrand + b_phi;

        double tht_avg = 0, phi_avg = 0;
        int count = 0;
        for (const auto& pt : layerMap[z3]) {
            tht_avg += pt.theta;
            phi_avg += pt.phi;
            count++;
        }
        if (count == 0) continue;
        tht_avg /= count;
        phi_avg /= count;

        double dtheta = tht_target - tht_avg;
        double dphi = phi_target - phi_avg;

        for (size_t j = i + 3; j < sortedZ.size(); ++j) {
            for (auto& pt : layerMap[sortedZ[j]]) {
                pt.theta += dtheta;
                pt.phi += dphi;
            }
        }
    }

    // Flatten result
    hits.clear();
    for (const auto& [z, pts] : layerMap)
        hits.insert(hits.end(), pts.begin(), pts.end());
}

// Assumes HitPoint has getZ(), theta, phi, and copyable behavior
void MultipleScatteringRotationCorrectionBackward(std::vector<HitPoint>& hits) {
    if (hits.empty()) return;

    std::unordered_map<double, std::vector<HitPoint>> layerMap;
    for (const auto& hit : hits)
        layerMap[hit.GetZ()].push_back(hit);

    std::vector<double> sortedZ;
    for (const auto& entry : layerMap)
        sortedZ.push_back(entry.first);
    std::sort(sortedZ.begin(), sortedZ.end());

    if (sortedZ.size() < 4) return;

    std::default_random_engine gen;
    std::uniform_real_distribution<double> smear(-0.015, 0.015); // ±150 µm

    for (int i = sortedZ.size() - 1; i >= 3; --i) {
        double z1 = sortedZ[i];
        double z2 = sortedZ[i - 1];
        double z3 = sortedZ[i - 2];

        const auto& hits1 = layerMap[z1];
        const auto& hits2 = layerMap[z2];
        const auto& hits3 = layerMap[z3];

        std::vector<double> thtVals, phiVals, zVals;
        for (const auto& h : hits1) { thtVals.push_back(h.theta); phiVals.push_back(h.phi); zVals.push_back(h.GetZ()); }
        for (const auto& h : hits2) { thtVals.push_back(h.theta); phiVals.push_back(h.phi); zVals.push_back(h.GetZ()); }
        for (const auto& h : hits3) { thtVals.push_back(h.theta); phiVals.push_back(h.phi); zVals.push_back(h.GetZ()); }

        auto linearFit = [](const std::vector<double>& x, const std::vector<double>& y) {
            double sx = 0, sy = 0, sxy = 0, sx2 = 0;
            for (size_t i = 0; i < x.size(); ++i) {
                sx += x[i];
                sy += y[i];
                sxy += x[i] * y[i];
                sx2 += x[i] * x[i];
            }
            double n = x.size();
            double denom = n * sx2 - sx * sx;
            if (std::abs(denom) < 1e-9) return std::make_pair(0.0, 0.0);
            double m = (n * sxy - sx * sy) / denom;
            double b = (sy - m * sx) / n;
            return std::make_pair(m, b);
        };

        auto [m_tht, b_tht] = linearFit(zVals, thtVals);
        auto [m_phi, b_phi] = linearFit(zVals, phiVals);

        for (int j = i - 3; j >= 0; --j) {
            double z = sortedZ[j];
            auto& hitsToModify = layerMap[z];
            for (auto& hit : hitsToModify) {
                double deltaZ = smear(gen);
                double targetTheta = m_tht * (hit.GetZ() + deltaZ) + b_tht;
                double targetPhi = m_phi * (hit.GetZ() + deltaZ) + b_phi;
                hit.theta = targetTheta;
                hit.phi = targetPhi;
            }
        }
    }

    hits.clear();
    for (const auto& z : sortedZ)
        for (const auto& h : layerMap[z])
            hits.push_back(h);
}


void MultipleScatteringSymmetricCorrection(std::vector<HitPoint>& hits) {
    std::vector<HitPoint> original = hits;
    std::vector<HitPoint> forward = hits;
    std::vector<HitPoint> backward = hits;

    MultipleScatteringRotationCorrection(forward);
    MultipleScatteringRotationCorrectionBackward(backward);

    // Combine
    hits.clear();
    for (size_t i = 0; i < original.size(); ++i) {
        HitPoint avg(
            0.5 * (forward[i].theta + backward[i].theta),
            0.5 * (forward[i].phi + backward[i].phi),
            original[i].GetZ()
        );
        hits.push_back(avg);
    }
}

std::vector<HitPointErrors> ComputeLayerwiseAverages(const std::vector<HitPoint>& hits) {
    std::map<double, std::vector<HitPoint>> layerMap;
    for (const auto& hit : hits) {
      layerMap[hit.GetZ()].push_back(hit);
    }

    std::vector<HitPointErrors> result;
    for (const auto& [zval, layerHits] : layerMap) {
        const size_t n = layerHits.size();
        if (n < 2){

	  // Need to catch case where there is one hit
	  result.emplace_back(layerHits[0].theta,layerHits[0].phi,zval,0.00002,0.001,0.015);
	}
	else
	  {
        double sum_theta = 0, sum_phi = 0;
        for (const auto& h : layerHits) {
            sum_theta += h.theta;
            sum_phi += h.phi;
        }

        double mean_theta = sum_theta / n;
        double mean_phi = sum_phi / n;

        double var_theta = 0, var_phi = 0;
        for (const auto& h : layerHits) {
            var_theta += (h.theta - mean_theta) * (h.theta - mean_theta);
            var_phi += (h.phi - mean_phi) * (h.phi - mean_phi);
        }
	  

        double sigma_theta = std::sqrt(var_theta / (n - 1)) / std::sqrt(n);
        double sigma_phi = std::sqrt(var_phi / (n - 1)) / std::sqrt(n);
        double sigma_z = 0.0; // all hits in layer have same Z

        result.emplace_back(mean_theta, mean_phi, zval, sigma_theta, sigma_phi, sigma_z);
	  }
    }

    return result;
}
