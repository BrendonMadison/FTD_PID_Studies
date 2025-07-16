#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

void LikelihoodTreeCovCalc2(const std::string& input_rootfile,
                                 const std::string& tree_name,
                                 const std::string& base_name) {
    // Open ROOT file and tree
    TFile* file = TFile::Open(input_rootfile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open input file: " << input_rootfile << std::endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
    if (!tree) {
        std::cerr << "Failed to retrieve tree: " << tree_name << std::endl;
        file->Close();
        return;
    }

    const std::vector<std::string> labels = {
        "eLhood_ene", "eLhood_helix", "eLhood_cnt", "eLhood_tof"
    };

    const int N = labels.size();
    std::vector<double> vals(N, 0.0);
    for (int i = 0; i < N; ++i)
        tree->SetBranchAddress(labels[i].c_str(), &vals[i]);

    std::vector<double> sum(N, 0.0), sum2(N, 0.0);
    std::vector<std::vector<double>> cov(N, std::vector<double>(N, 0.0));

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        for (int j = 0; j < N; ++j) {
            sum[j] += vals[j];
            sum2[j] += vals[j] * vals[j];
            for (int k = 0; k <= j; ++k)
                cov[j][k] += vals[j] * vals[k];
        }
    }

    std::vector<double> mean(N), var(N);
    for (int j = 0; j < N; ++j) {
        mean[j] = sum[j] / nentries;
        var[j] = sum2[j] / nentries - mean[j] * mean[j];
    }

    // Symmetrize covariance matrix
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k <= j; ++k) {
            cov[j][k] = cov[j][k] / nentries - mean[j] * mean[k];
            cov[k][j] = cov[j][k];
        }
    }

    // Write CSVs
    std::ofstream f_labels(base_name + "_labels.csv");
    std::ofstream f_means(base_name + "_means.csv");
    std::ofstream f_vars(base_name + "_var.csv");
    std::ofstream f_covar(base_name + "_covar.csv");

    f_labels << std::fixed;
    for (const auto& label : labels)
        f_labels << label << "\n";

    f_means << std::setprecision(10);
    for (const auto& m : mean)
        f_means << m << "\n";

    f_vars << std::setprecision(10);
    for (const auto& v : var)
        f_vars << v << "\n";

    f_covar << std::setprecision(10);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f_covar << cov[i][j];
            if (j < N - 1) f_covar << ",";
        }
        f_covar << "\n";
    }

    f_labels.close();
    f_means.close();
    f_vars.close();
    f_covar.close();

    file->Close();
    std::cout << "Exported likelihood model to: " << base_name << "_*.csv" << std::endl;
}

