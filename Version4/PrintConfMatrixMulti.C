//To print the confusion matrix:
void PrintConfMatrixMulti(const char *INFILE, const TString& methodName = "MLP_MC", float QUALITY = 0.0, bool DENSITY = false)
{
  const std::vector<TString> classLabels = {"electron","photon","muon","kMinus","piMinus","tauMinus","positron", "low quality"};
  const int nClasses = 8;

    TFile *f = TFile::Open(INFILE);
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open file." << std::endl;
        return;
    }

    TTree *testTree = (TTree*)f->Get("dataset/TestTree");
    if (!testTree) {
        std::cerr << "Failed to find TestTree in file." << std::endl;
        return;
    }

    float score_electron, score_photon, score_muon;
    int classID;
    float mlp_scores[10];
    testTree->SetBranchAddress(methodName.Data(), mlp_scores);

    //testTree->SetBranchAddress((methodName + "__electron").Data(), &score_electron);
    //testTree->SetBranchAddress((methodName + "__photon").Data(), &score_photon);
    //testTree->SetBranchAddress((methodName + "__muon").Data(), &score_muon);
    //testTree->SetBranchAddress(methodName + ".electron", &score_electron);
    //testTree->SetBranchAddress(methodName + ".photon", &score_photon);
    //testTree->SetBranchAddress(methodName + ".muon", &score_muon);
    testTree->SetBranchAddress("classID", &classID);

    // Confusion matrix: rows = true label, cols = predicted label
    int confMat[nClasses][nClasses] = {}; // last row/col is "low quality"

    const Long64_t nEvents = testTree->GetEntries();
    for (Long64_t i = 0; i < nEvents; ++i) {
        testTree->GetEntry(i);

	//float score_electron = mlp_scores[0];
	//float score_photon   = mlp_scores[1];
	//float score_muon     = mlp_scores[2];
        //float scores[3] = {score_electron, score_photon, score_muon};
        int pred = nClasses-1; // default is "low quality"
        float maxScore = -1.0;

        for (int j = 0; j < nClasses-1; ++j) {
            if (mlp_scores[j] > maxScore) {
                maxScore = mlp_scores[j];
                pred = j;
            }
        }

        if (maxScore < QUALITY)
            pred = nClasses-1; // assign to "low quality"

        if (classID < 0 || classID > 8) {
            std::cerr << "Invalid classID = " << classID << " at entry " << i << std::endl;
            continue;
        }

        confMat[classID][pred]++;
    }

    // Print header
    std::cout << "\n=== Confusion Matrix for method: " << methodName << " ===" << std::endl;
    std::cout << "Rows = True Class | Columns = Predicted Class" << std::endl;
    std::cout << "         ";
    for (const auto& label : classLabels)
        std::cout << Form("%14s", label.Data());
    std::cout << std::endl;

    for (int i = 0; i < nClasses; ++i) {
        std::cout << Form("%10s", classLabels[i].Data());
        int rowSum = 0;
        for (int j = 0; j < nClasses; ++j) rowSum += confMat[i][j];

        for (int j = 0; j < nClasses; ++j) {
            if (DENSITY && rowSum > 0)
                std::cout << Form("%14.2f", 100.0 * confMat[i][j] / rowSum);
            else
                std::cout << Form("%14d", confMat[i][j]);
        }
        std::cout << std::endl;
    }

    std::cout << (DENSITY ? "\nValues are row-normalized percentages.\n" : "\nRaw classification counts.\n") << std::endl;
}

