#include <TFile.h>
#include <TTree.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"

void TMVA_MultiClass_v3() {
    TMVA::Tools::Instance();

    TString outfileName = "TMVA_MultiClassOut_2.root";
    TFile* outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("TMVAMultiClass", outputFile,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Multiclass");

    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    // Load the ROOT files and trees
    TFile *eleFile = TFile::Open("BigEleFix.root");
    TFile *gamFile = TFile::Open("BigGamFix.root");
    TFile *muFile  = TFile::Open("BigMuMFix.root");
    TFile *posFile = TFile::Open("BigElePFix.root");    
    TFile *tauFile  = TFile::Open("BigTauMFix.root");
    TFile *piFile  = TFile::Open("BigPiMFix.root");
    TFile *kFile  = TFile::Open("BigKMinusFix.root");

    TTree *eleTree = (TTree*)eleFile->Get("ExTree");
    TTree *gamTree = (TTree*)gamFile->Get("ExTree");
    TTree *muTree  = (TTree*)muFile->Get("ExTree");
    TTree *kTree = (TTree*)kFile->Get("ExTree");
    TTree *piTree = (TTree*)piFile->Get("ExTree");
    TTree *tauTree  = (TTree*)tauFile->Get("ExTree");
    TTree *posTree = (TTree*)posFile->Get("ExTree");
    
    // Add training variables
    dataloader->AddVariable("eLhood_helix", 'F');
    dataloader->AddVariable("eLhood_cnt", 'F');
    dataloader->AddVariable("eLhood_tof", 'F');
    dataloader->AddVariable("eLhood_ene", 'F');
    dataloader->AddVariable("cal_LR_Asym", 'F');
    dataloader->AddVariable("cal_UD_Asym", 'F');
    dataloader->AddVariable("deltaR", 'F');
    dataloader->AddVariable("shower_ini", 'F');
    //to add each quadrant deltaR value
    //dataloader->AddVariable("v_deltaR[1]", 'F');
    //dataloader->AddVariable("v_deltaR[2]", 'F');
    //dataloader->AddVariable("v_deltaR[3]", 'F');
    //dataloader->AddVariable("v_deltaR[4]", 'F');
    dataloader->AddVariable("qene", 'F');
    dataloader->AddVariable("qpt", 'F');
    dataloader->AddVariable("qpz", 'F');
    
    // Add each tree with a class name
    dataloader->AddTree(eleTree, "electron", 1.0);
    dataloader->AddTree(gamTree, "photon", 1.0);
    dataloader->AddTree(muTree, "muon", 1.0);
    dataloader->AddTree(kTree, "kMinus", 1.0);
    dataloader->AddTree(piTree, "piMinus", 1.0);
    dataloader->AddTree(tauTree, "tauMinus", 1.0);
    dataloader->AddTree(posTree, "positron", 1.0);
    
    // Prepare dataset
    dataloader->PrepareTrainingAndTestTree("", "",
        "SplitMode=Random:NormMode=NumEvents:!V");

    // Book a multiclass BDT
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG_MC",
        "!H:!V:NTrees=300:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:GradBaggingFraction=0.5:MaxDepth=6:nCuts=20");

    // Book a multiclass neural network
    factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP_MC",
        "!H:!V:VarTransform=N:NCycles=500:HiddenLayers=N+5:NeuronType=tanh:TrainingMethod=BFGS");

    // Train, test, evaluate
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    std::cout << "=== Multiclass training complete ===" << std::endl;
    std::cout << "Run: TMVA::TMVAGui(\"" << outfileName << "\") to view results." << std::endl;
}
