
void OptimizeTrackerSetup(Int_t N_TRACKER_LAYERS, Double_t TRACKER_THICKNESS_UM, char *PLOT_NAME)
{
  gStyle->SetOptStat(0);
  // Number of simulation events and bins
  const int nEvents = 100000;
  const int nBins = 200;
  const double xMin = 0, xMax = 5000 * N_TRACKER_LAYERS*TRACKER_THICKNESS_UM/(50.0*30.0);

  // Output TGraph
  const int nSteps = 500;
  TGraphErrors* gOverlap = new TGraphErrors(nSteps);
  
  TRandom3 randGen(0);

  double ONE_MIP_MPV = 80.0*TRACKER_THICKNESS_UM/300.0; //Approximate scaling
  double ONE_MIP_WID = 15.0*sqrt(TRACKER_THICKNESS_UM/300.0); //also appx scaling

  for (int step = 0; step < nSteps; ++step) {
    double perc_cut = 0.01 + 5 * step * 0.0099;  // Range: 0.01 to 1.0

    // Histograms for this cut
    TH1D* hd  = new TH1D("hd",  "2 MIPs", nBins, xMin, xMax);
    TH1D* h2d = new TH1D("h2d", "1 MIP",  nBins, xMin, xMax);

    for (int i = 0; i < nEvents; ++i) {
      double sum1 = 0, sum2 = 0;
      for (int j = 0; j < N_TRACKER_LAYERS; ++j) {
        //double tmp1 = randGen.Landau(2*80.0/10.0, 2*15 * sqrt(0.1));
        //double tmp2 = randGen.Landau(80.0/10.0, 15 * sqrt(0.1));
        double tmp1 = randGen.Landau(2*ONE_MIP_MPV, 2*ONE_MIP_WID);
        double tmp2 = randGen.Landau(ONE_MIP_MPV, ONE_MIP_WID);
        if (fabs(tmp1 - 2.0*ONE_MIP_MPV) < perc_cut * 2.0*ONE_MIP_MPV) sum1 += tmp1;
        if (fabs(tmp2 - ONE_MIP_MPV)  < perc_cut * ONE_MIP_MPV) sum2 += tmp2;
      }
      hd->Fill(sum1);
      h2d->Fill(sum2);
    }

    // Normalize histograms
    hd->Scale(1.0 / hd->Integral());
    h2d->Scale(1.0 / h2d->Integral());

    // Compute overlap
    double overlap = 0;
    for (int i = 1; i <= nBins; ++i) {
      double c1 = hd->GetBinContent(i);
      double c2 = h2d->GetBinContent(i);
      overlap += std::min(c1, c2);
    }

    // Fill TGraph
    gOverlap->SetPoint(step, perc_cut, overlap);
    gOverlap->SetPointError(step, 0.0, overlap / sqrt(100000.0-1.0));

    delete hd;
    delete h2d;
  }

  // Draw the graph
  TCanvas* c1 = new TCanvas("c1", "Overlap vs Cut", 1200, 800);
  c1->Divide(2,1);
  c1->cd(1);
  gOverlap->SetTitle("Percent Overlap vs MPV-Centered Energy Window Cut;MPV Window Cut (%);Overlap Fraction (%)");
  gOverlap->SetLineWidth(2);
  gOverlap->SetLineColor(kBlue+1);
  gOverlap->Draw("ALPE");

  // Find and print the minimum overlap point, especially so that we can plot it on the second pad
  double minVal = 1e9, bestCut = -1;
  for (int i = 0; i < gOverlap->GetN(); ++i) {
    double x, y;
    gOverlap->GetPoint(i, x, y);
    if (y < minVal) {
      minVal = y;
      bestCut = x;
    }
  }
  std::cout << "Best perc_cut = " << bestCut << ", Min overlap = " << minVal << std::endl;

  TLegend *tl = new TLegend(0.4,0.55,0.7,0.7);
  tl->AddEntry("",Form("Best cut = %.4f (%)",bestCut),"");
  tl->AddEntry("",Form("Overlap = %.4f (%)",minVal),"");
  tl->SetTextSize(0.045);
  tl->SetBorderSize(0);
  tl->Draw("same");
  
  c1->cd(2);
  gPad->SetLeftMargin(0.15);  // To make more room for y axis title
  // Histograms for this cut
  TH1D* hd  = new TH1D("hd",  "2 MIPs", nBins, xMin, xMax);
  hd->SetLineColor(kRed);
  hd->SetTitle(Form("MAPS %.1f #mu m %i layer Tracker, Best Cut Value, Idealized 1 MIP vs 2 MIP;Tracker Energy Deposited Sum After Cuts (keV);Normalized Counts for 100k Events",TRACKER_THICKNESS_UM,N_TRACKER_LAYERS));
  TH1D* h2d = new TH1D("h2d", "1 MIP",  nBins, xMin, xMax);
  h2d->SetTitle(Form("MAPS %.1f #mu m %i layer Tracker, Best Cut Value, Idealized 1 MIP vs 2 MIP;Tracker Energy Deposited Sum After Cuts (keV);Normalized Counts for 100k Events",TRACKER_THICKNESS_UM,N_TRACKER_LAYERS));
  
  for (int i = 0; i < nEvents; ++i) {
    double sum1 = 0, sum2 = 0;
    for (int j = 0; j < N_TRACKER_LAYERS; ++j) {
      double tmp1 = randGen.Landau(2*ONE_MIP_MPV,2.0*ONE_MIP_WID);
      double tmp2 = randGen.Landau(ONE_MIP_MPV,ONE_MIP_WID);
      if (fabs(tmp1 - 2.0*ONE_MIP_MPV) < bestCut * 2.0*ONE_MIP_MPV) sum1 += tmp1;
      if (fabs(tmp2 - ONE_MIP_MPV)  < bestCut * ONE_MIP_MPV) sum2 += tmp2;
    }
    hd->Fill(sum1);
    h2d->Fill(sum2);
  }

  // Normalize histograms
  hd->Scale(1.0 / hd->Integral());
  h2d->Scale(1.0 / h2d->Integral());
  h2d->Draw("hist");
  hd->Draw("hist same");

  TLegend *tl2 = new TLegend(0.5,0.3,0.8,0.8);
  tl2->SetBorderSize(0.0);
  tl2->SetTextSize(0.045);
  tl2->AddEntry(h2d,"1 MIP","l");
  tl2->AddEntry("",Form("#mu = %.1f (keV)",h2d->GetMean()),"");
  tl2->AddEntry("",Form("#sigma = %.1f (keV)",h2d->GetStdDev()),"");
  tl2->AddEntry(hd,"2 MIP","l");
  tl2->AddEntry("",Form("#mu = %.1f (keV)",hd->GetMean()),"");
  tl2->AddEntry("",Form("#sigma = %.1f (keV)",hd->GetStdDev()),"");
  tl2->Draw("same");

  c1->cd(0);
  c1->SaveAs(Form("%s.pdf",PLOT_NAME));
}
