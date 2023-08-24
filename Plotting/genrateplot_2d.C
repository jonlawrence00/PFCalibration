#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"TGraph2D.h"

void genrateplot_2d(TString had="EH")
{
   TCanvas *canvas = new TCanvas("covariance", "Covariance between search region");
    canvas->cd();
    //    canvas->SetLeftMargin(0.08);
    canvas->SetRightMargin(0.14);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  int n=9;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.50, 0.50, 1.00, 1.00, 1.00 };
  Double_t green[NRGBs] = { 0.50, 1.00, 1.00, 0.60, 0.50 };
  Double_t blue[NRGBs]  = { 1.00, 1.00, 0.50, 0.40, 0.50 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f ");
  TFile* inputFile;
  if(had=="EH") inputFile= new TFile("output/target_trueE/Trained_2_350_eta3_epoch100/Output_eta2pt75_EH_2d.root", "READ");
  else inputFile= new TFile("output/target_trueE/Trained_2_350_eta3_epoch100/Output_eta2pt75_H_2d.root", "READ");
  TH2D *h_num;
  h_num = (TH2D*)inputFile->FindObjectAny("hist2D_norm_Eta_predTrue_Train");
  h_num->SetTitle(0);
  if(had=="EH")	 h_num->SetTitle("Response for EH hadrons");
  else h_num->SetTitle("Response for H hadrons");

  h_num->GetYaxis()->SetTitle("E_{true} (GeV)");
  h_num->GetYaxis()->SetTitleSize(0.045);
  h_num->GetYaxis()->SetLabelSize(0.035);
  h_num->GetYaxis()->SetTitleOffset(0.98);
  h_num->GetYaxis()->SetRangeUser(0,340);
  h_num->GetXaxis()->SetTitle("|#eta|");
  h_num->GetXaxis()->SetTitleSize(0.045);
  h_num->GetXaxis()->SetLabelSize(0.035);
  h_num->GetXaxis()->SetTitleOffset(0.98);
  h_num->GetXaxis()->SetRangeUser(0,2.94);
  h_num->GetZaxis()->SetRangeUser(-0.3,0.3);
  h_num->SetMarkerSize(1.25);

  h_num->Draw("colz text e");
  TLine *line = new TLine(0,280,2.94,280);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  TString pdf;
  if(had=="EH") pdf="output/target_trueE/Trained_2_350_eta3_epoch100/2d_EH";
  else pdf="output/target_trueE/Trained_2_350_eta3_epoch100/2d_H";
  canvas->SaveAs(pdf+".pdf");
  canvas->SaveAs(pdf+".png");

}
