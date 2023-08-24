#include <stdio.h>
// #include<conio.h>

void genrateplot_resovsE_separate(TString legendname, TString plot_name, TString path2, TString plot_name2, TString region)
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2

    TString vname,vname_,tname,tname_;
  if(region.Contains("barrel")){
    vname="Valid_response_0.000000_1.550000";
    vname_="Valid_reso_0.000000_1.550000";
    tname="Train_response_0.000000_1.550000";
    tname_="Train_reso_0.000000_1.550000";
  }
  else if(region.Contains("ec_in")){
    vname="Valid_response_1.550000_2.500000";
    vname_="Valid_reso_1.550000_2.500000";
    tname="Train_response_1.550000_2.500000";
    tname_="Train_reso_1.550000_2.500000";
  }
  else if(region.Contains("ec_out")){
    vname="Valid_response_2.500000_2.750000";
    vname_="Valid_reso_2.500000_2.750000";
    tname="Train_response_2.500000_2.750000";
    tname_="Train_reso_2.500000_2.750000";
  }
  else{
    vname="Valid_norm_pred_trueEn";
    vname_="Valid_norm_pred_trueEn_reso";
    tname="Train_norm_pred_trueEn";
    tname_="Train_norm_pred_trueEn_reso";
  }
  TString hname,hname1,pdf="plots/pdf/"+plot_name+".pdf",gif="plots/gif/"+plot_name+".gif",png="plots/png/"+plot_name+".png";
  bool xaxis_pt=false;
  TString path_root1="../../Run3/trial2/response/";
  path_root1="../../Run3_v2/noPU_v2/noPU_v2_trial/response/";
  path_root1="chi2_results/response/";
  TString path_root2=path2;
  pdf=path2+"/pdf/"+plot_name+"_reso.pdf";gif=path2+"/pdf/"+plot_name+"_reso.gif";png=path2+"/pdf/"+plot_name+"_reso.png";

  if(xaxis_pt==false){
    hname1= path_root2+plot_name2+".root";
    hname= path_root1+"resp_reso_"+plot_name+".root";
  }

  cout<<"file name --> "<<plot_name<<endl;  
  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,106,500,472);
  Canvas_1_n2->cd();
  /*
  Canvas_1_n2->SetFillColor(0);
  Canvas_1_n2->SetBorderMode(0);
  Canvas_1_n2->SetBorderSize(2);
  Canvas_1_n2->SetGridx(0);
  Canvas_1_n2->SetGridy(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetLeftMargin(0.13);
  */
  TFile * inputfile1_ = new TFile(hname,"READ");
  TFile * inputfile2_ = new TFile(hname1,"READ");
  TFile * inputfile3_ = new TFile(hname1,"READ");
  TGraph* graph1_ = (TGraph*) inputfile1_ -> Get("resolution;1");
  TGraph* graph2_ = (TGraph*) inputfile2_ -> Get(vname_);
  TGraph* graph3_ = (TGraph*) inputfile3_ -> Get(tname_);
  TH2F *respHisto__1 = (TH2F*) inputfile1_ -> Get("resoHisto");


  TPad *pad1 = new TPad("pad1","pad1",0.03,0,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetBottomMargin(0.355);
  pad1->SetRightMargin(0.03);
  pad1->Draw();pad1->SetGridx(0);

  pad1->cd();

    respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   respHisto__1->SetLineColor(ci);
   respHisto__1->GetXaxis()->SetRangeUser(4,340);
   respHisto__1->GetYaxis()->SetRangeUser(0.005,0.6);
   respHisto__1->Draw("");
   respHisto__1->GetXaxis()->SetTitleSize(0.045);
   respHisto__1->GetYaxis()->SetTitleSize(0.048);
   respHisto__1->GetYaxis()->SetTitleOffset(1.0);   
   graph1_->SetMarkerColor(kRed);
   graph2_->SetMarkerColor(kBlue);
   graph3_->SetMarkerColor(kGreen);
   graph1_->SetMarkerStyle(22);
   graph2_->SetMarkerStyle(22);
   graph3_->SetMarkerStyle(22);
   graph1_->SetMarkerSize(0.8);
   graph2_->SetMarkerSize(0.75);
   graph3_->SetMarkerSize(0.8);

   // graph3->Draw("p");
   // graph2->Draw("p");
   //   graph1->Draw("p");
   graph3_->Draw("p");
   graph2_->Draw("p");
   graph1_->Draw("p");
   TLine *line_ = new TLine(250,0.005,250,0.25);
   line_->SetLineColor(2);
   line_->SetLineStyle(2);
   line_->SetLineWidth(2);
   line_->Draw();

   TLegend* legends = new TLegend(0.42, 0.7, 0.97, 0.9,"","brNDC"); // the numbers determine the position of the box
   
   legends->SetFillColor(0);
   legends->SetHeader(legendname,"C");
   legends->AddEntry(graph1_,"using Run3 calib (#chi^{2})","P");
   legends->AddEntry(graph2_,"using DRN (Validation)","P");
   legends->AddEntry(graph3_,"using DRN (Training)","P");
   legends->SetTextSize(0.04);
   legends->Draw();
   TGraph *ratioGraph = new TGraph();
   TH1D *h_ratio = new TH1D("h_ratio","ratio",500,0,500);
   // Calculate the ratio for each point in the TGraphs                                                                                                                     
   int numPoints = 50;
   for (int i = 0; i <= numPoints; ++i) {
     double xNum, yNum;
     graph2_->GetPoint(i, xNum, yNum);
     double xDenom, yDenom;
     graph1_->GetPoint(i, xDenom, yDenom);
     double ratio = (yDenom) / yNum;
     ratioGraph->SetPoint(i, xNum, ratio);
     cout<<" xNum = "<<xNum<<" , yNum = "<<yNum<<" , yDenom = "<<yDenom<<endl;
     h_ratio->Fill(xNum,0);
   }
   for (int i = 0; i <= h_ratio->GetNbinsX(); ++i) {
     h_ratio->SetBinError(i,0);
   }
   ratioGraph->Print("all");

   TPad *pad2 = new TPad("pad2","pad2",0.001,0.0,1,0.35);
   pad2->Draw();
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.3);
   pad2->SetRightMargin(0.03);
   pad2->cd();
   pad2->SetGridx(0);
   pad2->SetGridy(0);
   h_ratio->SetTitle(0);
   ratioGraph->SetMarkerStyle(20);
   ratioGraph->SetMarkerSize(0.5);
   gStyle->SetOptStat(0);                                                                                                                                               
   h_ratio->Draw("p");
   h_ratio->GetXaxis()->SetRangeUser(4,340);
   if(plot_name2.Contains("_EH")){
     h_ratio->SetMaximum(2.9);
     h_ratio->SetMinimum(0.5);
   }
   else{
     h_ratio->SetMaximum(1.9);
     h_ratio->SetMinimum(0.1);
   }    
   h_ratio->GetYaxis()->SetNdivisions(5);
   h_ratio->GetYaxis()->SetTitle("#frac{chi2}{DRN (Validation)}");
   h_ratio->GetYaxis()->SetTitle("#chi^{2} / DRN (Valid.)");
   //   h_ratio->GetYaxis()->SetTitleCentered();
   h_ratio->GetYaxis()->SetTitleSize(0.125);
   h_ratio->GetYaxis()->SetLabelSize(0.10);
   h_ratio->GetXaxis()->SetTitle("E_{true} [GeV]");
   h_ratio->GetXaxis()->SetTitleSize(0.12);
   h_ratio->GetXaxis()->SetLabelSize(0.10);
   h_ratio->GetXaxis()->SetTitleOffset(1.0);
   h_ratio->GetYaxis()->SetTitleOffset(0.39);
   h_ratio->GetXaxis()->SetLabelOffset(0.016);

   ratioGraph->Draw("p");
   TLine *line = new TLine(4,1,340,1);
   if(xaxis_pt) line = new TLine(0,0,100,0);
   line->SetLineColor(1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();         
   TLine *line1 = new TLine(250,1,250,2.9);
   line1->SetLineColor(2);
   line1->SetLineStyle(2);
   line1->SetLineWidth(2);
   line1->Draw();

   Canvas_1_n2->SaveAs(pdf);
   Canvas_1_n2->SaveAs(png);
   Canvas_1_n2->SaveAs(gif);

}
