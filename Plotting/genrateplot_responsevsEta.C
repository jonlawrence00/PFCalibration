#include <stdio.h>
// #include<conio.h>

void genrateplot_responsevsEta(TString legendname, TString plot_name, TString path2, TString plot_name2, TString region)
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
  TString hname,hname1,pdf="plots/pdf/"+plot_name+"_eta.pdf",gif="plots/gif/"+plot_name+"_eta.gif",png="plots/png/"+plot_name+"_eta.png";
  bool xaxis_pt=false;
  TString path_root1="../../Run3/trial2/response/";
  path_root1="../../Run3_v2/noPU_v2/noPU_v2_trial/response/";
  path_root1="chi2_results/response/";
  TString path_root2=path2;
  pdf=path2+"/pdf/"+plot_name+".pdf";gif=path2+"/pdf/"+plot_name+".gif";png=path2+"/pdf/"+plot_name+".png";

  if(xaxis_pt==false){
    hname1= path_root2+plot_name2+".root";
    hname= path_root1+"resp_reso_"+plot_name+".root";
  }

  cout<<"file name --> "<<plot_name<<endl;  
  //  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,106,500,472);
  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,162,714,444);
  Canvas_1_n2->Range(-0.5064935,-0.5,3.38961,0.5);
  Canvas_1_n2->SetFillColor(0);
  Canvas_1_n2->SetBorderMode(0);
  Canvas_1_n2->SetBorderSize(2);
  Canvas_1_n2->SetGridx(0);
  Canvas_1_n2->SetGridy(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetLeftMargin(0.13);
  TFile * inputfile1_ = new TFile(hname,"READ");
  TFile * inputfile2_ = new TFile(hname1,"READ");
  TFile * inputfile3_ = new TFile(hname1,"READ");
  TGraph* graph1_ = (TGraph*) inputfile1_ -> Get("Graph;1");
  TH2F *respHisto__1 = (TH2F*) inputfile1_ -> Get("respHisto");
  TGraph* graph2_ = (TGraph*) inputfile2_ -> Get(vname);
  TGraph* graph3_ = (TGraph*) inputfile3_ -> Get(tname);



    respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   respHisto__1->SetLineColor(ci);
   respHisto__1->GetXaxis()->SetTitleSize(0.045);
   respHisto__1->GetXaxis()->SetRangeUser(0,2.75);
   respHisto__1->GetYaxis()->SetRangeUser(-0.25,0.25);
   respHisto__1->GetXaxis()->SetTitleSize(0.06);
   respHisto__1->GetXaxis()->SetTitleOffset(0.79);
   respHisto__1->GetXaxis()->SetLabelSize(0.05);
   respHisto__1->GetXaxis()->SetLabelOffset(0.005);

   respHisto__1->GetYaxis()->SetTitleSize(0.05);
   respHisto__1->GetYaxis()->SetTitleOffset(0.8);
   respHisto__1->GetYaxis()->SetLabelSize(0.05);
   respHisto__1->GetYaxis()->SetLabelOffset(0.005);

   respHisto__1->Draw("");
   
   graph1_->SetMarkerColor(kRed);
   graph2_->SetMarkerColor(kBlue);
   graph3_->SetMarkerColor(kGreen);
   graph1_->SetMarkerStyle(22);
   graph2_->SetMarkerStyle(22);
   graph3_->SetMarkerStyle(22);
   graph1_->SetMarkerSize(0.8);
   graph2_->SetMarkerSize(0.8);
   graph3_->SetMarkerSize(0.8);

   // graph3->Draw("p");
   // graph2->Draw("p");
   //   graph1->Draw("p");
   graph3_->Draw("p");
   graph2_->Draw("p");
   graph1_->Draw("p");

   // TLine *line = new TLine(1.2273,0,498.0851,0);
   TLine *line = new TLine(0,0,2.75,0);
   if(xaxis_pt) line = new TLine(0,0,100,0);
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);

   TLegend* legends = new TLegend(0.55, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box                                                        
   legends->SetFillColor(0);
   legends->SetHeader(legendname,"C");
   legends->AddEntry(graph1_,"using Run3 calib (#chi^{2})","P");
   legends->AddEntry(graph2_,"using DRN (Validation)","P");
   legends->AddEntry(graph3_,"using DRN (Training)","P");
   legends->SetTextSize(0.04);
   legends->Draw();
   Canvas_1_n2->SaveAs(pdf);
   Canvas_1_n2->SaveAs(png);
   Canvas_1_n2->SaveAs(gif);

}
