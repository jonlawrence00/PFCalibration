#include <stdio.h>
// #include<conio.h>
void genrateplot_responsevsE_UL2018wrtrun3(TString legendname, TString plot_name, TString path2, TString plot_name2, TString region)
{
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
  //  path_root1="../../Run3_v2/noPU_v2/usingrechit/response/";
  TString path_root2=path2;
  pdf=path2+"/pdf/"+plot_name+".pdf";gif=path2+"/pdf/"+plot_name+".gif";png=path2+"/pdf/"+plot_name+".png";
  //TString path_root2="usingUL2018/response/";
  //  pdf="usingUL2018/plots/pdf/"+plot_name+".pdf";gif="usingUL2018/plots/gif/"+plot_name+".gif";png="usingUL2018/plots/png/"+plot_name+".png";
  if(xaxis_pt==false){
    hname1= path_root2+plot_name2+".root";
    hname= path_root1+"resp_reso_"+plot_name+".root";

  }
  else if(xaxis_pt==true){
    hname= path_root1+"resp_reso_"+plot_name+".root";
    hname1= path_root2+"resp_reso_"+plot_name2+".root";
  }
        
  cout<<"file name --> "<<plot_name<<endl;



  //  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",2142,137);
  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",2142,137,1020,565);
  //   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
  Canvas_1_n2->Divide(2,1);
  Canvas_1_n2->SetFillColor(0);
  Canvas_1_n2->SetBorderMode(0);
  Canvas_1_n2->SetBorderSize(2);
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  TFile * inputfile1 = new TFile(hname,"READ");
   TFile * inputfile2 = new TFile(hname1,"READ");

   TGraph* graph1 = (TGraph*) inputfile1 -> Get("response;1");
   TH2F *respHisto_1 = (TH2F*) inputfile1 -> Get("respHisto");
   TGraph* graph1_ = (TGraph*) inputfile1 -> Get("resolution;1");
   TH2F *resoHisto_1 = (TH2F*) inputfile1 -> Get("resoHisto");
   TGraph* graph2 = (TGraph*) inputfile2 -> Get(vname);
   TGraph* graph2_ = (TGraph*) inputfile2 -> Get(vname_);
   TGraph* graph3 = (TGraph*) inputfile2 -> Get(tname);
   TGraph* graph3_ = (TGraph*) inputfile2 -> Get(tname_);

   /*
   TGraph* graph2 = (TGraph*) inputfile2 -> Get("Valid_norm_pred_trueEn");
   //   TH2F *respHisto_2 = (TH2F*) inputfile2 -> Get("respHisto");
   TGraph* graph2_ = (TGraph*) inputfile2 -> Get("Valid_norm_pred_trueEn_reso");
   // TH2F *resoHisto_2 = (TH2F*) inputfile2 -> Get("resoHisto");
   TGraph* graph3 = (TGraph*) inputfile2 -> Get("Train_norm_pred_trueEn");
   TGraph* graph3_ = (TGraph*) inputfile2 -> Get("Train_norm_pred_trueEn_reso");
   */
   // respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");

   Canvas_1_n2->cd(1);
   gPad->SetGrid();
   //   Canvas_1_n2->SetGridy();
   Canvas_1_n2->Range(-68.50026,-0.6269485,353.1421,0.6232721);

   respHisto_1->SetLineColor(ci);
   respHisto_1->GetXaxis()->SetRangeUser(2,350);
   respHisto_1->Draw("");
   
   graph1->SetName("Graph1");
   graph1->SetTitle("Graph");
   graph1->SetFillStyle(1000);
   graph1->SetMarkerColor(kRed);
   graph1->SetMarkerStyle(22);
   graph1->SetMarkerSize(0.8);   

   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","Graph",100,0,554.4);
   Graph_Graph01->SetMinimum(-540.4823);
   Graph_Graph01->SetMaximum(5940.616);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetLabelFont(42);
   if(xaxis_pt==false)
     Graph_Graph01->GetXaxis()->SetLabelSize(0.035);
   else if(xaxis_pt==true)
     Graph_Graph01->GetXaxis()->SetLabelSize(0.0);
   Graph_Graph01->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetXaxis()->SetTitleFont(42);
   Graph_Graph01->GetYaxis()->SetLabelFont(42);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph01->GetYaxis()->SetTitleFont(42);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph1->SetHistogram(Graph_Graph01);
   graph1->GetXaxis()->SetRangeUser(4,190);
   graph1->Draw("p");
   graph2->SetMarkerStyle(22);   
   graph2->SetMarkerColor(kBlue);
   graph2->SetMarkerSize(0.8);
   graph2->Draw("p");
   graph3->SetMarkerStyle(22);
   graph3->SetMarkerColor(kOrange);
   graph3->SetMarkerSize(0.8);
   graph3->Draw("p");

   Canvas_1_n2->cd(2);
   gPad->SetGrid();
   resoHisto_1->SetLineColor(ci);
   resoHisto_1->GetXaxis()->SetRangeUser(4,350);
   resoHisto_1->Draw("");

   graph1_->SetName("Graph1_");
   graph1_->SetTitle("Graph");
   graph1_->SetFillStyle(1000);
   graph1_->SetMarkerColor(kRed);
   graph1_->SetMarkerStyle(22);
   //   graph1_->SetMarkerSize(1.1);
   TH1F *Graph_Graph02 = new TH1F("Graph_Graph02","Graph",100,0,554.4);
   Graph_Graph02->SetMinimum(-540.4823);
   Graph_Graph02->SetMaximum(5940.616);
   Graph_Graph02->SetDirectory(0);
   Graph_Graph02->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph02->SetLineColor(ci);
   Graph_Graph02->GetXaxis()->SetLabelFont(42);
   if(xaxis_pt==false)
     Graph_Graph02->GetXaxis()->SetLabelSize(0.035);
   else if(xaxis_pt==true)
     Graph_Graph02->GetXaxis()->SetLabelSize(0.0);
   Graph_Graph02->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph02->GetXaxis()->SetTitleFont(42);
   Graph_Graph02->GetYaxis()->SetLabelFont(42);
   Graph_Graph02->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph02->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph02->GetYaxis()->SetTitleOffset(0);
   Graph_Graph02->GetYaxis()->SetTitleFont(42);
   Graph_Graph02->GetZaxis()->SetLabelFont(42);
   Graph_Graph02->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph02->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph02->GetZaxis()->SetTitleFont(42);
   graph1->SetHistogram(Graph_Graph02);
   graph1_->GetXaxis()->SetRangeUser(4,190);
   graph1_->Draw("p");
   graph2_->SetMarkerStyle(22);
   graph2_->SetMarkerColor(kBlue);
   graph2_->Draw("p");
   graph3_->SetMarkerStyle(22);
   graph3_->SetMarkerColor(kOrange);
   graph3_->SetMarkerSize(0.8);
   graph3_->Draw("p");


   
   // TLine *line = new TLine(1.2273,0,498.0851,0);
   Canvas_1_n2->cd(1);    
   TLine *line;
   if (xaxis_pt==false)
     line = new TLine(4,0,350,0);
   else if(xaxis_pt==true)
     line = new TLine(0,0,350,0);
     
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   //Canvas_1_n2->Modified();
   //   Canvas_1_n2->cd();
   //
   //Canvas_1_n2->SetSelected(Canvas_1_n2);


   // TLegend* legends = new TLegend(0.44, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   TLegend* legends = new TLegend(0.35, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   legends->SetFillColor(0);
   legends->SetHeader(legendname,"C"); 
   legends->AddEntry(graph1,"using Run3 calib (chi2)","P");
   legends->AddEntry(graph2,"using GNN (Validation)","P");
   legends->AddEntry(graph3,"using GNN (Training)","P");
   legends->SetTextSize(0.04);
   legends->Draw();
   Canvas_1_n2->SaveAs(pdf);
   Canvas_1_n2->SaveAs(png);
   Canvas_1_n2->SaveAs(gif);

}
