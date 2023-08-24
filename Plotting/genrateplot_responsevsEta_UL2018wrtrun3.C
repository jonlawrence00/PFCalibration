#include <stdio.h>
// #include<conio.h>
void genrateplot_responsevsEta_UL2018wrtrun3(TString legendname, TString plot_name,  TString plot_name2)
{
  TString hname,hname1,pdf="plots/pdf/"+plot_name+".pdf",gif="plots/gif/"+plot_name+".gif",png="plots/png/"+plot_name+".png";
  bool xaxis_pt=false;
  TString path_root1="../../Run3/trial2/response/";
  TString path_root2="./";
  pdf="trial2/plots/pdf/"+plot_name+".pdf";gif="trial2/plots/gif/"+plot_name+".gif";png="trial2/plots/png/"+plot_name+".png";
  //  pdf="usingUL2018/plots/pdf/"+plot_name+".pdf";gif="usingUL2018/plots/gif/"+plot_name+".gif";png="usingUL2018/plots/png/"+plot_name+".png";
  hname1= path_root2+plot_name2+".root";
  hname= path_root1+"resp_reso_"+plot_name+".root";

        
  cout<<"file name --> "<<hname<<endl;



  TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",700,400);
   //   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
  //Canvas_1_n2->Divide(2,1);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TFile * inputfile1 = new TFile(hname,"READ");
   TFile * inputfile2 = new TFile(hname1,"READ");

   
   TGraph* graph1 = (TGraph*) inputfile1 -> Get("Graph;1");
   TH2F *respHisto_1 = (TH2F*) inputfile1 -> Get("respHisto");
   TGraph* graph2 = (TGraph*) inputfile2 -> Get("Valid_norm_pred_trueEn");
   TGraph* graph3 = (TGraph*) inputfile2 -> Get("Train_norm_pred_trueEn");
 

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");

   Canvas_1_n2->cd();
   gPad->SetGrid();
   //   Canvas_1_n2->SetGridy();

   respHisto_1->SetLineColor(ci);
   respHisto_1->GetXaxis()->SetRangeUser(0,2.5);

   respHisto_1->GetXaxis()->SetTitleSize(0.06);
   respHisto_1->GetXaxis()->SetTitleOffset(0.79);
   respHisto_1->GetXaxis()->SetLabelSize(0.05);
   respHisto_1->GetXaxis()->SetLabelOffset(0.005);
   
   respHisto_1->GetYaxis()->SetTitleSize(0.05);
   respHisto_1->GetYaxis()->SetTitleOffset(0.8);
   respHisto_1->GetYaxis()->SetLabelSize(0.05);
   respHisto_1->GetYaxis()->SetLabelOffset(0.005);

   respHisto_1->Draw("");
   
   graph1->SetName("Graph1");
   graph1->SetTitle("Graph");
   graph1->SetFillStyle(1000);
   graph1->SetMarkerColor(kRed);
   graph1->SetMarkerStyle(22);
   graph1->SetMarkerSize(0.8);   
   graph1->GetXaxis()->SetRangeUser(0,2.5);
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
   graph1->Draw("p");
   graph2->SetMarkerStyle(22);
   graph2->SetMarkerSize(0.8);
   graph2->SetMarkerColor(kBlue);
   graph2->Draw("p");
   graph3->SetMarkerStyle(22);
   graph3->SetMarkerColor(kOrange);
   graph3->SetMarkerSize(0.8);
   graph3->Draw("p");


   TLine *line;
   line = new TLine(0,0,2.5,0);
     
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   //Canvas_1_n2->Modified();
   //   Canvas_1_n2->cd();
   //
   //Canvas_1_n2->SetSelected(Canvas_1_n2);


   // TLegend* legends = new TLegend(0.44, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   TLegend* legends = new TLegend(0.55, 0.72, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   legends->SetFillColor(0); 
   legends->SetHeader(legendname,""); 
   legends->AddEntry(graph1,"using Run3 calib (#chi^{2})","P");
   legends->AddEntry(graph2,"using DRN (Validation)","P");
   legends->AddEntry(graph3,"using DRN (Training)","P");

   legends->SetTextSize(0.04);
    //   legends->SetMarkerStyle(1);
   legends->Draw();
   Canvas_1_n2->SaveAs(pdf);
   Canvas_1_n2->SaveAs(png);
   Canvas_1_n2->SaveAs(gif);

}
