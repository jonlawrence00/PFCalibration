#include <stdio.h>
// #include<conio.h>
void overlay_script_responsevsEta()
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2
//=========  (Mon Apr  1 19:35:54 2019) by ROOT version 6.14/06
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* plot_name = new char[200];
  char* legendname = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  char* path1 = new char[200];
  char* path2 = new char[200];
  char* rootfile1 = new char[200];
  char* rootfile2 = new char[200];


  sprintf(legendname,"EH hadrons");
  sprintf(plot_name, "EH");//correta_H_barrel_10_6_0_pre2_pt.root

  sprintf(hname1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/response_plots_root/Eta_correta_%s.root",plot_name);
  sprintf(hname,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/response_plots_root/E_correta_%s.root",plot_name);
  sprintf(path,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/response_plots/pdf/Eta_correta_%s.pdf",plot_name);
  sprintf(path1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/response_plots/png/Eta_correta_%s.png",plot_name);
  sprintf(path2,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/response_plots/gif/Eta_correta_%s.gif",plot_name);
  sprintf(new_rel,"new calib (10_6_0)");
  sprintf(old_rel,"old calib (10_0_3)");


   TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,52,525,527);
   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetGridx();
   Canvas_1_n2->SetGridy();
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TFile * inputfile1 = new TFile(hname,"READ");
   TFile * inputfile2 = new TFile(hname1,"READ");
   TGraph* graph1 = (TGraph*) inputfile1 -> Get("Graph;1");
   TGraph* graph2 = (TGraph*) inputfile2 -> Get("Graph;1");
   TH2F *respHisto__1 = (TH2F*) inputfile1 -> Get("respHisto");

   // respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   respHisto__1->SetLineColor(ci);
   respHisto__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.4076218,0.9338535,0.5923782,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Response");
   pt->Draw();

   
   graph1->SetName("Graph1");
   graph1->SetTitle("Graph");
   graph1->SetFillStyle(1000);
   graph1->SetMarkerColor(2);
   graph1->SetMarkerStyle(22);
   graph1->SetMarkerSize(1.1);

   graph2->SetName("Graph2");
   graph2->SetTitle("Graph");
   graph2->SetFillStyle(1000);
   graph2->SetMarkerColor(4);
   graph2->SetMarkerStyle(22);
   graph2->SetMarkerSize(1.1);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","Graph",100,0,554.4);
   Graph_Graph01->SetMinimum(-540.4823);
   Graph_Graph01->SetMaximum(5940.616);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetLabelFont(42);
   Graph_Graph01->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetXaxis()->SetTitleFont(42);
   Graph_Graph01->GetYaxis()->SetLabelFont(42);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleOffset(0);
   Graph_Graph01->GetYaxis()->SetTitleFont(42);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph1->SetHistogram(Graph_Graph01);
   

   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,0,554);
   Graph_Graph12->SetMinimum(0.0);
   Graph_Graph12->SetMaximum(0.0);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleOffset(0);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   graph2->SetHistogram(Graph_Graph12);


   Double_t x[100], y[100];
   Int_t n = 20;
   // for (Int_t i=0;i<n;i++) {
   //   x[i] = i*0.1;
   //   y[i] = 10*sin(x[i]);
   // }
   x[19]=0.5;
   y[19]=11;

   // gr = new TGraph(n,x,y);
   // gr->Draw("AC*");   
   graph2->Draw("p");
   graph1->Draw("p");

   // TLine *line = new TLine(1.2273,0,498.0851,0);
   TLine *line = new TLine(0,0,3,0);
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);


   TLegend* legends = new TLegend(0.44, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   legends->SetFillColor(0); 
   legends->SetHeader(legendname); 
   // legends->AddEntry(graph1,"Raw (no correction)","P");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   // legends->AddEntry(graph2,"final correction","P");
   legends->AddEntry(graph1,old_rel,"P");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   legends->AddEntry(graph2,new_rel,"P");
   legends->SetTextSize(0.04);
   //   legends->SetMarkerStyle(1);
   legends->Draw();
   // Canvas_1_n2->SaveAs(path);
   // Canvas_1_n2->SaveAs(path1);
   // Canvas_1_n2->SaveAs(path2);

}
