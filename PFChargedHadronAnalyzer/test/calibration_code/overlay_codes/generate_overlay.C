#include <stdio.h>
// #include<conio.h>
void generate_overlay()
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2
//=========  (Mon Apr  1 19:35:54 2019) by ROOT version 6.14/06
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* plot_name = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  char* path_root1 = new char[200];
  char* path_root2 = new char[200];
  char* rootfile1 = new char[200];
  char* rootfile2 = new char[200];

  // sprintf(old_release, "10_0_3");
  // sprintf(new_release, "10_6_0_pre2");
  // sprintf(plot_name, "correta_H_ec_out");
  sprintf(hist, "");
  sprintf(path_root1,"");
  sprintf(path_root2,"");
  sprintf(rootfile1,"");
  sprintf(rootfile2,"");


  // sprintf(stat_name,"Etrue bin %d-%d GeV",n[i],n[i+1]);
  // cout<<"\n i "<<i<<" = "<<n[i];}
  // cout<<stat_name<<"\n";


  sprintf(hname,"%s.root",path_root1);//,rootfile1);
  sprintf(hname1,"%s.root",path_root2);//,rootfile2);
  // sprintf(plot_name,"%s",plot);
  sprintf(hist_name,"%s",hist);
  sprintf(path,"../overlay_response_plots/%s_vs_%s.gif",rootfile1,rootfile2);
  // sprintf(old_rel,"%s",old_release);
  // sprintf(new_rel,"%s",new_release);



      // sprintf(hname,"%s.root",plot);    
      // cout<<"1 \n";

      // sprintf(hname1,"%s.root",plot);
      // cout<<hname1;
      // sprintf(plot_name,"%s",plot);
      // sprintf(hist_name,"%s%d",hist,n[i]);
      // sprintf(path,"../1D_response_plots/%s/%s.gif",plot,hist_name);
      // sprintf(old_rel,"%s",old_release);
      // sprintf(new_rel,"%s",new_release);


   TCanvas *Canvas_1_n2 = new TCanvas(hist_name, hist_name,65,52,525,527);
   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetGridx();
   Canvas_1_n2->SetGridy();
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TFile * inputfile1 = new TFile(hname,"READ");
   TFile *inputfile2 = new TFile(hname1,"READ");
   TH1D *resp1 = (TH1D*) inputfile1 -> Get(hist_name);
   TH1D *resp2 = (TH1D*) inputfile2 -> Get(hist_name);
      // cout<<"2 \n";

   // respHisto__1->SetStats(0);

   // Int_t ci;      // for color index setting
   // TColor *color; // for color definition with alpha
   // ci = TColor::GetColor("#000099");
   // respHisto__1->SetLineColor(ci);
   // respHisto__1->Draw("");
   // TPaveStats *ptstats = new TPaveStats(0.65,0.50,0.9,0.75,"brNDC");
   // //   TPaveStats *ptstats = new TPaveStats(0.6260745,0.6433121,0.8810888,0.8832272,"brNDC");
   // //  ptstats->SetName(old_rel);
   // ptstats->SetBorderSize(1);
   // ptstats->SetFillColor(0);
   // ptstats->SetLineColor(2);
   // ptstats->SetTextAlign(12);
   // ptstats->SetTextColor(2);
   // ptstats->SetTextFont(42);
   // TText *ptstats_LaTex = ptstats->AddText(old_rel);
   // ptstats_LaTex->SetTextSize(0.0315317);
   // ptstats->SetOptStat(1111);
   // ptstats->SetOptFit(10001);
   // ptstats->Draw();
   // resp1->GetListOfFunctions()->Add(ptstats);
   // ptstats->SetParent(resp1);

   // resp1->SetLineColor(2);
   // resp1->SetLineWidth(2);
   // resp1->Scale(1.0/resp1->Integral());
   // resp1->GetXaxis()->SetRangeUser(xmin,xmax); 
   // resp1->GetYaxis()->SetRangeUser(ymin , ymax);                                                                                                                                                            
   // resp1->SetName(old_rel);
   // resp1->Draw("hists");

   TPaveStats *ptstats = new TPaveStats(0.590,0.7,0.9,0.9,"brNDC");
   // ptstats = new TPaveStats(0.6260745,0.403397,0.8810888,0.6433121,"brNDC");
   //   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetLineColor(4);
   ptstats->SetTextAlign(12);
   ptstats->SetTextColor(4);
   ptstats->SetTextFont(42);
   ptstats->AddText(hist_name);
   ptstats->SetTextSize(0.0315317);
   ptstats->SetOptStat(1111);
   // ptstats->SetOptFit(10001);
   ptstats->Draw();
   resp2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(resp2);

   resp2->SetLineColor(4);
   resp2->SetLineWidth(2);
   resp2->Scale(1.0/resp2->Integral());
   resp2->GetXaxis()->SetRangeUser(xmin , xmax);                                                                                                                                                               resp2->GetYaxis()->SetRangeUser(ymin , ymax);                                                                                                                                                               resp2->SetName(stat_name);
   // resp2->Fit("gaus");
   resp2->Draw(" hists");
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);
   resp1->Draw("sames");
   // TF1* f1 = (TF1*)resp1->GetFunction("gaus");
   // TF1* f2 = (TF1*)resp2->GetFunction("gaus");
   // // f1->SetLineColor(0);
   // // f2->SetLineColor(0);
   // f1->Delete();
   // f2->Delete();
   
   // TLegend* legends = new TLegend(0.65, 0.75, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   // legends->SetFillColor(0); 
   // legends->SetHeader(plot_name); 
   // // legends->AddEntry(resp1, "for old (10_0_3)","lep");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   // // legends->AddEntry(resp2, "for new (10_6_0_pre2)","lep");
   // legends->Draw();
   Canvas_1_n2->SaveAs(path);
   if (n[i]==494) continue;
    }
}
