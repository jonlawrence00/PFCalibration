#include <stdio.h>
// #include<conio.h>
void generate_1Dplots()
{
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* full_path = new char[200];
  char* path1 = new char[200];
  char* full_path1 = new char[200];
  char* path2 = new char[200];
  char* full_path2 = new char[200];
  char* plot_name = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  // sprintf(old_release, "corrEta");
  // sprintf(new_release, "corrEta");
  //  sprintf(plot_name, "EH_barrel");
   sprintf(plot_name, "EH_ec_out");
  sprintf(hist, "histcorhybrid");

  //comparision for old calib and new calib
  sprintf(hname1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/1D_response/correta_%s_10_6_0.root",plot_name);
  sprintf(hname,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/for_10_0_3_calib/1D_response/correta_%s_10_0_3.root",plot_name);
  sprintf(path,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/1D_response/gif/correta_%s",plot_name);
  sprintf(path1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/1D_response/pdf/correta_%s",plot_name);
  sprintf(path2,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_10_0_3/1D_response/png/correta_%s",plot_name);
  sprintf(new_rel,"new calib (10_6_0)");
  sprintf(old_rel,"old calib (10_0_3)");

  //comparision for old calib and new calib
  // sprintf(hname1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/1D_response/correta_%s_10_6_0.root",plot_name);
  // sprintf(hname,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/raw/1D_response/Raw_%s_10_6_0_pre2_response.root",plot_name);
  // sprintf(path,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_raw/1D_response/gif/correta_%s",plot_name);
  // sprintf(path1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_raw/1D_response/pdf/correta_%s",plot_name);
  // sprintf(path2,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/trial3/overlay_plots_with_raw/1D_response/png/correta_%s",plot_name);
  // sprintf(new_rel,"After Calibration");
  // sprintf(old_rel,"Without Calibration");

  //========= Raw E comparision for UL2016 and UL2017=====================
  // sprintf(hname1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/raw/1D_response/Raw_%s_10_6_0_pre2_response.root",plot_name);
  // sprintf(hname,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2017/raw/1D_response/Raw_%s_10_6_0_pre2.root",plot_name);
  // sprintf(path,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/raw/gif/overlay_plots_with_UL2017/1D_response/Raw_%s",plot_name);
  // sprintf(path1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/raw/pdf/overlay_plots_with_UL2017/1D_response/Raw_%s",plot_name);
  // sprintf(path2,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2016/raw/png/overlay_plots_with_UL2017/1D_response/Raw_%s",plot_name);
  // sprintf(new_rel,"for UL2016");
  // sprintf(old_rel,"for UL2017");

      
  double xmin= -1., xmax =1.;
  double ymin= 0, y1, y2, ymax=0.1;
  int n[69];
  n[0]=0;
  n[1]=2;
  for(int i=2; i<69; i++)
    {
      if (n[i-1]<12)
	n[i]=n[i-1]+2;
      else if (n[i-1]>=12 && n[i-1]<104)
	n[i]=n[i-1]+4;
      else
	n[i]=n[i-1]+10;
    }

  for(int i=2; i<69;i++)
    { 
      sprintf(stat_name,"Etrue bin %d-%d GeV",n[i],n[i+1]);

      sprintf(hist_name,"%s%d",hist,n[i]);
      sprintf(full_path,"%s/%s.gif",path,hist_name);
      sprintf(full_path1,"%s/%s.pdf",path1,hist_name);
      sprintf(full_path2,"%s/%s.png",path2,hist_name);
 
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
   TH1D *resp1 = (TH1D*)inputfile1->Get(hist_name);

   TFile *inputfile2 = new TFile(hname1,"READ");
   TH1D *resp3 = (TH1D*)inputfile2->Get(hist_name);

   // respHisto__1->SetStats(0);

   // Int_t ci;      // for color index setting
   // TColor *color; // for color definition with alpha
   // ci = TColor::GetColor("#000099");
   // respHisto__1->SetLineColor(ci);
   // respHisto__1->Draw("");
   TPaveStats *ptstats = new TPaveStats(0.65,0.67,0.9,0.9,"brNDC");
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
   ptstats->SetOptFit(10001);
   ptstats->Draw();
   resp3->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(resp3);

   TPaveStats *ptstats2 = new TPaveStats(0.65,0.47,0.9,0.67,"brNDC");
   //   TPaveStats *ptstats = new TPaveStats(0.6260745,0.6433121,0.8810888,0.8832272,"brNDC");
   //  ptstats->SetName(old_rel);
   ptstats2->SetBorderSize(1);
   ptstats2->SetFillColor(0);
   ptstats2->SetLineColor(kRed);
   ptstats2->SetTextAlign(12);
   ptstats2->SetTextColor(kRed);
   ptstats2->SetTextFont(42);
   TText *ptstats_LaTex = ptstats2->AddText(old_rel);
   ptstats_LaTex->SetTextSize(0.0315317);
   ptstats2->SetOptStat(1111);
   ptstats2->SetOptFit(10001);
   ptstats2->Draw();
   resp1->GetListOfFunctions()->Add(ptstats2);
   ptstats2->SetParent(resp1);

   resp3->SetTitle(stat_name);
   resp3->SetLineColor(4);
   resp3->SetLineWidth(2);
   //
   resp3->Scale(1.0/resp3->Integral());
   //
   resp3->GetXaxis()->SetTitle("Energy response");
   resp3->GetXaxis()->SetTitleSize(22);
   resp3->GetXaxis()->SetTitleFont(43);
   resp3->GetXaxis()->SetTitleOffset(0.8);  
   resp3->SetName(new_rel);

   //   resp3->Fit("gaus");
   TF1* ga1 = (TF1*)resp1->GetFunction("gaus");
   ga1->Delete();
      // cout<<"2 \n";
   TF1* ga2 = (TF1*)resp3->GetFunction("gaus");
   ga2->Delete();


   resp1->SetLineColor(kRed);
   resp1->SetLineWidth(2);
   //
   resp1->Scale(1.0/resp1->Integral());

   y1=resp1->GetMaximum();
   y2=resp3->GetMaximum();
   if (y1>y2)
     ymax=y1 + 0.01;
   else
     ymax=y2 + 0.01;

   cout<<"ymax    "<<ymax<<endl;

   resp3->GetXaxis()->SetRangeUser(xmin , xmax);                                                                                                                                                          
   resp3->GetYaxis()->SetRangeUser(ymin , ymax);

   resp1->GetXaxis()->SetRangeUser(xmin,xmax); 
   //
   resp1->GetYaxis()->SetRangeUser(ymin , ymax);                                                                                                                                                            
   resp1->SetName(old_rel);


   resp3->Draw();

     resp1->Draw("sames");



   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);


   Canvas_1_n2->SaveAs(full_path);
   Canvas_1_n2->SaveAs(full_path1);
   Canvas_1_n2->SaveAs(full_path2);
   if (n[i]==494) continue;
    }
}
