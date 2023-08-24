#include <stdio.h>
// #include<conio.h>
void genrateplot_1dresponsevsE_UL2018wrtrun3(TString legendname, TString plot_name, TString path1, TString plot_name2,TString region, int xmin=20, int xmax=24)
{
  gStyle->SetOptStat(0);
  TString hname,hname1,pdf,png,gif;
  bool xaxis_pt=false;
  TString path_root1="../../Run3/trial2/1d_response/";
  path_root1="../../Run3_v2/noPU_v2/noPU_v2_trial/1d_response/";
  TString path_root2=path1;
  //Output_eta1_100k.root
  hname1= path_root2+plot_name2+".root";
  hname= path_root1+"projections_"+plot_name+".root";

  cout<<"file name --> "<<hname<<endl;

  int n[69];
  n[0]=0;
  n[1]=2;
  n[2]=8;
  n[3]=12;
  for(int i=4; i<69; i++)
    {
      if (n[i-1]<6)
        n[i]=n[i-1]+6;
      else if (n[i-1]>=12 && n[i-1]<104)
        n[i]=n[i-1]+4;
      else
        n[i]=n[i-1]+10;
    }
  for(int i=1; i<69; i++)
    cout<<i<<", "<<n[i]<<", "<<endl;

  //  for(int i=3; i<50;i++)
    {
      //    if(legendname_=="EH_ec_out" && i<3) continue; 
      int i=18;
      xmin=n[i-1];
      xmax=n[i];
	char* hist_name = new char[1000];
      char* stat_name = new char[200];
      char* hist_name1 = new char[1000];
      char* hist_name2 = new char[1000];
      char* hist_save = new char[1000];

      //      TString hist_name,hist_name1,hist_name2,stat_name;
      
      //      TString hist_name="histcorhybrid";//+n[i];
      sprintf(hist_name1,"histcorhybrid%d",xmin);
      if(xmax<10) sprintf(hist_save,"histcorhybrid_00_%d",xmax);
      else if(xmax<100) sprintf(hist_save,"histcorhybrid_01_%d",xmax);
      else sprintf(hist_save,"histcorhybrid_02_%d",xmax);
      //      sprintf(hist_name,"Valid_norm_pred_trueEn_TrueEn_%d_to_%d",xmin,xmax);
      //     sprintf(hist_name2,"Train_norm_pred_trueEn_TrueEn_%d_to_%d",xmin,xmax);
      if(region.Contains("barrel")){
	//	hist_name="Valid_norm_pred_eta_0.000000_1.550000_E_"+xmin+"_to"
	sprintf(hist_name,"Valid_norm_pred_eta_0.000000_1.550000_E_%d_%d",xmin,xmax);
	sprintf(hist_name2,"Train_norm_pred_eta_0.000000_1.550000_E_%d_%d",xmin,xmax);
      }
      else if(region.Contains("ec_in")){
	sprintf(hist_name,"Valid_norm_pred_eta_1.550000_2.500000_E_%d_%d",xmin,xmax);
        sprintf(hist_name2,"Train_norm_pred_eta_1.550000_2.500000_E_%d_%d",xmin,xmax);

      }
      else if(region.Contains("ec_out")){
	sprintf(hist_name,"Valid_norm_pred_eta_2.500000_2.750000_E_%d_%d",xmin,xmax);
        sprintf(hist_name2,"Train_norm_pred_eta_2.500000_2.750000_E_%d_%d",xmin,xmax);

      }
      else{
	sprintf(hist_name,"Valid_norm_pred_trueEn_TrueEn_%d_to_%d",xmin,xmax);
	sprintf(hist_name2,"Train_norm_pred_trueEn_TrueEn_%d_to_%d",xmin,xmax);
      }

      sprintf(stat_name,"E_{true} =  {%d GeV - %d GeV}",xmin,xmax);
      cout<<"hist name --> "<<hist_name<<endl;
      pdf=path1+"/pdf/"+region+"/"+hist_save+".pdf";
      png=path1+"/pdf/"+region+"/"+hist_save+".png";

      pdf=path1+"/pdf/thesis/"+region+"_corr_68.pdf";
      png=path1+"/pdf/thesis/"+region+"_corr_68.png";
      /*
      pdf="UL2018wrtrun3/plots/pdf/"+legendname_+"/"+hist_name+"_v2.pdf";
      png="UL2018wrtrun3/plots/png/"+legendname_+"/"+hist_name+"_v2.png";
      gif="UL2018wrtrun3/plots/gif/"+legendname_+"/"+hist_name+"_v2.gif";
      */
      TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",600,600);

      Canvas_1_n2->SetFillColor(0);
      Canvas_1_n2->SetBorderMode(0);
      Canvas_1_n2->SetBorderSize(2);
      Canvas_1_n2->SetFrameBorderMode(0);
      Canvas_1_n2->SetFrameBorderMode(0);
      TFile * inputfile1 = new TFile(hname,"READ");
      TFile * inputfile2 = new TFile(hname1,"READ");
      TH1F* graph1 = (TH1F*)inputfile2->Get(hist_name);
      TH1D* graph2 = (TH1D*)inputfile1->Get(hist_name1);
      TH1D* graph3 = (TH1D*)inputfile2->Get(hist_name2);

      //     TF1* ga1 = (TF1*)graph1->GetFunction("gaus");
      //      ga1->Delete();
      //TF1* ga2 = (TF1*)graph2->GetFunction("gaus");
      // ga2->Delete();
      //TF1* ga3 = (TF1*)graph3->GetFunction("gaus");
      //      ga3->Delete();
      
      TPaveStats *ptstats1 = new TPaveStats(0.60,0.67,0.9,0.85,"brNDC");
      ptstats1->SetBorderSize(1);
      ptstats1->SetFillColor(0);
      ptstats1->SetLineColor(kGreen);
      ptstats1->AddText("Raw");
      ptstats1->SetTextAlign(12);
      ptstats1->SetTextColor(kGreen);
      ptstats1->SetTextFont(42);
      ptstats1->SetTextSize(0.0315317);
     ptstats1->SetOptStat(1110);
     //      ptstats1->SetOptStat(0000);
      ptstats1->SetOptFit(00000);
      //ptstats1->Draw();
      //graph3->GetListOfFunctions()->Add(ptstats1);
      // ptstats1->SetParent(graph3);

      TPaveStats *ptstats2 = new TPaveStats(0.60,0.5,0.9,0.67,"brNDC");
      ptstats2->SetBorderSize(1);
      ptstats2->SetFillColor(0);
      ptstats2->SetLineColor(kBlue);
      ptstats2->AddText("Raw");
      ptstats2->SetTextAlign(12);
      ptstats2->SetTextColor(kBlue);
      ptstats2->SetTextFont(42);
      ptstats2->SetTextSize(0.0315317);
      ptstats2->SetOptStat(0000);
      ptstats2->SetOptStat(1110);
      ptstats2->SetOptFit(00000);
      //ptstats2->Draw();                                                                                                                                                   
      //graph1->GetListOfFunctions()->Add(ptstats2);
      //ptstats2->SetParent(graph1);

      
      TPaveStats *ptstats3 = new TPaveStats(0.60,0.3,0.9,0.5,"brNDC");
      ptstats3->SetBorderSize(1);
      ptstats3->SetFillColor(0);
      ptstats3->SetLineColor(kRed);
      ptstats3->SetTextAlign(12);
      ptstats3->SetTextColor(kRed);
      ptstats3->SetTextFont(42);
      ptstats3->SetTextSize(0.0315317);
      ptstats3->SetOptStat(1110);
      //ptstats3->SetOptStat(0000);
      ptstats3->SetOptFit(00000);
      //   ptstats3->Draw();                                                                                                                                                   
      // graph2->GetListOfFunctions()->Add(ptstats3);
      //ptstats3->SetParent(graph2);
      if(region.Contains("ec_out")){
	graph1->Rebin(2);
	graph2->Rebin(2);
	graph3->Rebin(2);
      }
      graph3->Scale(1.0/graph3->Integral());
      graph2->Scale(1.0/graph2->Integral());
      graph1->Scale(1.0/graph1->Integral());
      double ymax;
      double y1=graph1->GetMaximum();
      double y2=graph2->GetMaximum();
      if (y1>y2)
	ymax=y1 + 0.01;
      else
	ymax=y2 + 0.01;
      
      cout<<"ymax    "<<ymax<<endl;

      graph3->GetXaxis()->SetRangeUser(-2 , 2);
      graph3->GetYaxis()->SetRangeUser(0 , ymax);
      graph3->GetYaxis()->SetRangeUser(0 , ymax);

      graph3->SetLineColor(kGreen);
      graph3->SetLineWidth(2);
      graph3->SetMarkerStyle(22);
      graph3->SetMarkerSize(0.8);
       graph3->Draw("hist");
      
      graph1->SetLineWidth(2);
      graph1->SetLineColor(kBlue);
      graph1->SetMarkerStyle(22);
      graph1->SetMarkerSize(0.8);   
      
      graph1->Draw("hist sames");
      
      graph2->SetLineWidth(2);
      graph2->SetLineColor(kRed);
      graph2->SetMarkerStyle(22);
      graph2->SetMarkerSize(0.8);


      graph2->Draw("hist sames");

      gPad->SetGrid();
      TLegend* legend1 = new TLegend(0.1, 0.82, 0.48, 0.9,"","brNDC");      
      legend1->SetFillColor(0);
      legend1->SetHeader(stat_name,"C");
      legend1->SetTextSize(0.038);
      legend1->Draw();

      TLegend* legend2 = new TLegend(0.52, 0.82, 0.9, 0.9,"","brNDC");
      legend2->SetFillColor(0);
      legend2->SetHeader(legendname,"C");
      legend2->SetTextSize(0.035);
      legend2->Draw();
      
      TLegend* legends = new TLegend(0.1, 0.65, 0.465, 0.8,"","brNDC"); // the numbers determine the position of the box 
      legends->SetFillColor(0); 
      //      legends->SetHeader(legendname,"C"); 
      legends->AddEntry(graph3,"using DRN (Training)","l");
      legends->AddEntry(graph1,"using DRN (Validation)","l");
      legends->AddEntry(graph2,"using Run3 calib (#chi^{2})","l");
      legends->SetTextSize(0.032);
      legends->Draw();

      
      Canvas_1_n2->SaveAs(pdf);
      Canvas_1_n2->SaveAs(png);
      //      Canvas_1_n2->SaveAs(gif);
      
    }
}
