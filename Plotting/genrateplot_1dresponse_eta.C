#include <stdio.h>
// #include<conio.h>
void genrateplot_1dresponse_eta(TString legendname="EH hadrons", TString plot_name="corrEtaDependenceEH", TString path1="output/target_trueE/Trained_2_350_eta3_epoch100/", TString plot_name2="output_eta2pt75_EH_eta",TString region="trans_EH")
{
  float xmin1,xmin;
  float xmax1,xmax;
  gStyle->SetOptStat(0);
  TString hname,hname1,pdf,png,gif;
  bool xaxis_pt=false;
  TString path_root1="../../Run3/trial2/1d_response/";
  path_root1="../../Run3_v2/noPU_v2/noPU_v2_trial/1d_response/";
  path_root1="chi2_results/1d_response/";
  TString path_root2=path1;
  //Output_eta1_100k.root
  hname1= path_root2+plot_name2+".root";
  //  hname= path_root1+"projections_"+plot_name+".root"; 
  hname= path_root1+"projections_Full_"+plot_name+"_eta.root";

  cout<<"file name --> "<<hname<<endl;
  double arr[74],brr[74];
  arr[0]=0.00;
  brr[0]=0.00;
  for (int i=1;i<74;i++)
    {
      arr[i]=arr[i-1]+0.04;
      brr[i]=brr[i-1]+0.04;

    }
  for (int i=0;i<74;i++)
    //    cout<<arr[i]<<", ";
    cout<<i<<" : "<< arr[i]<<" , "<<brr[i]<<endl;
  cout<<endl;
  cout<<"=================="<<endl;

  TString v_arr[74]={"0","0.04","0.08","0.12","0.16","0.2","0.24","0.28","0.32","0.36","0.4","0.44","0.48","0.52","0.56","0.6","0.64","0.68","0.72","0.76","0.8","0.84","0.88","0.92","0.96","1","1.04","1.08","1.12","1.16","1.2","1.24","1.28","1.32","1.36","1.4","1.44","1.48","1.52","1.56","1.6","1.64","1.68","1.72","1.76","1.8","1.84","1.88","1.92","1.96","2","2.04","2.08","2.12","2.16","2.2","2.24","2.28","2.32","2.36","2.4","2.44","2.48","2.52","2.56","2.6","2.64","2.68","2.72","2.76","2.8","2.84","2.88","2.92"};
  TString v_arr2[74]={"0","04","08","12","16","2","24","28","32","36","4","44","48","52","56","6","64","68","72","76","8","84","88","92","96","1","04","08","12","16","2","24","28","32","36","4","44","48","52","56","6","64","68","72","76","8","84","88","92","96","2","04","08","12","16","2","24","28","32","36","4","44","48","52","56","6","64","68","72","76","8","84","88","92"};

  //  for(int i=1; i<44;i++)
    {
      //    if(legendname_=="EH_ec_out" && i<3) continue; 
      int i=35;
      //      int i=38;
      xmin=arr[i-1];
      xmax=arr[i];
      xmin1=brr[i-1];
      xmax1=brr[i];

	char* hist_name = new char[1000];
      char* stat_name = new char[200];
      char* hist_name1 = new char[1000];
      char* hist_name2 = new char[1000];
      char* hist_save = new char[1000];
      TString statname;
      //      TString hist_name,hist_name1,hist_name2,stat_name;
      
      //      TString hist_name="histcorhybrid";//+n[i];
      sprintf(hist_name1,"histEta_%f",xmin);
      
      sprintf(hist_name,"Valid_norm_pred_trueEn_eta_%f_to_%f",xmin1,xmax1);
      sprintf(hist_name2,"Train_norm_pred_trueEn_eta_%f_to_%f",xmin1,xmax1);

      if(i==44){
	sprintf(hist_name,"Valid_norm_pred_trueEn_eta_1.740000_to_1.780000");//,brr[i-1],brr[i]);
	sprintf(hist_name2,"Train_norm_pred_trueEn_eta_1.740000_to_1.780000");
      }
      statname= v_arr[i-1]+" < |#eta| < "+v_arr[i];      
      sprintf(stat_name,"eta =  {%f - %f}",xmin,xmax);
      sprintf(hist_save,"histEta_%f",xmin);
      cout<<"hist name --> "<<hist_name<<endl;
      pdf=path1+"/pdf/"+region+"/"+hist_save+".pdf";
      png=path1+"/pdf/"+region+"/"+hist_save+".png";
      if(i<25){
	pdf=path1+"/pdf/thesis/"+region+"_corrEta_0pt"+v_arr2[i-1]+".pdf";
	png=path1+"/pdf/thesis/"+region+"_corrEta_0pt"+v_arr2[i-1]+".png";
      }
      if(i<50){
        pdf=path1+"/pdf/thesis/"+region+"_corrEta_1pt"+v_arr2[i-1]+".pdf";
        png=path1+"/pdf/thesis/"+region+"_corrEta_1pt"+v_arr2[i-1]+".png";
      }
      else{
        pdf=path1+"/pdf/thesis/"+region+"_corrEta_2pt"+v_arr2[i-1]+".pdf";
        png=path1+"/pdf/thesis/"+region+"_corrEta_2pt"+v_arr2[i-1]+".png";
      }


      /*
      pdf="UL2018wrtrun3/plots/pdf/"+legendname_+"/"+hist_name+"_v2.pdf";
      png="UL2018wrtrun3/plots/png/"+legendname_+"/"+hist_name+"_v2.png";
      gif="UL2018wrtrun3/plots/gif/"+legendname_+"/"+hist_name+"_v2.gif";
      */
      TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2","Canvas_1_n2",600,600);

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
      
      TPaveStats *ptstats1 = new TPaveStats(0.60,0.7,0.9,0.85,"brNDC");
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
      /*
      ptstats1->Draw();
      graph3->GetListOfFunctions()->Add(ptstats1);
      ptstats1->SetParent(graph3);
      */
      TPaveStats *ptstats2 = new TPaveStats(0.60,0.55,0.9,0.7,"brNDC");
      ptstats2->SetBorderSize(1);
      ptstats2->SetFillColor(0);
      ptstats2->SetLineColor(kBlue);
      ptstats2->AddText("Raw");
      ptstats2->SetTextAlign(12);
      ptstats2->SetTextColor(kBlue);
      ptstats2->SetTextFont(42);
      ptstats2->SetTextSize(0.031);
      ptstats2->SetOptStat(1110);
      ptstats2->SetOptFit(00000);
      /*
      ptstats2->Draw();                                                                                                                                                   
      graph1->GetListOfFunctions()->Add(ptstats2);
      ptstats2->SetParent(graph1);
      */
      TPaveStats *ptstats3 = new TPaveStats(0.60,0.4,0.9,0.55,"brNDC");
      ptstats3->SetBorderSize(1);
      ptstats3->SetFillColor(0);
      ptstats3->SetLineColor(kRed);
      ptstats3->AddText("Raw");
      ptstats3->SetTextAlign(12);
      ptstats3->SetTextColor(kRed);
      ptstats3->SetTextFont(42);
      ptstats3->SetTextSize(0.031);
      ptstats3->SetOptStat(1110);
      ptstats3->SetOptFit(00000);
      /*
      ptstats3->Draw();
      graph2->GetListOfFunctions()->Add(ptstats3);
      ptstats3->SetParent(graph2);
      */
      //      graph2->Rebin(2);      

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

      graph3->GetXaxis()->SetRangeUser(-1 ,1);
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
      TLegend* legend1 = new TLegend(0.1, 0.85, 0.45, 0.9,"","brNDC");      
      legend1->SetFillColor(0);
      //      legend1->SetHeader(stat_name,"C");
      legend1->SetHeader(statname,"C");
      legend1->SetTextSize(0.035);
      legend1->Draw();

      TLegend* legend2 = new TLegend(0.55, 0.85, 0.9, 0.9,"","brNDC");
      legend2->SetFillColor(0);
      legend2->SetHeader(legendname,"C");
      legend2->SetTextSize(0.035);
      legend2->Draw();
      
      TLegend* legends = new TLegend(0.1, 0.7, 0.45, 0.85,"","brNDC"); // the numbers determine the position of the box 
      legends->SetFillColor(0); 
      //      legends->SetHeader(legendname,"C"); 
      legends->AddEntry(graph3,"using DRN (Training)","l");
      legends->AddEntry(graph1,"using DRN (Validation)","l");
      legends->AddEntry(graph2,"using Run3 calib (#chi^{2})","l");
      legends->SetTextSize(0.03);
      legends->Draw();

      
      Canvas_1_n2->SaveAs(pdf);
      Canvas_1_n2->SaveAs(png);
      //      Canvas_1_n2->SaveAs(gif);
      
    }
}

