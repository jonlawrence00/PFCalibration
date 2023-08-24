#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TGraph.h>

void genrateplot_responsevsE_UL2018wrtrun3_v2(){
    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Three TPads in Canvas", 800, 600);

    // Divide the canvas into three TPads
    canvas->Divide(1, 3);

    // Create and draw a histogram in the first TPad
    canvas->cd(1);
    TH1F* histogram = new TH1F("histogram", "Example Histogram", 50, 0, 10);
    histogram->FillRandom("gaus", 1000);
    histogram->Draw();

    // Create and draw a TGraph in the second TPad
    canvas->cd(2);
    TGraph* graph = new TGraph();
    for (int i = 0; i < 10; ++i) {
        graph->SetPoint(i, i, i * i);
    }
    graph->SetMarkerStyle(20);
    graph->Draw("AP");

    // Create and draw some text in the third TPad
    canvas->cd(3);
    TPaveText* text = new TPaveText(0.1, 0.1, 0.9, 0.9);
    text->AddText("Example Text");
    text->AddText("in the third TPad");
    text->SetFillColor(18);
    text->SetTextAlign(22);
    text->SetTextSize(0.05);
    text->Draw();

    // Update the canvas to display the changes
    canvas->Update();

    // To save the canvas as an image file (optional)
    canvas->SaveAs("ThreeTPadsCanvas.png");
}
