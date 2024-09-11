#include "TRatioPlot.h"
#define FMT_HEADER_ONLY
#include "fmt-10.2.0/include/fmt/format.h"
#include <sstream>
using namespace std;
// useful resources for color choices
//  https://personal.sron.nl/~pault/
int cols[] = {TColor::GetColor("#000000"), TColor::GetColor("#0077bb"),
              TColor::GetColor("#ee3377"), TColor::GetColor("#aa4499"),
               TColor::GetColor("#ee99aa"), TColor::GetColor("#009988"),
                TColor::GetColor("#228833"),TColor::GetColor("#228833")};

double fontsize = 0.05;
std::string componets[6] =  {"Total"};

std::string bins[11] = {"[0 < P_{t} < 0.075]",	"[0.075 < P_{t} < 0.15]",	"[0.15  < P_{t} < 0.25]",	"[0.25 < P_{t} < 0.325]",
                        "[0.325 < P_{t} <0 .4]",	"[0.4 < P_{t} < 0.475]",  "[0.465 < P_{t} < 0.55]", "[0.55 <  P_{t} < 0.7]" ,	
                        "[0.7 < P_{t} < 0.85]" , "[0.85 < P_{t} < 1.0]"	, " [1.0 < P_{t} < 2.5]"};

std::string bins_pz[11] = { "[1.5 < P_{z} < 3.5]", "[3.5 <  P_{z} < 4.5]" ,	
                        "[4.5 < P_{z} < 7]" , "[7.0 < P_{z} < 8.0]"	, " [8.0 < P_{z} < 10]", " [10.0 < P_{z} < 20]"};


void DrawStatisticalUncertainty(TH1 *mc, TLegend *leg , std::string name) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(cols[0], 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(cols[0], 0);
  mc->SetFillColorAlpha(cols[0], 0.1);
  mc->SetStats(0);
  leg->AddEntry(mc, name.c_str() , "f");
  mc->Draw("ehist");
}

void DrawMCBand(TH1 *mc, TLegend *leg , std::string name) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(cols[1], 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(cols[1], 0);
  mc->SetFillColorAlpha(cols[1], 0.5);
  mc->SetStats(0);
  leg->AddEntry(mc, name.c_str() , "f");
  mc->Draw("E2SAME");
}

/*
void DrawStatisticalUncertaintyBand(TH1 *mc, TLegend *leg) {
  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(TColor::GetColor("#000000"), 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(TColor::GetColor("#000000"), 0);
  mc->SetFillColorAlpha(TColor::GetColor("#000000"), 0.5);

  leg->AddEntry(mc, "1 #sigma (stat)", "f");
  mc->Draw("E2SAME");
}
*/
int gcol = TColor::GetColor("#bbbbbb");
void DrawNonQE(TH1 *mc, TLegend *leg) {

  mc->SetLineWidth(0);
  mc->SetLineColorAlpha(gcol, 0);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(gcol, 0);
  mc->SetFillColorAlpha(gcol, 0.5);

  leg->AddEntry(mc, "Non-CCQE Contribution", "f");
  mc->Draw("HISTSAME");
}

void DrawMCCV_total(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname, double binmax) {
  mc->SetLineWidth(4);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);
  mc->SetStats(0);
  
  
  mc->SetName(newhnistname.c_str());
  mc->SetTitle(mc->GetTitle()); 
  mc->GetYaxis()->SetTitle("cm^{2}/GeV^{2}");
  mc->GetYaxis()->SetRangeUser(0,1.1* binmax);
  mc->GetXaxis()->SetTitle("Available Energy (GeV)");
  leg->AddEntry(mc, mc->GetName() , "l");
  /**
  TPaveText *myText = new TPaveText(0.2,0.7,0.4,0.85, "NDC");
  myText->SetTextSize(0.04);
  myText->SetFillColor(0); //white background myText->SetTextAlign(12);
  myText->AddText(nevents.c_str()); 
  myText->Draw();*/
  //leg->AddEntry(mc, nevents.c_str() , "l");
  mc->Draw("HISTSAME");

  

}

void DrawMCCV_total_TPad(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname, double binmax, std::string whichptslice ) {
  mc->SetLineWidth(2);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);
  mc->SetStats(0);
  mc->SetTitleSize(4*fontsize);
  TLatex Tl;
  mc->SetName(newhnistname.c_str());
  mc->SetTitle(whichptslice.c_str()); 
  mc->GetYaxis()->SetTitle("cm^{2}/GeV^{2}");
  mc->GetYaxis()->SetRangeUser(0,1.1 * binmax);
  leg->AddEntry(mc, mc->GetName(), "l");
  mc->Draw("EHISTSAME");
}


void DrawMCCV(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname, double binmax) {
  mc->SetLineWidth(4);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);
  mc->SetStats(0);
  mc->SetName(newhnistname.c_str());
  mc->SetTitle(mc->GetTitle()); 
  mc->GetYaxis()->SetTitle("cm^{2}/GeV^{2}");
  mc->GetYaxis()->SetRangeUser(0,1.1 * binmax);
  mc->GetXaxis()->SetTitle("Available Energy (GeV)");
  leg->AddEntry(mc, mc->GetName(), "l");
  mc->Draw("HISTSAME");
}

void DrawMCCV_CCQE(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname, double binmax) {
  mc->SetLineWidth(4);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);
  mc->SetStats(0);
  mc->SetName(newhnistname.c_str());
  mc->SetTitle(mc->GetTitle()); 
  mc->GetYaxis()->SetTitle("cm^{2}/GeV^{2}");
  mc->GetYaxis()->SetRangeUser(0,1.1 * binmax);
  mc->GetXaxis()->SetTitle("Available Energy (GeV)");
  leg->AddEntry(mc, mc->GetName(), "l");
  mc->Draw("EHISTSAME");
}


void DrawMCCV_TPad(TH1 *mc, TLegend *leg, int col ,  std::string newhnistname, double binmax, std::string whichptslice) {
  mc->SetLineWidth(2);
  mc->SetLineColorAlpha(col, 1);
  mc->SetLineStyle(1);
  mc->SetMarkerSize(0);
  mc->SetMarkerColorAlpha(col, 0);
  mc->SetFillColorAlpha(col, 0);
  mc->SetStats(0);
  mc->SetName(newhnistname.c_str());
  mc->SetTitleSize(4*fontsize); 
  mc->SetTitle(whichptslice.c_str());
  mc->GetYaxis()->SetRangeUser(0,1.1 * binmax);
  leg->AddEntry(mc, mc->GetName(), "l");
  mc->Draw("EHISTSAME");
}

void plotting_EAvailable() {
  //TFile *fdatar = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/EAvailable_bparams0_160224_finerbinning.root");
  TFile *fdatar = TFile::Open("/root/software/nuisance_version2/nuisance/dune3minervalike_v2/dune3minervalike/NuWro_10e6_scaled.root");
  //NuWro_10e6_scaled.root");
  //18thmarch_comparingdifferentbinning.root");
    auto fdata_3D =fdatar->Get<TH3>("DUNE_3D_MINERvALike_EAvail"); //get the 2D histogram
    int number_of_pt_bins = fdata_3D->GetXaxis()->GetNbins();
    int number_of_pz_bins = fdata_3D->GetYaxis()->GetNbins();
    int number_of_EAvail_bins = fdata_3D->GetZaxis()->GetNbins();
    std::cout<< "number_of_pt_bins = " << number_of_pt_bins <<std::endl;
    std::cout<< "number_of_pz_bins = " << number_of_pz_bins <<std::endl;

     int i=0;
    for(int yi = 0 ; yi< number_of_pz_bins ; yi++){ //slices in pt
        //now draw canvas in minerva style 
        std::stringstream pzslice;
        pzslice << fdata_3D->GetYaxis()->GetBinLowEdge(yi + 1) << " < p_z < "
           << fdata_3D->GetYaxis()->GetBinUpEdge(yi + 1) << " [GeV/c]";
        //std::stringstream pzslice = fdata_3D->GetYaxis()->GetBinLowEdge(yi + 1) << " < p_z < " << fdata_3D->GetYaxis()->GetBinUpEdge(yi + 1) << " [GeV/c]";
        TCanvas *c2 = new TCanvas("c2", "", 60000, 40000 );
        
        
        //lat->Draw();  // drawn on canvas
        int numberofslicestoplot = 0;
        double binmax_intotalhist3D = fdata_3D->GetMaximum(); 
        
        c2->SetLeftMargin(0.01);
        c2->SetRightMargin(0.01);
        c2->SetTopMargin(0.8);
        c2->SetBottomMargin(0.01);
        c2->Divide(3,3,0.001,0.001);

        for(int xi =0 ; xi < number_of_pt_bins; xi++){

          std::stringstream ptslice;
          ptslice << fdata_3D->GetXaxis()->GetBinLowEdge(xi + 1) << " < p_t < "
          << fdata_3D->GetXaxis()->GetBinUpEdge(xi + 1) << " [GeV/c]";
    
    
            std::string ptslicename = ptslice.str().c_str();
            std::string dist = fmt::format( "_x{}_y{}", xi,yi);
        
            std::string distribution = dist;
                auto hname = "DUNE_3D_MINERvALike_EAvail" + distribution;
                //hname += dist.c_str() , "_x0_y1";
                //get hist
                auto fdata =
                    fdatar->Get<TH1>(hname.c_str());

                    if(!fdata){ std::cout << "Is the name: " << hname << " definitely correct? I couldn't find that hist in the file." << std::endl;}
                    if(!fdata) { std::cout << "Failed to read hist <bla>" << std::endl;
                      fdatar -> ls();
                      abort();
                    } 
                  
            fdata->SetName(("DUNE_3D_MINERvALike_EAvail"+distribution).c_str());
            double numberofevents_total = fdata->GetEntries();
           
            double binmax_intotalhist = fdata->GetMaximum(); 

          //std::string distribution = dist;
          //auto hname = "f3DHist_f3DHist_EAvail" + distribution;
          //hname += dist.c_str() , "_x0_y1";
          //get hist
          
          /*
          auto fdata =
              fdatar->Get<TH1>(hname.c_str());

              if(!fdata){ std::cout << "Is the name: " << hname << " definitely correct? I couldn't find that hist in the file." << std::endl;}
              if(!fdata) { std::cout << "Failed to read hist <bla>" << std::endl; } 
          fdata->SetName(("DUNE_3D_MINERvALike_EAvail"+distribution).c_str());
          //fdata->Scale(cross_section_toeventrate , "width");
          std::cout<< " fdata->Integral(width) === " <<  fdata->Integral("width") <<std::endl;
          
          double binmax_intotalhist = fdata->GetMaximum(); 
          
          std::cout << " bin maximum  =  " << binmax_intotalhist << std::endl; 
          */

          
          std::string hname2 = "DUNE_3D_MINERvALike_EAvail_CCQE"+ distribution ;
          if(!fdata){ std::cout << "Is the name: " << hname2 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
          auto fdatar_mc_CCQE = fdatar->Get<TH1>(hname2.c_str());
          fdatar_mc_CCQE->SetName(hname2.c_str());

        /*
        std::string hname2_plus1sigma = "DUNE_3D_MINERvALike_CCQE"+ distribution;
        if(!fdatar_plus1sigma){ std::cout << "Is the name: " << hname2_plus1sigma << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCQE_plus1sigma = fdatar_plus1sigma ->Get<TH1>(hname2_plus1sigma .c_str());
        fdatar_mc_CCQE_plus1sigma ->SetName(hname2_plus1sigma.c_str());

        std::string hname2_minus1sigma = "DUNE_3D_MINERvALike_EAvail_CCQE"+ distribution;
        if(!fdatar_minus1sigma){ std::cout << "Is the name: " << hname2_minus1sigma << " definitely correct? I couldn't find that hist in the file." << std::endl;}
        auto fdatar_mc_CCQE_minus1sigma = fdatar_minus1sigma ->Get<TH1>(hname2_minus1sigma .c_str());
        fdatar_mc_CCQE_minus1sigma ->SetName(hname2_minus1sigma.c_str());
        */
        //fdatar_mc_CCQE->Scale(cross_section_toeventrate , "width");

          std::string hname3 = "DUNE_3D_MINERvALike_EAvail_CC2p2h" + distribution;
          if(!fdata){ std::cout << "Is the name: " << hname3 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
          auto fdatar_mc_CC2p2h = fdatar->Get<TH1>(hname3.c_str());
          fdatar_mc_CC2p2h->SetName(hname3.c_str());
          //fdatar_mc_CC2p2h->Scale(cross_section_toeventrate  , "width");

          std::string hname4 = "DUNE_3D_MINERvALike_EAvail_CC1p1pi"  + distribution ;
          if(!fdata){ std::cout << "Is the name: " << hname4 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
          auto fdatar_mc_CC1pi = fdatar->Get<TH1>(hname4.c_str());
          fdatar_mc_CC1pi->SetName(hname4.c_str());
        // fdatar_mc_CC1pi->Scale(cross_section_toeventrate , "width");

          std::string hname5 = "DUNE_3D_MINERvALike_EAvail_CCDIS" + distribution ;
          if(!fdata){ std::cout << "Is the name: " << hname5 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
          auto fdatar_mc_CCDIS = fdatar->Get<TH1>(hname5.c_str());
          fdatar_mc_CCDIS->SetName(hname5.c_str());
        // fdatar_mc_CCDIS->Scale(cross_section_toeventrate, "width");
          
          std::string hname6 = "DUNE_3D_MINERvALike_EAvail_Other" +  distribution ;
          if(!fdata){ std::cout << "Is the name: " << hname6 << " definitely correct? I couldn't find that hist in the file." << std::endl;}
          auto fdatar_mc_CCOther = fdatar->Get<TH1>(hname6.c_str());
          fdatar_mc_CCOther->SetName(hname6.c_str());
          //fdatar_mc_CCOther->Scale(cross_section_toeventrate  , "width");

          /*
          std::string throw_hist = "error_bands/DUNE_3D_MINERvALike_CCQE" + distribution  + "_prof" ;
          auto fmc_throws = fthrows->Get<TH1>((throw_hist).c_str());
          std::cout << throw_hist << std::endl;
          fmc_throws->GetYaxis()->SetTitle(fmc_throws->GetYaxis()->GetTitle());
          //fmc_throws->Scale(cross_section_toeventrate   , "width");

          */
        
        /*
        auto fdatar_mc_GiBUU =
            fdatar->Get<TH1>((std::string("GiBUU_") + dist + "Plot").c_str());
        fdatar_mc_GiBUU->SetName((std::string("GiBUU_") + dist + "Plot").c_str());
        auto fdatar_mc_neut =
            fdatar->Get<TH1>((std::string("NEUT_") + dist + "Plot").c_str());
        fdatar_mc_neut->SetName((std::string("NEUT_") + dist + "Plot").c_str());
        */
    /*
        auto fmc_throws = fthrows->Get<TH1>(
            (std::string("error_bands/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MC")
                .c_str());
        fdata->GetYaxis()->SetTitle(fmc_throws->GetYaxis()->GetTitle());
        fmc_throws->Scale(1E38);
        auto fmc_cv = fthrows->Get<TH1>(
            (std::string("nominal_throw/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MC")
                .c_str());
        fmc_cv->Scale(1E38);
        auto fmc_cv_qe = fthrows->Get<TH1>(
            (std::string("nominal_throw/MicroBooNE_CC1Mu1p_XSec_1D") + dist +
            "_nu_MODES_CCQE;"+mc_name_cycle[fdata->GetName()])
                .c_str());
        fmc_cv_qe->Scale(1E38);

        auto fmc_cv_nonqe = dynamic_cast<TH1 *>(fmc_cv->Clone(
            (std::string("MicroBooNE_CC1Mu1p_XSec_1D") + dist + "_nu_nonQE")
                .c_str()));
        fmc_cv_nonqe->Add(fmc_cv_qe, -1.0);
        */
        TLegend *legendl = new TLegend(0.3, 0.60, 0.8, 0.85);
        legendl->SetTextFont(132);
        legendl->SetTextSize(fontsize * 0.95);
        legendl->SetBorderSize(0);
        legendl->SetFillStyle(0);
        legendl->SetNColumns(1);

        /*
        TPaveText *pt = new TPaveText(0.01,0.1,0.2,0.1);
        pt->SetTextSize(0.04);
        pt->SetFillColor(0);
        */
        TCanvas c1("c1", "", 1200, 1200);
        c1.SetLeftMargin(0.1);
        c1.SetRightMargin(0.1);
        c1.SetTopMargin(0.1);
        c1.SetBottomMargin(0.1);
        TPad pad1("pad1", "", 0, 0.3, 1, 1);
        pad1.Draw();
        pad1.cd();
        /*
        auto ratio_plot_plus = new TRatioPlot(fdatar_mc_CCQE_plus1sigma, fdatar_mc_CCQE);
        auto ratio_plot_minus = new TRatioPlot(fdatar_mc_CCQE_minus1sigma, fdatar_mc_CCQE);
          
        ratio_plot_plus->Draw();
        ratio_plot_plus->GetLowYaxis()->SetNdivisions(505);

        ratio_plot_minus->Draw();
        ratio_plot_minus->GetLowYaxis()->SetNdivisions(505);
        */
    

        DrawStatisticalUncertainty(fdata, legendl, "Statistical Uncertainty");
        DrawMCCV_total(fdata, legendl, cols[0] ,componets[0], binmax_intotalhist );
      // DrawMCCV(fdatar_mc_CCQE_minus1sigma, legendl, cols[6],componets[6],binmax_intotalhist);
        //DrawMCCV(fdatar_mc_CCQE_plus1sigma, legendl, cols[7],componets[7],binmax_intotalhist);
        DrawMCCV_CCQE(fdatar_mc_CCQE, legendl, cols[1],componets[1],binmax_intotalhist);
        DrawMCCV_CCQE(fdatar_mc_CC1pi, legendl, cols[2],componets[2],binmax_intotalhist);
        DrawMCCV_CCQE(fdatar_mc_CC2p2h, legendl, cols[3],componets[3],binmax_intotalhist);
        DrawMCCV_CCQE(fdatar_mc_CCDIS, legendl, cols[4],componets[4],binmax_intotalhist);

      // std::cout << "binmax_intotalhist " << binmax_intotalhist <<std::endl;
        //ratio_plot_plus->Draw();
        //ratio_plot_minus->Draw();
        //auto 3d_hist_total =  DUNE_3D_MINERvALike
        //double_t integral =  fdata->Integral("width") ;
      //std::cout<<"integral= " << integral<<std::endl;
        
        
      
        //DrawMCCV(fdatar_mc_CCOther, legendl, cols[5], componets[5]);
        pad1.Update();
      
        //DrawData(fdata, legendl, "SAMES");
        //std::cout<<bins[i].c_str()<<std::endl;
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05); 

        legendl->Draw();
        latex.SetTextSize(0.03); 
        latex.DrawLatexNDC(0.1,0.005,"#it{A Peake}");
        c1.Update();
        int numberofevents = fdatar_mc_CCQE->GetEntries();
        std::string nevents = "CCQE Events = " +  std::to_string(numberofevents);
        latex.DrawLatexNDC(0.1, 0.03,nevents.c_str());

        c1.Update();
        c1.Draw();
        c1.SaveAs((std::string("/root/plots/EAvailable/finer_binning/multiplots/Nuwro/newNuWroscaled_") + dist + "_.pdf").c_str());
        i = i+ 1;
        /*
        for(int i = 0; i < fdatar_mc_CCQE_plus1sigma->GetXaxis()->GetNbins(); ++i){ 
            std::cout << "bin " << i << " error = " << fdatar_mc_CCQE_plus1sigma->GetBinError(i+1) << std::endl;
            std::cout << "bin " << i << " content = " << fdatar_mc_CCQE_plus1sigma->GetBinContent(i+1) << std::endl;
            std::cout << "bin " << i << "  content / bin error = " << (fdatar_mc_CCQE_plus1sigma->GetBinContent(i+1)) /(fdatar_mc_CCQE_plus1sigma->GetBinError(i+1)) << std::endl;


        }
        for(int i = 0; i < fdatar_mc_CCQE->GetXaxis()->GetNbins(); ++i){ 
            std::cout << "bin " << i << " error = " << fdatar_mc_CCQE->GetBinError(i+1) << std::endl;
            std::cout << "bin " << i << " content = " << fdatar_mc_CCQE->GetBinContent(i+1) << std::endl;
            std::cout << "bin " << i << "  content / bin error = " << (fdatar_mc_CCQE->GetBinContent(i+1)) /(fdatar_mc_CCQE->GetBinError(i+1)) << std::endl;


        }
        */  
        //c2->cd(xi+1);
        //gPad->SetTickx(2);
        //gPad->SetTicky(2);
        
        if(numberofevents > 500){
          numberofslicestoplot = numberofslicestoplot +1;
          c2->cd(xi+1);
          //c2->cd();
        //DrawStatisticalUncertainty(fdata, legendl, "Statistical Uncertainty");
        DrawMCCV_total_TPad(fdata, legendl, cols[0] ,componets[0], binmax_intotalhist, ptslicename);
        DrawMCCV_TPad(fdatar_mc_CCQE, legendl, cols[1],componets[1],binmax_intotalhist,ptslicename);
        DrawMCCV_TPad(fdatar_mc_CC1pi, legendl, cols[2],componets[2],binmax_intotalhist,ptslicename);
        DrawMCCV_TPad(fdatar_mc_CC2p2h, legendl, cols[3],componets[3],binmax_intotalhist,ptslicename);
        DrawMCCV_TPad(fdatar_mc_CCDIS, legendl, cols[4],componets[4],binmax_intotalhist,ptslicename);
        
        //c2->Update();
        }
        else{
          std::cout<<"only " << numberofevents_total << "in this slice of pt" << std::endl;
        }
        
        }
        c2->Update();
        //TLegend *legendlc2 = new TLegend(0.3, 0.60, 0.8, 0.85);
        //
        c2->Draw();
        TLatex lat;
        TPad *pad5 = new TPad("all","all",0,0,1,1);
          pad5->SetFillStyle(4000);  // transparent
          pad5->Draw();
          pad5->cd();
   


        lat.SetTextSize(0.02); 
        //lat.SetTextColor("kRed"); 
        //TLatex *lat = new DrawLatex(.3,.95,pzslice.str().c_str());
        lat.DrawLatexNDC(0.78,.3,pzslice.str().c_str());

         
        //legendl2->Draw();
        c2->Update();
        c2->Draw();
        std::string s = (std::to_string(yi).c_str());
        c2->SaveAs((std::string("/root/plots/EAvailable/finer_binning/multiplots/newNuWro") + s + "_EAvailable_maxevents_errorband_scaled.pdf").c_str());

            
      }
}

