#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TProfile.h"
#include "TF1.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
void getLowXHighX( TH1F * h1, double &xL, double &xH);
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);

double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {


	TApplication *myapp = new TApplication("myapp",0,0);

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	cout << "*************************************************************" << endl;
        cout << "                  ___                  ___" << endl;
        cout << "\\    /\\    / /\\   |  \\ |\\  | | |\\  |  /" << endl;
        cout << " \\  /  \\  / /--\\  |-\\/ | \\ | | | \\ | |  --\\" << endl;
        cout << "  \\/    \\/ /    \\ |  \\ |  \\| | |  \\|  \\___/" << endl << endl;
        cout << "Run this code on cosmic data" << endl << endl;
        cout << "*************************************************************" << endl;

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;
	double c = 29.9792;

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	TH1F ** h1_tdc_diff = new TH1F * [nHistos];
	TH1F ** h1_ftdc_diff = new TH1F * [nHistos];
	TH2F ** h2_empty = new TH2F * [nHistos];
	TH2F ** h2_meantime = new TH2F * [nHistos];
	for( int i = 0 ; i < nHistos ; i++ ){
		h1_tdc_diff [i] = new TH1F(Form("h1_tdc_diff_%i" ,i),"",400,-20,20);
		h1_ftdc_diff[i] = new TH1F(Form("h1_ftdc_diff_%i",i),"",400,-20,20);
		h2_empty    [i] = new TH2F(Form("h2_empty_%i"    ,i),";TDC L-R;FADC L-R",400,-20,20,400,-20,20);
		h2_meantime [i] = new TH2F(Form("h2_meantime_%i" ,i),";#sqrt{ADC_{L}*ADC_{R}};#frac{L+R}{2} TDC - FADC",500,0,25000,1000,1045,1075);
	
		PrettyTH2F(h2_empty    [i],"","");
		PrettyTH2F(h2_meantime [i],"","");
	}

	double par_pad[nHistos][2] = {0};
	double val_slope = 0;
	int ctr_slope = 0;
	TH1D * h1_slopes = new TH1D("h1_slopes",";slope",100,0.5,1.5);

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	//Read Dictionary of Hipo File  // new hipo4
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
        //factory.show();               // new hipo4

	BBand    band_hits   (factory.getSchema("BAND::hits" ));
	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));
        hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));

	//One also needs a hipo::event object which is called from the reader for each event to get
        //the information for each bank
        hipo::event readevent;  // new hipo4

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		//Reader has to load information about event in hipo::event class
                reader.read(readevent); // new hipo4

                //Load explicitly all information for each bank for the event
                readevent.getStructure(band_hits);	// new hipo4
		readevent.getStructure(BAND_ADC );	// new hipo4
		readevent.getStructure(BAND_TDC );	// new hipo4

                //Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
                //The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

                if(event_counter%1000==0) cout << "event: " << event_counter << endl;
                event_counter++;

                //BAND_ADC.show();
                //BAND_TDC.show();
                //RUN_config.show();

                int nADC = BAND_ADC.getRows();
                int nTDC = BAND_TDC.getRows();

                // Skip events with no entries
                if(nADC==0||nTDC==0) continue;

                //long timestamp = RUN_config.getLong(4,0);
                //double phaseCorr = getTriggerPhase(timestamp);

                for(int aIdx1 = 0 ; aIdx1 < nADC ; aIdx1++){

                        int   ADC1_sector    = BAND_ADC.getInt  (0,aIdx1);
                        int   ADC1_layer     = BAND_ADC.getInt  (1,aIdx1);
                        int   ADC1_component = BAND_ADC.getInt  (2,aIdx1);
                        int   ADC1_order     = BAND_ADC.getInt  (3,aIdx1);
                        float ADC1_adc       = (float)(BAND_ADC.getInt(4,aIdx1));
                        float ADC1_time      = BAND_ADC.getFloat(5,aIdx1);

                        for(int aIdx2 = 0 ; aIdx2 < nADC ; aIdx2++){

                                int   ADC2_sector    = BAND_ADC.getInt  (0,aIdx2);
                                int   ADC2_layer     = BAND_ADC.getInt  (1,aIdx2);
                                int   ADC2_component = BAND_ADC.getInt  (2,aIdx2);
                                int   ADC2_order     = BAND_ADC.getInt  (3,aIdx2);
                                float ADC2_adc       = (float)(BAND_ADC.getInt(4,aIdx2));
                                float ADC2_time      = BAND_ADC.getFloat(5,aIdx2);

                                // Matching each ADC to the corresponding ADC from the other side of the bar
                                if(             (ADC1_sector   ==ADC2_sector          )&&
                                                (ADC1_layer    ==ADC2_layer           )&&
                                                (ADC1_component==ADC2_component       )&&
                                                (ADC1_order+1==ADC2_order             )
                                  ){

                                        bool already_matched = false;
                                        for(int tIdx1 = 0 ; tIdx1 < nTDC ; tIdx1++){

                                                int   TDC1_sector    = BAND_TDC.getInt  (0,tIdx1);
                                                int   TDC1_layer     = BAND_TDC.getInt  (1,tIdx1);
                                                int   TDC1_component = BAND_TDC.getInt  (2,tIdx1);
                                                int   TDC1_order     = BAND_TDC.getInt  (3,tIdx1);
                                                float TDC1_tdc       = (float)(BAND_TDC.getInt(4,tIdx1));
                                                float TDC1_time      = TDC1_tdc*0.02345;

                                                //TDC1_time -= phaseCorr;

                                                // Matching these ADCs to a TDC
                                                if(             (ADC1_sector   ==TDC1_sector   )&&
                                                                (ADC1_layer    ==TDC1_layer    )&&
                                                                (ADC1_component==TDC1_component)&&
                                                                (ADC1_order+2  ==TDC1_order    )
                                                  ){

                                                        for(int tIdx2 = 0 ; tIdx2 < nTDC ; tIdx2++){

                                                                int   TDC2_sector    = BAND_TDC.getInt  (0,tIdx2);
                                                                int   TDC2_layer     = BAND_TDC.getInt  (1,tIdx2);
                                                                int   TDC2_component = BAND_TDC.getInt  (2,tIdx2);
                                                                int   TDC2_order     = BAND_TDC.getInt  (3,tIdx2);
                                                                float TDC2_tdc       = (float)(BAND_TDC.getInt(4,tIdx2));
                                                                float TDC2_time      = TDC2_tdc*0.02345;

                                                                //TDC2_time -= phaseCorr;

                                                                if(             (!already_matched                     )&& // Avoid multi-hits
                                                                                (TDC1_sector   ==TDC2_sector          )&&
                                                                                (TDC1_layer    ==TDC2_layer           )&&
                                                                                (TDC1_component==TDC2_component       )&&
                                                                                (TDC1_order+1==TDC2_order             )
                                                                  ){

                                                                        already_matched = true;

                                                                        int barKey = 100*ADC1_sector + 10*ADC1_layer + ADC1_component;

									if( ADC1_adc <  6000 || ADC2_adc <  6000 ) continue;
                        						if( ADC1_adc > 20000 || ADC2_adc > 20000 ) continue;

									double difftimeTdc  = TDC1_time-TDC2_time;
									double difftimeFadc = ADC1_time-ADC2_time;
									double meantimeTdc  = (TDC1_time+TDC2_time)/2.;
									double meantimeFadc = (ADC1_time+ADC2_time)/2.;

                        						h1_tdc_diff [barKey]->Fill( difftimeTdc );
                       							h1_ftdc_diff[barKey]->Fill( difftimeFadc );
                        						h2_empty    [barKey]->Fill( difftimeTdc,difftimeFadc);
                        						h2_meantime [barKey]->Fill( TMath::Sqrt(ADC1_adc*ADC2_adc) ,(meantimeTdc-meantimeFadc));

                                                                }
                                                        }
                                                }
                                        }
                                }
                        }

                }


	}// end file

	// Get the width for each histogram:
	ofstream lr_offsets,effective_velocity;
	lr_offsets.open("TEST_lr_offsets.txt");
	effective_velocity.open("TEST_effective_velocity.txt");
	for( int sector = 1 ; sector < 6 ; sector++){
		for( int layer = 1 ; layer < 7 ; layer++ ){
			for( int component = 1 ; component <= slc[layer-1][sector-1] ; component++ ){
				int i = sector*100 + layer*10 + component;

				double tdc_vEff,tdc_lr_off,ftdc_vEff,ftdc_lr_off;
				if( h1_tdc_diff[i]->Integral() == 0 ){
					tdc_vEff = 0.;
					tdc_lr_off = 0.;
					ftdc_vEff = 0.;
					ftdc_lr_off = 0.;
				}
				else{
					double barLength = bandlen[sector-1];

					// For FADC
					int ftdc_low_Bin = h1_ftdc_diff[i]->FindFirstBinAbove(h1_tdc_diff[i]->GetMaximum()*0.1);
					int ftdc_high_Bin = h1_ftdc_diff[i]->FindLastBinAbove(h1_tdc_diff[i]->GetMaximum()*0.1);
					double ftdc_left = h1_ftdc_diff[i]->GetXaxis()->GetBinCenter(ftdc_low_Bin);
					double ftdc_right = h1_ftdc_diff[i]->GetXaxis()->GetBinCenter(ftdc_high_Bin);
					double ftdc_width = ftdc_right - ftdc_left;
					ftdc_vEff = (2*barLength) / ftdc_width;
					ftdc_lr_off = (ftdc_left+ftdc_right)/2.;

					// For TDC
					int tdc_low_Bin = h1_tdc_diff[i]->FindFirstBinAbove(h1_tdc_diff[i]->GetMaximum()*0.1);
					int tdc_high_Bin = h1_tdc_diff[i]->FindLastBinAbove(h1_tdc_diff[i]->GetMaximum()*0.1);
					double tdc_left = h1_tdc_diff[i]->GetXaxis()->GetBinCenter(tdc_low_Bin);
					double tdc_right = h1_tdc_diff[i]->GetXaxis()->GetBinCenter(tdc_high_Bin);
					double tdc_width = tdc_right - tdc_left;
					tdc_vEff = (2*barLength) / tdc_width;
					tdc_lr_off = (tdc_left+tdc_right)/2.;

					if( isinf(tdc_vEff) ) tdc_vEff = 0.;
					if( isinf(ftdc_vEff) ) ftdc_vEff = 0.;
				}

				effective_velocity << sector << "\t" << layer << "\t" << component << "\t" << tdc_vEff << "\t" << ftdc_vEff << "\t" << 0.000 << "\t" << 0.000 << "\n";
				lr_offsets << sector << "\t" << layer << "\t" << component << "\t" << tdc_lr_off << "\t" << ftdc_lr_off << "\t" << 0.000 << "\t" << 0.000 << "\n";
			}
		}
	}
	lr_offsets.close();
	effective_velocity.close();

	// --------------------------------------------------------------------------------------------------------
	// Fitting lines to extract the slope and try to determine the needed TDC -> time conversion factor
        TProfile ** p_empty = new TProfile*[nHistos];
        // Fitting paddle-to-paddle parameters
        for(int i = 0 ; i < nHistos ; i++){
                int notEmpty = (h2_empty[i]->Integral());
                if(notEmpty){
                        int nBins = h2_empty[i] -> GetXaxis() -> GetNbins();
                        p_empty[i] = h2_empty[i]->ProfileX(Form("px_%i",i), 1, nBins );
                        TF1 * f_line = new TF1("f_line","pol1");
                        p_empty[i] -> Fit("f_line","Q");
                        par_pad[i][0] = f_line->GetParameter(1);
                        par_pad[i][1] = f_line->GetParError (1);
                
			val_slope += par_pad[i][0];
        		ctr_slope++;

			h1_slopes -> Fill(par_pad[i][0]);
		}
        }

	double avg_m = val_slope/((double)(ctr_slope));
	double std_TDC2ns  = 0.02345;
	double std_FADC2ns = 0.0625;
	cout << "Average slope = " << avg_m << endl;
	cout << "TDC  to ns conversion factor (assuming FADC t conversion factor is " << std_FADC2ns << "): " << avg_m*std_TDC2ns  << endl;
	cout << "FADC to ns conversion factor (assuming TDC  t conversion factor is " << std_TDC2ns  << "): " << std_FADC2ns/avg_m << endl;

	TCanvas * c1 = new TCanvas("c1","c1");
        h1_slopes -> Draw();
        c1 -> Modified();
        c1 -> Update();

	// --------------------------------------------------------------------------------------------------------
	// Create plots for all the L-R distributions of TDC and FADC
	TCanvas *** c_tdc_diff  = new TCanvas**[5];
	TCanvas *** c_empty     = new TCanvas**[5];
	TCanvas *** c_meantime  = new TCanvas**[5];

	int min, max;

	for(int is = 0 ; is < 5 ; is++){
		c_tdc_diff [is] = new TCanvas*[6];
		c_empty    [is] = new TCanvas*[6];
		c_meantime [is] = new TCanvas*[6];

		for(int il = 0 ; il < 5 ; il++){
			c_tdc_diff [is][il] = new TCanvas(Form("c_tdc_diff_S%iL%i" ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_tdc_diff [is][il] -> Divide(2,7);
			c_empty    [is][il] = new TCanvas(Form("c_empty_S%iL%i"    ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_empty    [is][il] -> Divide(2,4);
			c_meantime [is][il] = new TCanvas(Form("c_meantime_S%iL%i" ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_meantime [is][il] -> Divide(2,4);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h1_tdc_diff[identifier]->Integral()+h1_ftdc_diff[identifier]->Integral());
				if(notEmpty){
					double low_x,high_x;
					int low_Bin,high_Bin;
					double tdc_len,fadc_len;
					c_tdc_diff[is][il] -> cd(2*cIdx+1);
					gPad -> SetBottomMargin(0.26);
					// Draw TDC
					PrettyTH1F( h1_tdc_diff[identifier] , "TDC L-R [ns]","", 2);
					h1_tdc_diff[identifier]->Draw();
					// make two lines based on TDC height of 0.1 of max
					low_Bin  = h1_tdc_diff[identifier]->FindFirstBinAbove(h1_tdc_diff[identifier]->GetMaximum()*0.1);
					high_Bin = h1_tdc_diff[identifier]->FindLastBinAbove(h1_tdc_diff[identifier]->GetMaximum()*0.1);
					low_x = h1_tdc_diff[identifier]->GetXaxis()->GetBinCenter(low_Bin);
					high_x = h1_tdc_diff[identifier]->GetXaxis()->GetBinCenter(high_Bin);
					tdc_len = high_x-low_x;
					TLine * lLow = new TLine(low_x,0,low_x,h1_tdc_diff[identifier]->GetMaximum());
					TLine * lHigh = new TLine(high_x,0,high_x,h1_tdc_diff[identifier]->GetMaximum());
					TLine * across = new TLine(low_x,h1_tdc_diff[identifier]->GetMaximum()/2.,high_x,h1_tdc_diff[identifier]->GetMaximum()/2.);
					lLow->SetLineColor(4);
					lHigh->SetLineColor(4);
					across->SetLineColor(4);
					lLow  -> Draw("same");
					lHigh -> Draw("same");
					across->Draw("Same");

					// --------------------------------------------------------------
					c_tdc_diff[is][il] -> cd(2*cIdx+2);
					gPad -> SetBottomMargin(0.26);
					// Draw FADC
					PrettyTH1F( h1_ftdc_diff[identifier] , "FADC L-R [ns]","", 1);
					h1_ftdc_diff[identifier]->Draw();
					// make two lines
					low_Bin  = h1_ftdc_diff[identifier]->FindFirstBinAbove(h1_tdc_diff[identifier]->GetMaximum()*0.1);
					high_Bin = h1_ftdc_diff[identifier]->FindLastBinAbove(h1_tdc_diff[identifier]->GetMaximum()*0.1);
					low_x = h1_ftdc_diff[identifier]->GetXaxis()->GetBinCenter(low_Bin);
					high_x = h1_ftdc_diff[identifier]->GetXaxis()->GetBinCenter(high_Bin);
					fadc_len = high_x - low_x;
					TLine * lLow2 = new TLine(low_x,0,low_x,h1_ftdc_diff[identifier]->GetMaximum());
					TLine * lHigh2 = new TLine(high_x,0,high_x,h1_ftdc_diff[identifier]->GetMaximum());
					TLine * across2 = new TLine(low_x,h1_ftdc_diff[identifier]->GetMaximum()/2.,high_x,h1_ftdc_diff[identifier]->GetMaximum()/2.);

					lLow2->SetLineColor(8);
					lHigh2->SetLineColor(8);
					across2->SetLineColor(8);
					lLow2  -> Draw("same");
					lHigh2 -> Draw("same");
					across2->Draw("same");

					// --------------------------------------------------------------
					c_empty[is][il] -> cd(cIdx+1);
					gPad -> SetTopMargin   (0.00);
                                        gPad -> SetBottomMargin(0.26);
                                        gPad -> SetLeftMargin  (0.26);	
					h2_empty[identifier]->SetTitle(Form("FADC Length: %f / TDC Length: %f",fadc_len,tdc_len));
					h2_empty[identifier]->SetStats(0);
					h2_empty[identifier]->Draw("col");
					p_empty [identifier]->Draw("same");
					TLine * slope1 = new TLine(-20,-20,20,20);
					slope1->SetLineColor(1);
					slope1->Draw("same");

					// --------------------------------------------------------------
					c_meantime[is][il] -> cd(cIdx+1);
					gPad -> SetTopMargin   (0.00);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin  (0.26);
					h2_meantime[identifier]->SetStats(0);
					min = getMinNonEmptyBin(h2_meantime[identifier]);
					max = getMaxNonEmptyBin(h2_meantime[identifier]);
					h2_meantime[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_meantime[identifier] -> Draw("col");
					//h2_empty[identifier]->SetMaximum(h1_ftdc_diff[identifier]->GetMaximum());
					//h2_empty[identifier]->Draw();
					//lLow  -> Draw("same");
					//lHigh -> Draw("same");
					//lLow2  -> Draw("same");
					//lHigh2 -> Draw("same");
					//across->Draw("same");
					//across2->Draw("same");

				}
			}		
			c_tdc_diff [is][il] -> Modified();	c_tdc_diff [is][il] -> Update();	
			c_empty    [is][il] -> Modified();	c_empty    [is][il] -> Update();
			c_meantime [is][il] -> Modified();	c_meantime [is][il] -> Update();
		}
	}

	// Saving to pdf
	TCanvas * c0 = new TCanvas("c0","",900,900);
	c0 -> Print("results_offsets_veff_rawBanks.pdf(");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_tdc_diff[is][il] -> Print("results_offsets_veff_rawBanks.pdf");
		}
	}
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_empty[is][il] -> Print("results_offsets_veff_rawBanks.pdf");
		}
	}
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_meantime[is][il] -> Print("results_offsets_veff_rawBanks.pdf");
		}
	}
	c0 -> Print("results_offsets_veff_rawBanks.pdf)");

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetYaxis() -> CenterTitle();

	h1 -> GetXaxis() -> SetTitleSize(0.07);
	h1 -> GetXaxis() -> SetLabelSize(0.07);
	h1 -> GetXaxis() -> SetNdivisions(107);

	h1 -> GetYaxis() -> SetTitleSize(0.07);
	h1 -> GetYaxis() -> SetLabelSize(0.07);
	h1 -> GetYaxis() -> SetNdivisions(107);
	h1 -> GetYaxis() -> SetTitleOffset(0.63);
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2,TString titx,TString tity) {
	//h2 -> GetXaxis() -> SetTitle(titx);
	//h2 -> GetYaxis() -> SetTitle(tity);
	h2 -> GetXaxis() -> SetTitleSize(0.10);
        h2 -> GetYaxis() -> SetTitleSize(0.10);
	h2 -> GetXaxis() -> SetLabelSize(0.10);
        h2 -> GetYaxis() -> SetLabelSize(0.10);
	h2 -> GetXaxis() -> CenterTitle();
        h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetNdivisions(105);
        h2 -> GetYaxis() -> SetNdivisions(107);
}
// ========================================================================================================================================
void getLowXHighX( TH1F * h1, double &xL, double &xH){
	int low_Bin  = h1->FindFirstBinAbove(h1->GetMaximum()*0.1);
	int high_Bin = h1->FindLastBinAbove(h1->GetMaximum()*0.1);
	xL = h1->GetXaxis()->GetBinCenter(low_Bin);
	xH = h1->GetXaxis()->GetBinCenter(high_Bin);

}
// ========================================================================================================================================
int getMinNonEmptyBin(TH2F * h2){
        int nBinX = h2 -> GetXaxis() -> GetNbins();
        int nBinY = h2 -> GetYaxis() -> GetNbins();

        int binContent = 0;
	int nextBinContent = 0;

        for(int i = 1 ; i <= nBinY ; i++){
                for(int j = 1 ; j <= nBinX ; j++){
                        binContent     = h2 -> GetBinContent(j  ,i);
                        nextBinContent = h2 -> GetBinContent(j+1,i);
			if(binContent!=0&&nextBinContent!=0) return i;
			else{
				binContent     = 0;
				nextBinContent = 0;
			}
                }
        }
        return 1;
}
// ========================================================================================================================================
int getMaxNonEmptyBin(TH2F * h2){
        int nBinX = h2 -> GetXaxis() -> GetNbins();
        int nBinY = h2 -> GetYaxis() -> GetNbins();

        int binContent = 0;
	int nextBinContent = 0;

        for(int i = nBinY ; i > 0 ; i--){
                for(int j = 1 ; j <= nBinX ; j++){
                        binContent = h2 -> GetBinContent(j,i);
                        nextBinContent = h2 -> GetBinContent(j-1,i);
			if(binContent!=0) return i;
			 else{
                                binContent     = 0;
                                nextBinContent = 0;
                        }
                }
        }
        return nBinY;
}
// ========================================================================================================================================
