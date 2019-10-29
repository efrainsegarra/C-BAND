#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Forward-declaring functions
void PrettyTH2F(TH2F * h2);
double getTriggerPhase( long timeStamp );
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetTitleSize(0.3,"t");
	gStyle->SetOptStat(0);

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double parL[nHistos][4] = {{0}};
	double parR[nHistos][4] = {{0}};

	TH2F ** h2_tdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R = new TH2F * [nHistos];

	TH2F * h2_tdc_adc_tot_L = new TH2F("h2_tdc_adc_tot_L","Left PMTs;ADC;t_{TDC}-t_{FADC} [ns]" ,400,1000,14000,400,382,396);
	TH2F * h2_tdc_adc_tot_R = new TH2F("h2_tdc_adc_tot_R","Right PMTs;ADC;t_{TDC}-t_{FADC} [ns]",400,1000,14000,400,382,396);
	h2_tdc_adc_tot_L -> GetYaxis() -> SetTitleOffset(1.40);
	h2_tdc_adc_tot_L -> GetXaxis() -> SetNdivisions(509);
	h2_tdc_adc_tot_R -> GetYaxis() -> SetTitleOffset(1.40);
	h2_tdc_adc_tot_R -> GetXaxis() -> SetNdivisions(509);

	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC;t_{TDC}-t_{FADC} [ns]",400,1000,14000,400,382,396);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC;t_{TDC}-t_{FADC} [ns]",400,1000,14000,400,382,396);

		PrettyTH2F(h2_tdc_adc_L[i]);
		PrettyTH2F(h2_tdc_adc_R[i]);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	//Read Dictionary of Hipo File  // new hipo4
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
        //factory.show();               // new hipo4

	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));	// new hipo4
	hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));	// new hipo4
	//hipo::bank RUN_config(factory.getSchema("RUN::config"));	// new hipo4

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
                readevent.getStructure(BAND_ADC   );   // new hipo4
                readevent.getStructure(BAND_TDC   );   // new hipo4
                //readevent.getStructure(RUN_config );   // new hipo4 

                //Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
                //The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		//BAND_ADC.show();
		//BAND_TDC.show();
		//RUN_config.show();

		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();

		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;

		// Skip events with less than 200 entries (the laser lights all bars at the same time)
		if(nADC<200||nTDC<200) continue;

		//long timestamp = RUN_config.getLong(4,0);
		//double phaseCorr = getTriggerPhase(timestamp);
		double phaseCorr = 0;
		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){

			int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
			int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
			int   ADC_component = BAND_ADC.getInt  (2,aIdx);
			int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			float ADC_adc       = (float)(BAND_ADC.getInt(4,aIdx));
			float ADC_time      = BAND_ADC.getFloat(5,aIdx);

			for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
				int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
				int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
				int   TDC_component = BAND_TDC.getInt  (2,tIdx);
				int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
				float TDC_tdc       = (float)(BAND_TDC.getInt(4,tIdx));
				float TDC_time      = TDC_tdc*0.02345;

				if(		(ADC_sector   ==TDC_sector   )&&
						(ADC_layer    ==TDC_layer    )&&
						(ADC_component==TDC_component)&&
						(ADC_order +2 ==TDC_order    )
				  ){
					int barKey = 100*ADC_sector + 10*ADC_layer + ADC_component;

					double deltaT_phaseCorr = (TDC_time - phaseCorr) - ADC_time;

					// Left PMTs
					if(ADC_order==0){
						h2_tdc_adc_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_L[barKey] -> Fill(ADC_adc,deltaT_phaseCorr);
						h2_tdc_adc_tot_L     -> Fill(ADC_adc,deltaT_phaseCorr);
					}

					// Right PMTs
					if(ADC_order==1){
						h2_tdc_adc_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_R[barKey] -> Fill(ADC_adc,deltaT_phaseCorr);
						h2_tdc_adc_tot_R     -> Fill(ADC_adc,deltaT_phaseCorr);
					}

				}

			}
		}

	}// end file

	// -------------------------------------------------------------------------------------------------
	// Fitting functions
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_tdc_adc_L[i]->Integral()+h2_tdc_adc_R[i]->Integral());
		if(notEmpty){
			TF1 * f_timewalk_L = new TF1("f_timewalk_L","[0]/TMath::Sqrt(x)+[1]");
			TF1 * f_timewalk_R = new TF1("f_timewalk_R","[0]/TMath::Sqrt(x)+[1]");

			h2_tdc_adc_L[i] -> Fit("f_timewalk_L","Q");
			h2_tdc_adc_R[i] -> Fit("f_timewalk_R","Q");

			parL[i][0] = f_timewalk_L -> GetParameter(0);	//A
			parL[i][1] = f_timewalk_L -> GetParameter(1);	//B
			parL[i][2] = f_timewalk_L -> GetParError (0);	//eA
			parL[i][3] = f_timewalk_L -> GetParError (1);	//eB

			parR[i][0] = f_timewalk_R -> GetParameter(0);   //A
			parR[i][1] = f_timewalk_R -> GetParameter(1);   //B
			parR[i][2] = f_timewalk_R -> GetParError (0);   //eA
			parR[i][3] = f_timewalk_R -> GetParError (1);   //eB
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0");
	c0 -> Divide(2,1);
	c0 -> cd(1);	h2_tdc_adc_tot_L -> Draw("COLZ");
	c0 -> cd(2);	h2_tdc_adc_tot_R -> Draw("COLZ");
	c0 -> Modified();
	c0 -> Update();

	TCanvas *** cSLC = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC[is] = new TCanvas*[6];
		for(int il = 0 ; il < 5 ; il++){
			cSLC[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_tdc_adc_L[identifier]->Integral()+h2_tdc_adc_R[identifier]->Integral());
				if(notEmpty){
					int min, max;

					cSLC[is][il] -> cd(2*cIdx+1);
					gPad -> SetBottomMargin(0.26);
					min = getMinNonEmptyBin(h2_tdc_adc_L[identifier]);
					max = getMaxNonEmptyBin(h2_tdc_adc_L[identifier]);
					h2_tdc_adc_L[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_tdc_adc_L[identifier] -> Draw("COLZ");

					cSLC[is][il] -> cd(2*cIdx+2);
					gPad -> SetBottomMargin(0.26);
					min = getMinNonEmptyBin(h2_tdc_adc_R[identifier]);
                                        max = getMaxNonEmptyBin(h2_tdc_adc_R[identifier]);
                                        h2_tdc_adc_R[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_tdc_adc_R[identifier] -> Draw("COLZ");
				}
			}
			cSLC[is][il] -> Modified();	cSLC[is][il] -> Update();
		}
	}
	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	ofstream tabL, tabR;
	tabL.open("timeWalkPar_L.txt");
	tabR.open("timeWalkPar_R.txt");

	for(int is = 1 ; is <= 5 ; is++){	
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;

				//if(parL[idx][0]!=0){
					tabL << is << "\t" << il << "\t" << ic << "\t" << "\t";
					if(il==6) tabL << "0\t0\t0\t0" << endl;
					else tabL << parL[idx][0] << "\t" << parL[idx][1] << "\t" << parL[idx][2]  << "\t" << parL[idx][3] << endl;
				//}
				// ---
				//if(parR[idx][0]!=0){
					tabR << is << "\t" << il << "\t" << ic << "\t" << "\t";
					if(il==6) tabR << "0\t0\t0\t0" << endl;
					else tabR << parR[idx][0] << "\t" << parR[idx][1] << "\t" << parR[idx][2]  << "\t" << parR[idx][3] << endl;
				//}
			}
		}
	}
	tabL.close();
	tabR.close();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_tdc_adc_L[i]->Integral()+h2_tdc_adc_R[i]->Integral());
		if(!notEmpty){
			delete h2_tdc_adc_L[i];
			delete h2_tdc_adc_R[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print("results_timeWalk.pdf(");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			cSLC[is][il] -> Print("results_timeWalk.pdf");
		}
	}
	c0 -> Print("results_timeWalk.pdf)");

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2) {
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetYaxis() -> CenterTitle();

	h2 -> SetTitleSize(0.5);

	h2 -> GetXaxis() -> SetLabelSize(0.13);
	h2 -> GetYaxis() -> SetLabelSize(0.13);
	h2 -> GetXaxis() -> SetTitleSize(0.13);
	h2 -> GetYaxis() -> SetTitleSize(0.13);

	h2 -> GetYaxis() -> SetTitleOffset(0.38);
	h2 -> GetYaxis() -> SetNdivisions(509);

	h2 -> GetXaxis() -> SetNdivisions(509);
}
// ========================================================================================================================================
double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}
// ========================================================================================================================================
int getMinNonEmptyBin(TH2F * h2){
	int nBinX = h2 -> GetXaxis() -> GetNbins();
	int nBinY = h2 -> GetYaxis() -> GetNbins();

	int binContent = 0;

	for(int i = 1 ; i <= nBinY ; i++){
		for(int j = 1 ; j <= nBinX ; j++){
			binContent = h2 -> GetBinContent(j,i);
			if(binContent!=0) return i;
		}
	}
	return 1;
}
// ========================================================================================================================================
int getMaxNonEmptyBin(TH2F * h2){
        int nBinX = h2 -> GetXaxis() -> GetNbins();
        int nBinY = h2 -> GetYaxis() -> GetNbins();

        int binContent = 0;

        for(int i = nBinY ; i > 0 ; i--){
                for(int j = 1 ; j <= nBinX ; j++){
                        binContent = h2 -> GetBinContent(j,i);
                        if(binContent!=0) return i;
                }
        }
        return nBinY;
}

