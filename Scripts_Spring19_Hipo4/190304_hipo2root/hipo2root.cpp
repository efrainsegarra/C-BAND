#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

// Loading files from include directory
#include "constants.h"
#include "colors.h"

using namespace std;
double TDC_GLOBSHIFT[600] = {0.};
double FADC_GLOBSHIFT[600] = {0.};
double parA_L[600] = {0.};		// loaded
double parB_L[600] = {0.};		// loaded
double parA_R[600] = {0.};		// loaded
double parB_R[600] = {0.};		// loaded
double TDC_TDIFF[600] = {0.};		// loaded
double FADC_TDIFF[600] = {0.};		// loaded
double TDC_P2P[600] = {0.};		// loaded
double TDC_L2L[600] = {0.};		// loaded
double FADC_P2P[600] = {0.};		// loaded
double FADC_L2L[600] = {0.};		// loaded
double TDC_VEFF[600] = {0.};		// loaded
double FADC_VEFF[600] = {0.};		// loaded
double FADC_ATTEN_LENGTH[600] = {0,};	// loaded
double globPos[600][3] = {0.};		
double BARLENGTHS[]  = {163.7,201.9,51.2,51.2,201.9};

double getTriggerPhase( long timeStamp );
void LoadGlobalShift();
void LoadTimeWalk();
void LoadLROffsets();
void LoadPaddleOffsets();
void LoadLayerOffsets();
void LoadVelocityMap();
void LoadAttenuation();
void CreateGeo();
char* getRunNumber( char* parse );
double getBeamEnergy( int runNum );

int main(int argc, char** argv) {


	// check number of arguments
	if( argc < 4 ){
		cerr << BOLD(FRED("Incorrect number of arguments. Instead use:\n\t./hipo2root [outputFile] [byHand (0,1)] [intputFiles]\n\n"));
		return -1;
	}
	cout << BOLD(FYEL("Using inputs: ")) << argv[1] << " " << argv[2] << " " << argv[3] << "\n";
	//cout << BOLD(FYEL("Loading global bar shifts...\n"));
	//LoadGlobalShift();
	//cout << BOLD(FGRN("...Done!\n"));

	bool doByHand = (atoi(argv[2]) == 1 ? true : false );
	if( doByHand ){
		// Load calibration constants from include DIR:
		cout << BOLD(FYEL("Doing byHand -- Loading constants...\n"));
		LoadTimeWalk();
		LoadLROffsets();
		LoadPaddleOffsets();
		LoadLayerOffsets();
		LoadVelocityMap();
		LoadAttenuation();
		CreateGeo();
		cout << BOLD(FGRN("...Done!\n"));
	}


	// Event selection cuts for electron
	double cut_ep      =     2; //GeV
	double cut_chi2pid =     5;
	double cut_min_vz  =   -15; //cm
	double cut_max_vz  =    10; //cm
	double cut_W       =     0; //GeV
	double cut_uvw     =    15; //cm
	double cut_Epcal   = 0.060; //GeV (60 MeV)
	double cut_tof_e   =    10; //ns

	// Create output tree
	cout << BOLD(FYEL("Creating output...\n"));
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");
	// BAND variables
	int nHits;
	int sector, layer, component;
	double adcLcorr, adcRcorr;
	double meantimeFadc, meantimeTdc;
	double difftimeFadc, difftimeTdc;
	double x,y,z;
	double dL, theta_n, phi_n, ToF;
	double beta, p_n, phi_nq, theta_nq;
	double E_n, phi_en, CosTheta_nq;
	double Xp, Wp, As;
	// Raw BAND variables 
	int nADC, nTDC;
	double phaseCorr;
	double adcLraw, adcRraw;
	double tTdcLraw, tTdcRraw;
	double tFadcLraw, tFadcRraw;
	// By hand variables using tables
	double byHand_adcL, byHand_adcR;
	double byHand_meantimeFadc, byHand_meantimeTdc;
	double byHand_difftimeFadc, byHand_difftimeTdc;
	double byHand_x, byHand_y, byHand_z;
	double byHand_dL;
	// CLAS variables
	double Ebeam, p_e, theta_e, phi_e;
	double q, theta_q, phi_q;
	double Q2, nu, xB, W2, EoP;
	double start_time;
	int sector_e;
	double vrt_x_e, vrt_y_e, vrt_z_e, t_e, beta_e, chi2pid_e;

	// Branches for CLAS
	outTree->Branch("Ebeam",		&Ebeam		,	"Ebeam/D");
	outTree->Branch("p_e",			&p_e		,	"p_e/D");
	outTree->Branch("theta_e",		&theta_e	,	"theta_e/D");
	outTree->Branch("phi_e",		&phi_e		,	"phi_e/D");
	outTree->Branch("q",			&q		,	"q/D");
	outTree->Branch("theta_q",		&theta_q	,	"theta_q/D");
	outTree->Branch("phi_q",		&phi_q		,	"phi_q/D");
	outTree->Branch("Q2",			&Q2		,	"Q2/D");
	outTree->Branch("nu",			&nu		,	"nu/D");
	outTree->Branch("xB",			&xB		,	"xB/D");
	outTree->Branch("W2",			&W2		,	"W2/D");
	outTree->Branch("EoP",			&EoP		,	"EoP/D");
	outTree->Branch("STTime",		&start_time	,	"STTime/D");
	outTree->Branch("sector_e",		&sector_e	,	"sector_e/I");
	outTree->Branch("vrt_x_e",		&vrt_x_e	,	"vrt_x_e/D");
	outTree->Branch("vrt_y_e",		&vrt_y_e	,	"vrt_y_e/D");
	outTree->Branch("vrt_z_e",		&vrt_z_e	,	"vrt_z_e/D");
	outTree->Branch("t_e",			&t_e		,	"t_e/D");
	outTree->Branch("beta_e",		&beta_e		,	"beta_e/D");
	outTree->Branch("chi2pid_e",		&chi2pid_e	,	"chi2pid_e/D");
	// Branches for BAND
	outTree->Branch("nHits",		&nHits		,	"nHits/I");
	outTree->Branch("sector",		&sector		,	"sector/I");
	outTree->Branch("layer",		&layer		,	"layer/I");
	outTree->Branch("component",		&component	,	"component/I");
	outTree->Branch("adcLcorr",		&adcLcorr	,	"adcLcorr/D");
	outTree->Branch("adcRcorr",		&adcRcorr	,	"adcRcorr/D");
	outTree->Branch("meantimeFadc",		&meantimeFadc	,	"meantimeFadc/D");
	outTree->Branch("meantimeTdc",		&meantimeTdc	,	"meantimeTdc/D");
	outTree->Branch("difftimeFadc",		&difftimeFadc	,	"difftimeFadc/D");
	outTree->Branch("difftimeTdc",		&difftimeTdc	,	"difftimeTdc/D");
	outTree->Branch("x",			&x		,	"x/D");
	outTree->Branch("y",			&y		,	"y/D");
	outTree->Branch("z",			&z		,	"z/D");
	outTree->Branch("dL",			&dL		,	"dL/D");
	outTree->Branch("theta_n",		&theta_n	,	"theta_n/D");
	outTree->Branch("phi_n",		&phi_n		,	"phi_n/D");
	outTree->Branch("ToF",			&ToF		,	"ToF/D");
	outTree->Branch("beta",			&beta		,	"beta/D");
	outTree->Branch("p_n",			&p_n		,	"p_n/D");
	outTree->Branch("phi_nq",		&phi_nq		,	"phi_nq/D");
	outTree->Branch("theta_nq",		&theta_nq	,	"theta_nq/D");
	outTree->Branch("E_n",			&E_n		,	"E_n/D");
	outTree->Branch("phi_en",		&phi_en		,	"phi_en/D");
	outTree->Branch("CosTheta_nq",		&CosTheta_nq	,	"CosTheta_nq/D");
	outTree->Branch("Xp",			&Xp		,	"Xp/D");
	outTree->Branch("Wp",			&Wp		,	"Wp/D");
	outTree->Branch("As",			&As		,	"As/D");
	// Raw branches for BAND
	outTree->Branch("nADC",			&nADC		,	"nADC/I");
	outTree->Branch("nTDC",			&nTDC		,	"nTDC/I");
	outTree->Branch("phaseCorr",		&phaseCorr	,	"phaseCorr/D");
	outTree->Branch("adcLraw",		&adcLraw	,	"adcLraw/D");
	outTree->Branch("adcRraw",		&adcRraw	,	"adcRraw/D");
	outTree->Branch("tTdcLraw",		&tTdcLraw	,	"tTdcLraw/D");
	outTree->Branch("tTdcRraw",		&tTdcRraw	,	"tTdcRraw/D");
	outTree->Branch("tFadcLraw",		&tFadcLraw	, 	"tFadcLraw/D");
	outTree->Branch("tFadcRraw",		&tFadcRraw	,	"tFadcRraw/D");
	// By-hand branches for BAND
	outTree->Branch("byHand_adcL",		&byHand_adcL	,	"byHand_adcL/D");
	outTree->Branch("byHand_adcR",		&byHand_adcR	, 	"byHand_adcR/D");
	outTree->Branch("byHand_meantimeFadc",	&byHand_meantimeFadc	,	"byHand_meantimeFadc/D");
	outTree->Branch("byHand_meantimeTdc",	&byHand_meantimeTdc	,	"byHand_meantimeTdc/D");
	outTree->Branch("byHand_difftimeFadc",	&byHand_difftimeFadc	,	"byHand_difftimeFadc/D");
	outTree->Branch("byHand_difftimeTdc",	&byHand_difftimeTdc	,	"byHand_difftimeTdc/D");
	outTree->Branch("byHand_x",		&byHand_x	,	"byHand_x/D");
	outTree->Branch("byHand_y",		&byHand_y	,	"byHand_y/D");
	outTree->Branch("byHand_z",		&byHand_z	,	"byHand_z/D");
	outTree->Branch("byHand_dL",		&byHand_dL	,	"byHand_dL/D");

	cout << BOLD(FGRN("...Done!\n"));
	
	// Load input file
	for( int i = 3 ; i < argc ; i++ ){
		int thisRun = atoi(getRunNumber(argv[i]));
		double fixed_Ebeam = getBeamEnergy(thisRun);
		// Setup initial vector for beam
		TVector3 beamVec(0,0,fixed_Ebeam);

		TString inputFile = argv[i];
		cout << BOLD(FBLU("Now working on file: ")) << inputFile << endl;
		hipo::reader reader;
		reader.open(inputFile);
		
		//Read Dictionary of Hipo File  // new hipo4
		hipo::dictionary  factory;      // new hipo4
		reader.readDictionary(factory); // new hipo4
		//factory.show();               // new hipo4

		// Banks for EMC-SRC physics
		BEvent		event		(factory.getSchema("REC::Event"       ));
		BParticle	particles	(factory.getSchema("REC::Particle"    ));
		BCalorimeter	calo		(factory.getSchema("REC::Calorimeter" ));
		BScintillator	scintillator	(factory.getSchema("REC::Scintillator"));
		BBand		band_hits	(factory.getSchema("BAND::hits"       ));
		hipo::bank 	BAND_ADC	(factory.getSchema("BAND::adc"));
		hipo::bank 	BAND_TDC	(factory.getSchema("BAND::tdc"));
		hipo::bank 	RUN_config	(factory.getSchema("RUN::config"));

		//One also needs a hipo::event object which is called from the reader for each event to get
		//the information for each bank
		hipo::event readevent;  // new hipo4

		// Loop over events in hipo fil
		int event_counter = -1;
		while(reader.next()==true){
			event_counter++;
			//Reader has to load information about event in hipo::event class
			reader.read(readevent); // new hipo4
			//if( event_counter == 100000 ) break;

			//Load explicitly all information for each bank for the event
			readevent.getStructure(event       );   // new hipo4
			readevent.getStructure(particles   );   // new hipo4
			readevent.getStructure(calo        );   // new hipo4
			readevent.getStructure(scintillator);   // new hipo4
			readevent.getStructure(band_hits   );   // new hipo4
			readevent.getStructure(BAND_ADC	   );	// new hipo4
			readevent.getStructure(BAND_TDC	   );	// new hipo4
			readevent.getStructure(RUN_config  );	// new hipo4
		

			//Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
			//The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.
			// Clear electron quantities:
			Ebeam,p_e,theta_e,phi_e,q,theta_q,phi_q,Q2,nu,xB,W2,EoP		= 0.;
			start_time, sector_e, vrt_x_e, vrt_y_e, vrt_z_e, t_e		= 0.;
			beta_e,chi2pid_e						= 0.;
	
			// Clear BAND REC quantities:
			nHits,sector,layer,component,adcLcorr,adcRcorr			= 0.;
			meantimeFadc,meantimeTdc,difftimeFadc,difftimeTdc 		= 0.;
			x,y,z,dL,theta_n,phi_n,ToF,beta,p_n,phi_nq,theta_nq		= 0.;
			E_n,phi_en,CosTheta_nq,Xp,Wp,As					= 0.;

			// Clear BAND RAW quantities:
			nADC, nTDC 							= 0.;
			phaseCorr							= 0.;
			adcLraw, adcRraw, tTdcLraw, tTdcRraw, tFadcLraw, tFadcRraw 	= 0.;
			byHand_adcL, byHand_adcR, byHand_meantimeFadc			= 0.;
			byHand_meantimeTdc, byHand_difftimeFadc, byHand_difftimeTdc	= 0.;
			byHand_x, byHand_y, byHand_z, byHand_dL				= 0.;

			if(event_counter%1000000==0) cout << BOLD(FBLU("\t on event: ")) << event_counter << endl;

			// Debugging print option
			//band_hits.show();
			//event.show();
			//scintillator.show();	

			// Grab the electron information:
			Ebeam = fixed_Ebeam;
			int pid0		= particles.getPid    (0);       // electron candidate id assigned by clas
			TVector3 eVert		= particles.getV3v    (0);       // electron candidate vertex vector
			TVector3 eVec 		= particles.getV3P    (0);       // electron candidate momentum vector
			double chr0     	= particles.getCharge (0);       // electron candidate charge
			beta_e		   	= particles.getBeta   (0);       // electron candidate beta = v/c
			chi2pid_e 	 	= particles.getChi2pid(0);       // electron candidate goodness of pid fit
			int eStatus    		= particles.getStatus (0);       // electron candidate status
			// Calorimeter bank
			double Epcal 		= calo.getPcalE(0); 
			double Ee    		= calo.getTotE (0);
			//float lU    = calo.getLU   (0);	// electron candidate distance on U-side [cm?]
			//float lV    = calo.getLV   (0);	// electron candidate distance on V-side [cm?]
			//float lW    = calo.getLW   (0);	// electron candidate distance on W-side [cm?]

			// Event vertex time calibrated from FTOF
			double t_vtx   		= event.getSTTime(0);
			double tof_e  		= t_e - t_vtx;		// electron candidate time-of-flight [ns]

			// Define the q vector
			TVector3 qVec 	= beamVec - eVec;

			// Only keep events for which the first particle is an electron
			if(             (pid0!=11              )||
					(chr0!=-1              ) 
					//(chi2pid>=cut_chi2pid  )||
					//(eP.Mag()<=cut_ep      )||
					//(eP.Mag()>=fixed_Ebeam       )||
					//(eP_vertex.Z()>cut_max_vz  )||
					//(eP_vertex.Z()<cut_min_vz  )||
					//(lU<cut_uvw            )||
					//(lV<cut_uvw            )||
					//(lW<cut_uvw            )||
					//(Epcal<cut_Epcal       )||
					//(TMath::Sqrt(W2)<=cut_W)||
					//(tof_e<cut_tof_e       )
			  ) continue;


			// Calculate electron kinematic variables
			p_e 		= eVec.Mag();
			theta_e 	= eVec.Theta();
			phi_e		= eVec.Phi();
			q		= qVec.Mag();
			theta_q		= qVec.Theta();
			phi_q		= qVec.Phi();
			nu 		= fixed_Ebeam - sqrt( p_e*p_e + mE*mE );
			Q2 		= q*q - nu*nu;
			xB 		= Q2 / (2.*mP*nu);
			W2     		= mP*mP - Q2 + 2*nu*mP;	
			EoP		= Ee / p_e;
			start_time 	= t_vtx;					
			sector_e	= calo.getSector(0);	// technically this is a bit wrong, but it should all have same sector...
			vrt_x_e		= eVert.X();
			vrt_y_e		= eVert.Y();
			vrt_z_e		= eVert.Z();
			t_e     	= scintillator.getTime(0);

			// Looking in BAND
			nHits = band_hits.getRows();
			if( nHits == 1){
				for(int hit = 0; hit < nHits; hit++) {

					int    barKey           = band_hits.getBarKey  		(hit); 
					sector          	= band_hits.getSector      	(hit);
					layer			= band_hits.getLayer		(hit);
					component		= band_hits.getComponent	(hit);

					adcLcorr         	= band_hits.getAdcLcorr    	(hit);
					adcRcorr         	= band_hits.getAdcRcorr    	(hit);
					meantimeFadc    	= band_hits.getMeantimeFadc	(hit);
					meantimeTdc		= band_hits.getMeantimeTdc	(hit);	

					difftimeFadc		= band_hits.getDifftimeFadc	(hit);
					difftimeTdc		= band_hits.getDifftimeTdc	(hit);

					x			= band_hits.getX		(hit);
					y			= band_hits.getY		(hit);
					z			= band_hits.getZ		(hit);

					// Form path length and angles
					dL  = sqrt( x*x + y*y + z*z );	// [cm]
					theta_n = acos( z/dL );
					phi_n = atan2( y, x );
				
					// Calculate ToF based on meantimeFadc - STTime - global offset
					ToF = meantimeFadc - t_vtx - FADC_GLOBSHIFT[barKey];	// [ns]
					// Calculate beta
					beta = dL/(cAir*ToF);	
					// Calculate momentum
					p_n = 1./sqrt( 1./(beta*beta) - 1.) * mN; 	// [GeV]

					// Create 3 vector
					TVector3 nVec; nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
	
					// Solve for theta_nq, phi_nq
					TVector3 norm_scatter = qVec.Cross( beamVec );
        				norm_scatter    = norm_scatter.Unit();
					TVector3 norm_reaction = qVec.Cross( nVec );
        				norm_reaction   = norm_reaction.Unit();
					phi_nq   = norm_scatter.Angle( norm_reaction );
        				theta_nq = nVec.Angle( qVec );
					TVector3 direction = norm_scatter.Cross(norm_reaction);
					if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
					}
					else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
						phi_nq *= (-1);
        				}

					// Solve for last variables:
					E_n = sqrt( p_n*p_n + mN*mN );		// [GeV]

					// Now look at electron-neutron quantities to build the physics
					phi_en = phi_n - phi_e;
					if (phi_en < -M_PI) phi_en += 2.*M_PI;
					if (phi_en >  M_PI) phi_en -= 2.*M_PI;
					CosTheta_nq = cos(theta_nq);
					Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
					Wp = sqrt((mD*mD) - Q2 + (mP*mP) + 2.*mD*(nu-E_n) -2.* nu * E_n + 2.*q*p_n*CosTheta_nq);
					As = (E_n - p_n*CosTheta_nq)/mN;


				}// end loop over band hits in an event
				
				
				// Now check the corresponding ADC banks -- we should only have 2 ADCs, 2 TDCs:
				if( doByHand ){
					nADC = BAND_ADC.getRows();
					nTDC = BAND_TDC.getRows();
					if( nADC == 2 && nTDC == 2){
						//RUN_config.show();
						long timestamp = RUN_config.getLong(4,0);
						phaseCorr = getTriggerPhase(timestamp);
						int adc_barKey, tdc_barKey;

						// Get the raw ADC information, uncorrected
						//BAND_ADC.show();
						for(int aIdx = 0 ; aIdx < nADC ; aIdx++){
							int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
							int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
							int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
							int   ADC_component = BAND_ADC.getInt  (2,aIdx);
							if( ADC_order == 0 ){
								adcLraw = (float)(BAND_ADC.getInt(4,aIdx));
								tFadcLraw = BAND_ADC.getFloat(5,aIdx);
								adc_barKey = ADC_sector*100 + ADC_layer*10 + ADC_component;
							}
							if( ADC_order == 1 ){
								adcRraw = (float)(BAND_ADC.getInt(4,aIdx));
								tFadcRraw = BAND_ADC.getFloat(5,aIdx);
								adc_barKey = ADC_sector*100 + ADC_layer*10 + ADC_component;
							}
						}

						// Get the raw TDC information, uncorrected
						//BAND_TDC.show();
						for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
							int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
							int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
							int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
							int   TDC_component = BAND_TDC.getInt  (2,tIdx);
							TDC_order -= 2;
							if( TDC_order == 0 ){
								tTdcLraw = (float)(BAND_TDC.getInt(4,tIdx))*0.02345;
								tdc_barKey = TDC_sector*100 + TDC_layer*10 + TDC_component;
							}
							if( TDC_order == 1 ){
								tTdcRraw = (float)(BAND_TDC.getInt(4,tIdx))*0.02345;
								tdc_barKey = TDC_sector*100 + TDC_layer*10 + TDC_component;
							}

						}

						// Correct everything by hand using tables in include DIR
						if( 	adcLraw != 0. && adcRraw != 0. && 			// ADC non zero
								tFadcLraw != 0. && tFadcRraw != 0. &&			// ADC time non zero
								tTdcLraw != 0. && tTdcRraw != 0. && 			// TDC time non zero
								adc_barKey == sector*100+layer*10+component && 		// matching bar ID for ADC
								tdc_barKey == sector*100+layer*10+component ){		/// matching bar ID for TDC

							int barID = sector*100+layer*10+component;					

							// TDC phase correction:
							tTdcLraw -= phaseCorr;
							tTdcRraw -= phaseCorr;
							// TDC time walk
							tTdcLraw = tTdcLraw - (parA_L[barID]/sqrt(adcLraw) + parB_L[barID]);
							tTdcRraw = tTdcRraw - (parA_R[barID]/sqrt(adcRraw) + parB_R[barID]);
							// TDiff:
							byHand_difftimeTdc = (tTdcLraw-tTdcRraw) - TDC_TDIFF[barID];
							byHand_difftimeFadc = (tFadcLraw-tFadcRraw) - FADC_TDIFF[barID];
							// Meantime:
							byHand_meantimeTdc = (tTdcLraw+tTdcRraw)/2. - fabs(TDC_TDIFF[barID]/2.) - TDC_P2P[barID] - TDC_L2L[barID];
							byHand_meantimeFadc = (tFadcLraw+tFadcRraw)/2. - fabs(FADC_TDIFF[barID]/2.) - FADC_P2P[barID] - FADC_L2L[barID];
							//cout << "Raw, corrected TDC mean time: " << (tTdcLraw+tTdcRraw)/2. << " " << byHand_meantimeTdc << "\n";
							//cout << sector << " " << layer << " " << component << " " 
								//<< parA_L[barID] << " " << parB_L[barID] << " "
								//<< parA_R[barID] << " " << parB_R[barID] << " "
								//<< TDC_TDIFF[barID] << " " << FADC_TDIFF[barID] << " "
								//<< TDC_TDIFF[barID] << " " 
								//<< TDC_P2P[barID] << " " << TDC_L2L[barID] << "\n"; 
								//<< FADC_P2P[barID] << " " << FADC_L2L[barID] << "\n";
								//<< TDC_VEFF[barID] << " " << FADC_VEFF[barID] << " "
								//<< FADC_ATTEN_LENGTH[barID] << " "
								//<< globPos[barID][0] << " " << globPos[barID][1] << " " << globPos[barID][2] << "\n";
							//cout << "HitFinder mean time: " << meantimeTdc << "\n";
							//cout << "difference : " << fabs(meantimeTdc-byHand_meantimeTdc) << "\n\n";
							//if( fabs(meantimeTdc-byHand_meantimeTdc) > 1 ) cout << "********BIG DIFF********\n\n";
															
							// ADC attenuation corr:
							double xpos_tdc = (-1./2) * byHand_difftimeTdc * TDC_VEFF[barID];
							double xpos_fadc = (-1./2) * byHand_difftimeFadc * FADC_VEFF[barID];
							double sectorLen = BARLENGTHS[sector-1];
							double mu_cm = FADC_ATTEN_LENGTH[barID];
							byHand_adcL = adcLraw * exp( (sectorLen/2.-xpos_fadc) / mu_cm );
							byHand_adcR = adcRraw * exp( (sectorLen/2.+xpos_fadc) / mu_cm );
							// Path length:
							byHand_x = (xpos_tdc+xpos_fadc)/2.;
							byHand_x += globPos[barID][0];
							byHand_y = globPos[barID][1];
							byHand_z = globPos[barID][2];
							byHand_dL = sqrt( byHand_x*byHand_x + byHand_y*byHand_y + byHand_z*byHand_z );

						}
					} // end if for raw ADC==2 and TDC =2
				}// end if doByHand
				
	
			} // end if for nHits == 1 for BAND

			outTree->Fill();


		}// end file
		
		cout << BOLD(FGRN("\t...done with file\n"));
	}// end loop over files
	

	cout << BOLD(FGRN("...writing and closing!")) << endl;
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}


double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

void LoadGlobalShift(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;
	
	f.open("../include/global_offset_fadc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	//f.open("../include/global_offset_tdc.txt");
	//while(!f.eof()){
	//	f >> sector;
	//	f >> layer;
	//	f >> component;
	//	barId = 100*sector + 10*layer + component;
	//	f >> pol0;
	//	f >> height;
	//	f >> mean;
	//	f >> sig;
	//	TDC_GLOBSHIFT[barId] = mean;
	//}
}

void LoadTimeWalk(){
	ifstream f;
	int sector, layer, component, barId;
	double parA, parB, temp;

	f.open("../include/time_walk_corr_left.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_L[barId] = parA;
		parB_L[barId] = parB;
	}
	f.close();

	f.open("../include/time_walk_corr_right.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_R[barId] = parA;
		parB_R[barId] = parB;
	}
	f.close();
	return;
}
void LoadLROffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_off, fadc_off, temp;

	f.open("../include/lr_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> tdc_off;
		f >> fadc_off;
		f >> temp;
		f >> temp;
		TDC_TDIFF[barId] = tdc_off;
		FADC_TDIFF[barId] = fadc_off;
	}
	f.close();
	return;
}
void LoadPaddleOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double offset_fadc, offset_tdc, temp;

	f.open("../include/paddle_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_fadc;
		f >> temp;
		FADC_P2P[barId] = offset_fadc;
	}
	f.close();
	f.open("../include/paddle_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_tdc;
		f >> temp;
		TDC_P2P[barId] = offset_tdc;
	}
	f.close();

	return;
}
void LoadLayerOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double offset_fadc, offset_tdc, temp;

	f.open("../include/layer_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_fadc;
		f >> temp;
		FADC_L2L[barId] = offset_fadc;
	}
	f.close();
	f.open("../include/layer_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_tdc;
		f >> temp;
		TDC_L2L[barId] = offset_tdc;
	}
	f.close();

	return;
}
void LoadVelocityMap(){
	ifstream f;
	int sector, layer, component, barId;
	double veff_tdc, veff_fadc, temp;

	f.open("../include/effective_velocity.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> veff_tdc;
		f >> veff_fadc;
		f >> temp;
		f >> temp;
		TDC_VEFF[barId] = veff_tdc;
		FADC_VEFF[barId] = veff_fadc;
	}
	f.close();

	return;
}
void LoadAttenuation(){
	ifstream f;
	int sector, layer, component, barId;
	double atten_len, temp;

	f.open("../include/attenuation_lengths.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> atten_len;
		f >> temp;
		FADC_ATTEN_LENGTH[barId] = atten_len;
	}
	f.close();

	return;
}

void CreateGeo(){
	// GEOMETRY PARAMETERS
	int sectNum = 5;	// Number of sectors (blocks)
	int layNum = 6; 	// Number of layer
	int compNum = 7;    // Maximum Number of components in a sector
	// Number of components per layer and sector [layer][sector] index
	int compNumSecLay[6][5] = { {3,7,6,6,2}, {3,7,6,6,2}, {3,7,6,6,2}, {3,7,6,6,2}, {3,7,5,5,0}, {3,7,6,6,2} };

	double distVetoLead = 17.463;			// distance from veto to lead wall (cm)

	double thickness = 7.2;				// thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};  // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6
	double zOffset = 100;                               // distance from center first layer to target.
	double surveyBox[4][3] = {  	{-24.05,-21.10,-302.69},
					{-24.05, 22.81,-302.69},
					{ 24.10,-21.06,-302.57},
					{ 24.37, 22.81,-302.64}  	};
	double lenLG = 8.9; // [cm] -- length of LG
	double lenET = 16.; // [cm] -- length of PMT tube for ET PMTs
	double lenHam = 13.3; // [cm] -- length of PMT tube for Hamamatsu PMTs

	double avgX = ( (surveyBox[0][0] + surveyBox[2][0]) + (surveyBox[1][0] + surveyBox[3][0]) )/2.;
	double avgY = ( (surveyBox[0][1] + thickness*3.) + (surveyBox[1][1] - thickness*3.) 
			+ (surveyBox[2][1] + thickness*3.) + (surveyBox[3][1] - thickness*3.) )/4.;
	double avgZ = ( surveyBox[0][2] + surveyBox[1][2] + surveyBox[2][2] + surveyBox[3][2] )/4.;

	double globPt[] = {avgX,avgY,avgZ}; // single global position

	double barLengthSector[] = {164, 202, 51, 51, 202} ;           // Bar length in each layer (cm)


	for( int layer = 1 ; layer < layNum + 1 ; layer++){
		double localZ = 0.;
		localZ += (layerGap[1])/2.; // taking this thickness because wrapping material
		// isn't 'squeezed' by the weight of the detector
		for( int i = 1 ; i < layer ; i++ ){
			localZ += layerGap[i-1];
		}

		for( int sector = 1 ; sector < sectNum + 1 ; sector++){
			int nBars = compNumSecLay[layer-1][sector-1];
			for( int bar = 1 ; bar < nBars + 1 ; bar++){
				int key = sector*100+layer*10+bar;
				double localY = 666666.;				
				double localX = 666666.;

				double secYOff = 666666.;

				if( sector == 1){
					secYOff = 10.;
					localX = 0.;
				}
				else if( sector == 2){
					secYOff = 3.;
					localX = 0.;
				}
				else if( sector == 3 || sector == 4){
					secYOff = -3.;
					if( sector == 3){
						localX = (surveyBox[2][0]+surveyBox[3][0])/2. + lenLG + lenET + barLengthSector[sector-1]/2.;
					}
					else if( sector == 4){
						localX = (surveyBox[0][0]+surveyBox[1][0])/2. - lenLG - lenET - barLengthSector[sector-1]/2.;
					}
				}
				else if( sector == 5){
					secYOff = -5.;
					localX = 0.;
				}
				localY = secYOff*thickness + (nBars - (bar-1) )*thickness - thickness/2.;

				globPos[key][0] =  localX + globPt[0];
				globPos[key][1] =  localY + globPt[1];
				globPos[key][2] =  localZ + globPt[2];
				
			}


		}


	}
	return;
}
char* getRunNumber( char* parse ){
        char * parse_copy = (char*) malloc( strlen(parse)+1 );
        char * parsed;

        strcpy( parse_copy, parse );
        char * loop = strtok( parse_copy, ".");
        while(loop){
                char* equals_sign = strchr(loop, '_');
                if (equals_sign){
                        *equals_sign = 0;
                        equals_sign++;
                        parsed = (char*)malloc(strlen(equals_sign) + 1);
                        strcpy(parsed, equals_sign);
                }
                loop = strtok(NULL, ".");
        }
        free(parse_copy);
        return parsed;
}

double getBeamEnergy( int runNum ){
        double thisEn = 0;

        if( runNum <= 6399 ) thisEn = 10.6;
        else{ thisEn = 10.2; }
        if( runNum == 6523 || runNum == 6524 || runNum == 6525 ) thisEn = 10.;

        return thisEn;
}
