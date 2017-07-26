#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "Riostream.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <sstream>
#include <cmath>
#include <TProfile.h>
#include <TXMLNode.h>
#include <signal.h>

#include "BeamGe_AIDA_mod.h"
//#include "BeamGe_AIDA.h"

using namespace std;

int main(int argc, char* argv[]){
	
	int i,test;
	char beamfile[128];
	char AIDAfile[128];
	char outfile[128];	
	char gamfile[128];
	char calfile[128];

	int runB, runA, runA_sub, runG;

	if(argc==5){
	

	runA = atoi(argv[3]);
	runA_sub = atoi(argv[4]);
	runG = atoi(argv[2]);
	runB = atoi(argv[1]);

    sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
	sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
	sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
	sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_%d.root",runA,runA_sub);
	sprintf(outfile,"$NP1306DIR/BeamGe_AIDA/TSCheck_B%04d_G%04d_A%04d_%d.root",runB,runG,runA,runA_sub);
	}else if (argc==4){

		runB = atoi(argv[1]);
        runA = atoi(argv[3]);
        runG = atoi(argv[2]);
    
    	sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
        sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
        sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
        sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d.root",runA);
        sprintf(outfile,"$NP1306DIR/BeamGe_AIDA/TSCheck_B%04d_G%04d_A%04d.root",runB,runG,runA);

	}else { cout<<"USAGE: ./Beta_gamma #BigRIPSRun #EURICARun #AIDARUN (#AIDASUBRUN)"<<endl;
                        return 0;
	}

	BigRIPS_reco beam(beamfile);
	AIDAraw aidaraw(AIDAfile);
	eurica gamma(gamfile);
	Calib calib(calfile);

	map<Long64_t,Long64_t> mtsb,mtsg,mgoodi,mgoode; //time stamp table for bigrips and aida ion and electron
	map<Long64_t,Long64_t>::iterator imtsb, imtsg, imgoodi,imgoode;
	multimap<Long64_t,Long64_t> mtsi,mtse;
	multimap<Long64_t,Long64_t>::iterator imtsi, imtse;

	multimap<Long64_t,Long64_t> mbi, mge;	//time stamp correlation between ion,beta and bigrips
	multimap<Long64_t,Long64_t>::iterator imbi ,imge;
	
	cout<<"start building time map for BigRIPS, beta, ion, EURICA"<<endl;

	//aidaraw.GetIonTsEntry(mtsi);
	aidaraw.GetTsEntry(mtsi,mtse);
	gamma.GetTsEntry(mtsg);
	beam.GetTsEntry(mtsb);
	//aidaion.GetGoodness(mgoodi);
	//aidabeta.GetGoodness(mgoode);
	
	Long64_t bEntry, iEntry, eEntry, gEntry;
	Long64_t betaraw_ts,eurica_ts,ionraw_ts,bigrips_ts,goodness;

	// The timestamps are in timestamp unit.	10ns/1stamp
	
	Long64_t ts_check_windowL = -50000; 
	Long64_t ts_check_windowU = 50000*2; // 
	
	int toffset =0;
	TString soutfile = Form("TS_check_raw_%04d_%04d_%04d_%d.root",runB,runG,runA,runA_sub);
	TFile* rootfile = new TFile(soutfile,"RECREATE");
//	TTree* tree = new TTree("tree","tree");
//	TTree* tree2 = new TTree("beta_gam","beta_gam");
	TTree* tree = new TTree("TS_beta","TS_beta");
	TTree* tree2 = new TTree("TS_ion","TS_ion");

	tree->Branch("betaraw_ts",&betaraw_ts,"beta_ts/L");
	tree->Branch("eurica_ts",&eurica_ts,"eurica_ts/L");
	tree->Branch("goodness",&goodness,"goodness/I");
	tree->Branch("toffset",&toffset,"toffset/I");

	tree2->Branch("ionraw_ts",&ionraw_ts,"ion_ts/L");
	tree2->Branch("bigrips_ts",&bigrips_ts,"bigrips_ts/L");
	tree2->Branch("goodness",&goodness,"goodness/I");
	tree2->Branch("toffset",&toffset,"toffset/I");


	BeamGe_AIDA outtree;
	
	std::cout<<mtse.size()<<" events for beta events"<<std::endl;
	std::cout<<mtsi.size()<<" events for ion events"<<std::endl;
	std::cout<<mtsg.size()<<" events for gamma events"<<std::endl;
	std::cout<<mtsb.size()<<" events for bigrips events"<<std::endl;

	

	Int_t n_bins = 1000;
	//TH1F* h_TS_beta_Eurica= new TH1F("h_TS_beta_Eurica","",n_bins, ts_check_windowL+toffset, ts_check_windowU+toffset);
	//TH1F* h_TS_ion_BigRIPS= new TH1F("h_TS_ion_BigRIPS","",n_bins, ts_check_windowL+toffset, ts_check_windowU+toffset);
	//TH1F* h_TS_beta_Eurica_g= new TH1F("h_TS_beta_Eurica_g","",n_bins, ts_check_windowL+toffset, ts_check_windowU+toffset);
	//TH1F* h_TS_ion_BigRIPS_g= new TH1F("h_TS_ion_BigRIPS_g","",n_bins, ts_check_windowL+toffset, ts_check_windowU+toffset);

	Long64_t entry = 0;
	Long64_t temp;

	Long64_t ts1;
	Long64_t ts2;

	//int goodness;
	
	std::cout<<"Scan coincidence timestamps beta with Gamma, beam with ion" << endl;
		for(imtse=mtse.begin();imtse!=mtse.end();imtse++){
			if(imtse->first >0){
			ts1 = imtse->first + ts_check_windowL+toffset;
			ts2 = ts1 + ts_check_windowU;
			imtsg = mtsg.lower_bound(ts1);
			//imgoode= mgoode.find(imtse->second);
			//goodness= imgoode->second;
			betaraw_ts=imtse->first;
				while(imtsg!=mtsg.end() && (imtsg->first) <ts2){
				//mge.insert(std::make_pair(imtse->second,imtsg->second));
				eurica_ts=imtsg->first;
				//h_TS_beta_Eurica->Fill(imtse->first-imtsg->first);
				//	if(goodness==1)
				tree->Fill();
				imtsg++;
				
				}
			}
		}
	std::cout<<"Beta and gamma Coincidence tree entry now"<<tree->GetEntries()<<std::endl;
	//std::cout<<"Beta and gamma Coincidence "<<h_TS_beta_Eurica->GetEntries()<<std::endl;
	//std::cout<<"Good Beta and gamma Coincidence "<<h_TS_beta_Eurica_g->GetEntries()<<std::endl;

		for(imtsi=mtsi.begin();imtsi!=mtsi.end();imtsi++){
		if(imtsi->first >0){
			ts1 = imtsi->first + ts_check_windowL+toffset;
			ts2 = ts1 + ts_check_windowU;
			imtsb = mtsb.lower_bound(ts1);
			//imgoodi= mgoodi.find(imtse->second);
			//goodness= imgoodi->second;
			ionraw_ts=imtsi->first;
			while(imtsb!=mtsb.end() && (imtsb->first) <ts2){
				//mge.insert(std::make_pair(imtse->second,imtsg->second));
				bigrips_ts=imtsb->first;
				//if(goodness==1) 
				tree2->Fill();
				imtsb++;
				
				}
			}
		}	
	std::cout<<"Ion and beam Coincidence tree entry now"<<tree2->GetEntries()<<std::endl;
	//std::cout<<"Ion and beam Coincidence "<<h_TS_ion_BigRIPS->GetEntries()<<std::endl;
	//std::cout<<"Good Ion and beam Coincidence "<<h_TS_ion_BigRIPS_g->GetEntries()<<std::endl;

	
	rootfile->Write();

	/*TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->Divide(2,2);
	c1->cd(1); h_TS_beta_Eurica->Draw();
	c1->cd(2); h_TS_beta_Eurica_g->Draw();
	c1->cd(3); h_TS_ion_BigRIPS->Draw();
	c1->cd(4); h_TS_ion_BigRIPS_g->Draw();
	c1->SaveAs(Form("./TSCheck/TSCheck_B%04d_G%04d_A%04d_%d_toff_%d.eps",runB,runG,runA,runA_sub,toffset));
//	rootfile->Close();
	h_TS_beta_Eurica->Delete();
	h_TS_beta_Eurica_g->Delete();
	h_TS_ion_BigRIPS->Delete();
	h_TS_ion_BigRIPS_g->Delete();*/
	std::cout<<std::endl<<"Job finished "<<soutfile<<" is made"<<std::endl;
	
}
