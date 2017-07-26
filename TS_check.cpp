#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "Riostream.h"

#include <time.h>
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

//#include "AIDA_betaevt.h"
#include "AIDA_ionevt.h"
#include "eurica.h"
#include "BigRIPS_reco.h"
#include "Calib.h"
//#include "BeamGe_AIDA_mod.h"
//#include "BeamGe_AIDA.h"

using namespace std;

int main(int argc, char* argv[]){
	
	int i,test;
	char beamfile[128];
	char AIDAfile[128];
	char outfile[128];	
	char gamfile[128];
	char calfile[128];
	time_t tm;
	int time0;
	int runB, runA, runA_sub, runG;

	if(argc==5){
	

	runA = atoi(argv[3]);
	runA_sub = atoi(argv[4]);
	runG = atoi(argv[2]);
	runB = atoi(argv[1]);

 	sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
	sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
	sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
	sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_%d_packed_eventbuild.root",runA,runA_sub);
	sprintf(outfile,"$NP1306DIR/newmerger/TSCheck_B%04d_G%04d_A%04d_%d.root",runB,runG,runA,runA_sub);
	}else if (argc==4){

		runB = atoi(argv[1]);
        runA = atoi(argv[3]);
        runG = atoi(argv[2]);
    
    	sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
        sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
        sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
        sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_packed_eventbuild.root",runA);
        sprintf(outfile,"$NP1306DIR/newmerger/TSCheck_B%04d_G%04d_A%04d.root",runB,runG,runA);

	}else { cout<<"USAGE: ./Beta_gamma #BigRIPSRun #EURICARun #AIDARUN (#AIDASUBRUN)"<<endl;
                        return 0;
	}

	BigRIPS_reco beam(beamfile);
	AIDA_betaevt aidabeta(AIDAfile);
	AIDA_ionevt aidaion(AIDAfile);
	eurica gamma(gamfile);
	Calib calib(calfile);

	map<Long64_t,Long64_t> mtsb,mtsg,mgoodi,mgoode; //time stamp table for bigrips and aida ion and electron
	map<Long64_t,Long64_t>::iterator imtsb, imtsg, imgoodi,imgoode;
	map<Long64_t,Long64_t> mtsi,mtse;
	map<Long64_t,Long64_t>::iterator imtsi, imtse;

	multimap<Long64_t,Long64_t> mbi, mge;	//time stamp correlation between ion,beta and bigrips
	multimap<Long64_t,Long64_t>::iterator imbi ,imge;
	
	cout<<"start building time map for BigRIPS, beta, ion, EURICA"<<endl;
	
	aidaion.GetTsEntry(mtsi) ;
	aidabeta.GetTsEntry(mtse);
	gamma.GetTsEntry(mtsg);
	beam.GetTsEntry(mtsb);
	aidaion.GetGoodness(mgoodi);
	aidabeta.GetGoodness(mgoode);
	
	Long64_t bEntry, iEntry, eEntry, gEntry;
	Long64_t beta_ts,eurica_ts[10000],ion_ts,bigrips_ts[10000],goodness,beta_cor_ts[10000],eurica_firstts,bigrips_firstts;
	Long64_t beta_ts_b[10000],ion_ts_b[10000],bigrips_ts_b;
	Int_t mult,mult_i,mult_b;
	Double_t ion_z,ion_aoq;

	
	// The timestamps are in timestamp unit.	10ns/1stamp
	
	Long64_t ts_check_windowL = -50000; 
	Long64_t ts_check_windowU = 50000*2; //  +-500us
	
	int toffset=0;
	TString soutfile = Form("TS_check_%04d_%04d_%04d_%d.root",runB,runG,runA,runA_sub);
	TFile* rootfile = new TFile(soutfile,"RECREATE");
//	TTree* tree = new TTree("tree","tree");
//	TTree* tree2 = new TTree("beta_gam","beta_gam");
	TTree* tree = new TTree("TS_beta","TS_beta");
	TTree* tree2 = new TTree("TS_ion","TS_ion");
	TTree* tree3 = new TTree("TS_B","TS_B");
	TTree* tree4 = new TTree("TS_IB","TS_IB");

	tree->Branch("beta_ts",&beta_ts,"beta_ts/L");
	tree->Branch("mult",&mult,"mult/I");
	tree->Branch("eurica_ts",eurica_ts,"eurica_ts[mult]/L");
	tree->Branch("eurica_firstts",&eurica_firstts,"eurica_firstts/L");
	tree->Branch("goodness",&goodness,"goodness/I");
//	tree->Branch("toffset",&toffset,"toffset/I");

	tree2->Branch("ion_ts",&ion_ts,"ion_ts/L");
	tree2->Branch("mult",&mult,"mult/I");
	tree2->Branch("bigrips_ts",bigrips_ts,"bigrips_ts[mult]/L");
	tree2->Branch("bigrips_firstts",&bigrips_firstts,"bigrips_firstts/L");
	tree2->Branch("goodness",&goodness,"goodness/I");
//	tree2->Branch("toffset",&toffset,"toffset/I");

	tree3->Branch("mult_i",&mult_i,"mult_i/I");
	tree3->Branch("mult_b",&mult_b,"mult_b/I");
	tree3->Branch("ion_ts",ion_ts_b,"ion_ts[mult_i]/L");
	tree3->Branch("beta_ts",beta_ts_b,"beta_ts[mult_b]/L");
	tree3->Branch("bigrips_ts",&bigrips_ts_b,"bigrips_ts/L");
//	tree3->Branch("ion_z",&ion_z,"ion_z/D");
//	tree3->Branch("ion_aoq",&ion_aoq,"ion_aoq/D");

	tree4->Branch("ion_ts",&ion_ts,"ion_ts/L");
	tree4->Branch("mult",&mult,"mult/I");
	tree4->Branch("beta_ts",beta_ts_b,"beta_ts[mult]/L");
//	tree4->Branch("eurica_ts",eurica_ts,"eurica_ts[mult]/L");
//	tree3->Branch("goodness",&goodness,"goodness/I");
//	tree3->Branch("toffset",&toffset,"toffset/I");

//	BeamGe_AIDA outtree;
	
	std::cout<<mtse.size()<<" events for beta events"<<std::endl;
	std::cout<<mtsi.size()<<" events for ion events"<<std::endl;
	std::cout<<mtsg.size()<<" events for gamma events"<<std::endl;
	std::cout<<mtsb.size()<<" events for bigrips events"<<std::endl;

	

	Int_t n_bins = 1000;
	TH1F* h_TS_beta_Eurica= new TH1F("h_TS_beta_Eurica","",n_bins, ts_check_windowL+toffset, ts_check_windowU/2+toffset);
	TH1F* h_TS_ion_BigRIPS= new TH1F("h_TS_ion_BigRIPS","",n_bins, ts_check_windowL+toffset, ts_check_windowU/2+toffset);
	TH1F* h_TS_beta_Eurica_g= new TH1F("h_TS_beta_Eurica_g","",n_bins, ts_check_windowL+toffset, ts_check_windowU/2+toffset);
	TH1F* h_TS_ion_BigRIPS_g= new TH1F("h_TS_ion_BigRIPS_g","",n_bins, ts_check_windowL+toffset, ts_check_windowU/2+toffset);

	Long64_t entry = 0;
	Long64_t temp;

	Long64_t ts1;
	Long64_t ts2;

	//int goodness;
	int tester=0;
	std::cout<<"Scan coincidence timestamps beta with Gamma, beam with ion" << std::endl;
		for(imtse=mtse.begin();imtse!=mtse.end();imtse++){
			if(imtse->first >0){
			i=0;
			ts1 = imtse->first + ts_check_windowL+toffset;
			ts2 = ts1 + ts_check_windowU;
	//		tester+=1;
	//		if(tester%100000==0) {std::cout<<imtse->first<<" tsceneter //"<<ts1<<" ts1 //"<<ts2<<"ts2"<<std::endl;}
			imtsg = mtsg.lower_bound(ts1);
			imgoode= mgoode.find(imtse->second);
			goodness= imgoode->second;
			beta_ts=imtse->first;
				while(imtsg!=mtsg.end() && (imtsg->first) <ts2){
				//mge.insert(std::make_pair(imtse->second,imtsg->second));
				eurica_ts[i]=imtsg->first;
				h_TS_beta_Eurica->Fill(imtse->first-imtsg->first);
				if(goodness==1) { h_TS_beta_Eurica_g->Fill(imtse->first-imtsg->first);}
				imtsg++;
				i++;
				}
				mult=i;
			imtsg = mtsg.lower_bound(imtse->first-6000);
			eurica_firstts=imtsg->first;
			tree->Fill();
			}
		}
	std::cout<<"Beta and gamma Coincidence tree entry now"<<tree->GetEntries()<<std::endl;
	//std::cout<<"Beta and gamma Coincidence "<<h_TS_beta_Eurica->GetEntries()<<std::endl;
	std::cout<<"Good Beta and gamma Coincidence "<<h_TS_beta_Eurica_g->GetEntries()<<std::endl;

	tester=0;
	i=0;	
		for(imtsi=mtsi.begin();imtsi!=mtsi.end();imtsi++){
                        if(imtsi->first >0){
	//		tester+=1;
	//		if(imtsi->first >0){
			i=0;
			ts1 = imtsi->first + ts_check_windowL;
			ts2 = ts1 + ts_check_windowU;
			

	//		if(tester%10000==0) {std::cout<<imtsi->first<<" tsceneter //"<<ts1<<" ts1 //"<<ts2<<"ts2"<<std::endl;}
			imtsb = mtsb.lower_bound(ts1);
			imgoodi= mgoodi.find(imtsi->second);
			goodness= imgoodi->second;
			ion_ts=imtsi->first;
			while(imtsb!=mtsb.end() && (imtsb->first) <ts2){
				h_TS_ion_BigRIPS->Fill(imtse->first-imtsg->first);
				//mge.insert(std::make_pair(imtse->second,imtsg->second));
				bigrips_ts[i]=imtsb->first;
				if(goodness==1) { h_TS_ion_BigRIPS_g->Fill(imtse->first-imtsg->first);}
				imtsb++;
				i++;
				}
			mult=i;
			imtsb = mtsb.lower_bound(imtsi->first-6000);
			bigrips_firstts=imtsb->first;
			tree2->Fill();


			i=0;
			ts1 = imtsi->first - 10000;
			ts2 = imtsi->first + 50000;
			
			imtse = mtse.lower_bound(ts1);
			while(imtse!=mtse.end() && imtse->first<ts2){
			beta_ts_b[i] = imtse->first;
			imtse++;
			i++;
			}
			mult=i;
			tree4->Fill();

			}
		}	

	std::cout<<"Ion and beam Coincidence tree entry now"<<tree2->GetEntries()<<std::endl;
	std::cout<<"Good Ion and beam Coincidence "<<h_TS_ion_BigRIPS_g->GetEntries()<<std::endl;

        std::cout<<"Scan coincidence timestamps beta with ion" << std::endl;

	Long64_t betaion_check_windowL = 100000; // 1000us before the implantation
	Long64_t betaion_check_windowU = 100000; // 1000us after the implantation
	time0=time(&tm);
//	int nlist = beam.GetEntriesFast();
	int nlist=1;
	for(int ii=0;ii<nlist;ii++){
		int iion=0;
		beam.GetEntry(ii);
                ts1 = beam.timestamp - 100000;
                ts2 = beam.timestamp + betaion_check_windowU;
		imtsi = mtsi.lower_bound(ts1);
	        imtse = mtse.lower_bound(ts1);
		bigrips_ts_b=beam.timestamp;
		if(bigrips_ts_b>0){
			while(imtsi!=mtsi.end() && (imtsi->first) < ts2+1900000){


			aidaion.GetEntry(imtsi->second);
			
			ion_ts_b[iion]= aidaion.extTstart*4;
			iion++;

			imtsi++;
			}

		int ibeta=0;
   			while(imtse!=mtse.end() && (imtse->first) <ts2+1900000){
//			imgoode  = mgoode.find(imtse->second);
//			goodness = imgoode->second;
			aidabeta.GetEntry(imtse->second);
//			if(aidabeta.start_z>-1){
			//if(aidaion.CorrPosition(aidabeta,aidaion)==1)
			//	{
				beta_ts_b[ibeta]=aidabeta.extTstart*4;
				ibeta++;
//				if(i>9999){ cout<<"The multiplicity exceed limit"<<endl<<
//		"The multiplicity emerged at ion: "<<aidaion.stop_z<<","<<aidaion.stop_x<<","<<aidaion.stop_y<<endl<<"The multiplicity emerged at beta: "<<aidabeta.start_z<<","<<aidabeta.start_x<<","<<aidabeta.start_y<<endl;
				//return -1;}

			//				}
			//}
                        imtse++;
//                    }else imtse++;
			}

//		if(i>0)
//		{
		mult_i=iion;
		mult_b=ibeta;
		if(iion!=0 || ibeta!=0) {tree3->Fill();}
//		}
		}
		if(ii%100000==1) cout<<"jentry "<<ii<<"done "<<time(&tm)-time0<<"s spent"<<endl;
	}
                
        std::cout<<"ion and beta Coincidence tree entry now"<<tree3->GetEntries()<<std::endl;
        //std::cout<<"Beta and gamma Coincidence "<<h_TS_beta_Eurica->GetEntries()<<std::endl;
  //      std::cout<<"Good Beta and gamma Coincidence "<<h_TS_beta_Eurica_g->GetEntries()<<std::endl;



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
