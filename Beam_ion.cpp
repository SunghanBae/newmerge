/*
This code is for bigrips and AIDAion merge 20170519. Bae.
*/


#include "TFile.h"
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

//#include "BeamGe_AIDA.h"
#include "BeamGe_AIDA_mod.h"

using namespace std;

int main(int argc, char* argv[]){
	
	int i,test;
	char beamfile[128];
	char AIDAfile[128];
	char outfile[128];	
	char gamfile[128];
	char calfile[128];

	int runB, runA, runA_sub;

	if(argc==4){
	runB = atoi(argv[1]);
	runA = atoi(argv[2]);
	runA_sub = atoi(argv[3]);
	

	//sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
	//sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
	sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
	sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_%d_packed_eventbuild.root",runA,runA_sub);
	sprintf(outfile,"$NP1306DIR/newmerge/BA_B%04d_A%04d_%d.root",runB,runA,runA_sub);
	}else if (argc==3){

        runB = atoi(argv[1]);
        runA = atoi(argv[2]);
    
    //sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
    //sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
        sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d_recoPID.root",runB);
        sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_packed_eventbuild.root",runA);
        sprintf(outfile,"$NP1306DIR/newmerge/BA_B%04d_A%04d.root",runB, runA);

	}else { cout<<"USAGE: ./Beam_Ion #BigRIPSRUN #AIDARUN (#AIDASUBRUN)"<<endl;
                        return 0;
	}


	BigRIPS_reco beam(beamfile);
	AIDA_ionevt aidaion(AIDAfile);
	
	map<Long64_t,Long64_t> mtsb,mtsg,mtsi; //time stamp table for bigrips and aida ion and electron
	map<Long64_t,Long64_t>::iterator imtsb, imtsg,imtsi;
	multimap<Long64_t,Long64_t> mtse,mtsi2;
	multimap<Long64_t,Long64_t>::iterator imtse,imtsi2;

	multimap<Long64_t,Long64_t> mbi, mge, mei;	//time stamp correlation between ion,beta and bigrips
	multimap<Long64_t,Long64_t>::iterator imbi ,imge ,imei;
	
	cout<<"start building time map for BigRIPS, ion"<<endl;
	aidaion.GetTsEntry(mtsi);
//	aidaion.GetTsEntry2(mtsi2);
	beam.GetTsEntry(mtsb);
		
	Long64_t bEntry, iEntry, eEntry, gEntry;
	
	TFile* rootfile = new TFile(outfile,"RECREATE");
	//TTree* tree = new TTree("tree","tree");
	TTree* tree2 = new TTree("ion","ion");
	
	BeamGe_AIDA outtree;
	//outtree.TreeBranch(tree);
	outtree.TreeBranchIon(tree2);

	std::cout<<mtsb.size()<<" events for BigRIPS"<<std::endl;
	std::cout<<mtsi.size()<<" events for ion"<<std::endl;
//	std::cout<<mtsi2.size()<<" events for ion 2"<<std::endl;

	// The timestamps are in 10ns/1stamp unit
	
//	Long64_t beam_ion_offset =1550;
	Long64_t beam_ion_tsL = -1800; // 
	Long64_t beam_ion_tsU = -1300; // 
	
	Long64_t entry = 0;
	Long64_t temp;

	Long64_t ts1;
	Long64_t ts2;
	
	std::cout<<"Scan coincidence timestamps BigRIPS with Ion" << endl;
	for(imtsi=mtsi.begin();imtsi!=mtsi.end();imtsi++){
		if(imtsi->first >0){ 		//If the timestamp > 0

			ts1 = imtsi->first + beam_ion_tsL;
			ts2 = imtsi->first + beam_ion_tsU;
			imtsb = mtsb.lower_bound(ts1);

			if(imtsb!=mtsb.end() && imtsb->first <ts2){
				mbi.insert(std::make_pair(imtsi->second, imtsb->second));
				if(mbi.count(imtsi->second)>1)
				{std::cout<<"Beam ion doubly assinged to AIDA ion."<<std::endl; imbi=mbi.lower_bound(imtsi->second);
			 	std::cout<<"First of the doubled AIDA ion entry :"<< imbi->first << " Second of the doubled AIDA ion entry"<<imtsi->second<<std::endl; 
				 return -1;}
			}
		}
	}
	cout<<"nentries = " <<mbi.size()<<" for BigRIPS x ion"<<endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//building beam_ion merge tree//

	//for(imtsi=mtsi.begin();imtsi!=mtsi.end();imtsi++){
	for(imbi=mbi.begin();imbi!=mbi.end();imbi++){
		iEntry = imbi->first;

		if(entry%1000==0) { std::cout<<"\rThe merged ion entry with beam :"<<iEntry; std::cout.flush();}
		entry+=1;

		temp = aidaion.LoadTree(iEntry);
		if(temp<0){
			break;
		}
//		if(imtsi->first !=0){ 		//If the timestamp != 0
//			imbi=mbi.find(imtsi->second);

			if(imbi!=mbi.end()){

				aidaion.GetEntry(iEntry);
				beam.GetEntry(imbi->second);
				outtree.SyncBeamIon(beam,aidaion);
				tree2->Fill();
			}
//		}
	}
		
rootfile->WriteTObject(tree2);
rootfile->Close();
std::cout<<"Job finished "<<outfile<<" is made"<<std::endl;
}
