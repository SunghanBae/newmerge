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

	if(argc==4){
	
	runA = atoi(argv[2]);
	runA_sub = atoi(argv[3]);
	runG = atoi(argv[1]);


	sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
	sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
	sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_%d_packed_eventbuild.root",runA,runA_sub);
	sprintf(outfile,"$NP1306DIR/newmerge/GA_G%04d_A%04d_%d.root",runG,runA,runA_sub);
	}else if (argc==3){

        runA = atoi(argv[2]);
        runG = atoi(argv[1]);

        sprintf(gamfile,"$NP1306DIR/EURICA_exp/run_%04d.root",runG);
        sprintf(calfile,"$NP1306DIR/Beam_Gamma/152Eu_eff_1000Hz_0001_calibdata.root");
        sprintf(AIDAfile,"$NP1306DIR/AIDA/R%04d_packed_eventbuild.root",runA);
        sprintf(outfile,"$NP1306DIR/newmerge/GA_G%04d_A%04d.root",runG,runA);

	}else { cout<<"USAGE: ./Beta_gamma #EURICARun #AIDARUN (#AIDASUBRUN)"<<endl;
                        return 0;
	}


	AIDA_betaevt aidabeta(AIDAfile);
	eurica gamma(gamfile);
	Calib calib(calfile);

	map<Long64_t,Long64_t> mtsb,mtsg; //time stamp table for bigrips and aida ion and electron
	map<Long64_t,Long64_t>::iterator imtsb, imtsg;
	map<Long64_t,Long64_t> mtsi,mtse;
	map<Long64_t,Long64_t>::iterator imtsi, imtse;

	multimap<Long64_t,Long64_t> mbi, mge, mei;	//time stamp correlation between ion,beta and bigrips
	multimap<Long64_t,Long64_t>::iterator imbi ,imge ,imei;
	
	cout<<"start building time map for BigRIPS, beta, ion, EURICA"<<endl;

	aidabeta.GetTsEntry(mtse);
	gamma.GetTsEntry(mtsg);
	
	Long64_t bEntry, iEntry, eEntry, gEntry;
	
	TFile* rootfile = new TFile(outfile,"RECREATE");
//	TTree* tree = new TTree("tree","tree");
	TTree* tree2 = new TTree("beta_gam","beta_gam");
	
	BeamGe_AIDA outtree;
	outtree.TreeBranchBetaGam(tree2);

	std::cout<<mtse.size()<<" events for beta events"<<std::endl;
	std::cout<<mtsg.size()<<" events for gamma events"<<std::endl;


	// The timestamps are in timestamp unit.	10ns/1stamp
	
	Long64_t beta_gam_tsL = -1800; 
	Long64_t beta_gam_tsU = -1300; // 5us window
	
	Long64_t entry = 0;
	Long64_t temp;

	Long64_t ts1;
	Long64_t ts2;
	
	std::cout<<"Scan coincidence timestamps beta with Gamma, beta with ion" << endl;
	for(imtse=mtse.begin();imtse!=mtse.end();imtse++){
		if(imtse->first >0){
			ts1 = imtse->first + beta_gam_tsL;
			ts2 = imtse->first + beta_gam_tsU;
			imtsg = mtsg.lower_bound(ts1);
			while(imtsg!=mtsg.end() && (imtsg->first) <ts2){
				mge.insert(std::make_pair(imtse->second,imtsg->second));
				imtsg++;
			}
		}
	}
	cout<<"nentries = " <<mge.size()<<" for beta x Gamma"<<endl;
	
		size_t gamma_n;
		outtree.LoadCalib(calib);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	for(imtse=mtse.begin();imtse!=mtse.end();imtse++){
	imge=mge.begin();
	while(imge!=mge.end()){
//		eEntry = imtse->second;
		eEntry = imge->first;
		
		if(entry%10000==0) { std::cout<<"\rThe merged beta entry with gamma :"<<eEntry; std::cout.flush();}
		entry+=1;

		temp = aidabeta.LoadTree(eEntry);
		if(temp<0){
			break;
		}
		
		//if(imtse->first >0){

			gamma_n = mge.count(eEntry);
		//if(mge.count(eEntry)>0) {std::cout<<mge.count(eEntry)<<"ngamma"<<std::endl;}
			outtree.SyncMulti(gamma_n,0);
		//	imge=mge.find(eEntry);

			i=0;
			while(i<gamma_n && imge!=mge.end() && imge->first==eEntry){
				gamma.GetEntry(imge->second);
				outtree.CalGamma(gamma,i);
				outtree.SyncGammaBeta(gamma,i);
				imge++;
				i++;
			}

			aidabeta.GetEntry(eEntry);
			outtree.SyncBeta(aidabeta);
			tree2->Fill();

			outtree.ResetGe();
	//	}
	}
rootfile->WriteTObject(tree2);
rootfile->Close();
std::cout<<std::endl<<"Job finished "<<outfile<<" is made"<<std::endl;
}
