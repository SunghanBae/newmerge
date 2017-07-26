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
#include <string.h>

//#include "BeamGe_AIDA.h"
#include "newmerge.h"

using namespace std;

int main(int argc, char* argv[]){
	Bool_t filltree;	
	int i,test;
	char BAfile[128];
	char GAfile[128];
	char outfile[128];	
//	char gamfile[128];
//	char calfile[128];

	int runB, runA, runA_sub, runG;

	if(argc==5){
	runB = atoi(argv[1]);
	runA = atoi(argv[3]);
	runA_sub = atoi(argv[4]);
	runG = atoi(argv[2]);

		if(runB == 1){
		string baname = argv[1];
		string ganame = argv[2];
		string outname = argv[3];

		strcpy(BAfile,baname.c_str());
		strcpy(GAfile,ganame.c_str());
		strcpy(outfile,outname.c_str());
		}else{

		sprintf(BAfile,"$NP1306DIR/newmerge/BA_B%04d_A%04d_%d.root",runB,runA,runA_sub);
		sprintf(GAfile,"$NP1306DIR/newmerge/GA_G%04d_A%04d_%d.root",runG,runA,runA_sub);
		sprintf(outfile,"$NP1306DIR/newmerge/newmerge_B%04d_G%04d_A%04d_%d.root",runB,runG,runA,runA_sub);
		}

	}else if (argc==4){

        runB = atoi(argv[1]);
        runA = atoi(argv[3]);
        runG = atoi(argv[2]);

	sprintf(BAfile,"$NP1306DIR/newmerge/BA_B%04d_A%04d.root",runB,runA);
	sprintf(GAfile,"$NP1306DIR/newmerge/GA_G%04d_A%04d.root",runG,runA);
	sprintf(outfile,"$NP1306DIR/newmerge/newmerge_B%04d_G%04d_A%04d.root",runB,runG,runA);

	}else { cout<<"USAGE: ./newmerge #BigRIPSRUN #EURICARun #AIDARUN (#AIDASUBRUN)"<<endl<<"Or : ./newmerge 1 BAfile GAfile outfile"<<endl;
                        return 0;
	}


//	BigRIPS_reco beam(beamfile);
	BA beamion(BAfile);
	GA betagam(GAfile);
//	eurica gamma(gamfile);
//	Calib calib(calfile);

	map<Long64_t,Long64_t> mtsb,mtsg; //time stamp table for bigrips and aida ion and electron
	map<Long64_t,Long64_t>::iterator imtsb, imtsg;
	map<Long64_t,Long64_t> mtsi,mtse;
	map<Long64_t,Long64_t>::iterator imtsi, imtse;

	map<Long64_t,Long64_t> mbi, mge;	//time stamp correlation between ion,beta and bigrips
	map<Long64_t,Long64_t>::iterator imbi ,imge;

	multimap<pair<Long64_t,Long64_t>, pair<Long64_t,Long64_t>> mei;
	multimap<pair<Long64_t,Long64_t>, pair<Long64_t,Long64_t>>::iterator imei;
	
	cout<<"start building time map for BigRIPS, beta, ion, EURICA"<<endl;
//	aidaion.GetTsEntryGood(mtsi);
	beamion.GetTsEntry(mtsi);
	betagam.GetTsEntry(mtse);
//	beam.GetTsEntry(mtsb);
//	aidabeta.GetTsEntryGood(mtse);
//	gamma.GetTsEntry(mtsg);
	Long64_t bEntry, iEntry, eEntry, gEntry;
	
	TFile* rootfile = new TFile(outfile,"RECREATE");
	TTree* tree = new TTree("tree","tree");
//	TTree* tree2 = new TTree("ion","ion");
	
	newmerge outtree;
	outtree.TreeBranch(tree);
//	outtree.TreeBranchIon(tree2);

//	std::cout<<mtsb.size()<<" events for BigRIPS"<<std::endl;
	std::cout<<mtsi.size()<<" events for ion"<<std::endl;
	std::cout<<mtse.size()<<" events for beta events"<<std::endl;
//	std::cout<<mtsg.size()<<" events for gamma events"<<std::endl;


	// The timestamps are in 10ns/1stamp
	
//	Long64_t beam_ion_tsL = -1800; // window
//	Long64_t beam_ion_tsU = -1300; 
//	Long64_t beta_gam_tsL = -1800; 
//	Long64_t beta_gam_tsU = -1300; // window
	Long64_t beta_ion_tsU = 300000000; // 3000ms window
	Long64_t beta_ion_tsL = 20000; // 50us from ion. To exclude bad event
	
	Long64_t entry = 0;
	Long64_t temp;

	Long64_t ts1;
	Long64_t ts2;
	
	std::cout<<"Scan correlated beta with ion using timestamps and positions"<<std::endl;

	Long64_t tempi=0, temptot;
	temptot=mtsi.size();
	double time0;
	time_t tm;
	int nb;
	time0=time(&tm);
	for(imtsi=mtsi.begin();imtsi!=mtsi.end();imtsi++){
			
		if(imtsi->first>0){
		ts1 = imtsi->first +beta_ion_tsL;
		ts2 = imtsi->first +beta_ion_tsU;
		imtse= mtse.lower_bound(ts1);
		nb=0;
			beamion.GetEntry(imtsi->second);
			while(imtse!=mtse.end() && imtse->first <ts2){
				betagam.GetEntry(imtse->second);
//				aidaion.GetEntry(imtsi->second);
//				test=outtree.CorrPosition(betagam, beamion);
				test=outtree.CorrPixel(betagam, beamion);
				if (test==1)
				{
					outtree.SyncGA(betagam,nb);	
					nb++;
				}
				imtse++;
			}
		if(nb>0){
		outtree.SyncBA(beamion,nb);
		tree->Fill();
		}
		outtree.ResetGA();		
		outtree.ResetBA();		
		}

	if(tempi%1000==0) {std::cout<<fixed<<setprecision(4)<<"\r scanning "<<100.0*tempi/temptot<<" \% done. Timespent "<<time(&tm)-time0<<"s"; std::cout.flush();}
	tempi++;
	}
	
rootfile->WriteTObject(tree);
//rootfile->WriteTObject(tree2);
rootfile->Close();
std::cout<<"Job finished "<<outfile<<" is made"<<std::endl;
}

