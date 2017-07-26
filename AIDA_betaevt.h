//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 21 12:47:54 2017 by ROOT version 6.10/00
// from TTree beta/beta
// found on file: ../AIDA/R1208_0_packed_eventbuild.root
//////////////////////////////////////////////////////////

#ifndef AIDA_betaevt_h
#define AIDA_betaevt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class AIDA_betaevt {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		Long64_t        tstart;
		Long64_t        tstop;
		Long64_t        extTstart;
		Long64_t        extTstop;
		Long64_t        t_del;
		Long64_t        extT_del;
		Int_t           mulasic;
		Int_t           multitot;
		Double_t        start_x;
		Double_t        start_y;
		Int_t           start_z;
		Int_t           flag_pulser;
		Int_t           goodness;
		Int_t           fire_flag[6][2];
		Double_t        EX[6];
		Double_t        EY[6];
		Double_t        X[6];
		Double_t        Y[6];
		Int_t           hitX[6];
		Int_t           hitY[6];
		Int_t           mulpix;
		Int_t           pixelx[900];   //[mulpix]
		Int_t           pixely[900];   //[mulpix]
		Int_t           pixelz[900];   //[mulpix]
		Double_t        pixelEX[900];   //[mulpix]
		Double_t        pixelEY[900];   //[mulpix]
		Int_t           asicNo[360];   //[mulasic]
		Int_t           dssdNo[360];   //[mulasic]
		Int_t           multi[360];   //[mulasic]
		Long64_t        timestamp[360][16];   //[mulasic]
		Long64_t        extTimestamp[360][16];   //[mulasic]
		Int_t           stripNo[360][16];   //[mulasic]
		Int_t           adcData[360][16];   //[mulasic]
		Int_t           rangeType[360][16];   //[mulasic]
		Double_t        adcE[360][16];   //[mulasic]

		// List of branches
		TBranch        *b_tstart;   //!
		TBranch        *b_tstop;   //!
		TBranch        *b_extTstart;   //!
		TBranch        *b_extTstop;   //!
		TBranch        *b_t_del;   //!
		TBranch        *b_extT_del;   //!
		TBranch        *b_mulasic;   //!
		TBranch        *b_multitot;   //!
		TBranch        *b_start_x;   //!
		TBranch        *b_start_y;   //!
		TBranch        *b_start_z;   //!
		TBranch        *b_flag_pulser;   //!
		TBranch        *b_goodness;   //!
		TBranch        *b_fire_flag;   //!
		TBranch        *b_EX;   //!
		TBranch        *b_EY;   //!
		TBranch        *b_X;   //!
		TBranch        *b_Y;   //!
		TBranch        *b_hitX;   //!
		TBranch        *b_hitY;   //!
		TBranch        *b_mulpix;   //!
		TBranch        *b_pixelx;   //!
		TBranch        *b_pixely;   //!
		TBranch        *b_pixelz;   //!
		TBranch        *b_pixelEX;   //!
		TBranch        *b_pixelEY;   //!
		TBranch        *b_asicNo;   //!
		TBranch        *b_dssdNo;   //!
		TBranch        *b_multi;   //!
		TBranch        *b_timestamp;   //!
		TBranch        *b_extTimestamp;   //!
		TBranch        *b_stripNo;   //!
		TBranch        *b_adcData;   //!
		TBranch        *b_rangeType;   //!
		TBranch        *b_adcE;   //!

		AIDA_betaevt(char* filename, TTree *tree=0);
		virtual ~AIDA_betaevt();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);
		virtual void     GetTsEntry(std::map <Long64_t, Long64_t> &mts);
		virtual void     GetTree(char *filename, TTree *tree);
		virtual void     GetGoodness(std::map<Long64_t,Long64_t> &mts);
		virtual int      GetEntriesFast();
		//   virtual int      CorrPosition(AIDA_betaevt &beta, AIDA_betaevt &ion);

};

#endif

AIDA_betaevt::AIDA_betaevt(char * filename, TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
		if (!f || !f->IsOpen()) {
			f = new TFile(filename);
		}
		f->GetObject("beta",tree);

	}
	Init(tree);
}

AIDA_betaevt::~AIDA_betaevt()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

void AIDA_betaevt::GetTsEntry(std::map<Long64_t,Long64_t> &mts)
{  Long64_t ts; 
	Long64_t nenties = fChain->GetEntriesFast();
	std::cout<<nenties<<"nenties in AIDAbeta "<<std::endl;
	Long64_t nbytes = 0, nb =0; 
	for(Long64_t jentry = 0;jentry < nenties ;jentry++)
	{   
		Long64_t ientry = LoadTree(jentry);
		if(ientry<0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb; 
		ts = extTstart*4;
		if(jentry%10000 ==0 && ts > 0) std::cout<<"jentry = "<<jentry<<" beta ts = "<<ts<<std::endl;
		if(ts != 0) mts.insert(std::pair<Long64_t,Long64_t> (ts,jentry));
	}   
}

int AIDA_betaevt::GetEntriesFast()
{
	return fChain->GetEntriesFast();
}

void AIDA_betaevt::GetGoodness(std::map<Long64_t,Long64_t> &mts)
{  Long64_t good;
	Long64_t nenties = fChain->GetEntriesFast();
	//std::cout<<nenties<<"nenties in AIDAion "<<std::endl;
	Long64_t nbytes = 0, nb =0; 
	for(Long64_t jentry = 0;jentry < nenties ;jentry++)
	{   
		Long64_t ientry = LoadTree(jentry);
		if(ientry<0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb; 
		good = goodness;
		//    if(jentry%10000 ==0 && ts > 0) std::cout<<"jentry = "<<jentry<<" ion ts = "<<ts<<std::endl;
		mts.insert(std::pair<Long64_t,Long64_t> (jentry,good));
	}
}

void AIDA_betaevt::GetTree(char* filename, TTree *tree)
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
		if (!f || !f->IsOpen()) {
			f = new TFile(filename);
		}
		f->GetObject("beta",tree);

	}
	Init(tree);
}


Int_t AIDA_betaevt::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t AIDA_betaevt::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void AIDA_betaevt::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("tstart", &tstart, &b_tstart);
	fChain->SetBranchAddress("tstop", &tstop, &b_tstop);
	fChain->SetBranchAddress("extTstart", &extTstart, &b_extTstart);
	fChain->SetBranchAddress("extTstop", &extTstop, &b_extTstop);
	fChain->SetBranchAddress("t_del", &t_del, &b_t_del);
	fChain->SetBranchAddress("extT_del", &extT_del, &b_extT_del);
	fChain->SetBranchAddress("mulasic", &mulasic, &b_mulasic);
	fChain->SetBranchAddress("multitot", &multitot, &b_multitot);
	fChain->SetBranchAddress("start_x", &start_x, &b_start_x);
	fChain->SetBranchAddress("start_y", &start_y, &b_start_y);
	fChain->SetBranchAddress("start_z", &start_z, &b_start_z);
	fChain->SetBranchAddress("flag_pulser", &flag_pulser, &b_flag_pulser);
	fChain->SetBranchAddress("goodness", &goodness, &b_goodness);
	fChain->SetBranchAddress("fire_flag", fire_flag, &b_fire_flag);
	fChain->SetBranchAddress("EX", EX, &b_EX);
	fChain->SetBranchAddress("EY", EY, &b_EY);
	fChain->SetBranchAddress("X", X, &b_X);
	fChain->SetBranchAddress("Y", Y, &b_Y);
	fChain->SetBranchAddress("hitX", hitX, &b_hitX);
	fChain->SetBranchAddress("hitY", hitY, &b_hitY);
	fChain->SetBranchAddress("mulpix", &mulpix, &b_mulpix);
	fChain->SetBranchAddress("pixelx", pixelx, &b_pixelx);
	fChain->SetBranchAddress("pixely", pixely, &b_pixely);
	fChain->SetBranchAddress("pixelz", pixelz, &b_pixelz);
	fChain->SetBranchAddress("pixelEX", pixelEX, &b_pixelEX);
	fChain->SetBranchAddress("pixelEY", pixelEY, &b_pixelEY);
	fChain->SetBranchAddress("asicNo", asicNo, &b_asicNo);
	fChain->SetBranchAddress("dssdNo", dssdNo, &b_dssdNo);
	fChain->SetBranchAddress("multi", multi, &b_multi);
	fChain->SetBranchAddress("timestamp", timestamp, &b_timestamp);
	fChain->SetBranchAddress("extTimestamp", extTimestamp, &b_extTimestamp);
	fChain->SetBranchAddress("stripNo", stripNo, &b_stripNo);
	fChain->SetBranchAddress("adcData", adcData, &b_adcData);
	fChain->SetBranchAddress("rangeType", rangeType, &b_rangeType);
	fChain->SetBranchAddress("adcE", adcE, &b_adcE);
	Notify();
}

Bool_t AIDA_betaevt::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void AIDA_betaevt::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t AIDA_betaevt::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
