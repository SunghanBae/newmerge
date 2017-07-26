//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 23 20:19:56 2017 by ROOT version 6.04/10
// from TTree ion/ion
// found on file: R1208_1_packed_eventbuild.root
//////////////////////////////////////////////////////////

#ifndef AIDA_ionevt_h
#define AIDA_ionevt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "AIDA_betaevt.h"
// Header file for the classes stored in the TTree if any.

class AIDA_ionevt {
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
   Double_t        stop_x;
   Double_t        stop_y;
   Int_t           stop_z;
   Int_t           goodness;
   Int_t           fire_flag[6][2];
   Double_t        EX[6];
   Double_t        EY[6];
   Double_t        X[6];
   Double_t        Y[6];
   Int_t           hitX[6];
   Int_t           hitY[6];
   Int_t           asicNo[96];   //[mulasic]
   Int_t           dssdNo[96];   //[mulasic]
   Int_t           multi[96];   //[mulasic]
   Long64_t        timestamp[96][16];   //[mulasic]
   Long64_t        extTimestamp[96][16];   //[mulasic]
   Int_t           stripNo[96][16];   //[mulasic]
   Int_t           adcData[96][16];   //[mulasic]
   Int_t           rangeType[96][16];   //[mulasic]
   Double_t        adcE[96][16];   //[mulasic]

   // List of branches
   TBranch        *b_tstart;   //!
   TBranch        *b_tstop;   //!
   TBranch        *b_extTstart;   //!
   TBranch        *b_extTstop;   //!
   TBranch        *b_t_del;   //!
   TBranch        *b_extT_del;   //!
   TBranch        *b_mulasic;   //!
   TBranch        *b_multitot;   //!
   TBranch        *b_stop_x;   //!
   TBranch        *b_stop_y;   //!
   TBranch        *b_stop_z;   //!
   TBranch        *b_goodness;   //!
   TBranch        *b_fire_flag;   //!
   TBranch        *b_EX;   //!
   TBranch        *b_EY;   //!
   TBranch        *b_X;   //!
   TBranch        *b_Y;   //!
   TBranch        *b_hitX;   //!
   TBranch        *b_hitY;   //!
   TBranch        *b_asicNo;   //!
   TBranch        *b_dssdNo;   //!
   TBranch        *b_multi;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_extTimestamp;   //!
   TBranch        *b_stripNo;   //!
   TBranch        *b_adcData;   //!
   TBranch        *b_rangeType;   //!
   TBranch        *b_adcE;   //!

   AIDA_ionevt(char* filename, TTree *tree=0);
   virtual ~AIDA_ionevt();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     GetTsEntry(std::map <Long64_t, Long64_t> &mts);
   virtual void     GetTree(char *filename, TTree *tree);
   virtual void     GetGoodness(std::map<Long64_t,Long64_t> &mts);
   virtual int	    GetEntriesFast();
   virtual int   CorrPosition(AIDA_betaevt &beta, AIDA_ionevt &ion);
};

#endif


AIDA_ionevt::AIDA_ionevt(char* filename, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      f->GetObject("ion",tree);

   }
   Init(tree);
}

AIDA_ionevt::~AIDA_ionevt()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void AIDA_ionevt::GetTsEntry(std::map<Long64_t,Long64_t> &mts)
{  Long64_t ts;
   Long64_t nenties = fChain->GetEntriesFast();
   std::cout<<nenties<<"nenties in AIDAion "<<std::endl;
   Long64_t nbytes = 0, nb =0;
   for(Long64_t jentry = 0;jentry < nenties ;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if(ientry<0) break;
      nb = fChain->GetEntry(jentry); nbytes += nb;
      ts = extTstop*4;
      if(jentry%10000 ==0 && ts > 0) std::cout<<"jentry = "<<jentry<<" ion ts = "<<ts<<std::endl;
      if(ts != 0) mts.insert(std::pair<Long64_t,Long64_t> (ts,jentry));
   }
}

int AIDA_ionevt::GetEntriesFast()
{
 return fChain->GetEntriesFast();
}


void AIDA_ionevt::GetGoodness(std::map<Long64_t,Long64_t> &mts)
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

void AIDA_ionevt::GetTree(char* filename, TTree *tree)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      f->GetObject("ion",tree);

   }
   Init(tree);
}

Int_t AIDA_ionevt::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AIDA_ionevt::LoadTree(Long64_t entry)
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

void AIDA_ionevt::Init(TTree *tree)
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
   fChain->SetBranchAddress("stop_x", &stop_x, &b_stop_x);
   fChain->SetBranchAddress("stop_y", &stop_y, &b_stop_y);
   fChain->SetBranchAddress("stop_z", &stop_z, &b_stop_z);
   fChain->SetBranchAddress("goodness", &goodness, &b_goodness);
   fChain->SetBranchAddress("fire_flag", fire_flag, &b_fire_flag);
   fChain->SetBranchAddress("EX", EX, &b_EX);
   fChain->SetBranchAddress("EY", EY, &b_EY);
   fChain->SetBranchAddress("X", X, &b_X);
   fChain->SetBranchAddress("Y", Y, &b_Y);
   fChain->SetBranchAddress("hitX", hitX, &b_hitX);
   fChain->SetBranchAddress("hitY", hitY, &b_hitY);
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

Bool_t AIDA_ionevt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AIDA_ionevt::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AIDA_ionevt::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

int AIDA_ionevt::CorrPosition(AIDA_betaevt &beta,AIDA_ionevt &ion)
{
//     if(beta.start_z>-1){
     if(abs(beta.start_y-ion.stop_y)< 3 && ion.stop_y >-1 && beta.start_y>-1){
	 if( abs(beta.start_x-ion.stop_x) <3 && ion.stop_x>-1 && beta.start_x>-1){
	     if( abs(beta.start_z-ion.stop_z) < 2 && ion.stop_z>-1){
	              return 1;
              }else{ return 0;}
     }else{return 0;}
    }else{return 0;}
//	}else{return 0;}
}




