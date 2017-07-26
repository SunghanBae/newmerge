//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  3 11:37:42 2017 by ROOT version 6.10/00
// from TTree ion/ion
// found on file: BA_B1193_A1208_1.root
//////////////////////////////////////////////////////////

#ifndef BA_h
#define BA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        ion_ts;
   Long64_t        ion_bigrips_ts;
   Double_t        bigrips_BeamZ;
   Double_t        bigrips_BeamAOQ;
   Double_t        ion_x;
   Double_t        ion_y;
   Double_t        ion_z;
   Double_t        ion_EX;
   Double_t        ion_EY;
   Int_t	   ion_goodness;
   Int_t           bigrips_F11PPAC1X1;
   Int_t           bigrips_F11PPAC1X2;
   Int_t           bigrips_F11PPAC2X1;
   Int_t           bigrips_F11PPAC2X2;
   Int_t           bigrips_F5PPAC1X1;
   Int_t           bigrips_F5PPAC1X2;
   Int_t           bigrips_F5PPAC2X1;
   Int_t           bigrips_F5PPAC2X2;
   Int_t           bigrips_F3PPAC1X1;
   Int_t           bigrips_F3PPAC1X2;
   Int_t           bigrips_F3PPAC2X1;
   Int_t           bigrips_F3PPAC2X2;

   // List of branches
   TBranch        *b_ion_ts;   //!
   TBranch        *b_ion_bigrips_ts;   //!
   TBranch        *b_bigrips_BeamZ;   //!
   TBranch        *b_bigrips_BeamAOQ;   //!
   TBranch        *b_ion_x;   //!
   TBranch        *b_ion_y;   //!
   TBranch        *b_ion_z;   //!
   TBranch        *b_ion_EX;   //!
   TBranch        *b_ion_EY;   //!
   TBranch        *b_ion_goodness;   //!
   TBranch        *b_bigrips_F11PPAC1X1;   //!
   TBranch        *b_bigrips_F11PPAC1X2;   //!
   TBranch        *b_bigrips_F11PPAC2X1;   //!
   TBranch        *b_bigrips_F11PPAC2X2;   //!
   TBranch        *b_bigrips_F5PPAC1X1;   //!
   TBranch        *b_bigrips_F5PPAC1X2;   //!
   TBranch        *b_bigrips_F5PPAC2X1;   //!
   TBranch        *b_bigrips_F5PPAC2X2;   //!
   TBranch        *b_bigrips_F3PPAC1X1;   //!
   TBranch        *b_bigrips_F3PPAC1X2;   //!
   TBranch        *b_bigrips_F3PPAC2X1;   //!
   TBranch        *b_bigrips_F3PPAC2X2;   //!

   BA(char* filename, TTree *tree=0);
   virtual ~BA();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     GetTsEntry(std::map <Long64_t, Long64_t> &mts);

};

#endif

//#ifdef BA_cxx
BA::BA(char* filename, TTree *tree) : fChain(0) 
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

BA::~BA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void BA::GetTsEntry(std::map<Long64_t,Long64_t> &mts)
{  Long64_t ts;
   Long64_t nenties = fChain->GetEntriesFast();
   std::cout<<nenties<<"nenties in BA "<<std::endl;
   Long64_t nbytes = 0, nb =0;
   for(Long64_t jentry = 0;jentry < nenties ;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if(ientry<0) break;
      nb = fChain->GetEntry(jentry); nbytes += nb;
      ts = ion_ts;
      if(jentry%10000 ==0 && ts > 0) std::cout<<"jentry = "<<jentry<<" ion ts = "<<ts<<std::endl;
      if(ts != 0) mts.insert(std::pair<Long64_t,Long64_t> (ts,jentry));
   }
}

Int_t BA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BA::LoadTree(Long64_t entry)
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

void BA::Init(TTree *tree)
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

   fChain->SetBranchAddress("ion_ts", &ion_ts, &b_ion_ts);
   fChain->SetBranchAddress("ion_bigrips_ts", &ion_bigrips_ts, &b_ion_bigrips_ts);
   fChain->SetBranchAddress("bigrips_BeamZ", &bigrips_BeamZ, &b_bigrips_BeamZ);
   fChain->SetBranchAddress("bigrips_BeamAOQ", &bigrips_BeamAOQ, &b_bigrips_BeamAOQ);
   fChain->SetBranchAddress("ion_x", &ion_x, &b_ion_x);
   fChain->SetBranchAddress("ion_y", &ion_y, &b_ion_y);
   fChain->SetBranchAddress("ion_z", &ion_z, &b_ion_z);
   fChain->SetBranchAddress("ion_EX", &ion_EX, &b_ion_EX);
   fChain->SetBranchAddress("ion_EY", &ion_EY, &b_ion_EY);
   fChain->SetBranchAddress("ion_EY", &ion_goodness, &b_ion_goodness);
   fChain->SetBranchAddress("bigrips_F11PPAC1X1", &bigrips_F11PPAC1X1, &b_bigrips_F11PPAC1X1);
   fChain->SetBranchAddress("bigrips_F11PPAC1X2", &bigrips_F11PPAC1X2, &b_bigrips_F11PPAC1X2);
   fChain->SetBranchAddress("bigrips_F11PPAC2X1", &bigrips_F11PPAC2X1, &b_bigrips_F11PPAC2X1);
   fChain->SetBranchAddress("bigrips_F11PPAC2X2", &bigrips_F11PPAC2X2, &b_bigrips_F11PPAC2X2);
   fChain->SetBranchAddress("bigrips_F5PPAC1X1", &bigrips_F5PPAC1X1, &b_bigrips_F5PPAC1X1);
   fChain->SetBranchAddress("bigrips_F5PPAC1X2", &bigrips_F5PPAC1X2, &b_bigrips_F5PPAC1X2);
   fChain->SetBranchAddress("bigrips_F5PPAC2X1", &bigrips_F5PPAC2X1, &b_bigrips_F5PPAC2X1);
   fChain->SetBranchAddress("bigrips_F5PPAC2X2", &bigrips_F5PPAC2X2, &b_bigrips_F5PPAC2X2);
   fChain->SetBranchAddress("bigrips_F3PPAC1X1", &bigrips_F3PPAC1X1, &b_bigrips_F3PPAC1X1);
   fChain->SetBranchAddress("bigrips_F3PPAC1X2", &bigrips_F3PPAC1X2, &b_bigrips_F3PPAC1X2);
   fChain->SetBranchAddress("bigrips_F3PPAC2X1", &bigrips_F3PPAC2X1, &b_bigrips_F3PPAC2X1);
   fChain->SetBranchAddress("bigrips_F3PPAC2X2", &bigrips_F3PPAC2X2, &b_bigrips_F3PPAC2X2);
   Notify();
}

Bool_t BA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BA::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef BA_cxx
