//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 21 16:42:24 2017 by ROOT version 6.10/00
// from TTree beta_gam/beta_gam
// found on file: GA_G1408_A1208_1.root
//////////////////////////////////////////////////////////

#ifndef GA_h
#define GA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any
const size_t kMaxPix = 1000;

class GA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        beta_ts;
   Double_t        beta_x;
   Double_t        beta_y;
   Double_t        beta_z;
   Double_t        beta_EX;
   Double_t        beta_EY;
   Int_t           beta_goodness;
   Int_t           beta_mulpix;
   Int_t           beta_pixelx[kMaxPix];   //[beta_mulpix]
   Int_t           beta_pixely[kMaxPix];   //[beta_mulpix]
   Int_t           beta_pixelz[kMaxPix];   //[beta_mulpix]
   Double_t        beta_pixelEX[kMaxPix];   //[beta_mulpix]
   Double_t        beta_pixelEY[kMaxPix];   //[beta_mulpix]
   Int_t           n_gamma;
   Long64_t        gc_ts[1];   //[n_gamma]
   Int_t           gc_hit[1];   //[n_gamma]
   Int_t           gc_ch[1][46];   //[n_gamma]
   Int_t           gc_DGFt[1][46];   //[n_gamma]
   Double_t        gc_Ecal[1][46];   //[n_gamma]
   Int_t           ab_hit[1];   //[n_gamma]
   Int_t           ab_DGFt[1][27];   //[n_gamma]
   Double_t        ab_E[1][27];   //[n_gamma]

   // List of branches
   TBranch        *b_beta_ts;   //!
   TBranch        *b_beta_x;   //!
   TBranch        *b_beta_y;   //!
   TBranch        *b_beta_z;   //!
   TBranch        *b_beta_EX;   //!
   TBranch        *b_beta_EY;   //!
   TBranch        *b_beta_goodness;   //!
   TBranch        *b_beta_mulpix;   //!
   TBranch        *b_beta_pixelx;   //!
   TBranch        *b_beta_pixely;   //!
   TBranch        *b_beta_pixelz;   //!
   TBranch        *b_beta_pixelEX;   //!
   TBranch        *b_beta_pixelEY;   //!
   TBranch        *b_n_gamma;   //!
   TBranch        *b_gc_ts;   //!
   TBranch        *b_gc_hit;   //!
   TBranch        *b_gc_ch;   //!
   TBranch        *b_gc_DGFt;   //!
   TBranch        *b_gc_Ecal;   //!
   TBranch        *b_ab_hit;   //!
   TBranch        *b_ab_DGFt;   //!
   TBranch        *b_ab_E;   //!

   GA(char* filename, TTree *tree=0);
   virtual ~GA();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     GetTsEntry(std::map <Long64_t, Long64_t> &mts);
};

#endif

//#ifdef GA_cxx
GA::GA(char* filename, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      f->GetObject("beta_gam",tree);

   }
   Init(tree);
}

GA::~GA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void GA::GetTsEntry(std::map<Long64_t,Long64_t> &mts)
{  Long64_t ts;
   Long64_t nenties = fChain->GetEntriesFast();
   std::cout<<nenties<<"nenties in GA` "<<std::endl;
   Long64_t nbytes = 0, nb =0;
   for(Long64_t jentry = 0;jentry < nenties ;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if(ientry<0) break;
      nb = fChain->GetEntry(jentry); nbytes += nb;
   //   if(start_z>-1){
      ts = beta_ts;
      if(ts != 0) mts.insert(std::pair<Long64_t,Long64_t> (ts,jentry));
     // }
      if(jentry%100000 ==0 && ts > 0) std::cout<<"jentry = "<<jentry<<" beta ts = "<<ts<<std::endl;
   }
}

Int_t GA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GA::LoadTree(Long64_t entry)
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

void GA::Init(TTree *tree)
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

   fChain->SetBranchAddress("beta_ts", &beta_ts, &b_beta_ts);
   fChain->SetBranchAddress("beta_x", &beta_x, &b_beta_x);
   fChain->SetBranchAddress("beta_y", &beta_y, &b_beta_y);
   fChain->SetBranchAddress("beta_z", &beta_z, &b_beta_z);
   fChain->SetBranchAddress("beta_EX", &beta_EX, &b_beta_EX);
   fChain->SetBranchAddress("beta_EY", &beta_EY, &b_beta_EY);
   fChain->SetBranchAddress("beta_goodness", &beta_goodness, &b_beta_goodness);
   fChain->SetBranchAddress("beta_mulpix", &beta_mulpix, &b_beta_mulpix);
   fChain->SetBranchAddress("beta_pixelx", beta_pixelx, &b_beta_pixelx);
   fChain->SetBranchAddress("beta_pixely", beta_pixely, &b_beta_pixely);
   fChain->SetBranchAddress("beta_pixelz", beta_pixelz, &b_beta_pixelz);
   fChain->SetBranchAddress("beta_pixelEX", beta_pixelEX, &b_beta_pixelEX);
   fChain->SetBranchAddress("beta_pixelEY", beta_pixelEY, &b_beta_pixelEY);
   fChain->SetBranchAddress("n_gamma", &n_gamma, &b_n_gamma);
   fChain->SetBranchAddress("gc_ts", gc_ts, &b_gc_ts);
   fChain->SetBranchAddress("gc_hit", gc_hit, &b_gc_hit);
   fChain->SetBranchAddress("gc_ch", gc_ch, &b_gc_ch);
   fChain->SetBranchAddress("gc_DGFt", gc_DGFt, &b_gc_DGFt);
   fChain->SetBranchAddress("gc_Ecal", gc_Ecal, &b_gc_Ecal);
   fChain->SetBranchAddress("ab_hit", ab_hit, &b_ab_hit);
   fChain->SetBranchAddress("ab_DGFt", ab_DGFt, &b_ab_DGFt);
   fChain->SetBranchAddress("ab_E", ab_E, &b_ab_E);
   Notify();
}

Bool_t GA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GA::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef GA_cxx
