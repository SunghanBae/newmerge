#ifndef BeamGe_AIDA_h
#define BeamGe_AIDA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"
//#include "eurica.h"
#include "BigRIPS_reco.h"
#include "eurica.h"
#include "AIDA_ionevt.h"
//#include "AIDA_betaevt.h"
#include "AIDAraw.h"
//#include "fixedpars.h"
#include "Calib.h"

class BeamGe_AIDA {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain
	
		size_t		n_ion_bi;
		Long64_t	ion_ts_bi;
		Long64_t	ion_bigrips_ts_bi;
		Double_t	bigrips_BeamZ_bi;
		Double_t	bigrips_BeamAOQ_bi;
		Double_t	ion_x_bi;
		Double_t	ion_y_bi;
		Double_t	ion_z_bi;
		Double_t	ion_EX_bi;
		Double_t	ion_EY_bi;
		Int_t		ion_goodness_bi;

                Int_t bigrips_F3PPAC1X1_bi;
                Int_t bigrips_F3PPAC1X2_bi;
                Int_t bigrips_F3PPAC2X1_bi;
                Int_t bigrips_F3PPAC2X2_bi;
                Int_t bigrips_F5PPAC1X1_bi;
                Int_t bigrips_F5PPAC1X2_bi;
                Int_t bigrips_F5PPAC2X1_bi;
                Int_t bigrips_F5PPAC2X2_bi;
                Int_t bigrips_F11PPAC1X1_bi;
                Int_t bigrips_F11PPAC1X2_bi;
                Int_t bigrips_F11PPAC2X1_bi;
                Int_t bigrips_F11PPAC2X2_bi;


		// Ion related parameters (BigRIPS + AIDA)
		size_t		n_ion;
		Long64_t	ion_ts[kMax_n_ion];
		Long64_t	ion_bigrips_ts[kMax_n_ion];
		Double_t	bigrips_BeamZ[kMax_n_ion];
		Double_t	bigrips_BeamAOQ[kMax_n_ion];
		Double_t	ion_x[kMax_n_ion];
		Double_t	ion_y[kMax_n_ion];
		Double_t	ion_z[kMax_n_ion];
		Double_t	ion_EX[kMax_n_ion];
		Double_t	ion_EY[kMax_n_ion];
		Int_t		ion_goodness[kMax_n_ion];

                Int_t bigrips_F3PPAC1X1[kMax_n_ion];
                Int_t bigrips_F3PPAC1X2[kMax_n_ion];
                Int_t bigrips_F3PPAC2X1[kMax_n_ion];
                Int_t bigrips_F3PPAC2X2[kMax_n_ion];
                Int_t bigrips_F5PPAC1X1[kMax_n_ion];
                Int_t bigrips_F5PPAC1X2[kMax_n_ion];
                Int_t bigrips_F5PPAC2X1[kMax_n_ion];
                Int_t bigrips_F5PPAC2X2[kMax_n_ion];
                Int_t bigrips_F11PPAC1X1[kMax_n_ion];
                Int_t bigrips_F11PPAC1X2[kMax_n_ion];
                Int_t bigrips_F11PPAC2X1[kMax_n_ion];
                Int_t bigrips_F11PPAC2X2[kMax_n_ion];


		// AIDA beta parameters
		Long64_t	beta_ts;
		Double_t	beta_EX;
		Double_t	beta_EY;
		Double_t	beta_x;
		Double_t	beta_y;
		Double_t	beta_z;
		Int_t		beta_goodness;
		
		Int_t           beta_mulpix;
		Int_t           beta_pixelx[900];   //[mulpix]
		Int_t           beta_pixely[900];   //[mulpix]
		Int_t           beta_pixelz[900];   //[mulpix]
		Double_t        beta_pixelEX[900];   //[mulpix]
		Double_t        beta_pixelEY[900];   //[mulpix]


//std::vector<ROOT::Math::XYZVector>* beta_hitpixel;
		
		// Eurica parameters
		size_t 		n_gamma;
		Long64_t	gc_ts[kMax_n_gamma];
		Int_t 		gc_hit[kMax_n_gamma];
		Int_t       gc_ch[kMax_n_gamma][kMaxGeCluster];   //[GeCluster_]
		UInt_t      gc_DGFt[kMax_n_gamma][kMaxGeCluster];   //[GeCluster_]
		UInt_t      gc_DGFe[kMaxGeCluster];   //[GeCluster_]
		Double_t	gc_Ecal[kMax_n_gamma][kMaxGeCluster];
		Int_t       ab_hit[kMax_n_gamma];
		Int_t       ab_ch[kMax_n_gamma][kMaxGeAddback];   //[GeAddback_]
		UInt_t      ab_DGFt[kMax_n_gamma][kMaxGeAddback];   //[GeAddback_]
		Double_t    ab_E[kMax_n_gamma][kMaxGeAddback];   //[GeAddback_]

		Double_t	slope[85];
		Double_t	offset[85];



		// Ion branches
		TBranch		*b_n_ion;
		TBranch		*b_ion_ts;
		TBranch		*b_ion_bigrips_ts;
		TBranch		*b_bigrips_BeamZ;
		TBranch		*b_bigrips_BeamAOQ;
		TBranch		*b_ion_x;
		TBranch		*b_ion_y;
		TBranch		*b_ion_z;
		TBranch		*b_ion_EX;
		TBranch		*b_ion_EY;
		TBranch		*b_ion_goodness;
                TBranch         *b_bigrips_F3PPAC1X1;
                TBranch         *b_bigrips_F3PPAC1X2;
                TBranch         *b_bigrips_F3PPAC2X1;
                TBranch         *b_bigrips_F3PPAC2X2;
                TBranch         *b_bigrips_F5PPAC1X1;
                TBranch         *b_bigrips_F5PPAC1X2;
                TBranch         *b_bigrips_F5PPAC2X1;
                TBranch         *b_bigrips_F5PPAC2X2;
                TBranch         *b_bigrips_F11PPAC1X1;
                TBranch         *b_bigrips_F11PPAC1X2;
                TBranch         *b_bigrips_F11PPAC2X1;
                TBranch         *b_bigrips_F11PPAC2X2;


		// AIDA beta branches
		TBranch		*b_beta_ts;
		TBranch		*b_beta_EX;
		TBranch		*b_beta_EY;
		TBranch		*b_beta_x;
		TBranch		*b_beta_y;
		TBranch		*b_beta_z;
		TBranch		*b_beta_goodness;
		TBranch         *b_beta_mulpix;   //!
		TBranch         *b_beta_pixelx;   //!
		TBranch         *b_beta_pixely;   //!
		TBranch         *b_beta_pixelz;   //!
		TBranch         *b_beta_pixelEX;   //!
		TBranch         *b_beta_pixelEY;   //!


		// Eurica branches
		TBranch		*b_gc_n_gamma;
		TBranch 	*b_gc_ts; //  [GeCluster_]
		TBranch 	*b_gc_hit;
		TBranch     *b_gc_ch;   //[GeCluster_]
		TBranch     *b_gc_DGFt;   //[GeCluster_]
		TBranch		*b_gc_Ecal;

		TBranch     *b_ab_hit;
		TBranch     *b_ab_ch;   //[GeAddback_]
		TBranch     *b_ab_DGFt;   //[GeAddback_]
		TBranch     *b_ab_E;   //[GeAddback_]

		BeamGe_AIDA();
		~BeamGe_AIDA();
		virtual void	 GetTree(char* filename, TTree *tree=0);
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
//		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);
//		virtual void	 GetTsEntry(std::map <Long64_t, Long64_t> &mts);

		//tree structure
		virtual void     TreeBranch(TTree *tree);
		virtual void     TreeBranchIon(TTree *tree);
		//virtual void     TreeBranchBeta(TTree *tree);
		virtual void  	 TreeBranchBetaGam(TTree *tree);

		//Sync functions
//		virtual void     SyncTSbeta(BigRIPS_reco &bigrips);
//		virtual void     SyncBigRIPS(BigRIPS_reco &bigrips);
		virtual void	 ResetGe();
		virtual void	 ResetBeam();
		virtual void     ResetBeta();
		//Calib functions
		virtual void	 LoadCalib(Calib &calib);
		virtual void	 CalGamma(eurica &gamma, Int_t gc_n);
		virtual void     SyncGammaBeta(eurica &gamma, Int_t gc_n);
		virtual void     SyncIonBeta(BigRIPS_reco &bigrips, AIDA_ionevt &ion, Int_t ion_n);
		virtual void     SyncBeamIon(BigRIPS_reco &bigrips, AIDA_ionevt &ion);
		virtual void 	 SyncMulti(size_t gc_n, size_t ion_n);
		virtual void 	 SyncBeta(AIDA_betaevt &beta);
		virtual int 	 CorrPosition(AIDA_betaevt &beta, AIDA_ionevt& ion);

};

#endif


BeamGe_AIDA::BeamGe_AIDA(){
}

BeamGe_AIDA::~BeamGe_AIDA()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

void BeamGe_AIDA::GetTree(char *filename, TTree *tree)
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
		if (!f || !f->IsOpen()) {
			f = new TFile(filename);
		}
		f->GetObject("tree",tree);

	}
	Init(tree);
}

Int_t BeamGe_AIDA::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}

Long64_t BeamGe_AIDA::LoadTree(Long64_t entry)
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

void BeamGe_AIDA::Init(TTree *tree)
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
	/******* ion brips and aida ******/	
	fChain->SetBranchAddress("n_ion",&n_ion,&b_n_ion);
	fChain->SetBranchAddress("ion_ts",ion_ts,&b_ion_ts);
	fChain->SetBranchAddress("ion_bigrips_ts",ion_bigrips_ts,&b_ion_bigrips_ts);
	fChain->SetBranchAddress("bigrips_BeamZ",bigrips_BeamZ,&b_bigrips_BeamZ);
	fChain->SetBranchAddress("bigrips_BeamAOQ",bigrips_BeamAOQ,&b_bigrips_BeamAOQ);
	fChain->SetBranchAddress("ion_x",ion_x,&b_ion_x);
	fChain->SetBranchAddress("ion_y",ion_y,&b_ion_y);
	fChain->SetBranchAddress("ion_z",ion_z,&b_ion_z);
	fChain->SetBranchAddress("ion_EX",ion_EX,&b_ion_EX);
	fChain->SetBranchAddress("ion_EY",ion_EY,&b_ion_EY);
	fChain->SetBranchAddress("ion_goodness",ion_goodness,&b_ion_goodness);
        fChain->SetBranchAddress("bigrips_F11PPAC1X1",bigrips_F11PPAC1X1,&b_bigrips_F11PPAC1X1);
        fChain->SetBranchAddress("bigrips_F11PPAC1X2",bigrips_F11PPAC1X2,&b_bigrips_F11PPAC1X2);
        fChain->SetBranchAddress("bigrips_F11PPAC2X1",bigrips_F11PPAC2X1,&b_bigrips_F11PPAC2X1);
        fChain->SetBranchAddress("bigrips_F11PPAC2X2",bigrips_F11PPAC2X2,&b_bigrips_F11PPAC2X2);
        fChain->SetBranchAddress("bigrips_F5PPAC1X1",bigrips_F5PPAC1X1,&b_bigrips_F5PPAC1X1);
        fChain->SetBranchAddress("bigrips_F5PPAC1X2",bigrips_F5PPAC1X2,&b_bigrips_F5PPAC1X2);
        fChain->SetBranchAddress("bigrips_F5PPAC2X1",bigrips_F5PPAC2X1,&b_bigrips_F5PPAC2X1);
        fChain->SetBranchAddress("bigrips_F5PPAC2X2",bigrips_F5PPAC2X2,&b_bigrips_F5PPAC2X2);
        fChain->SetBranchAddress("bigrips_F3PPAC1X1",bigrips_F3PPAC1X1,&b_bigrips_F3PPAC1X1);
        fChain->SetBranchAddress("bigrips_F3PPAC1X2",bigrips_F3PPAC1X2,&b_bigrips_F3PPAC1X2);
        fChain->SetBranchAddress("bigrips_F3PPAC2X1",bigrips_F3PPAC2X1,&b_bigrips_F3PPAC2X1);
        fChain->SetBranchAddress("bigrips_F3PPAC2X2",bigrips_F3PPAC2X2,&b_bigrips_F3PPAC2X2);


	/******* beta aida ******/	
	fChain->SetBranchAddress("beta_ts",&beta_ts,&b_beta_ts);
	fChain->SetBranchAddress("beta_x",&beta_x,&b_beta_x);
	fChain->SetBranchAddress("beta_y",&beta_y,&b_beta_y);
	fChain->SetBranchAddress("beta_z",&beta_z,&b_beta_z);
	fChain->SetBranchAddress("beta_EX",&beta_EX,&b_beta_EX);
	fChain->SetBranchAddress("beta_EY",&beta_EY,&b_beta_EY);
	fChain->SetBranchAddress("beta_goodness",&beta_goodness,&b_beta_goodness);
	fChain->SetBranchAddress("beta_mulpix",&beta_mulpix,&b_beta_mulpix);
	fChain->SetBranchAddress("beta_pixelx",beta_pixelx,&b_beta_pixelx);   //[mulpix]
	fChain->SetBranchAddress("beta_pixely",beta_pixely,&b_beta_pixely);   //[mulpix]
	fChain->SetBranchAddress("beta_pixelz",beta_pixelz,&b_beta_pixelz);   //[mulpix]
	fChain->SetBranchAddress("beta_pixelEX",beta_pixelEX,&b_beta_pixelEX);   //[mulpix]
	fChain->SetBranchAddress("beta_pixelEY",beta_pixelEY,&b_beta_pixelEY);  //[mulpix]

	/****** without addback ******/
	fChain->SetBranchAddress("n_gamma",&n_gamma,&b_gc_n_gamma);
	fChain->SetBranchAddress("gc_ts",gc_ts,&b_gc_ts);
	fChain->SetBranchAddress("gc_hit",gc_hit,&b_gc_hit);
	fChain->SetBranchAddress("gc_ch",gc_ch,&b_gc_ch);
	fChain->SetBranchAddress("gc_DGFt",gc_DGFt,&b_gc_DGFt);
	fChain->SetBranchAddress("gc_Ecal",gc_Ecal,&b_gc_Ecal);

	/****** with addback ******/
	fChain->SetBranchAddress("ab_hit",&ab_hit,&b_ab_hit);
	fChain->SetBranchAddress("ab_DGFt",ab_DGFt,&b_ab_DGFt);
	fChain->SetBranchAddress("ab_E",ab_E,&b_ab_E);

	Notify();
}

Bool_t BeamGe_AIDA::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void BeamGe_AIDA::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t BeamGe_AIDA::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

void BeamGe_AIDA::TreeBranch(TTree *tree)
{
	tree->Branch("n_ion",&n_ion,"n_ion/I");
	tree->Branch("ion_ts",ion_ts,"ion_ts[n_ion]/L");
	tree->Branch("ion_bigrips_ts",ion_bigrips_ts,"ion_bigrips_ts[n_ion]/L");
	tree->Branch("bigrips_BeamZ",bigrips_BeamZ,"bigrips_BeamZ[n_ion]/D");
	tree->Branch("bigrips_BeamAOQ",bigrips_BeamAOQ,"bigrips_BeamAOQ[n_ion]/D");
        tree->Branch("bigrips_F11PPAC1X1",bigrips_F11PPAC1X1,"bigrips_F11PPAC1X1[n_ion]/I");
        tree->Branch("bigrips_F11PPAC1X2",bigrips_F11PPAC1X2,"bigrips_F11PPAC1X2[n_ion]/I");
        tree->Branch("bigrips_F11PPAC2X1",bigrips_F11PPAC2X1,"bigrips_F11PPAC2X1[n_ion]/I");
        tree->Branch("bigrips_F11PPAC2X2",bigrips_F11PPAC2X2,"bigrips_F11PPAC2X2[n_ion]/I");
        tree->Branch("bigrips_F5PPAC1X1",bigrips_F5PPAC1X1,"bigrips_F5PPAC1X1[n_ion]/I");
        tree->Branch("bigrips_F5PPAC1X2",bigrips_F5PPAC1X2,"bigrips_F5PPAC1X2[n_ion]/I");
        tree->Branch("bigrips_F5PPAC2X1",bigrips_F5PPAC2X1,"bigrips_F5PPAC2X1[n_ion]/I");
        tree->Branch("bigrips_F5PPAC2X2",bigrips_F5PPAC2X2,"bigrips_F5PPAC2X2[n_ion]/I");
        tree->Branch("bigrips_F3PPAC1X1",bigrips_F3PPAC1X1,"bigrips_F3PPAC1X1[n_ion]/I");
        tree->Branch("bigrips_F3PPAC1X2",bigrips_F3PPAC1X2,"bigrips_F3PPAC1X2[n_ion]/I");
        tree->Branch("bigrips_F3PPAC2X1",bigrips_F3PPAC2X1,"bigrips_F3PPAC2X1[n_ion]/I");
        tree->Branch("bigrips_F3PPAC2X2",bigrips_F3PPAC2X2,"bigrips_F3PPAC2X2[n_ion]/I");


	tree->Branch("ion_x",ion_x,"ion_x[n_ion]/D");
	tree->Branch("ion_y",ion_y,"ion_y[n_ion]/D");
	tree->Branch("ion_z",ion_z,"ion_z[n_ion]/D");
	tree->Branch("ion_EX",ion_EX,"ion_EX[n_ion]/D");
	tree->Branch("ion_EY",ion_EY,"ion_EY[n_ion]/D");
	tree->Branch("ion_goodness",ion_goodness,"ion_goodness[n_ion]/I");


	/****** beta aida ******/
	tree->Branch("beta_ts",&beta_ts,"beta_ts/L");
	tree->Branch("beta_x",&beta_x,"beta_x/D");
	tree->Branch("beta_y",&beta_y,"beta_y/D");
	tree->Branch("beta_z",&beta_z,"beta_z/D");
	tree->Branch("beta_EX",&beta_EX,"beta_EX/D");
	tree->Branch("beta_EY",&beta_EY,"beta_EY/D");
	tree->Branch("beta_goodness",&beta_goodness,"beta_goodness/I");
	tree->Branch("beta_mulpix",&beta_mulpix,"beta_mulpix/I");
	tree->Branch("beta_pixelx",beta_pixelx,"beta_pixelx[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixely",beta_pixely,"beta_pixely[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixelz",beta_pixelz,"beta_pixelz[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixelEX",beta_pixelEX,"beta_pixelEX[beta_mulpix]/D");   //[mulpix]
	tree->Branch("beta_pixelEY",beta_pixelEY,"beta_pixelEY[beta_mulpix]/D");  //[mulpix]

	/****** without addback ******/
	tree->Branch("n_gamma",&n_gamma,"n_gamma/I");
	tree->Branch("gc_ts",gc_ts,"gc_ts[n_gamma]/L");
	tree->Branch("gc_hit",gc_hit,"gc_hit[n_gamma]/I");
	tree->Branch("gc_ch",gc_ch,"gc_ch[n_gamma][46]/I");
	tree->Branch("gc_DGFt",gc_DGFt,"gc_DGFt[n_gamma][46]/I");
	tree->Branch("gc_Ecal",gc_Ecal,"gc_Ecal[n_gamma][46]/D");
	/****** with addback ******/
	tree->Branch("ab_hit",ab_hit,"ab_hit[n_gamma]/I");
	tree->Branch("ab_DGFt",ab_DGFt,"ab_DGFt[n_gamma][27]/I");
	tree->Branch("ab_E",ab_E,"ab_E[n_gamma][27]/D");

}

void BeamGe_AIDA::TreeBranchBetaGam(TTree *tree)
{
	/****** beta aida ******/
	tree->Branch("beta_ts",&beta_ts,"beta_ts/L");
	tree->Branch("beta_x",&beta_x,"beta_x/D");
	tree->Branch("beta_y",&beta_y,"beta_y/D");
	tree->Branch("beta_z",&beta_z,"beta_z/D");
	tree->Branch("beta_EX",&beta_EX,"beta_EX/D");
	tree->Branch("beta_EY",&beta_EY,"beta_EY/D");
	tree->Branch("beta_goodness",&beta_goodness,"beta_goodness/I");
//	tree->Branch("beta_hitpixel",&beta_hitpixel);
	tree->Branch("beta_mulpix",&beta_mulpix,"beta_mulpix/I");
	tree->Branch("beta_pixelx",beta_pixelx,"beta_pixelx[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixely",beta_pixely,"beta_pixely[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixelz",beta_pixelz,"beta_pixelz[beta_mulpix]/I");   //[mulpix]
	tree->Branch("beta_pixelEX",beta_pixelEX,"beta_pixelEX[beta_mulpix]/D");   //[mulpix]
	tree->Branch("beta_pixelEY",beta_pixelEY,"beta_pixelEY[beta_mulpix]/D");  //[mulpix]

	/****** without addback ******/
	tree->Branch("n_gamma",&n_gamma,"n_gamma/I");
	tree->Branch("gc_ts",gc_ts,"gc_ts[n_gamma]/L");
	tree->Branch("gc_hit",gc_hit,"gc_hit[n_gamma]/I");
	tree->Branch("gc_ch",gc_ch,"gc_ch[n_gamma][46]/I");
	tree->Branch("gc_DGFt",gc_DGFt,"gc_DGFt[n_gamma][46]/I");
	tree->Branch("gc_Ecal",gc_Ecal,"gc_Ecal[n_gamma][46]/D");

	/****** with addback ******/
	tree->Branch("ab_hit",ab_hit,"ab_hit[n_gamma]/I");
	tree->Branch("ab_DGFt",ab_DGFt,"ab_DGFt[n_gamma][27]/I");
	tree->Branch("ab_E",ab_E,"ab_E[n_gamma][27]/D");

}



void BeamGe_AIDA::TreeBranchIon(TTree *tree)
{
	//tree->Branch("n_ion",&n_ion_bi,"n_ion/I");
	tree->Branch("ion_ts",&ion_ts_bi,"ion_ts/L");
	tree->Branch("ion_bigrips_ts",&ion_bigrips_ts_bi,"ion_bigrips_ts/L");
	tree->Branch("bigrips_BeamZ",&bigrips_BeamZ_bi,"bigrips_BeamZ/D");
	tree->Branch("bigrips_BeamAOQ",&bigrips_BeamAOQ_bi,"bigrips_BeamAOQ/D");
	tree->Branch("ion_x",&ion_x_bi,"ion_x/D");
	tree->Branch("ion_y",&ion_y_bi,"ion_y/D");
	tree->Branch("ion_z",&ion_z_bi,"ion_z/D");
	tree->Branch("ion_EX",&ion_EX_bi,"ion_EX/D");
	tree->Branch("ion_EY",&ion_EY_bi,"ion_EY/D");
	tree->Branch("ion_goodness",&ion_goodness_bi,"ion_goodness/I");
        tree->Branch("bigrips_F11PPAC1X1",&bigrips_F11PPAC1X1_bi,"bigrips_F11PPAC1X1/I");
        tree->Branch("bigrips_F11PPAC1X2",&bigrips_F11PPAC1X2_bi,"bigrips_F11PPAC1X2/I");
        tree->Branch("bigrips_F11PPAC2X1",&bigrips_F11PPAC2X1_bi,"bigrips_F11PPAC2X1/I");
        tree->Branch("bigrips_F11PPAC2X2",&bigrips_F11PPAC2X2_bi,"bigrips_F11PPAC2X2/I");
        tree->Branch("bigrips_F5PPAC1X1",&bigrips_F5PPAC1X1_bi,"bigrips_F5PPAC1X1/I");
        tree->Branch("bigrips_F5PPAC1X2",&bigrips_F5PPAC1X2_bi,"bigrips_F5PPAC1X2/I");
        tree->Branch("bigrips_F5PPAC2X1",&bigrips_F5PPAC2X1_bi,"bigrips_F5PPAC2X1/I");
        tree->Branch("bigrips_F5PPAC2X2",&bigrips_F5PPAC2X2_bi,"bigrips_F5PPAC2X2/I");
        tree->Branch("bigrips_F3PPAC1X1",&bigrips_F3PPAC1X1_bi,"bigrips_F3PPAC1X1/I");
        tree->Branch("bigrips_F3PPAC1X2",&bigrips_F3PPAC1X2_bi,"bigrips_F3PPAC1X2/I");
        tree->Branch("bigrips_F3PPAC2X1",&bigrips_F3PPAC2X1_bi,"bigrips_F3PPAC2X1/I");
        tree->Branch("bigrips_F3PPAC2X2",&bigrips_F3PPAC2X2_bi,"bigrips_F3PPAC2X2/I");



}
/*
void BeamGe_AIDA::SyncTS(BigRIPS_reco &bigrips){
		ts = bigrips.timestamp;
}
*/
void BeamGe_AIDA::SyncBeamIon(BigRIPS_reco &bigrips, AIDA_ionevt &ion){
	
		ion_ts_bi=ion.extTstop*4;
		ion_x_bi=ion.stop_x	;
		ion_y_bi=ion.stop_y	;
		ion_z_bi=ion.stop_z	;
		ion_goodness_bi=ion.goodness	;
		ion_EX_bi=ion.EX[(int)ion.stop_z]	;
		ion_EY_bi=ion.EY[(int)ion.stop_z]	;

    		ion_bigrips_ts_bi = bigrips.timestamp;	 
		bigrips_BeamZ_bi = bigrips.zet;
		bigrips_BeamAOQ_bi = bigrips.aoq[0];
                 bigrips_F3PPAC1X1_bi = 2.*bigrips.PPAC3X[0]-bigrips.PPAC3Xdiff[0];
                 bigrips_F3PPAC1X2_bi = 2.*bigrips.PPAC3X[0]+bigrips.PPAC3Xdiff[0];
                 bigrips_F3PPAC2X1_bi = 2.*bigrips.PPAC3X[1]-bigrips.PPAC3Xdiff[1];
                 bigrips_F3PPAC2X2_bi = 2.*bigrips.PPAC3X[1]+bigrips.PPAC3Xdiff[1];
                 bigrips_F5PPAC1X1_bi = 2.*bigrips.PPAC5X[0]-bigrips.PPAC5Xdiff[0];
                 bigrips_F5PPAC1X2_bi = 2.*bigrips.PPAC5X[0]+bigrips.PPAC5Xdiff[0];
                 bigrips_F5PPAC2X1_bi = 2.*bigrips.PPAC5X[1]-bigrips.PPAC5Xdiff[1];
                 bigrips_F5PPAC2X2_bi = 2.*bigrips.PPAC5X[1]+bigrips.PPAC5Xdiff[1];
                 bigrips_F11PPAC1X1_bi = 2.*bigrips.PPAC11X[0]-bigrips.PPAC11Xdiff[0];
                 bigrips_F11PPAC1X2_bi = 2.*bigrips.PPAC11X[0]+bigrips.PPAC11Xdiff[0];
                 bigrips_F11PPAC2X1_bi = 2.*bigrips.PPAC11X[1]-bigrips.PPAC11Xdiff[1];
                 bigrips_F11PPAC2X2_bi = 2.*bigrips.PPAC11X[1]+bigrips.PPAC11Xdiff[1];
	
}

void BeamGe_AIDA::LoadCalib(Calib &calib){

//	double sl, off;
	slope[0]=0;
	offset[0]=0;

	for(int i=0; i<84; i++){
		calib.GetEntry(i);
		slope[i+1]=calib.slope;
		offset[i+1]=calib.offset;
		}
}


void BeamGe_AIDA::SyncIonBeta(BigRIPS_reco &bigrips, AIDA_ionevt &ion, Int_t ion_n){
		
		ion_ts[ion_n]=ion.extTstop*4;
		ion_x[ion_n]=ion.stop_x	;
		ion_y[ion_n]=ion.stop_y	;
		ion_z[ion_n]=ion.stop_z ;
		ion_goodness[ion_n]=ion.goodness ;
		ion_EX[ion_n]=ion.EX[(int)ion.stop_z]	;
		ion_EY[ion_n]=ion.EY[(int)ion.stop_z]	;

	    	ion_bigrips_ts[ion_n] = bigrips.timestamp;	 
		bigrips_BeamZ[ion_n] = bigrips.zet;
		bigrips_BeamAOQ[ion_n] = bigrips.aoq[0];

                 bigrips_F3PPAC1X1[ion_n] = 2.*bigrips.PPAC3X[0]-bigrips.PPAC3Xdiff[0];
                 bigrips_F3PPAC1X2[ion_n] = 2.*bigrips.PPAC3X[0]+bigrips.PPAC3Xdiff[0];
                 bigrips_F3PPAC2X1[ion_n] = 2.*bigrips.PPAC3X[1]-bigrips.PPAC3Xdiff[1];
                 bigrips_F3PPAC2X2[ion_n] = 2.*bigrips.PPAC3X[1]+bigrips.PPAC3Xdiff[1];
                 bigrips_F5PPAC1X1[ion_n] = 2.*bigrips.PPAC5X[0]-bigrips.PPAC5Xdiff[0];
                 bigrips_F5PPAC1X2[ion_n] = 2.*bigrips.PPAC5X[0]+bigrips.PPAC5Xdiff[0];
                 bigrips_F5PPAC2X1[ion_n] = 2.*bigrips.PPAC5X[1]-bigrips.PPAC5Xdiff[1];
                 bigrips_F5PPAC2X2[ion_n] = 2.*bigrips.PPAC5X[1]+bigrips.PPAC5Xdiff[1];
                 bigrips_F11PPAC1X1[ion_n] = 2.*bigrips.PPAC11X[0]-bigrips.PPAC11Xdiff[0];
		 bigrips_F11PPAC1X2[ion_n] = 2.*bigrips.PPAC11X[0]+bigrips.PPAC11Xdiff[0];
                 bigrips_F11PPAC2X1[ion_n] = 2.*bigrips.PPAC11X[1]-bigrips.PPAC11Xdiff[1];
                 bigrips_F11PPAC2X2[ion_n] = 2.*bigrips.PPAC11X[1]+bigrips.PPAC11Xdiff[1];

		

}
void BeamGe_AIDA::SyncMulti(size_t gc_n, size_t ion_n){
	n_ion = ion_n;
	n_gamma = gc_n;
}


void BeamGe_AIDA::CalGamma(eurica &gamma, Int_t gc_n){

   gc_hit[gc_n] = gamma.GeCluster_;
  // memcpy(gc_ch[gc_n]   ,gamma.GeCluster_channel,sizeof(gamma.GeCluster_channel));
  // memcpy(gc_Ecal[gc_n] ,gamma.GeCluster_fEnergy, sizeof(gamma.GeCluster_fEnergy));
   memcpy(gc_DGFe, gamma.GeCluster_fADCe, sizeof(gamma.GeCluster_fADCe));

   for(int i=0; i<gamma.GeCluster_; i++){
   	gc_ch[gc_n][i]=gamma.GeCluster_channel[i];
   	gc_Ecal[gc_n][i]= gc_DGFe[i]*slope[gc_ch[gc_n][i]]+offset[gc_ch[gc_n][i]];
   }
}

void BeamGe_AIDA::SyncGammaBeta(eurica &gamma, Int_t gc_n){

   gc_ts[gc_n]  = gamma.EventInfo_timestamp[0];
   //memcpy(gc_DGFt[gc_n], gamma.GeCluster_fADCt, sizeof(gamma.GeCluster_fADCt));
    for(int i=0; i<gc_hit[gc_n]; i++){
   	gc_DGFt[gc_n][i]=gamma.GeCluster_fADCt[i];
   }

   ab_hit[gc_n] = gamma.GeAddback_;
   for(int i=0; i<ab_hit[gc_n]; i++){
   ab_ch[gc_n][i]= gamma.GeAddback_channel[i];
   ab_DGFt[gc_n][i]= gamma.GeAddback_fADCt[i];
   ab_E[gc_n][i]= gamma.GeAddback_fEnergy[i];
   }

}

void BeamGe_AIDA::ResetBeta(){
        beta_ts=0;
        beta_EX=0;
        beta_EY=0;
        beta_x=0;
        beta_y=0;
        beta_z=0;
	beta_goodness=0;
//	beta_hitpixel->clear();
	beta_mulpix= 0;
	for(int i = 0; i<beta_mulpix;i++){
	beta_pixelx[i]=-1;   //[mulpix]
	beta_pixely[i]=-1;
	beta_pixelz[i]=-1;
	beta_pixelEX[i]=0;
	beta_pixelEY[i]=0;
	}
}

void BeamGe_AIDA::SyncBeta(AIDA_betaevt &beta){
	beta_ts=beta.extTstart*4;
	beta_EX=beta.EX[beta.start_z]	;
	beta_EY=beta.EY[beta.start_z]	;
	beta_x=beta.start_x ;
	beta_y=beta.start_y ;
	beta_z=beta.start_z ;
	beta_goodness=beta.goodness;
//	beta_hitpixel=beta.hitpixel;
	beta_mulpix= beta.mulpix;
	for(int i = 0; i<beta_mulpix;i++){
	beta_pixelx[i]=beta.pixelx[i];   //[mulpix]
	beta_pixely[i]=beta.pixely[i];   //[mulpix]
	beta_pixelz[i]=beta.pixelz[i];   //[mulpix]
	beta_pixelEX[i]=beta.pixelEX[i];   //[mulpix]
	beta_pixelEY[i]=beta.pixelEY[i];   //[mulpix]
	}
}
void BeamGe_AIDA::ResetGe(){
	for(int resetgei=0; resetgei<kMax_n_gamma;resetgei++){
	   	for(int j=0;j<kMaxGeCluster;j++){
	   		gc_ch[resetgei][j]=0;
	   		gc_Ecal[resetgei][j]= 0;
	   		gc_DGFt[resetgei][j]= 0;
		   	ab_ch[resetgei][j] = 0;
		   	ab_DGFt[resetgei][j]= 0;
   			ab_E[resetgei][j]= 0;
   		}
		gc_hit[resetgei] = 0;
		ab_hit[resetgei] = 0;
	}
	n_gamma=0;
}
void BeamGe_AIDA::ResetBeam(){
	n_ion=0;
	for(int resetbeami=0; resetbeami<kMax_n_ion;resetbeami++){
		ion_ts[resetbeami]=0;
		ion_x[resetbeami]=0;
		ion_y[resetbeami]=0	;
		ion_z[resetbeami]=0;
		ion_EX[resetbeami]=0;
		ion_EY[resetbeami]=0;
		ion_goodness[resetbeami]=0;

    	ion_bigrips_ts[resetbeami] = 0;	 
		bigrips_BeamZ[resetbeami] = 0;
		bigrips_BeamAOQ[resetbeami] = 0;	
	}
}

int BeamGe_AIDA::CorrPosition(AIDA_betaevt &beta, AIDA_ionevt &ion){

//	if(ion.goodEvent==1){
//	    if(beta.goodEvent==1){
		if( abs(beta.start_z-ion.stop_z) < -2 ){
//			if(beta.position[0][2]-ion.stop_xyz[2] <2){
				if((abs(beta.start_y-ion.stop_z)< 3 && abs(beta.start_x-ion.stop_z) <3)){
			return 1;
			}else{ return -1;}
		}else{return -1;}
//	}else {return -1;}

}


//#endif // #ifdef BigRIPS_cxx
