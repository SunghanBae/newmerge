#ifndef newmerge_h
#define newmerge_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"
//#include "eurica.h"
//#include "BigRIPS_reco.h"
//#include "eurica.h"
//#include "AIDA_ionevt.h"
//#include "AIDA_betaevt.h"
#include "AIDAraw.h"
#include "GA.h"
#include "BA.h"
//#include "fixedpars.h"
#include "Calib.h"


const size_t		kMaxbeta = 100;
//const size_t		kMaxPix = 1000;
const int		kMaxGeCluster =46;
const int		kMaxGeAddback =27;

class newmerge {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain
	
		// Ion related parameters (BigRIPS + AIDA)
                size_t		n_beta;
		Long64_t	ion_ts;
		Long64_t	ion_bigrips_ts;
		Double_t	bigrips_BeamZ;
		Double_t	bigrips_BeamAOQ;
		Double_t	ion_x;
		Double_t	ion_y;
		Double_t	ion_z;
		Double_t	ion_EX;
		Double_t	ion_EY;
		Int_t		ion_good;
                Int_t 		bigrips_F3PPAC1X1;
                Int_t 		bigrips_F3PPAC1X2;
                Int_t 		bigrips_F3PPAC2X1;
                Int_t 		bigrips_F3PPAC2X2;
                Int_t 		bigrips_F5PPAC1X1;
                Int_t		bigrips_F5PPAC1X2;
                Int_t		bigrips_F5PPAC2X1;
                Int_t 		bigrips_F5PPAC2X2;
                Int_t 		bigrips_F11PPAC1X1;
                Int_t 		bigrips_F11PPAC1X2;
                Int_t 		bigrips_F11PPAC2X1;
                Int_t 		bigrips_F11PPAC2X2;


		// AIDA beta parameters
		Long64_t	beta_ts[kMaxbeta];
		Double_t	beta_EX[kMaxbeta];
		Double_t	beta_EY[kMaxbeta];
		Double_t	beta_x[kMaxbeta];
		Double_t	beta_y[kMaxbeta];
		Double_t	beta_z[kMaxbeta];
		Int_t		beta_good[kMaxbeta];

		Int_t		n_pixel;
		Int_t           pixelx[kMaxPix];   //[mulpix]
		Int_t           pixely[kMaxPix];   //[mulpix]
		Int_t           pixelz[kMaxPix];   //[mulpix]
		Double_t        pixelEX[kMaxPix];   //[mulpix]
		Double_t        pixelEY[kMaxPix];   //[mulpix]
		Long64_t	pixelts[kMaxPix];



//		std::vector<ROOT::Math::XYZVector> beta_hitpixel[kMaxbeta]; 

		// Eurica parameters
		// size_t 		n_gamma;
		Long64_t	gc_ts[kMaxbeta];
		Int_t 		gc_hit[kMaxbeta];
		Int_t           gc_ch[kMaxbeta][kMaxGeCluster];   //[GeCluster_]
		UInt_t          gc_DGFt[kMaxbeta][kMaxGeCluster];   //[GeCluster_]
		UInt_t          gc_DGFe[kMaxbeta][kMaxGeCluster];   //[GeCluster_]
		Double_t	gc_Ecal[kMaxbeta][kMaxGeCluster];
		Int_t           ab_hit[kMaxbeta];
		Int_t           ab_ch[kMaxbeta][kMaxGeAddback];   //[GeAddback_]
		UInt_t          ab_DGFt[kMaxbeta][kMaxGeAddback];   //[GeAddback_]
		Double_t        ab_E[kMaxbeta][kMaxGeAddback];   //[GeAddback_]

		Double_t	slope[85];
		Double_t	offset[85];



		// Ion branches
		TBranch		*b_n_beta;
		TBranch		*b_ion_ts;
		TBranch		*b_ion_bigrips_ts;
		TBranch		*b_bigrips_BeamZ;
		TBranch		*b_bigrips_BeamAOQ;
		TBranch		*b_ion_x;
		TBranch		*b_ion_y;
		TBranch		*b_ion_z;
		TBranch		*b_ion_EX;
		TBranch		*b_ion_EY;
		TBranch		*b_ion_good;
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
		TBranch		*b_beta_good;
		//TBranch		*b_beta_hitpixel;
		TBranch		*b_n_pixel;
		TBranch         *b_pixelx;   //[mulpix]
		TBranch         *b_pixely;   //[mulpix]
		TBranch         *b_pixelz;   //[mulpix]
		TBranch         *b_pixelEX;   //[mulpix]
		TBranch         *b_pixelEY;   //[mulpix]
		TBranch		*b_pixelts;

		// Eurica branches
		TBranch		*b_gc_n_gamma;
		TBranch 	*b_gc_ts; //  [GeCluster_]
		TBranch 	*b_gc_hit;
		TBranch		*b_gc_ch;   //[GeCluster_]
		TBranch         *b_gc_DGFt;   //[GeCluster_]
		TBranch		*b_gc_Ecal;

		TBranch         *b_ab_hit;
		TBranch         *b_ab_ch;   //[GeAddback_]
		TBranch         *b_ab_DGFt;   //[GeAddback_]
		TBranch         *b_ab_E;   //[GeAddback_]

		newmerge();
		~newmerge();
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
//		virtual void     TreeBranchIon(TTree *tree);
		//virtual void     TreeBranchBeta(TTree *tree);
//		virtual void  	 TreeBranchBetaGam(TTree *tree);

		//Sync functions
//		virtual void     SyncTSbeta(BigRIPS_reco &bigrips);
//		virtual void     SyncBigRIPS(BigRIPS_reco &bigrips);
		virtual void	 ResetGA();
		virtual void	 ResetBA();
//		virtual void     ResetBeta();
		//Calib functions
//		virtual void	 LoadCalib(Calib &calib);
//		virtual void	 CalGamma(eurica &gamma, Int_t gc_n);
		virtual void     SyncGA(GA &beta, Int_t beta_n);
		virtual void     SyncBA(BA &ion, Int_t beta_n);
//		virtual void     SyncBeamIon(BigRIPS_reco &bigrips, AIDA_ionevt &ion);
//		virtual void 	 SyncMulti(size_t gc_n, size_t ion_n);
//		virtual void 	 SyncBeta(AIDA_betaevt &beta);
		virtual int 	 CorrPosition(GA &beta, BA &ion);
		virtual int	 CorrPixel(GA &beta, BA &ion);

};

#endif


newmerge::newmerge(){
}

newmerge::~newmerge()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

void newmerge::GetTree(char *filename, TTree *tree)
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

Int_t newmerge::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}

Long64_t newmerge::LoadTree(Long64_t entry)
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

void newmerge::Init(TTree *tree)
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
	fChain->SetBranchAddress("ion_ts",&ion_ts,&b_ion_ts);
	fChain->SetBranchAddress("ion_bigrips_ts",&ion_bigrips_ts,&b_ion_bigrips_ts);
	fChain->SetBranchAddress("bigrips_BeamZ",&bigrips_BeamZ,&b_bigrips_BeamZ);
	fChain->SetBranchAddress("bigrips_BeamAOQ",&bigrips_BeamAOQ,&b_bigrips_BeamAOQ);
	fChain->SetBranchAddress("ion_x",&ion_x,&b_ion_x);
	fChain->SetBranchAddress("ion_y",&ion_y,&b_ion_y);
	fChain->SetBranchAddress("ion_z",&ion_z,&b_ion_z);
	fChain->SetBranchAddress("ion_EX",&ion_EX,&b_ion_EX);
	fChain->SetBranchAddress("ion_EY",&ion_EY,&b_ion_EY);
	fChain->SetBranchAddress("ion_good",&ion_good,&b_ion_good);
        fChain->SetBranchAddress("bigrips_F11PPAC1X1",&bigrips_F11PPAC1X1,&b_bigrips_F11PPAC1X1);
        fChain->SetBranchAddress("bigrips_F11PPAC1X2",&bigrips_F11PPAC1X2,&b_bigrips_F11PPAC1X2);
        fChain->SetBranchAddress("bigrips_F11PPAC2X1",&bigrips_F11PPAC2X1,&b_bigrips_F11PPAC2X1);
        fChain->SetBranchAddress("bigrips_F11PPAC2X2",&bigrips_F11PPAC2X2,&b_bigrips_F11PPAC2X2);
        fChain->SetBranchAddress("bigrips_F5PPAC1X1",&bigrips_F5PPAC1X1,&b_bigrips_F5PPAC1X1);
        fChain->SetBranchAddress("bigrips_F5PPAC1X2",&bigrips_F5PPAC1X2,&b_bigrips_F5PPAC1X2);
        fChain->SetBranchAddress("bigrips_F5PPAC2X1",&bigrips_F5PPAC2X1,&b_bigrips_F5PPAC2X1);
        fChain->SetBranchAddress("bigrips_F5PPAC2X2",&bigrips_F5PPAC2X2,&b_bigrips_F5PPAC2X2);
        fChain->SetBranchAddress("bigrips_F3PPAC1X1",&bigrips_F3PPAC1X1,&b_bigrips_F3PPAC1X1);
        fChain->SetBranchAddress("bigrips_F3PPAC1X2",&bigrips_F3PPAC1X2,&b_bigrips_F3PPAC1X2);
        fChain->SetBranchAddress("bigrips_F3PPAC2X1",&bigrips_F3PPAC2X1,&b_bigrips_F3PPAC2X1);
        fChain->SetBranchAddress("bigrips_F3PPAC2X2",&bigrips_F3PPAC2X2,&b_bigrips_F3PPAC2X2);


	/******* beta aida ******/	
	fChain->SetBranchAddress("n_beta",&n_beta,&b_n_beta);
	fChain->SetBranchAddress("beta_ts",beta_ts,&b_beta_ts);
	fChain->SetBranchAddress("beta_x",beta_x,&b_beta_x);
	fChain->SetBranchAddress("beta_y",beta_y,&b_beta_y);
	fChain->SetBranchAddress("beta_z",beta_z,&b_beta_z);
	fChain->SetBranchAddress("beta_EX",beta_EX,&b_beta_EX);
	fChain->SetBranchAddress("beta_EY",beta_EY,&b_beta_EY);
	fChain->SetBranchAddress("beta_good",beta_good,&b_beta_good);
//	fChain->SetBranchAddress("beta_hitpixel",beta_hitpixel,&b_beta_hitpixel);
	fChain->SetBranchAddress("n_pixel",&n_pixel,&b_n_pixel);
	fChain->SetBranchAddress("pixelx",pixelx,&b_pixelx);   //[mulpix]
	fChain->SetBranchAddress("pixely",pixely,&b_pixely);   //[mulpix]
	fChain->SetBranchAddress("pixelz",pixelz,&b_pixelz);   //[mulpix]
	fChain->SetBranchAddress("pixelEX",pixelEX,&b_pixelEX);   //[mulpix]
	fChain->SetBranchAddress("pixelEY",pixelEY,&b_pixelEY);  //[mulpix]
	fChain->SetBranchAddress("pixelts",pixelts,&b_pixelts);  //[mulpix]

	/****** without addback ******/
	//fChain->SetBranchAddress("n_gamma",&n_gamma,&b_gc_n_gamma);
	fChain->SetBranchAddress("gc_ts",gc_ts,&b_gc_ts);
	fChain->SetBranchAddress("gc_hit",gc_hit,&b_gc_hit);
	fChain->SetBranchAddress("gc_ch",gc_ch,&b_gc_ch);
	fChain->SetBranchAddress("gc_DGFt",gc_DGFt,&b_gc_DGFt);
	fChain->SetBranchAddress("gc_Ecal",gc_Ecal,&b_gc_Ecal);

	/****** with addback ******/
	fChain->SetBranchAddress("ab_hit",ab_hit,&b_ab_hit);
	fChain->SetBranchAddress("ab_DGFt",ab_DGFt,&b_ab_DGFt);
	fChain->SetBranchAddress("ab_E",ab_E,&b_ab_E);

	Notify();
}

Bool_t newmerge::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void newmerge::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t newmerge::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

void newmerge::TreeBranch(TTree *tree)
{
	tree->Branch("ion_ts",&ion_ts,"ion_ts/L");
	tree->Branch("ion_bigrips_ts",&ion_bigrips_ts,"ion_bigrips_ts/L");
	tree->Branch("bigrips_BeamZ",&bigrips_BeamZ,"bigrips_BeamZ/D");
	tree->Branch("bigrips_BeamAOQ",&bigrips_BeamAOQ,"bigrips_BeamAOQ/D");
        tree->Branch("bigrips_F11PPAC1X1",&bigrips_F11PPAC1X1,"bigrips_F11PPAC1X1/I");
        tree->Branch("bigrips_F11PPAC1X2",&bigrips_F11PPAC1X2,"bigrips_F11PPAC1X2/I");
        tree->Branch("bigrips_F11PPAC2X1",&bigrips_F11PPAC2X1,"bigrips_F11PPAC2X1/I");
        tree->Branch("bigrips_F11PPAC2X2",&bigrips_F11PPAC2X2,"bigrips_F11PPAC2X2/I");
        tree->Branch("bigrips_F5PPAC1X1",&bigrips_F5PPAC1X1,"bigrips_F5PPAC1X1/I");
        tree->Branch("bigrips_F5PPAC1X2",&bigrips_F5PPAC1X2,"bigrips_F5PPAC1X2/I");
        tree->Branch("bigrips_F5PPAC2X1",&bigrips_F5PPAC2X1,"bigrips_F5PPAC2X1/I");
        tree->Branch("bigrips_F5PPAC2X2",&bigrips_F5PPAC2X2,"bigrips_F5PPAC2X2/I");
        tree->Branch("bigrips_F3PPAC1X1",&bigrips_F3PPAC1X1,"bigrips_F3PPAC1X1/I");
        tree->Branch("bigrips_F3PPAC1X2",&bigrips_F3PPAC1X2,"bigrips_F3PPAC1X2/I");
        tree->Branch("bigrips_F3PPAC2X1",&bigrips_F3PPAC2X1,"bigrips_F3PPAC2X1/I");
        tree->Branch("bigrips_F3PPAC2X2",&bigrips_F3PPAC2X2,"bigrips_F3PPAC2X2/I");


	tree->Branch("ion_x",&ion_x,"ion_x/D");
	tree->Branch("ion_y",&ion_y,"ion_y/D");
	tree->Branch("ion_z",&ion_z,"ion_z/D");
	tree->Branch("ion_EX",&ion_EX,"ion_EX/D");
	tree->Branch("ion_EY",&ion_EY,"ion_EY/D");
	tree->Branch("ion_good",&ion_good,"ion_good/I");




	/****** beta aida ******/
	tree->Branch("n_beta",&n_beta,"n_beta/I");
	tree->Branch("beta_ts",beta_ts,"beta_ts[n_beta]/L");
	tree->Branch("beta_x",beta_x,"beta_x[n_beta]/D");
	tree->Branch("beta_y",beta_y,"beta_y[n_beta]/D");
	tree->Branch("beta_z",beta_z,"beta_z[n_beta]/D");
	tree->Branch("beta_EX",beta_EX,"beta_EX[n_beta]/D");
	tree->Branch("beta_EY",beta_EY,"beta_EY[n_beta]/D");
	tree->Branch("beta_good",beta_good,"beta_good[n_beta]/I");
//	tree->Branch("beta_hitpixel",beta_hitpixel,"beta_hitpixel[n_beta]");

	tree->Branch("n_pixel",&n_pixel,"n_pixel/I");
	tree->Branch("pixelx",pixelx,"pixelx[n_pixel]/I");   //[mulpix]
	tree->Branch("pixely",pixely,"pixely[n_pixel]/I");   //[mulpix]
	tree->Branch("pixelz",pixelz,"pixelz[n_pixel]/I");   //[mulpix]
	tree->Branch("pixelEX",pixelEX,"pixelEX[n_pixel]/D");   //[mulpix]
	tree->Branch("pixelEY",pixelEY,"pixelEY[n_pixel]/D");  //[mulpix]
	tree->Branch("pixelts",pixelts,"pixelts[n_pixel]/L");  //[mulpix]


	/****** without addback ******/
//	tree->Branch("n_gamma",n_gamma,"n_gamma/I");
	tree->Branch("gc_ts",gc_ts,"gc_ts[n_beta]/L");
	tree->Branch("gc_hit",gc_hit,"gc_hit[n_beta]/I");
	tree->Branch("gc_ch",gc_ch,"gc_ch[n_beta][46]/I");
	tree->Branch("gc_DGFt",gc_DGFt,"gc_DGFt[n_beta][46]/I");
	tree->Branch("gc_Ecal",gc_Ecal,"gc_Ecal[n_beta][46]/D");
	/****** with addback ******/
	tree->Branch("ab_hit",ab_hit,"ab_hit[n_beta]/I");
	tree->Branch("ab_DGFt",ab_DGFt,"ab_DGFt[n_beta][27]/I");
	tree->Branch("ab_E",ab_E,"ab_E[n_beta][27]/D");

}

/*
void newmerge::SyncTS(BigRIPS_reco &bigrips){
		ts = bigrips.timestamp;
}
*/
/*
void newmerge::LoadCalib(Calib &calib){

//	double sl, off;
	slope[0]=0;
	offset[0]=0;

	for(int i=0; i<84; i++){
		calib.GetEntry(i);
		slope[i+1]=calib.slope;
		offset[i+1]=calib.offset;
		}
}
*/

void newmerge::SyncBA(BA &ion, Int_t beta_n){
		
		ion_ts=ion.ion_ts;
		ion_x=ion.ion_x	;
		ion_y=ion.ion_y	;
		ion_z=ion.ion_z ;
		ion_EX=ion.ion_EX	;
		ion_EY=ion.ion_EY	;
		ion_good=ion.ion_goodness;

	    	ion_bigrips_ts = ion.ion_bigrips_ts;	 
		bigrips_BeamZ = ion.bigrips_BeamZ;
		bigrips_BeamAOQ = ion.bigrips_BeamAOQ;

                bigrips_F3PPAC1X1 = ion.bigrips_F3PPAC1X1;
                bigrips_F3PPAC1X2 = ion.bigrips_F3PPAC1X2;
                bigrips_F3PPAC2X1 = ion.bigrips_F3PPAC2X1;
                bigrips_F3PPAC2X2 = ion.bigrips_F3PPAC2X2;
                bigrips_F5PPAC1X1 = ion.bigrips_F5PPAC1X1;
                bigrips_F5PPAC1X2 = ion.bigrips_F5PPAC1X2;
                bigrips_F5PPAC2X1 = ion.bigrips_F5PPAC2X1;
                bigrips_F5PPAC2X2 = ion.bigrips_F5PPAC2X2;
                bigrips_F11PPAC1X1 = ion.bigrips_F11PPAC1X1;
		bigrips_F11PPAC1X2 = ion.bigrips_F11PPAC1X2;
	        bigrips_F11PPAC2X1 = ion.bigrips_F11PPAC2X1;
                bigrips_F11PPAC2X2 = ion.bigrips_F11PPAC2X2;

		n_beta=beta_n;		

}

void newmerge::SyncGA(GA &beta, Int_t beta_n)
{
	beta_ts[beta_n] = beta.beta_ts;
	beta_x[beta_n] = beta.beta_x;
	beta_y[beta_n] = beta.beta_y;
	beta_z[beta_n] = beta.beta_z;
	beta_EX[beta_n] = beta.beta_EX;
	beta_EY[beta_n] = beta.beta_EY;
	beta_good[beta_n] = beta.beta_goodness;

	for(int i=n_pixel; i<n_pixel+beta.beta_mulpix;i++){
	pixelts[i]	= beta.beta_ts;
	pixelx[i]	= beta.beta_pixelx[i];
	pixely[i]	= beta.beta_pixely[i];
	pixelz[i]	= beta.beta_pixelz[i];
	pixelEX[i]	= beta.beta_pixelEX[i];
	pixelEY[i]	= beta.beta_pixelEY[i];
	}

	n_pixel		+=beta.beta_mulpix;

	gc_ts[beta_n] = beta.gc_ts[0] ;
	gc_hit[beta_n] = beta.gc_hit[0] ;
	for(int iga=0;iga<gc_hit[beta_n];iga++)
	{
	gc_ch[beta_n][iga] = beta.gc_ch[0][iga] ;
	gc_DGFt[beta_n][iga] = beta.gc_DGFt[0][iga] ;
	gc_Ecal[beta_n][iga] = beta.gc_Ecal[0][iga] ;
	}
	ab_hit[beta_n] = beta.ab_hit[0] ;
	for(int iab=0;iab<ab_hit[beta_n];iab++)
	{
	ab_DGFt[beta_n][iab] = beta.ab_DGFt[0][iab] ;
	ab_E[beta_n][iab] = beta.ab_E[0][iab] ;
	}

}

/*void newmerge::CalGamma(GA &gamma, Int_t beta_n){

   gc_hit[beta_n] = gamma.GeCluster_;
  // memcpy(gc_ch[gc_n]   ,gamma.GeCluster_channel,sizeof(gamma.GeCluster_channel));
  // memcpy(gc_Ecal[gc_n] ,gamma.GeCluster_fEnergy, sizeof(gamma.GeCluster_fEnergy));
   memcpy(gc_DGFe, gamma.GeCluster_fADCe, sizeof(gamma.GeCluster_fADCe));

   for(int i=0; i<gamma.GeCluster_; i++){
   	gc_ch[gc_n][i]=gamma.GeCluster_channel[i];
   	gc_Ecal[gc_n][i]= gc_DGFe[i]*slope[gc_ch[gc_n][i]]+offset[gc_ch[gc_n][i]];
   }
}*/

/*void newmerge::SyncGammaBeta(eurica &gamma, Int_t gc_n){

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

}*/

void newmerge::ResetGA(){
   for(int kk=0;kk<kMaxbeta;kk++){

	beta_ts[kk] = 0;
	beta_x[kk] =-999;
	beta_y[kk] =-999;
	beta_z[kk] =-999;
	beta_EX[kk]=-999;
	beta_EY[kk]= -999;
	beta_good[kk]= 0;
	gc_ts[kk] = 0;
	gc_hit[kk] = -999;
	for(int iga=-999;iga<kMaxGeCluster;iga++)
	{
	gc_ch[kk][iga] = -999;
	gc_DGFt[kk][iga] =-999 ;
	gc_Ecal[kk][iga] = -999;
	}
	ab_hit[kk] = -999;
	for(int iab=-999;iab<kMaxGeAddback;iab++)
	{
	ab_DGFt[kk][iab] = -999;
	ab_E[kk][iab] =-999;
	}
   }
	n_pixel=0;
}

/*void newmerge::SyncBeta(AIDA_betaevt &beta){
	beta_ts=beta.extTstart*4;
	beta_EX=beta.EX[beta.start_z]	;
	beta_EY=beta.EY[beta.start_z]	;
	beta_x=beta.start_x ;
	beta_y=beta.start_y ;
	beta_z=beta.start_z ;
}*/


void newmerge::ResetBA(){
		ion_ts=0;
		ion_x= -999;
		ion_y=-999;
		ion_z=-999;
		ion_EX=-999;
		ion_EY=-999;
		ion_good=0;

	    	ion_bigrips_ts =0;
		bigrips_BeamZ = -999;
		bigrips_BeamAOQ = -999;

                bigrips_F3PPAC1X1 = -999;
                bigrips_F3PPAC1X2 = -999;
                bigrips_F3PPAC2X1 = -999;
                bigrips_F3PPAC2X2 = -999;
                bigrips_F5PPAC1X1 = -999;
                bigrips_F5PPAC1X2 = -999;
                bigrips_F5PPAC2X1 = -999;
                bigrips_F5PPAC2X2 = -999;
                bigrips_F11PPAC1X1 =-999;
		bigrips_F11PPAC1X2 =-999;
	        bigrips_F11PPAC2X1 =-999;
                bigrips_F11PPAC2X2 =-999;

		n_beta=0;
}
/*
void newmerge::ResetBeam(){
	for(int resetbeami=0; resetbeami<kMax_n_ion;resetbeami++){
		ion_ts[resetbeami]=-999;
		ion_x[resetbeami]=-999;
		ion_y[resetbeami]=-999	;
		ion_z[resetbeami]=-999;
		ion_EX[resetbeami]=-999;
		ion_EY[resetbeami]=-999;

	    	ion_bigrips_ts[resetbeami] = -999;	 
		bigrips_BeamZ[resetbeami] =-999;
		bigrips_BeamAOQ[resetbeami] = -999;	
	}
}
*/
int newmerge::CorrPosition(GA &beta, BA &ion){

//	if(ion.goodEvent==1){
//	    if(beta.goodEvent==1){
	if(beta.beta_z>-1){
		if(abs(beta.beta_EX-beta.beta_EY)<600 && beta.beta_EX>0 && beta.beta_EY>0){
//		if(beta.beta_EX>0 && beta.beta_EY>0){
//		if(beta.beta_x>-1 && beta.beta_y>-1 && ion.ion_x>-1 && ion.ion_y>-1){
		if( abs(beta.beta_z-ion.ion_z) < 2 ){
//			if(beta.position[0][2]-ion.stop_xyz[2] <2){
				if((abs(beta.beta_y-ion.ion_y)< 4 && abs(beta.beta_x-ion.ion_x) <4)){
			return 1;
			}else{ return -1;}
		}else{return -1;}
//	}else {return -1;}
	}else{return -1;}
	}else{return -1;}
}

int newmerge::CorrPixel(GA &beta, BA &ion){
//	if(beta.beta_goodness==1 && ion.ion_goodness==1){
	if(ion.ion_goodness==1){
	for(int i=0; i<beta.beta_mulpix;i++){
		if( abs(beta.beta_pixelz[i]-ion.ion_z) < 2 ){
				if((abs(beta.beta_pixely[i]-ion.ion_y)< 4 && abs(beta.beta_pixelx[i]-ion.ion_x) <4)){
			return 1;
			}else{ return -1;}
		}else{return -1;}
	}
//	}else return -1;
}else return -1;
}

//#endif // #ifdef BigRIPS_cxx
