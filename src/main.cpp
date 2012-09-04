#include "reweight.h"
#include <fstream>
#include <iostream>
#include "photonpreselection.h"
#include <sstream>
#include <string>
#include "TEventList.h"
#include <RooVoigtian.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include "TFile.h"
#include "TH2D.h"
#include "toolbox.h"
#include "TROOT.h"
#include "TTree.h"
#include <vector>
#include "fits.h"
#include <event.h>

using namespace std;

//-----------------------------------------------
//const char * eleDY = "samplesDYToEE.txt";
//const char * muDY = "samplesDYToMuMu.txt";

const char * eleDY = "samplesDoubleElectron.txt";
const char * muDY = "samplesDoubleMu.txt";
//-----------------------------------------------

//-----------------------------------------------
const char * gElePresel = "samplesPhoton_ele.txt";
const char * gMuPresel = "samplesPhoton_mu.txt";

const char * gPtRewElePresel = "samplesPhotonPtRew_ele.txt";
const char * gPtRewMuPresel = "samplesPhotonPtRew_mu.txt";

const char * gNjetRewElePresel = "samplesPhotonJetRew_ele.txt";
const char * gNjetRewMuPresel = "samplesPhotonJetRew_mu.txt";

const char * gFinalElePresel = "samplesPhotonRew_ele.txt";
const char * gFinalMuPresel = "samplesPhotonRew_mu.txt";

const char * zptEleFileName = "zptEle_MC.ps";
const char * gptEleFileName = "gptEle_MC.ps";
const char * zptMuFileName = "zptMu_MC.ps";
const char * gptMuFileName = "gptMu_MC.ps";
//-----------------------------------------------

/*
//-----------------------------------------------
const char * gElePresel = "samplesPhotonDATA_ele.txt";
const char * gMuPresel = "samplesPhotonDATA_mu.txt";

const char * gPtRewElePresel = "samplesPhotonDATAPtRew_ele.txt";
const char * gPtRewMuPresel = "samplesPhotonDATAPtRew_mu.txt";

const char * gNjetRewElePresel = "samplesPhotonDATAJetRew_ele.txt";
const char * gNjetRewMuPresel = "samplesPhotonDATAJetRew_mu.txt";

const char * gFinalElePresel = "samplesPhotonDATARew_ele.txt";
const char * gFinalMuPresel = "samplesPhotonDATARew_mu.txt";

const char * zptEleFileName = "zptEle_DATA.ps";
const char * gptEleFileName = "gptEle_DATA.ps";
const char * zptMuFileName = "zptMu_DATA.ps";
const char * gptMuFileName = "gptMu_DATA.ps";
//-----------------------------------------------
*/

int main(int argc, const char * argv[]) {
	gROOT->ProcessLine("#include <vector>");

	fitMass("SamplesPreselected/MC7TeV_DYJetsToLL_muPresel.root", "muMassFit.root");
	readFitMass("muMassFit.root", "zMass_mu.ps");

	fitMass("SamplesPreselected/MC7TeV_DYJetsToLL_elePresel.root", "eleMassFit.root");
	readFitMass("eleMassFit.root", "zMass_ele.ps");

//	vector<Fit> zfitsEle;
//	fitPt(zfitsEle, eleDY, false, zptEleFileName);
//	vector<Fit> gfitsEle;
//	fitPt(gfitsEle, gElePresel, true, gptEleFileName);
//	reweightPt( gElePresel, gPtRewElePresel, zfitsEle, gfitsEle );
//
//	vector<Fit> zfitsMu;
//	fitPt(zfitsMu, muDY, false, zptMuFileName);
//	vector<Fit> gfitsMu;
//	fitPt(gfitsMu, gMuPresel, true, gptMuFileName);
//	reweightPt( gMuPresel, gPtRewMuPresel, zfitsMu, gfitsMu );
//
//	vector<string> zweights;
//	vector<string> gweights;
//	gweights.push_back("Weight");
//	gweights.push_back("ptWeight");
//
//	TH1D * zNJETele = new TH1D("zNJETele", "zNJETele", 15, 0, 15);
//	fillHistogram( eleDY, "NJET", zweights, zNJETele);
//	TH1D * gNJETele = new TH1D("gNJETele", "gNJETele", 15, 0, 15);
//	fillHistogram( gPtRewElePresel, "NJET", gweights, gNJETele);
//	reweight(gPtRewElePresel, gNjetRewElePresel, gNJETele, zNJETele, "NJET", "njetWeight");
//
//	TH1D * zNJETmu = new TH1D("zNJETmu", "zNJETmu", 15, 0, 15);
//	fillHistogram( muDY, "NJET", zweights, zNJETmu);
//	TH1D * gNJETmu = new TH1D("gNJETmu", "gNJETmu", 15, 0, 15);
//	fillHistogram( gPtRewMuPresel, "NJET", gweights, gNJETmu);
//	reweight(gPtRewMuPresel, gNjetRewMuPresel, gNJETmu, zNJETmu, "NJET", "njetWeight");
//
//	gweights.push_back("njetWeight");
//
//	TH1D * zNVTXele = new TH1D("zNVTXele", "zNVTXele", 50, 0, 50);
//	fillHistogram( eleDY, "NVTX", zweights, zNVTXele);
//	TH1D * gNVTXele = new TH1D("gNVTXele", "gNVTXele", 50, 0, 50);
//	fillHistogram( gNjetRewElePresel, "NVTX", gweights, gNVTXele);
//	reweight(gNjetRewElePresel, gFinalElePresel, gNVTXele, zNVTXele, "NVTX", "nvtxWeight");
//
//	TH1D * zNVTXmu = new TH1D("zNVTXmu", "zNVTXmu", 50, 0, 50);
//	fillHistogram( muDY, "NVTX", zweights, zNVTXmu);
//	TH1D * gNVTXmu = new TH1D("gNVTXmu", "gNVTXmu", 50, 0, 50);
//	fillHistogram( gNjetRewMuPresel, "NVTX", gweights, gNVTXmu);
//	reweight(gNjetRewMuPresel, gFinalMuPresel, gNVTXmu, zNVTXmu, "NVTX", "nvtxWeight");
}
