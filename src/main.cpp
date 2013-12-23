#include <event.h>
#include "fits.h"
#include <fstream>
#include <iostream>
#include <options.h>
#include <preselectionCMG.h>
#include "reweight.h"
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooVoigtian.h>
#include <RooWorkspace.h>
#include "samples.h"
#include <sstream>
#include <string>
#include <TCanvas.h>
#include "TEventList.h"
#include "TFile.h"
#include "TH2D.h"
#include "toolbox.h"
#include "TROOT.h"
#include "TTree.h"
#include <vector>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/reader.h>

using namespace std;

int main(int argc, const char * argv[]) {
	if (argc != 2) {
		cout << "USAGE: " << argv[0] << " CONFIG_FILE" << endl;
		return -1;
	}

	gROOT->ProcessLine("#include <vector>");
	try {
		Options::getInstance().readInOptions(argv[1]);

		bool isData = Options::getInstance().checkBoolOption("DATA");
		bool isMC = Options::getInstance().checkBoolOption("MC");
		if (isData == isMC)
			throw string("ERROR: Problem with isData and isMC flags!");

		string inputSamples = Options::getInstance().checkStringOption("INPUT_SAMPLES");

		ifstream input(inputSamples);
		if (!input.is_open())
			throw std::string("ERROR: Can't open JSON file!");

		Json::Value root;
		Json::Reader reader;
		bool parsingSuccessful = reader.parse( input, root );
		if ( !parsingSuccessful ) {
			cout  << "Failed to parse configuration : " << reader.getFormattedErrorMessages() << endl;
			return -1;
		}

		vector<PhotonSample> photonSamplesEle;
		vector<PhotonSample> photonSamplesMu;
		vector<PhotonSample> DE2012;
		vector<PhotonSample> DM2012;
		vector<PhotonSample> EWKEle;
		vector<PhotonSample> EWKMu;
		for (unsigned i = 0; i < root.size(); ++i) {
			string type = root[i].get("type", "").asString();
			if (type == "PhotonEle")
				photonSamplesEle.push_back( PhotonSample(root[i]) );
			else if (type == "PhotonMu")
				photonSamplesMu.push_back( PhotonSample(root[i]) );
			else if (type == "ZJetsEle")
				DE2012.push_back( PhotonSample(root[i]) );
			else if (type == "ZJetsMu")
				DM2012.push_back( PhotonSample(root[i]) );
			else if (type == "EWKEle")
				EWKEle.push_back( PhotonSample(root[i]) );
			else if (type == "EWKMu")
				EWKMu.push_back( PhotonSample(root[i]) );
			else {
				throw string("ERROR: Can't recognize sample type!");
			}
		}

		// Electron DY Mass Fit

		//string depath = DE2012[0].getInputDir() + "/" + DE2012[0].getSampleName() + ".root";
		//fitMass( depath, "eleMassFit2012Data.root" );
		RooWorkspace * eleMassFit = readMassFit("eleMassFit2012Data.root", "zMass_ele2012.ps");

		// Muon DY Mass Fit

		//string dmpath = DM2012[0].getInputDir() + "/" + DM2012[0].getSampleName() + ".root";
		//fitMass( dmpath, "muMassFit2012Data.root" );
		RooWorkspace * muMassFit = readMassFit("muMassFit2012Data.root", "zMass_mu2012.ps" );


		// Electron Preselection

		//for (unsigned i = 0; i < photonSamplesEle.size(); ++i) {
		//	if (isData)
		//		Options::getInstance().addBoolOption("DATA", true);
		//	Options::getInstance().addStringOption( "SAMPLE_NAME", photonSamplesEle[i].getSampleName() );
		//	Options::getInstance().addStringOption( "INPUT_DIR", photonSamplesEle[i].getInputDir() );
		//	Options::getInstance().addStringOption( "OUTPUT_DIR", photonSamplesEle[i].getOutputDir() );
		//	LeptonPreselectionCMG(PHOT, eleMassFit);
		//}
		for (unsigned i = 0; i < EWKEle.size(); ++i) {
			if (isData)
				Options::getInstance().addBoolOption("DATA", true);
			Options::getInstance().addStringOption( "SAMPLE_NAME", EWKEle[i].getSampleName() );
			Options::getInstance().addStringOption( "INPUT_DIR", EWKEle[i].getInputDir() );
			Options::getInstance().addStringOption( "OUTPUT_DIR", EWKEle[i].getOutputDir() );
			LeptonPreselectionCMG(PHOT, eleMassFit);
		}

		// Muon Preselection

		//for (unsigned i = 0; i < photonSamplesMu.size(); ++i) {
		//	if (isData)
		//		Options::getInstance().addBoolOption("DATA", true);
		//	Options::getInstance().addStringOption( "SAMPLE_NAME", photonSamplesMu[i].getSampleName() );
		//	Options::getInstance().addStringOption( "INPUT_DIR", photonSamplesMu[i].getInputDir() );
		//	Options::getInstance().addStringOption( "OUTPUT_DIR", photonSamplesMu[i].getOutputDir() );
		//	LeptonPreselectionCMG(PHOT, muMassFit);
		//}
		for (unsigned i = 0; i < EWKMu.size(); ++i) {
			if (isData)
				Options::getInstance().addBoolOption("DATA", true);
			Options::getInstance().addStringOption( "SAMPLE_NAME", EWKMu[i].getSampleName() );
			Options::getInstance().addStringOption( "INPUT_DIR", EWKMu[i].getInputDir() );
			Options::getInstance().addStringOption( "OUTPUT_DIR", EWKMu[i].getOutputDir() );
			LeptonPreselectionCMG(PHOT, muMassFit);
		}

		// PT Fitting

		//while(fitPt(DE2012, false, true, 1, "zEE2012PtFit"));
		//while(fitPt(DM2012, false, true, 1, "zMM2012PtFit"));
		//while (fitPt(photonSamplesEle, true, true, 1, "g2012PtFit"));

		//while(fitPt(DE2012, false, true, 2, "zEE2012PtFit"));
		//while(fitPt(DM2012, false, true, 2, "zMM2012PtFit"));
		//while (fitPt(photonSamplesEle, true, true, 2, "g2012PtFit"));

		//while(fitPt(DE2012, false, true, 3, "zEE2012PtFit"));
		//while(fitPt(DM2012, false, true, 3, "zMM2012PtFit"));
		//while (fitPt(photonSamplesEle, true, true, 3, "g2012PtFit"));

		// Pt Reweighting

		vector<string> gFits;
		gFits.push_back("FitFunction");
		gFits.push_back("FitFunctionDown");
		gFits.push_back("FitFunctionUp");
		vector<string> zFits;
		zFits.push_back("FitFunction");
		zFits.push_back("FitFunctionUp");
		zFits.push_back("FitFunctionDown");
		vector<string> systName;
		systName.push_back("");
		systName.push_back("_ZJETS_DOWN");
		systName.push_back("_ZJETS_UP");

		for (unsigned k = 0; k < gFits.size(); ++k) {
			vector<string> fitsNames;
			fitsNames.push_back(gFits[k]);
			fitsNames.push_back(zFits[k]);
			fitsNames.push_back(systName[k]);

//			for (unsigned i = 0; i < photonSamplesEle.size(); ++i) {
//				reweightPt( photonSamplesEle[i], "zEE2012PtFit", "g2012PtFit", fitsNames );
//			}
			for (unsigned i = 0; i < EWKEle.size(); ++i) {
				reweightPt( EWKEle[i], "zEE2012PtFit", "g2012PtFit", fitsNames );
			}
//			for (unsigned i = 0; i < photonSamplesMu.size(); ++i) {
//				reweightPt( photonSamplesMu[i], "zMM2012PtFit", "g2012PtFit", fitsNames );
//			}
			for (unsigned i = 0; i < EWKMu.size(); ++i) {
				reweightPt( EWKMu[i], "zMM2012PtFit", "g2012PtFit", fitsNames );
			}
//
//			vector<string> zweights;
//			vector<string> gweights;
//			gweights.push_back("Weight");
//			gweights.push_back("ptWeight");
//
//			TH2D * zNJETele = new TH2D( ("zNJETele" + fitsNames[2]).c_str(), ("zNJETele" + fitsNames[2]).c_str(), 15, 0, 15, 3, 1, 4);
//			fillHistogram<int>( DE2012, &PhotonSample::getPreselectedSampleName, "NHARDJET", zweights, zNJETele, true);
//			TH2D * gNJETele = new TH2D(("gNJETele" + fitsNames[2]).c_str(), ("gNJETele" + fitsNames[2]).c_str(), 15, 0, 15, 3, 1, 4);
//			fillHistogram<int>( photonSamplesEle, &PhotonSample::getPtRewSampleName, "NHARDJET", gweights, gNJETele, true);
//
//			TH2D * zNJETmu = new TH2D(("zNJETmu" + fitsNames[2]).c_str(), ("zNJETmu" + fitsNames[2]).c_str(), 15, 0, 15, 3, 1, 4);
//			fillHistogram<int>( DM2012, &PhotonSample::getPreselectedSampleName, "NHARDJET", zweights, zNJETmu, true);
//			TH2D * gNJETmu = new TH2D(("gNJETmu" + fitsNames[2]).c_str(), ("gNJETmu" + fitsNames[2]).c_str(), 15, 0, 15, 3, 1, 4);
//			fillHistogram<int>( photonSamplesMu, &PhotonSample::getPtRewSampleName, "NHARDJET", gweights, gNJETmu, true);
//
//			reweight(photonSamplesEle, &PhotonSample::getPtRewSampleName, &PhotonSample::getNjRewSampleName, gNJETele, zNJETele, "NHARDJET", "njetWeight", fitsNames[2]);
//			reweight(photonSamplesMu, &PhotonSample::getPtRewSampleName, &PhotonSample::getNjRewSampleName, gNJETmu, zNJETmu, "NHARDJET", "njetWeight", fitsNames[2]);
//
//			gweights.push_back("njetWeight");
//
//			TH2D * zNVTXele = new TH2D(("zNVTXele" + fitsNames[2]).c_str(), ("zNVTXele" + fitsNames[2]).c_str(), 50, 0, 50, 3, 1, 4);
//			fillHistogram<int>( DE2012, &PhotonSample::getPreselectedSampleName, "NVTX", zweights, zNVTXele, true);
//			TH2D * gNVTXele = new TH2D(("gNVTXele" + fitsNames[2]).c_str(), ("gNVTXele" + fitsNames[2]).c_str(), 50, 0, 50, 3, 1, 4);
//			fillHistogram<int>( photonSamplesEle, &PhotonSample::getNjRewSampleName, "NVTX", gweights, gNVTXele, true);
//
//			TH2D * zNVTXmu = new TH2D(("zNVTXmu" + fitsNames[2]).c_str(), ("zNVTXmu" + fitsNames[2]).c_str(), 50, 0, 50, 3, 1, 4);
//			fillHistogram<int>( DM2012, &PhotonSample::getPreselectedSampleName, "NVTX", zweights, zNVTXmu, true);
//			TH2D * gNVTXmu = new TH2D(("gNVTXmu" + fitsNames[2]).c_str(), ("gNVTXmu" + fitsNames[2]).c_str(), 50, 0, 50, 3, 1, 4);
//			fillHistogram<int>( photonSamplesMu, &PhotonSample::getNjRewSampleName, "NVTX", gweights, gNVTXmu, true);
//
//			reweight(photonSamplesEle, &PhotonSample::getNjRewSampleName, &PhotonSample::getNvRewSampleName, gNVTXele, zNVTXele, "NVTX", "nvtxWeight", fitsNames[2]);
//			reweight(photonSamplesMu, &PhotonSample::getNjRewSampleName, &PhotonSample::getNvRewSampleName, gNVTXmu, zNVTXmu, "NVTX", "nvtxWeight", fitsNames[2]);
//
//			gweights.push_back("nvtxWeight");
//
//			TH2D * zNormele = new TH2D(("zNormele" + fitsNames[2]).c_str(), ("zNormele" + fitsNames[2]).c_str(), 1, 0, 50, 3, 1, 4);
//			fillHistogram<double>( DE2012, &PhotonSample::getPreselectedSampleName, "PFMET", zweights, zNormele, false);
//			TH2D * gNormele = new TH2D(("gNormele" + fitsNames[2]).c_str(), ("gNormele" + fitsNames[2]).c_str(), 1, 0, 50, 3, 1, 4);
//			fillHistogram<double>( photonSamplesEle, &PhotonSample::getNvRewSampleName, "PFMET", gweights, gNormele, false);
//
//			TH2D * zNormmu = new TH2D(("zNormmu" + fitsNames[2]).c_str(), ("zNormmu" + fitsNames[2]).c_str(), 1, 0, 50, 3, 1, 4);
//			fillHistogram<double>( DM2012, &PhotonSample::getPreselectedSampleName, "PFMET", zweights, zNormmu, false);
//			TH2D * gNormmu = new TH2D(("gNormmu" + fitsNames[2]).c_str(), ("gNormmu" + fitsNames[2]).c_str(), 1, 0, 50, 3, 1, 4);
//			fillHistogram<double>( photonSamplesMu, &PhotonSample::getNvRewSampleName, "PFMET", gweights, gNormmu, false);
//
//			normalize(photonSamplesEle, &PhotonSample::getNvRewSampleName, &PhotonSample::getNormRewSampleName, gNormele, zNormele, "normWeight", fitsNames[2]);
//			normalize(photonSamplesMu, &PhotonSample::getNvRewSampleName, &PhotonSample::getNormRewSampleName, gNormmu, zNormmu, "normWeight", fitsNames[2]);
		}

	} catch (const string & e) {
		cout << e << endl;
	}
}
