#include "reweight.h"
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include "TCanvas.h"
#include "event.h"
#include "TH1D.h"
#include <fstream>
#include <cstdlib>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <sstream>
#include <RooAddPdf.h>
#include "samples.h"

using std::stringstream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

void reweightPt(const PhotonSample & sample, const string & zfits, const string & gfits, const vector<string> & fitNames) {
	cout << "PT reweighting: " << sample.getSampleName() << endl;
	string inputFile = sample.getOutputDir() + "/" + sample.getSampleName() + sample.getPreselectedSampleName() + ".root";
	cout << "\tInputSample: " << inputFile << endl;
	string outputFile = sample.getOutputDir() + "/" + sample.getSampleName() + sample.getPtRewSampleName() + fitNames[2] + ".root";
	cout << "\tOutputSample: " << outputFile << endl;
	TFile in( inputFile.c_str() );

	TH1D * nEvHisto = (TH1D *)  in.Get("nevt");
	TH1D * puHisto = (TH1D *) in.Get("pileup");
	TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );

	TFile out( outputFile.c_str(), "recreate" );

	TH1D outNEvHisto( *nEvHisto );
	outNEvHisto.Write("nevt");
	TH1D outPuHisto( *puHisto );
	outPuHisto.Write("pileup");

	Event ev(tree);
	double  * pt = ev.getSVA<double>("ZPT");
	TTree * smallTree = tree->CopyTree( "ZPT>55 && ZPT<755" );
	double ptWeight;
	TBranch * wb = smallTree->Branch("ptWeight", &ptWeight, "ptWeight/D");

	TGraph zfits1 = readPtFit(zfits, 1, 55, 755, fitNames[1]);
	TGraph zfits2 = readPtFit(zfits, 2, 55, 755, fitNames[1]);
	TGraph zfits3 = readPtFit(zfits, 3, 55, 755, fitNames[1]);
	TGraph gfits1 = readPtFit(gfits, 1, 55, 755, fitNames[0]);
	TGraph gfits2 = readPtFit(gfits, 2, 55, 755, fitNames[0]);
	TGraph gfits3 = readPtFit(gfits, 3, 55, 755, fitNames[0]);

	unsigned long entries = smallTree->GetEntries();
	int * nhardjet = ev.getSVA<int>("NHARDJET");
	int * nsoftjet = ev.getSVA<int>("NSOFTJET");
	double * delEta = ev.getSVA<double>("DETAJJ");
	double * mjj = ev.getSVA<double>("MJJ");
	for ( unsigned long i = 0; i < entries; i++ ) {
		smallTree->GetEntry( i );
		unsigned cat = evCategory(*nhardjet, *nsoftjet, *delEta, *mjj, false);
		switch (cat) {
			case 1:
				ptWeight = zfits1.Eval(*pt) / gfits1.Eval(*pt);
				break;
			case 2:
				ptWeight = zfits2.Eval(*pt) / gfits2.Eval(*pt);
				break;
			case 3:
				ptWeight = zfits3.Eval(*pt) / gfits3.Eval(*pt);
				break;
		}

		wb->Fill();
	}
	out.cd();
	smallTree->Write("", TObject::kOverwrite);
}

void reweight(const vector<PhotonSample> samples, NameMethod inName, NameMethod outName, const TH2D * ghisto, const TH2D * zhisto, const char * varName, const char * branchName, const std::string & syst) {
	string outFileName(branchName);
	TCanvas canv("canv", "canv", 800, 600);
	zhisto->Clone()->Draw("LEGO2Z 0");
	canv.SaveAs( (outFileName + "_zevents.ps").c_str() );
	ghisto->Clone()->Draw("LEGO2Z 0");
	canv.SaveAs( (outFileName + "_gevents.ps").c_str() );

	TH2D weights( *zhisto );
	weights.Divide( ghisto );
	weights.Draw("LEGO2Z 0");
	canv.SaveAs( (outFileName + "_weights.ps").c_str() );

	for ( unsigned i = 0; i < samples.size(); ++i ) {
		cout << varName << " reweighting:" << endl;
		std::string inputFile = samples[i].getOutputDir() + "/" + samples[i].getSampleName() + (samples[i].*inName)() + syst + ".root";
		cout << "Input sample: " << inputFile << endl;
		std::string outputFile = samples[i].getOutputDir() + "/" + samples[i].getSampleName() + (samples[i].*outName)() + syst + ".root";
		cout << "Output sample: " << outputFile << endl;
		TFile in( inputFile.c_str() );
		TH1D * nEvHisto = (TH1D *)  in.Get("nevt");
		TH1D * puHisto = (TH1D *) in.Get("pileup");
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		TFile out( outputFile.c_str(), "recreate" );
		TH1D outNEvHisto( *nEvHisto );
		outNEvHisto.Write("nevt");
		TH1D outPuHisto( *puHisto );
		outPuHisto.Write("pileup");
		TTree * smallTree = tree->CopyTree( "" );
		double newWeight;
		string branchType(branchName);
		branchType += "/D";
		TBranch * wb = smallTree->Branch( branchName, &newWeight, branchType.c_str() );
		int varValue;
		smallTree->SetBranchAddress(varName, &varValue);

		unsigned long entries = smallTree->GetEntries();
		int * nhardjet = ev.getSVA<int>("NHARDJET");
		int * nsoftjet = ev.getSVA<int>("NSOFTJET");
		double * delEta = ev.getSVA<double>("DETAJJ");
		double * mjj = ev.getSVA<double>("MJJ");
		for ( unsigned long j = 0; j < entries; j++ ) {
			tree->GetEntry( j );
			smallTree->GetEntry( j );
			unsigned cat = evCategory(*nhardjet, *nsoftjet, *delEta, *mjj, false);
			if (varValue < weights.GetNbinsX() || varValue >= 0) {
				newWeight = weights.GetBinContent(varValue + 1, cat);
			} else
				newWeight = 1;
			wb->Fill();
		}
		smallTree->Write("", TObject::kOverwrite);
	}
}

void normalize(
		const vector<PhotonSample> samples, NameMethod inName, NameMethod outName,
		const TH2D * ghisto, const TH2D * zhisto,
		const char * branchName, const std::string & syst
		) {

	TH2D weights( *zhisto );
	weights.Divide( ghisto );

	for ( unsigned i = 0; i < samples.size(); ++i ) {
		cout << "Normalization reweighting:" << endl;
		std::string inputFile = samples[i].getOutputDir() + "/" + samples[i].getSampleName() + (samples[i].*inName)() + syst + ".root";
		cout << "Input sample: " << inputFile << endl;
		std::string outputFile = samples[i].getOutputDir() + "/" + samples[i].getSampleName() + (samples[i].*outName)() + syst + ".root";
		cout << "Output sample: " << outputFile << endl;
		TFile in( inputFile.c_str() );
		TH1D * nEvHisto = (TH1D *)  in.Get("nevt");
		TH1D * puHisto = (TH1D *) in.Get("pileup");
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		TFile out( outputFile.c_str(), "recreate" );
		TH1D outNEvHisto( *nEvHisto );
		outNEvHisto.Write("nevt");
		TH1D outPuHisto( *puHisto );
		outPuHisto.Write("pileup");
		TTree * smallTree = tree->CopyTree( "" );
		double newWeight;
		string branchType(branchName);
		branchType += "/D";
		TBranch * wb = smallTree->Branch( branchName, &newWeight, branchType.c_str() );

		unsigned long entries = smallTree->GetEntries();
		int * nhardjet = ev.getSVA<int>("NHARDJET");
		int * nsoftjet = ev.getSVA<int>("NSOFTJET");
		double * delEta = ev.getSVA<double>("DETAJJ");
		double * mjj = ev.getSVA<double>("MJJ");
		for ( unsigned long j = 0; j < entries; j++ ) {
			tree->GetEntry( j );
			smallTree->GetEntry( j );
			unsigned cat = evCategory(*nhardjet, *nsoftjet, *delEta, *mjj, false);
			newWeight = weights.GetBinContent(1, cat);
			wb->Fill();
		}
		smallTree->Write("", TObject::kOverwrite);
	}
}
