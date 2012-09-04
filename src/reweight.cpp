#include "reweight.h"
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include "TCanvas.h"
#include "event.h"
#include "TH1D.h"
#include <fstream>
#include <cstdlib>

using std::vector;
using std::string;
using std::cout;
using std::endl;

double calculatePtWeight(const vector<Fit> & zfits, const vector<Fit> & gfits, double pt) {
	int zIdx = -1;
	int gIdx = -1;
	for (unsigned i = 0; i < zfits.size(); ++i) {
		if (pt >= zfits[i].low && pt < zfits[i].high) {
			zIdx = i;
			break;
		}
	}
	for (unsigned i = 0; i < gfits.size(); ++i) {
		if (pt >= gfits[i].low && pt < gfits[i].high) {
			gIdx = i;
			break;
		}
	}
	if (zIdx < 0 || gIdx < 0) {
		return 0;
	}
	return zfits[zIdx].ptFun->Eval(pt) / gfits[gIdx].ptFun->Eval(pt);
}

void reweightPt(const char * inFile, const char * outFile, const vector<Fit> & zfits, const vector<Fit> & gfits) {
	ifstream inputSamples( inFile );
	vector<string> inSamples;
	if( !inputSamples.is_open() ) {
		cout << "ERROR: Can't open file with samples (" << inFile << ")!" << endl;
		exit( EXIT_FAILURE );
	}
	while( !inputSamples.eof() ) {
		string temp;
		getline( inputSamples, temp );
		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
			continue;
		}
		inSamples.push_back( temp );
	}
	ifstream outputSamples( outFile );
	vector<string> outSamples;
	if( !outputSamples.is_open() ) {
		cout << "ERROR: Can't open file with samples (" << outFile << ")!" << endl;
		exit( EXIT_FAILURE );
	}
	while( !outputSamples.eof() ) {
		string temp;
		getline( outputSamples, temp );
		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
			continue;
		}
		outSamples.push_back( temp );
	}
	if (inSamples.size() != outSamples.size() ) {
		cout << "ERROR: Number of input samples different from number of output samples!" << endl;
		exit( EXIT_FAILURE );
	}

	for ( unsigned i = 0; i < inSamples.size(); ++i ) {
		cout << "PT reweighting: " << inSamples[i] << endl;
		TFile in( inSamples[i].c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		TFile out( outSamples[i].c_str(), "recreate" );
		TTree * smallTree = tree->CopyTree( "ZPT<650" );
		double ptWeight;
		TBranch * wb = smallTree->Branch("ptWeight", &ptWeight, "ptWeight/D");
		double pt;
		smallTree->SetBranchAddress("ZPT", &pt);

		unsigned long entries = smallTree->GetEntries();
		for ( unsigned long i = 0; i < entries; i++ ) {
			smallTree->GetEntry( i );
			ptWeight = calculatePtWeight(zfits, gfits, pt);
			wb->Fill();
		}
		smallTree->Write("", TObject::kOverwrite);
	}
}

void fillHistogram(const string & fileName, const string & varName, const vector<string> & weights, TH1D * histo ) {
	ifstream input( fileName.c_str() );
	vector<string> samples;
	if( !input.is_open() ) {
		cout << "ERROR: Can't open file with samples (" << fileName << ")!" << endl;
		exit( EXIT_FAILURE );
	}
	while( !input.eof() ) {
		string temp;
		getline( input, temp );
		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
			continue;
		}
		samples.push_back( temp );
	}

	histo->Sumw2();
	for ( unsigned i = 0; i < samples.size(); ++i ) {
		cout << samples[i] << endl;
		TFile in(  samples[i].c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		unsigned long entries = tree->GetEntries();
		for ( unsigned long j = 0; j < entries; j++ ) {
			tree->GetEntry( j );
			int tempVarValue = ev.getSVV<int>(varName);
			double tempWeight = 1;
			for ( unsigned k = 0; k < weights.size(); ++k )
				tempWeight *= ev.getSVV<double>( weights[k] );
			if (histo)
				histo->Fill( tempVarValue, tempWeight );
		}
	}
	histo->Scale( 1.0 / histo->Integral() );
}

void reweight(const char * inFile, const char * outFile, const TH1D * ghisto, const TH1D * zhisto, const char * varName, const char * branchName) {
	ifstream inputSamples( inFile );
	vector<string> inSamples;
	if( !inputSamples.is_open() ) {
		cout << "ERROR: Can't open file with samples (" << inFile << ")!" << endl;
		exit( EXIT_FAILURE );
	}
	while( !inputSamples.eof() ) {
		string temp;
		getline( inputSamples, temp );
		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
			continue;
		}
		inSamples.push_back( temp );
	}
	ifstream outputSamples( outFile );
	vector<string> outSamples;
	if( !outputSamples.is_open() ) {
		cout << "ERROR: Can't open file with samples (" << outFile << ")!" << endl;
		exit( EXIT_FAILURE );
	}
	while( !outputSamples.eof() ) {
		string temp;
		getline( outputSamples, temp );
		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
			continue;
		}
		outSamples.push_back( temp );
	}
	if (inSamples.size() != outSamples.size() ) {
		cout << "ERROR: Number of input samples different from number of output samples!" << endl;
		exit( EXIT_FAILURE );
	}

	TCanvas canv("canv", "canv", 800, 600);

	TH1D weights( *zhisto );
	weights.Divide( ghisto );
	weights.Draw();
	string outFileName(branchName);
	outFileName += ".ps";
	canv.SaveAs( outFileName.c_str() );

	for ( unsigned i = 0; i < inSamples.size(); ++i ) {
		cout << varName << " reweighting: " << inSamples[i] << endl;
		TFile in( inSamples[i].c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		TFile out( outSamples[i].c_str(), "recreate" );
		TTree * smallTree = tree->CopyTree( "" );
		double newWeight;
		string branchType(branchName);
		branchType += "/D";
		TBranch * wb = smallTree->Branch( branchName, &newWeight, branchType.c_str() );
		int varValue;
		smallTree->SetBranchAddress(varName, &varValue);

		unsigned long entries = smallTree->GetEntries();
		for ( unsigned long j = 0; j < entries; j++ ) {
			smallTree->GetEntry( j );
			if (varValue < weights.GetNbinsX() || varValue >= 0)
				newWeight = weights.GetBinContent(varValue + 1);
			else
				newWeight = 1;
			wb->Fill();
		}
		smallTree->Write("", TObject::kOverwrite);
	}
}
