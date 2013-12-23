#ifndef REWEIGHT_H
#define REWEIGHT_H

#include <vector>
#include <string>
#include "fits.h"
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include "samples.h"
#include <event.h>
#include <toolbox.h>
#include <iostream>

typedef std::string (PhotonSample::*NameMethod)() const;
class TH1D;

void reweightPt(const PhotonSample & sample, const std::string & zfits, const std::string & gfits, const std::vector<std::string> & fitNames);
void reweight(const std::vector<PhotonSample> samples, NameMethod inName, NameMethod outName, const TH2D * ghisto, const TH2D * zhisto, const char * varName, const char * branchName, const std::string & syst);
void normalize(const std::vector<PhotonSample> samples, NameMethod inName, NameMethod outName, const TH2D * ghisto, const TH2D * zhisto, const char * branchName, const std::string & syst);
template <typename T>
void fillHistogram(
		const std::vector<PhotonSample> & samples, NameMethod name, const std::string & varName,
		const std::vector<std::string> & weights, TH2D * histo, bool normalize
	) {
	histo->Sumw2();
	for ( unsigned i = 0; i < samples.size(); ++i ) {
		std::string inputFile = samples[i].getOutputDir() + "/" + samples[i].getSampleName() + (samples[i].*name)() + ".root";
		TFile in(  inputFile.c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		double xsweight = samples[i].getCrossSection() / samples[i].getNumberOfEvents();
		unsigned long entries = tree->GetEntries();
		int * nhardjet = ev.getSVA<int>("NHARDJET");
		int * nsoftjet = ev.getSVA<int>("NSOFTJET");
		double * delEta = ev.getSVA<double>("DETAJJ");
		double * mjj = ev.getSVA<double>("MJJ");
		for ( unsigned long j = 0; j < entries; j++ ) {
			tree->GetEntry( j );
			unsigned cat = evCategory(*nhardjet, *nsoftjet, *delEta, *mjj, false);
			int tempVarValue = ev.getSVV<T>(varName);
			double tempWeight = xsweight;
			for ( unsigned k = 0; k < weights.size(); ++k )
				tempWeight *= ev.getSVV<double>( weights[k] );
			if (histo)
				histo->Fill( tempVarValue, cat, tempWeight );
		}
	}
	if (normalize)
		histo->Scale(1.0 / histo->Integral());
}

#endif
