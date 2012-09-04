#include "event.h"
#include "fits.h"
#include <fstream>
#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include "RooZPtPdf.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <vector>

using std::vector;
using std::string;
using std::cout;
using std::endl;

Fit::Fit( double l, double h ) : low(l), high(h), a(0), b(0), c(0), N(0) {
	ptFun = new TF1("ptFun", "[0]*TMath::Exp([1]*x*x+[2]*TMath::Sqrt(x*x*x)+[3]*x)", low, high);
	alow = -1;
	ahigh = 0.01;
	blow = -1;
	bhigh = 0.01;
	clow = -10;
	chigh = 0.01;
	histo = 0;
}

Fit::Fit( const Fit & fit ) {
	low = fit.low;
	high = fit.high;
	a = fit.a;
	b = fit.b;
	c = fit.c;
	N = fit.N;
	ptFun = new TF1("ptFun", "[0]*TMath::Exp([1]*x*x+[2]*TMath::Sqrt(x*x*x)+[3]*x)", low, high);
	ptFun->SetParameters(fit.N, fit.a, fit.b, fit.c);
	alow = fit.alow;
	blow = fit.blow;
	clow = fit.clow;
	ahigh = fit.ahigh;
	bhigh = fit.bhigh;
	chigh = fit.chigh;
	histo = 0;
}

Fit & Fit::operator=( const Fit & fit ) {
	low = fit.low;
	high = fit.high;
	a = fit.a;
	b = fit.b;
	c = fit.c;
	N = fit.N;
	if (ptFun)
		delete ptFun;
	ptFun = new TF1("ptFun", "[0]*TMath::Exp([1]*x*x+[2]*TMath::Sqrt(x*x*x)+[3]*x)", low, high);
	ptFun->SetParameters(fit.N, fit.a, fit.b, fit.c);
	alow = fit.alow;
	blow = fit.blow;
	clow = fit.clow;
	ahigh = fit.ahigh;
	bhigh = fit.bhigh;
	chigh = fit.chigh;
	if (histo)
		delete histo;
	histo = 0;
	return *this;
}

Fit::~Fit() {
	delete ptFun;
	if (histo)
		delete histo;
}

TH1D * Fit::buildHisto(int nBins, double lowE, double highE, Color_t col) {
	histo = new TH1D("fitHisto", "fitHisto", nBins, lowE, highE);
	histo->SetDirectory(0);
	for (int i = 1; i <= nBins; ++i) {
		double width = histo->GetBinWidth(i);
		double integ = ptFun->Integral(lowE + (i - 1) * width, lowE + i * width);
		if (lowE + i * width > low && lowE + (i - 1) * width < high)
			histo->SetBinContent(i, integ);
	}
	histo->SetLineColor( col );
	histo->SetLineWidth(1);
	return histo;
}

void fitRange(const string & fileName, const string & varName, Fit & fit, bool useWeights, TH1D * histo ) {
	std::ifstream input( fileName.c_str() );
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

	RooRealVar pt("pt", "pt", fit.low, fit.high, "GeV");
	RooRealVar weight("weight", "weight", 1, 0, 1000);
	RooRealVar a("a", "a", 0, fit.alow, fit.ahigh);
	RooRealVar b("b", "b", -0.01, fit.blow, fit.bhigh);
	RooRealVar c("c", "c", 0, fit.clow, fit.chigh);
	RooZPtPdf ptPDF("ptPDF", "ptPDF", pt, a, b, c);
	RooDataSet ptEvents("ptEvents", "ptEvents", RooArgSet(pt, weight), RooFit::WeightVar("weight"));

	double nEvInRange = 0;
	for ( unsigned i = 0; i < samples.size(); ++i ) {
		TFile in(  samples[i].c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		unsigned long entries = tree->GetEntries();
		for ( unsigned long i = 0; i < entries; i++ ) {
			tree->GetEntry( i );
			double tempPt = ev.getSVV<double>(varName);
			double tempWeight = 1;
			if (useWeights)
				tempWeight = ev.getSVV<double>("Weight");
			if (tempWeight > 1000) {
				cout << "ERROR: Weight larger than 100!" << endl;
				exit( EXIT_FAILURE );
			}
			if ( tempPt >= fit.low && tempPt < fit.high ) {
				pt = tempPt;
				ptEvents.add(RooArgSet(pt), tempWeight);
				if (histo)
					histo->Fill( pt.getVal(), tempWeight );
				nEvInRange += tempWeight;
			}
		}
	}

	RooFitResult * fitR = ptPDF.fitTo(ptEvents, RooFit::Range(fit.low, fit.high), RooFit::Save(true), RooFit::SumW2Error(true) );

	RooRealVar * par1 = (RooRealVar *) fitR->floatParsFinal().find("a");
	RooRealVar * par2 = (RooRealVar *) fitR->floatParsFinal().find("b");
	RooRealVar * par3 = (RooRealVar *) fitR->floatParsFinal().find("c");
	fit.a = par1->getVal();
	fit.b = par2->getVal();
	fit.c = par3->getVal();

	fit.ptFun->SetParameters( 1, fit.a, fit.b, fit.c );
	double integral = fit.ptFun->Integral(fit.low, fit.high);
	fit.N = nEvInRange / integral;
	fit.ptFun->SetParameters( fit.N, fit.a, fit.b, fit.c );
	fit.ptFun->SetRange(fit.low, fit.high);
}

void fitPt(vector<Fit> & fits, const char * fileName, bool useWeights, const char * outName) {
	int nBins = 130;
	double maxPt = 650;
	TH1D histo("ptHisto", "ptHisto", nBins, 0, maxPt);
	histo.Sumw2();

//	Fit fit1(25, 35);
//	fitRange(fileName, "ZPT", fit1, useWeights, &histo);
//	fits.push_back( fit1 );

//	Fit fit2(35, 55);
//	fitRange(fileName, "ZPT", fit2, useWeights, &histo);
//	fits.push_back( fit2 );

	Fit fit3(55, 80);
	fitRange(fileName, "ZPT", fit3, useWeights, &histo);
	fits.push_back( fit3 );

	Fit fit4(80, 95);
	fitRange(fileName, "ZPT", fit4, useWeights, &histo);
	fits.push_back( fit4 );

	Fit fit5(95,170);
	fitRange(fileName, "ZPT", fit5, useWeights, &histo);
	fits.push_back( fit5 );

	Fit fit6(170, maxPt);
	fit6.ahigh = 0.000001;
	fitRange(fileName, "ZPT", fit6, useWeights, &histo);
	fits.push_back( fit6 );

	TCanvas canv("canv", "canv", 800, 600);
	canv.SetLogy();

	histo.SetMarkerSize(0.5);
	histo.SetMarkerStyle(8);
	histo.Draw("PE");

//	fit1.buildHisto(nBins, 0, maxPt, kRed)->Draw("SAME");
//	fit2.buildHisto(nBins, 0, maxPt, kBlue)->Draw("SAME");
	fit3.buildHisto(nBins, 0, maxPt, kOrange)->Draw("SAME");
	fit4.buildHisto(nBins, 0, maxPt, kViolet)->Draw("SAME");
	fit5.buildHisto(nBins, 0, maxPt, kTeal)->Draw("SAME");
	fit6.buildHisto(nBins, 0, maxPt, kRed)->Draw("SAME");
	histo.Draw("SAMEPE");
	
	canv.SaveAs(outName);
}

void fitMass(const std::string & inputFileName, const std::string & outputFileName) {
	double mMin = 60;
	double mMax = 120;
	RooRealVar mass("mass", "mass", mMin, mMax, "GeV");

	RooRealVar meanBW("meanBW", "meanBW", 91.1876);
	RooRealVar widthBW("widthBW", "widthBW", 2.4952);
	RooBreitWigner signal("signal", "signal", mass, meanBW, widthBW);

	RooRealVar meanCB("meanCB", "meanCB", -1.5, -5, 5);
	RooRealVar sigmaCB("sigmaCB", "sigmaCB", 1.5, 0.0, 5.0);
	RooRealVar alphaCB("alphaCB", "alphaCB", 1.5, 0.1, 5.0);
	RooRealVar nCB("nCB", "nCB", 1.5, 0.0, 10.0);
	RooCBShape resolution("resolution", "resolution", mass, meanCB, sigmaCB, alphaCB, nCB);

	RooRealVar c("c", "c", -0.05, -0.1, 0.0);
	RooExponential background("background", "background", mass, c);

	mass.setBins(10000, "cache");
	RooFFTConvPdf sigRes("sigRes", "sigRes", mass, signal, resolution);

	RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.05, 0.0, 0.2);
	RooAddPdf massPDF("massPDF", "massPDF", background, sigRes, bkgfrac);

	RooDataSet massEvents("massEvents", "massEvents", RooArgSet(mass));

	TFile input( inputFileName.c_str() );
	TTree * tree = ( TTree * ) input.Get( "HZZ2l2nuAnalysis" );
	Event ev(tree);
	unsigned long entries = tree->GetEntries();
	for ( unsigned long i = 0; i < entries; i++ ) {
		tree->GetEntry( i );
		double tmpMass = ev.getSVV<double>("ZMASS");
		if ( tmpMass >= mMin && tmpMass < mMax ) {
			mass = tmpMass;
			massEvents.add(RooArgSet(mass));
		}
	}
	RooFitResult * fitR = massPDF.fitTo(massEvents, RooFit::Range(mMin, mMax), RooFit::Save(true), RooFit::SumW2Error(true) );

	RooWorkspace * w = new RooWorkspace("zMassFit", "zMassFit");
	w->import(massPDF);
	w->import(massEvents);
	w->import(*fitR);

	w->writeToFile( outputFileName.c_str() );
}

RooAbsPdf * readFitMass(const string & inputFileName, const string & plotFileName) {
	TFile f( inputFileName.c_str() );
	RooWorkspace * w = (RooWorkspace*) f.Get("zMassFit");

	RooRealVar * zmass = w->var("mass");
	RooAbsPdf * pdf = w->pdf("massPDF");
	RooAbsData * data = w->data("massEvents");

	int nentries = data->numEntries();
	RooDataSet * events = pdf->generate(*zmass, nentries);

	RooPlot * xframe = zmass->frame(RooFit::Title("Z mass [GeV]")) ;
	data->plotOn(xframe) ;
	events->plotOn(xframe, RooFit::MarkerColor(kRed)) ;
	pdf->plotOn(xframe) ;
	TCanvas canv("canv","canv", 1600, 1200);
	gPad->SetLeftMargin(0.15);
	xframe->GetYaxis()->SetTitleOffset(1.4);
	xframe->Draw();
	canv.SetLogy();
	canv.SaveAs(plotFileName.c_str());

	return pdf;
}
