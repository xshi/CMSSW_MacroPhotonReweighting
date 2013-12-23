#include <iomanip>
#include "event.h"
#include "fits.h"
#include <fstream>
#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooCMSShape.h>
#include <RooGaussian.h>
#include <RooParametricStepFunction.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFermi.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include "RooZPtPdf.h"
#include "samples.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include "toolbox.h"
#include <TTree.h>
#include <vector>
#include <options.h>
#include <RooNumIntConfig.h>
#include <RooMinuit.h>
#include "Math/MinimizerOptions.h"
#include <TLine.h>
#include "TGraph.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cin;

bool fitPt(const vector<PhotonSample> & photonSamples, bool useWeights, bool isData,
		unsigned cat, const string & fitName) {

	RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-16);
	RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-16);
//	RooAbsReal::RooMinuit::setEps(1e-16);

	Options::getInstance().addDoubleOption("minNLL", 1e15);

	string workspaceName = fitName + "_c" + double2string(cat);
	string command = "mkdir -p " + workspaceName;
	if (system( command.c_str() ))
		throw string("ERROR: Can't create directory!");
	ofstream outLog( workspaceName + "/log.txt", std::ios_base::app );
	outLog.precision(10);
	time_t rawtime;
	struct tm * timeinfo;
	time( &rawtime );
	timeinfo = localtime(&rawtime);
	outLog << endl;
	outLog << asctime(timeinfo);

	ifstream input( workspaceName + "/fitOptions.txt" );
//	bool optionsExist = false;
	if (input.is_open()) {
//		optionsExist = true;
		Options::getInstance().readInOptions( workspaceName + "/fitOptions.txt" );
	}
	input.close();

	double xmin = 55;
	double xmax = 755;

	RooRealVar pt("pt", "pt", xmin, xmax, "GeV");
	RooRealVar weight("weight", "weight", 1, 0, 1e4);

	vector<RooRealVar *> parameters;
	parameters.push_back(new RooRealVar("p1", "p1", 8.4, 0, 1e3));
//	parameters.push_back(new RooRealVar("p1", "p1", 8.4));

	parameters.push_back(new RooRealVar("p2", "p2", 3.6, 0, 1e2));
//	parameters.push_back(new RooRealVar("p2", "p2", 3.6));

	parameters.push_back(new RooRealVar("p3", "p3", 10, 0, 1e3));
//	parameters.push_back(new RooRealVar("p3", "p3", 0.0));

	parameters.push_back(new RooRealVar("p4", "p4", 10, 0, 1e2));
//	parameters.push_back(new RooRealVar("p4", "p4", 11.27));

	parameters.push_back(new RooRealVar("p5", "p5", 1.0, 0.01, 10));
//	parameters.push_back(new RooRealVar("p5", "p5", 0.0));

	Options::getInstance().addDoubleOption("p1_var", 10.0);
	Options::getInstance().addDoubleOption("p2_var", 1.0);
	Options::getInstance().addDoubleOption("p3_var", 10.0);
	Options::getInstance().addDoubleOption("p4_var", 1.0);
	Options::getInstance().addDoubleOption("p5_var", 0.1);

//	if (optionsExist) {
//		for (unsigned i = 0; i < 2; ++i) {
//			string parName = parameters[i]->GetName();
//			double centralV = 0.0;
//			double parVar = Options::getInstance().checkDoubleOption(parName + "_var");
//			try {
//				centralV = Options::getInstance().checkDoubleOption(parName);
//			} catch (const std::string & exc) {
//				centralV = parameters[i]->getVal();
//			}
//			parameters[i]->setVal(centralV);
//			if (centralV - parVar > parameters[i]->getMin())
//				parameters[i]->setMin(centralV - parVar);
//			else if (centralV - parVar / 10 < parameters[i]->getMin())
//				outLog << "Parameter is close to limit: " << parName << endl;
//			if (centralV + parVar < parameters[i]->getMax())
//				parameters[i]->setMax(centralV + parVar);
//			else if (centralV + parVar / 10 > parameters[i]->getMax())
//				outLog << "Parameter is close to limit: " << parName << endl;
//		}
//	}

	RooZPtPdf * ptPDF = new RooZPtPdf("ptPDFRaw", "ptPDFRaw", pt, *parameters[0], *parameters[1], *parameters[2], *parameters[3], *parameters[4]);

	RooDataSet ptEvents("ptEvents", "ptEvents", RooArgSet(pt, weight), RooFit::WeightVar("weight"));
	TH1F hist("ratio", "ratio", 100, xmin, xmax);
	hist.Sumw2();
	
	
	for ( unsigned i = 0; i < photonSamples.size(); ++i ) {
		string inFileName = photonSamples[i].getOutputDir() + "/" + photonSamples[i].getSampleName() + photonSamples[i].getPreselectedSampleName() + ".root";
		TFile in(  inFileName.c_str() );
		TTree * tree = ( TTree * ) in.Get( "HZZ2l2nuAnalysis" );
		Event ev(tree);
		double weight = photonSamples[i].getCrossSection() / photonSamples[i].getNumberOfEvents();
		unsigned long entries = tree->GetEntries();
		int * category = ev.getSVA<int>("CATEGORY");
		for ( unsigned long i = 0; i < entries; i++ ) {
			tree->GetEntry( i );
			double tmpPt = ev.getSVV<double>("ZPT");
			double tmpWeight = 1;
			if (useWeights) {
				if (isData)
					tmpWeight = ev.getSVV<double>("Weight");
					//tmpWeight = ev.getSVV<double>("Weight") * 1.0e-3;
				else
					tmpWeight = weight;
			}
			if (tmpWeight > 1000)
				throw string("ERROR: Weight larger than 1!");
			if ( tmpWeight != 0 && tmpPt >= xmin && tmpPt < xmax && *category == int(cat)) {
				pt = tmpPt;
				ptEvents.add(RooArgSet(pt), tmpWeight);
				hist.Fill(tmpPt, tmpWeight);
			}
		}
	}

	RooPlot * xframe = pt.frame(RooFit::Title("G P_T [GeV]"));
	ptEvents.plotOn(xframe, RooFit::DataError(RooAbsData::SumW2));

	pt.setRange("R", xmin, xmax);
	//RooFitResult * fitR = ptPDF->fitTo(*(ptEvents.binnedClone()), RooFit::Range("R"), RooFit::Save(true), RooFit::SumW2Error(true) );
	RooFitResult * fitR = ptPDF->fitTo(ptEvents, RooFit::Range("R"), RooFit::Save(true), RooFit::SumW2Error(true), RooFit::NumCPU(8), RooFit::Minimizer("Minuit"), RooFit::Strategy(2) );

	for (unsigned i = 0; i < parameters.size(); ++i) {
		string parName = parameters[i]->GetName();
		double centralV = parameters[i]->getVal();
		double error = parameters[i]->getError();
		parameters[i]->setMin(centralV - 10 * error);
		parameters[i]->setMax(centralV + 10 * error);
	}

	TCanvas canv("canv","canv", 1600, 1200);
	for (unsigned k = 0; k < parameters.size(); ++k) {
		for (unsigned l = 0; l < k; ++l) {
			if (!fitR->floatParsFinal().find(parameters[k]->GetName()) || !fitR->floatParsFinal().find(parameters[l]->GetName()))
				continue;
			double xc = parameters[k]->getVal();
			double xe = parameters[k]->getError();
			double yc = parameters[l]->getVal();
			double ye = parameters[l]->getError();
			RooPlot * corr = new RooPlot(xc-2*xe, xc+2*xe, yc-2*ye, yc+2*ye);
			fitR->plotOn(corr, *parameters[k], *parameters[l], "MEABHV12");
			corr->Draw();
			canv.SaveAs( (workspaceName + "/correlations_" + string(parameters[k]->GetName()) + "_" + string(parameters[l]->GetName()) + ".ps").c_str() );
		}
	}

	TF1 fun = readInFitParameters(fitR, workspaceName, xmin, xmax, parameters);
	cout << "Calculating integral ..." << endl;
	double integral = fun.Integral(xmin, xmax);
	cout << "Done!" << endl;
	string cut = "pt>" + double2string(xmin) + " && pt<" + double2string(xmax);
	double nEv = ptEvents.sumEntries(cut.c_str());
	double normalization = nEv / integral * 7.0;
	fun.SetParameter(5, normalization);

	TF1 part1("part1", "[0] * TMath::Power(1 + x / [1], -[2])", xmin, xmax);
	part1.SetParameter(0, fun.GetParameter(5));
	part1.SetParameter(1, fun.GetParameter(0));
	part1.SetParameter(2, fun.GetParameter(1));

	TF1 part2("part2", "[0] * TMath::Power(1 + x / [1], -[2])", xmin, xmax);
	part2.SetParameter(0, fun.GetParameter(5) * fun.GetParameter(4));
	part2.SetParameter(1, fun.GetParameter(2));
	part2.SetParameter(2, fun.GetParameter(3));

	for (int i = 1; i <= hist.GetNbinsX(); ++i) {
		double ratio = hist.GetBinContent(i) / fun.Eval( hist.GetBinCenter(i) );
		hist.SetBinContent(i, ratio);
	}
	hist.SetMaximum(2);
	hist.SetMinimum(0);
	hist.SetMarkerStyle(9);	
	hist.SetMarkerSize(2);
	hist.Draw("HISTP");
	canv.SaveAs( (workspaceName + "/ratio.ps").c_str() );
	bool improved = true;
	bool same = false;
	try {
		double oldMinNLL = Options::getInstance().checkDoubleOption("minNLL");
		double minNLL = fitR->minNll();
		if (fabs((minNLL - oldMinNLL) / oldMinNLL) < 1e-9) {
			same = true;
			improved = false;
		} else if (minNLL > oldMinNLL)
			improved = false;
	} catch (const string & exc) {
	}

	if (improved || same) {
		ptPDF->plotOn(xframe, RooFit::LineColor(kBlue));
		xframe->addObject(fun.Clone(), "SAME");
		part1.SetLineColor(kOrange);
		xframe->addObject(part1.Clone(), "SAME");
		part2.SetLineColor(kViolet);
		xframe->addObject(part2.Clone(), "SAME");
		ptEvents.plotOn(xframe, RooFit::DataError(RooAbsData::SumW2));

		gPad->SetLeftMargin(0.15);
		xframe->GetYaxis()->SetTitleOffset(1.4);
		xframe->SetMaximum(1e8);
		xframe->SetMinimum(1e-1);
		xframe->Draw();
		canv.SetLogy();
		//canv.SetLogx();
		canv.SaveAs((workspaceName + "/fit.ps").c_str());

		ofstream output( workspaceName + "/fitOptions.txt" );
		output.precision(10);
		output << "minNLL=" << fitR->minNll() << endl;
		for (unsigned i = 0; i < parameters.size(); ++i) {
			string parName = parameters[i]->GetName();
			RooRealVar * par = (RooRealVar *) fitR->floatParsFinal().find(parName.c_str());
			if (!par)
				continue;
			output << parName << "=" << par->getVal() << endl;
		}

		TGraph lowBand;
		TGraph highBand;
		if (same) {
			vector<double> minErrV;
			vector<double> maxErrV;
			vector<double> xV;
			for (double i = 55; i <= 755; i += 10) {
				xV.push_back(i);
				double funVal = fun.Eval(i);
				maxErrV.push_back(funVal);
				minErrV.push_back(funVal);
			}

			ofstream test("test.txt");

			RooPlot * xframe1 = pt.frame(RooFit::Title("G P_T [GeV]"));
			ptEvents.plotOn(xframe1);
			fun.SetLineColor(kBlue);
			xframe1->addObject(fun.Clone(), "SAME");
			for (int i = 0; i < 100; ++i) {
				RooAbsData * tmpData = bootstrap(ptEvents, "bootstrap" + double2string(i));
				RooAbsPdf * tmpPdf = (RooAbsPdf *) ptPDF->clone( ("bootstrap" + double2string(i)).c_str() );
				RooFitResult * tmpFitR = tmpPdf->fitTo(*tmpData, RooFit::Range("R"), RooFit::Save(true), RooFit::SumW2Error(true), RooFit::NumCPU(8), RooFit::Minimizer("Minuit"), RooFit::Strategy(2) );
				TF1 tmpFun = readInFitParameters(tmpFitR, "bootstrap" + double2string(i), xmin, xmax, parameters);
				double tmpIntegral = tmpFun.Integral(xmin, xmax);
				tmpFun.SetParameter(5, nEv / tmpIntegral * 7.0);
				for (unsigned i = 0; i < xV.size(); i++) {
					double tmpPt = xV[i];
					double funVal = tmpFun.Eval(tmpPt);
					if (funVal > maxErrV[i])
						maxErrV[i] = funVal;
					if (funVal < minErrV[i])
						minErrV[i] = funVal;
				}
				test << std::setw(5) << i << std::setw(15) << tmpIntegral << std::setw(10) << tmpData->sumEntries() << endl;

				//ptPDF->plotOn(xframe1, RooFit::LineColor(kOrange));

				RooPlot * xframe2 = pt.frame(RooFit::Title("G P_T [GeV]"));
				tmpData->plotOn(xframe2);
				ptPDF->plotOn(xframe2, RooFit::LineColor(kBlue));
				xframe2->SetMaximum(1e8);
				xframe2->SetMinimum(1e-1);
				xframe2->Draw();
				canv.SetLogy();
				canv.SetLogx();
				canv.SaveAs( (workspaceName + "/bootstrap" + double2string(i) + ".ps").c_str() );
				delete xframe2;
			}
			double xA[xV.size()];
			double minErrA[xV.size()];
			double maxErrA[xV.size()];
			for (unsigned j = 0; j < xV.size(); j++) {
				xA[j] = xV[j];
				minErrA[j] = minErrV[j];
				maxErrA[j] = maxErrV[j];
			}
			lowBand = TGraph(xV.size(), xA, minErrA);
			highBand = TGraph(xV.size(), xA, maxErrA);
			lowBand.SetLineColor(kRed);
			xframe1->addObject(lowBand.Clone(), "SAME");
			highBand.SetLineColor(kRed);
			xframe1->addObject(highBand.Clone(), "SAME");

			xframe1->SetMaximum(1e8);
			xframe1->SetMinimum(1e-1);
			xframe1->Draw();
			canv.SetLogy();
			canv.SetLogx();
			canv.SaveAs( (workspaceName + "/errors.ps").c_str() );
		}

		RooWorkspace * w = new RooWorkspace(workspaceName.c_str(), workspaceName.c_str());
		w->import( (TObject &) (*ptPDF), "FitPDF", kFALSE);
		w->import(*fitR, "FitResults", kFALSE);
		w->import(fun, "FitFunction", kFALSE);
		w->import(lowBand, "FitFunctionDown", kFALSE);
		w->import(highBand, "FitFunctionUp", kFALSE);
		w->writeToFile( (workspaceName + "/workspace.root").c_str() );

	}

	outLog << "minNLL: " << fitR->minNll() << endl;
	fun.Print();
	cout << "status: " << fitR->status() << endl;
	cout << "covQual: " << fitR->covQual() << endl;
	cout << "minNLL: " << fitR->minNll() << endl;
	cout << "improved: " << improved << endl;
	cout << "same: " << same << endl;
	fitR->Print("v");
	fitR->correlationMatrix().Print("v");

	return improved;
}

RooAbsData * bootstrap(const RooDataSet & dataset, const string & newName) {
	RooAbsData * newData = dataset.emptyClone(newName.c_str());
	int nEntries = dataset.numEntries();
	srand(time(NULL));
	for (int i = 0; i < nEntries; ++i) {
		int j = rand() % nEntries;
		const RooArgSet * tmpRow = dataset.get(j);
		double tmpWeight = dataset.weight();
		newData->add(*tmpRow, tmpWeight);
	}
	return newData;
}

Double_t ptSpectrum(Double_t * x, Double_t * par) {
	Double_t zpt = x[0];
	Double_t N = par[5];
	return N*ptFunc(zpt, 5, par);
}

TF1 readInFitParameters(const RooFitResult * fitR, const string & name, double xmin, double xmax, const vector<RooRealVar *> & params ) {
	Double_t (*funPtr)(Double_t * x, Double_t * par) = ptSpectrum;
	const int nPar = params.size();
	Double_t parameters[nPar + 1];
	TF1 ptFun(name.c_str(), funPtr, xmin, xmax, nPar + 1);
	for (int i = 0; i < nPar; i++) {
		std::stringstream ss;
		ss << 'p' << i + 1;
		RooRealVar * tmpPar = (RooRealVar *) fitR->floatParsFinal().find(ss.str().c_str());
		if (tmpPar)
			parameters[i] = tmpPar->getVal();
		else
			parameters[i] = params[i]->getVal();
	}
	parameters[nPar] = 1.0;
	ptFun.SetParameters( parameters );
	return ptFun;
}

void fitMass(const std::string & inputFileName, const std::string & outputFileName) {
	double mMin = 76.0;
	double mMax = 106.0;
	RooRealVar mass("mass", "mass", mMin, mMax, "GeV");

	RooRealVar meanBW("meanBW", "meanBW", 91.1876);
	RooRealVar widthBW("widthBW", "widthBW", 2.4952);
	RooBreitWigner signal("signal", "signal", mass, meanBW, widthBW);

	RooRealVar meanCB("meanCB", "meanCB", -1.5, -5, 5);
	RooRealVar sigmaCB("sigmaCB", "sigmaCB", 1.5, 0.0, 5.0);
	RooRealVar alphaCB("alphaCB", "alphaCB", 1.5, 0.1, 5.0);
	RooRealVar nCB("nCB", "nCB", 1.5, 0.0, 10.0);
	RooCBShape resolution("resolution", "resolution", mass, meanCB, sigmaCB, alphaCB, nCB);

	RooRealVar alpha("alpha", "alpha", 60.0, 30.0, 150.0);
	RooRealVar beta("beta", "beta", 0.1, 0.0, 0.4);
	RooRealVar gamma("gamma", "gamma", 0.1, 0.0, 0.4);
	RooRealVar peak("peak", "peak", 90.0);
	RooCMSShape background("background", "background", mass, alpha, beta, gamma, peak);

	mass.setBins(10000, "cache");
	RooFFTConvPdf sigRes("sigRes", "sigRes", mass, signal, resolution);

	RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.05, 0.0, 0.3);
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
	RooFitResult * fitR = massPDF.fitTo(massEvents, RooFit::Range(mMin, mMax), RooFit::Save(true), RooFit::SumW2Error(true), RooFit::NumCPU(8) );

	RooWorkspace w("zMassFit", "zMassFit");
	w.import(massPDF);
	w.import(massEvents);
	w.import(*fitR);

	w.writeToFile( outputFileName.c_str() );
}

RooWorkspace * readMassFit(const string & inputFileName, const string & plotFileName) {
	TFile f( inputFileName.c_str() );
	RooWorkspace * w = (RooWorkspace*) f.Get("zMassFit");

	RooRealVar * zmass = w->var("mass");
	RooAbsPdf * pdf = w->pdf("massPDF");
	RooAbsData * data = w->data("massEvents");

	RooPlot * xframe = zmass->frame(RooFit::Title("Z mass [GeV]")) ;
	data->plotOn(xframe) ;
	pdf->plotOn(xframe) ;
	TCanvas canv("canv","canv", 1600, 1200);
	gPad->SetLeftMargin(0.15);
	xframe->GetYaxis()->SetTitleOffset(1.4);
	xframe->Draw();
	canv.SetLogy();
	canv.SaveAs(plotFileName.c_str());

	return w;
}

TGraph readPtFit(const string & inputFileName, unsigned category, double low, double high, const std::string & fitName) {
	string name = inputFileName + "_c" + double2string(category);
	string filename = name + "/workspace.root";
	TFile f( filename.c_str() );
	RooWorkspace * w = (RooWorkspace*) f.Get(name.c_str());
//	string funName = inputFileName + "_" + double2string(low) + "_" + double2string(high)
//		+ "_c" + double2string(category);
//	w->Print();
	TObject * objPtr = w->obj(fitName.c_str());
	TF1 * funPtr = dynamic_cast<TF1 *>(objPtr);
	if (funPtr) {
		return TGraph(funPtr);
	}
	TGraph * graphPtr = dynamic_cast<TGraph *>(objPtr);
	if (graphPtr) {
		return TGraph(*graphPtr);
	}
	throw std::string("Can't interpret a fit!");
}
