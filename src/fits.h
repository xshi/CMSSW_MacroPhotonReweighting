#ifndef FITS_H
#define FITS_H

#include <TROOT.h>
#include <string>
#include <RooAbsPdf.h>

class TF1;
class TH1D;
class RooVoigtian;
class RooRealVar;

struct Fit {
	double low;
	double high;
	double a;
	double b;
	double c;
	double N;
	TF1 * ptFun;
	double alow;
	double ahigh;
	double blow;
	double bhigh;
	double clow;
	double chigh;
	TH1D * histo;

	Fit( double l, double h );
	Fit( const Fit & fit );
	~Fit();
	Fit & operator=( const Fit & fit );
	TH1D * buildHisto(int nBins, double lowE, double highE, Color_t col = kBlack);
};

void fitRange(const std::string & fileName, const std::string & varName, Fit & fit, bool useWeights = false, TH1D * histo = 0 );
void fitPt(std::vector<Fit> & fits, const char * fileName, bool useWeights, const char * outName);

void fitMass(const std::string & inputFileName, const std::string & outputFileName);
RooAbsPdf * readFitMass(const std::string & inputFileName, const std::string & plotFileName);

#endif
