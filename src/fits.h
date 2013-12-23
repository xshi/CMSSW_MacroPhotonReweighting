#ifndef FITS_H
#define FITS_H

#include <string>
#include <TF1.h>

class RooWorkspace;
class RooRealVar;
class RooAbsPdf;
class RooDataSet;
class RooAbsData;
class PhotonSample;
class RooFitResult;
class Event;
class TF1;
class TGraph;

TGraph readPtFit(const std::string & inputFileName, unsigned cat, double xmin, double xmax, const std::string & fitName);
RooWorkspace * readMassFit(const std::string & inputFileName, const std::string & plotFileName);
void fitMass(const std::string & inputFileName, const std::string & outputFileName);
bool fitPt(const std::vector<PhotonSample> & photonSamples, bool useWeights, bool isData, unsigned category, const std::string & fitName);
TF1 fitPtRange(RooRealVar & pt, RooAbsPdf & ptPDF, RooDataSet & ptEvents, unsigned category, const std::string & fitName, double xmin, double xmax);
TF1 readInFitParameters(const RooFitResult * fitR, const std::string & name, double xmin, double xmax, const std::vector<RooRealVar *> & params);
Double_t ptSpectrum(Double_t * x, Double_t * par);
RooAbsData * bootstrap(const RooDataSet & dataset, const std::string &);
#endif
