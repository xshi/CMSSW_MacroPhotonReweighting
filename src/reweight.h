#ifndef REWEIGHT_H
#define REWEIGHT_H

#include <vector>
#include "fits.h"

double calculatePtWeight(const std::vector<Fit> & zfits, const std::vector<Fit> & gfits, double pt);
void reweightPt(const char * inFile, const char * outFile, const std::vector<Fit> & zfits, const std::vector<Fit> & gfits);
void fillHistogram(const std::string & fileName, const std::string & varName, const std::vector<std::string> & weights, TH1D * histo );
void reweight(const char * inFile, const char * outFile, const TH1D * ghisto, const TH1D * zhisto, const char * varName, const char * branchName);

#endif
