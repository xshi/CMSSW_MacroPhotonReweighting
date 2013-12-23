#ifndef SAMPLES_H
#define SAMPLES_H

#include <string>
#include <vector>
#include <jsoncpp/json/value.h>

class PhotonSample {
	public :
//		PhotonSample( const std::string & sampleName_, const std::string & inputDir_, const std::string & outputDir_,
//				const std::string & preselName_, const std::string & ptRewName_, const std::string & njRewName_,
//				const std::string & nvRewName_, unsigned n, double xs) :
//			sampleName(sampleName_), inputDir(inputDir_), outputDir(outputDir_), preselName(preselName_),
//			ptRewName(ptRewName_), njRewName(njRewName_), nvRewName(nvRewName_), nEvt(n), xsection(xs) {}
		PhotonSample( const Json::Value & val );
		double getCrossSection() const {return xsection;}
		std::string getSampleName() const {return sampleName;}
		std::string getInputDir() const {return inputDir;}
		std::string getOutputDir() const {return outputDir;}
		std::string getPreselectedSampleName() const {return preselName;}
		std::string getPtRewSampleName() const {return ptRewName;}
		std::string getNjRewSampleName() const {return njRewName;}
		std::string getNvRewSampleName() const {return nvRewName;}
		std::string getNormRewSampleName() const {return normRewName;}
		unsigned getNumberOfEvents() const {return nEvt;}
	private :
		std::string sampleName;
		std::string inputDir;
		std::string outputDir;
		std::string preselName;
		std::string ptRewName;
		std::string njRewName;
		std::string nvRewName;
		std::string normRewName;
		unsigned nEvt;
		double xsection;
};

//std::vector<PhotonSample> readInPhotonSamples(const std::string & fileName);

#endif
