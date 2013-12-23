#include "samples.h"
#include <fstream>
#include <sstream>

using std::string;
using std::vector;
using std::ifstream;
using std::stringstream;

//vector<PhotonSample> readInPhotonSamples(const string & fileName) {
//	ifstream input( fileName.c_str() );
//	if( !input.is_open() )
//		throw string("ERROR: Can't open file with samples (" + fileName + ")!");
//
//	vector<PhotonSample> samples;
//	while( !input.eof() ) {
//		string temp;
//		getline( input, temp );
//		if( !temp.size() || temp[0] == '\n' || temp[0] == '#' ) {
//			continue;
//		}
//
//		stringstream in( temp );
//		string sampleName;
//		string inputDir;
//		string outputDir;
//		string preselName;
//		string ptRewName;
//		string njRewName;
//		string nvRewName;
//		string normRewName;
//		double xsection;
//		unsigned nEvt;
//		in >> sampleName >> inputDir >> outputDir >> preselName >> ptRewName >> njRewName >> nvRewName >> xsection >> nEvt;
//		samples.push_back( PhotonSample(sampleName, inputDir, outputDir, preselName, ptRewName, njRewName, nvRewName, nEvt, xsection) );
//	}
//	return samples;
//}

PhotonSample::PhotonSample( const Json::Value & value ) :
		sampleName(value.get("sampleName", "").asString()),
		inputDir(value.get("inputDir", "").asString()),
		outputDir(value.get("outputDir", "").asString()),
		preselName(value.get("preselName", "").asString()),
		ptRewName(value.get("ptRewName", "").asString()),
		njRewName(value.get("njRewName", "").asString()),
		nvRewName(value.get("nvRewName", "").asString()),
		normRewName(value.get("normRewName", "").asString()),
		nEvt(value.get("nEvt", 1).asInt()),
		xsection(value.get("xsection", 1).asDouble())
{}
