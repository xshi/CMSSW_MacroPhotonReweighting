/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooZPtPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

//ClassImp(RooZPtPdf);

RooZPtPdf::RooZPtPdf(const char *name, const char *title, RooAbsReal & _zpt, RooAbsReal & _a, RooAbsReal & _b, RooAbsReal & _c) :
						RooAbsPdf(name,title), zpt("zpt", "zpt", this, _zpt), a("a", "a", this, _a),
						b("b", "b", this, _b), c("c", "c", this, _c) {} 

RooZPtPdf::RooZPtPdf(const RooZPtPdf& other, const char* name) : 
						RooAbsPdf(other,name), zpt("zpt",this,other.zpt), a("a",this,other.a),
						b("b",this,other.b), c("c",this,other.c) {}



Double_t RooZPtPdf::evaluate() const { 
	// ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
	return TMath::Exp(a*zpt*zpt+b*TMath::Sqrt(zpt*zpt*zpt)+c*zpt);
}
