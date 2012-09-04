/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOZPTPDF
#define ROOZPTPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooZPtPdf : public RooAbsPdf {
	public:
		RooZPtPdf() {} ; 
		RooZPtPdf(const char *name, const char *title,
				RooAbsReal & _zpt,
				RooAbsReal & _a,
				RooAbsReal & _b,
				RooAbsReal & _c);
		RooZPtPdf(const RooZPtPdf& other, const char* name = 0) ;
		virtual TObject* clone(const char* newname) const { return new RooZPtPdf(*this, newname); }
		inline virtual ~RooZPtPdf() { }
	protected:
		RooRealProxy zpt ;
		RooRealProxy a ;
		RooRealProxy b ;
		RooRealProxy c ;

		Double_t evaluate() const ;
	private:
//		ClassDef(RooZPtPdf, 1) // Your description goes here...
};

#endif
