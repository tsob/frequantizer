//-----------------------------------------------------
// name: "Frequency Domain Quantizer"
// version: "0.1"
// author: "Tim O'Brien"
//
// Code generated with Faust 0.9.62 (http://faust.grame.fr)
//-----------------------------------------------------
#ifndef FAUSTFLOAT
#define FAUSTFLOAT float
#endif

typedef long double quad;
/* link with  */

#include <math.h>

#ifndef FAUSTCLASS
#define FAUSTCLASS mydsp
#endif

class mydsp : public dsp {
  private:
	FAUSTFLOAT 	fslider0;
	float 	fRec0[2];
	FAUSTFLOAT 	fslider1;
	FAUSTFLOAT 	fslider2;
	FAUSTFLOAT 	fcheckbox0;
  public:
	static void metadata(Meta* m) 	{
		m->declare("name", "Frequency Domain Quantizer");
		m->declare("version", "0.1");
		m->declare("author", "Tim O'Brien");
		m->declare("music.lib/name", "Music Library");
		m->declare("music.lib/author", "GRAME");
		m->declare("music.lib/copyright", "GRAME");
		m->declare("music.lib/version", "1.0");
		m->declare("music.lib/license", "LGPL with exception");
		m->declare("math.lib/name", "Math Library");
		m->declare("math.lib/author", "GRAME");
		m->declare("math.lib/copyright", "GRAME");
		m->declare("math.lib/version", "1.0");
		m->declare("math.lib/license", "LGPL with exception");
	}

	virtual int getNumInputs() 	{ return 1; }
	virtual int getNumOutputs() 	{ return 1; }
	static void classInit(int samplingFreq) {
	}
	virtual void instanceInit(int samplingFreq) {
		fSamplingFreq = samplingFreq;
		fslider0 = 0.0f;
		for (int i=0; i<2; i++) fRec0[i] = 0;
		fslider1 = 4.0f;
		fslider2 = 6e+01f;
		fcheckbox0 = 0.0;
	}
	virtual void init(int samplingFreq) {
		classInit(samplingFreq);
		instanceInit(samplingFreq);
	}
	virtual void buildUserInterface(UI* interface) {
		interface->openVerticalBox("dummy");
		interface->addHorizontalSlider("FFT block size (power of two)", &fslider2, 6e+01f, 6.0f, 18.0f, 1.0f);
		interface->addHorizontalSlider("Level (db)", &fslider0, 0.0f, -96.0f, 4.0f, 0.1f);
		interface->addCheckButton("On", &fcheckbox0);
		interface->addHorizontalSlider("Quantization bits", &fslider1, 4.0f, 2.0f, 16.0f, 1.0f);
		interface->closeBox();
	}
	virtual void compute (int count, FAUSTFLOAT** input, FAUSTFLOAT** output) {
		float 	fSlow0 = (0.0010000000000000009f * powf(10,(0.05f * float(fslider0))));
		float 	fSlow1 = ((float(fcheckbox0) * float(fslider2)) * float(fslider1));
		FAUSTFLOAT* input0 = input[0];
		FAUSTFLOAT* output0 = output[0];
		for (int i=0; i<count; i++) {
			fRec0[0] = (fSlow0 + (0.999f * fRec0[1]));
			output0[i] = (FAUSTFLOAT)(fSlow1 * ((float)input0[i] * fRec0[0]));
			// post processing
			fRec0[1] = fRec0[0];
		}
	}
};


