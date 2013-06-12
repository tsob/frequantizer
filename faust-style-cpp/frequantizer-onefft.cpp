//---------------------------------------------------------------
// name: "Frequency Domain Quantizer"
// version: "0.1"
// author: "Tim O'Brien"
//
// Code generated with Faust 0.9.58 (http://faust.grame.fr)
// and then modified by Tim O'Brien with FFT-based functionality
//---------------------------------------------------------------
#ifndef FAUSTFLOAT
#define FAUSTFLOAT float
#endif

typedef long double quad;
// remember to link with -lfftw3 -lm

#ifndef FAUSTCLASS
#define FAUSTCLASS mydsp
#endif

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include "quantize.h"
#include "window.h"

#define TWOPOW(x) (1 << (x))

float quantize(int bits, float input, float maxval)
{
      unsigned long quantval = 0;
      int sign;
      float output;

      sign = (input < 0);
      input = fabs(input);

      //quantize to bit code
      if (input >= maxval) {
	  quantval = (long)(TWOPOW(bits - 1) - 1);
      } else {
	  quantval = (long)(
		(
		    (unsigned long)(TWOPOW(bits-1) +
		    (TWOPOW(bits-1)-1)) * input/maxval + 1.0
		) / 2.0
	  );
      }
      if (sign) { quantval |= TWOPOW(bits - 1); }

      //dequantize back to float
      if (quantval & TWOPOW(bits - 1)) {
	  sign = -1;
	  quantval &= TWOPOW(bits - 1) - 1;
      } else {
	  sign = 1;
      }
      output = sign * 2.0 * maxval * ((float)quantval / (TWOPOW(bits) - 1));

      return output;
}

void sinewindow(double* inputdata, int length)
{
      int i;
      for (i = 0; i < length; ++i)
      {
	  inputdata[i] = inputdata[i] * sin(M_PI*(i+0.5)/(double)length);
      }
}

class mydsp : public dsp {
  private:
	FAUSTFLOAT 	fslider0;   //Level (dB)
	float 	        fRec0[2];
	FAUSTFLOAT 	fslider1;   //Quantization bits
	FAUSTFLOAT 	fslider2;   //FFT block size (power of two)
	FAUSTFLOAT 	fcheckbox0; //mute button

	// My variables
	fftw_complex *fftout; //output fft array
	double *in, *out; //input & output sample arrays - real-valued signals
	fftw_plan pf, pb;
	int n = 8192; //fft window size - samples
	int nc = (n/2) + 1; //buffer size for complex fft bins (using fftw r2c)
	int halfn = n/2; //50% overlap
	int siglength = n * 4;
	double *insig, *outsig, *halfframe;
	int i, j, samp=0, fftcountin=0, fftcountout=0;
	int bits = 8;

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
		fslider0      = 0.0f;   //Level (dB)
		for (int i=0; i<2; i++) fRec0[i] = 0;
		fslider1      = 4.0f;   //Quantization bits
		fslider2      = 6e+01f; //FFT block size (power of two)
		fcheckbox0    = 0.0;    //mute button
	}

	// FFT init function
	virtual void FFTInit( void ) {
	  //data arrays for use with fftw
	  in     = (double*)       fftw_malloc( sizeof( double ) * n);
	  out    = (double*)       fftw_malloc( sizeof( double ) * n);
	  fftout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc);
	  //forward & backward transform plans
	  pf = fftw_plan_dft_r2c_1d(n, in, fftout, FFTW_MEASURE);
	  pb = fftw_plan_dft_c2r_1d(n, fftout, out, FFTW_MEASURE);
	  //overlap-add frame...
	  halfframe = (double*) fftw_malloc( sizeof(double) * halfn );
	  for (i = 0; i < halfn; ++i) halfframe[i] = 0.0;
	  //buffers for queueing input and output
	  insig  = (double*) fftw_malloc( sizeof(double) * n);
	  outsig = (double*) fftw_malloc( sizeof(double) * n);
	  for (i = 0; i < n; ++i) {
	    insig[i]  = 0.0;
	    outsig[i] = 0.0;
	  }
	}
	// FFT destructor function
	virtual void FFTDestruct( void ) {
	  printf("Called destructor function...\n");
	  fftw_destroy_plan(pf);
	  fftw_destroy_plan(pb);
	  fftw_free(in);
	  fftw_free(out);
	  fftw_free(fftout);
	  fftw_free(insig);
	  fftw_free(outsig);
	  fftw_free(halfframe);
	}
	virtual void checksize( int count ) {
		//Hack: check that smallest fft buffer size is >= COUNT
		//      AND n must also be an integer multiple of count...
		if (2 * count > n) {
		  printf("FFT window size must be at least twice the"
			 " control block size.\n");
		  printf("FFT window: %d, control block size: %d.\n",
			 n, count);
		  exit(1);
		}
		else if (n % count) {
		  printf("FFT window size must be an integer multiple of the"
			 " control block size.\n");
		  printf("FFT window: %d, control block size: %d.\n",
			 n, count);
		  exit(1);
		}
	}

	~mydsp() {
	  printf("Calling destructor\n");
	  FFTDestruct();
	}
  virtual void init(int samplingFreq) {
		classInit(samplingFreq);
		instanceInit(samplingFreq);
		FFTInit(); //one fft for now
	}
	virtual void buildUserInterface(UI* interface) {
		interface->openVerticalBox("Frequency Domain Quantizer");
		interface->addHorizontalSlider(
		    "FFT block size (power of two)",
		    &fslider2, 6e+01f, 6.0f, 18.0f, 1.0f
		);
		interface->addHorizontalSlider(
		    "Level (db)",
		    &fslider0, 0.0f, -60.0f, 40.0f, 0.1f
		);
		interface->addCheckButton("On", &fcheckbox0);
		interface->addHorizontalSlider(
		    "Quantization bits",
		    &fslider1, 4.0f, 2.0f, 16.0f, 1.0f
		);
		interface->closeBox();
	}
	virtual void compute (
                        int count,
                        FAUSTFLOAT** input,
                        FAUSTFLOAT** output
                       )
	{
		checksize(count); //check ctrl block & FFT buffer sizes behave

		//linear amplitude factor from dB volume input variable
		float fSlow0 = (
		    0.0010000000000000009f * powf(10,(0.05f * fslider0))
		);

		//mute (really "on") button
		float fSlow1 = fcheckbox0;

		// incoming audio for this ctrl block
		FAUSTFLOAT* input0  = input[0];

		// outgoing audio for this ctrl block
		FAUSTFLOAT* output0 = output[0];


		//---------------
		for (i = 0; i < count; i++) insig[i+fftcountin] = input0[i];

		fftcountin += count;

		if (fftcountin == n) //perform the fft/quantization/ifft)
		{
		      for (i = 0; i < n; i++) in[i] = insig[i];
		      sinewindow(in, n); //window input
		      fftw_execute(pf);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc; i++){
			fftout[i][0] = quantize(bits, fftout[i][0], (int) n/2.0);
			fftout[i][1] = quantize(bits, fftout[i][1], (int) n/2.0);
		      }

		      fftw_execute(pb); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n; i++) out[i] /= n;

		      sinewindow(out, n); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn; i++){
			outsig[fftcountout+i] = fmin(
			    1.0, fmax(-1.0, out[i] + halfframe[i])
			);
			insig[i]     = insig[i + halfn];
			halfframe[i] = out[i + halfn];
		      }

		      fftcountin  -= halfn;
		      fftcountout += halfn;
		}

		if (fftcountout>0){
		      for (i = 0; i < count; ++i){
			  output0[i] = outsig[i];
		      }
		      for (i = 0; i<(fftcountout-count); i++){
			  outsig[i] = outsig[i + count];
		      }
		      fftcountout -= count;
		}
		else {
		      for (i = 0; i < count; ++i){
			  output0[i] = 0.0;
		      }
		}

		for (int i=0; i<count; i++) {
			fRec0[0] = (fSlow0 + (0.999f * fRec0[1]));
			output0[i] = (FAUSTFLOAT)(fSlow1 * ((float)output0[i] * fRec0[0]));
			// post processing
			fRec0[1] = fRec0[0];
		}

	}
};


