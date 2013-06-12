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
	// ------------
	//output fft arrays
	fftw_complex *fftout0, *fftout1, *fftout2, *fftout3, *fftout4, *fftout5;
	//input & output sample arrays - real-valued signals
	double *in0,  *in1,  *in2,  *in3,  *in4,  *in5,
	       *out0, *out1, *out2, *out3, *out4, *out5;
	// plans for fftw
	fftw_plan pf0, pf1, pf2, pf3, pf4, pf5,
		  pb0, pb1, pb2, pb3, pb4, pb5;
	int numwindows = 6; //number of different FFT window sizes to use
	int smalln = 2048;  //smalles FFT window size. Goes up by powers of 2.
	int *n; //array for fft window sizes
	int *nc; //array for complex fft window sizes
	int *halfn; //for 50% overlap
	double *insig0,  *insig1,  *insig2,  *insig3,  *insig4,  *insig5,
	       *outsig0, *outsig1, *outsig2, *outsig3, *outsig4, *outsig5,
	       *halfframe0, *halfframe1, *halfframe2, *halfframe3, *halfframe4,
	       *halfframe5;
	double *fftamount;
	int *fftcountin, *fftcountout;
	int i, j;
	int bits;

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

	virtual void FFTInit( void ) {
		//FFT window sizes
		n = (int*) malloc( sizeof( int ) * numwindows);
		for (i = 0; i < numwindows; i++) n[i] = ( smalln << i );
		//buffer sizes for complex fft bins (using fftw r2c)
		nc = (int*) malloc( sizeof( int ) * numwindows);
		for (i = 0; i < numwindows; i++) nc[i] = (n[i]/2) + 1;
		//for 50% overlap of windows
		halfn = (int*) malloc( sizeof( int ) * numwindows);
		for (i = 0; i < numwindows; i++) halfn[i] = n[i]/2;

		//data arrays for use with fftw
		in0 = (double*) fftw_malloc( sizeof( double ) * n[0]);
		in1 = (double*) fftw_malloc( sizeof( double ) * n[1]);
		in2 = (double*) fftw_malloc( sizeof( double ) * n[2]);
		in3 = (double*) fftw_malloc( sizeof( double ) * n[3]);
		in4 = (double*) fftw_malloc( sizeof( double ) * n[4]);
		in5 = (double*) fftw_malloc( sizeof( double ) * n[5]);
		out0 = (double*) fftw_malloc( sizeof( double ) * n[0]);
		out1 = (double*) fftw_malloc( sizeof( double ) * n[1]);
		out2 = (double*) fftw_malloc( sizeof( double ) * n[2]);
		out3 = (double*) fftw_malloc( sizeof( double ) * n[3]);
		out4 = (double*) fftw_malloc( sizeof( double ) * n[4]);
		out5 = (double*) fftw_malloc( sizeof( double ) * n[5]);
		fftout0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[0]);
		fftout1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[1]);
		fftout2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[2]);
		fftout3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[3]);
		fftout4 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[4]);
		fftout5 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc[5]);
		//forward & backward transform plans
		pf0 = fftw_plan_dft_r2c_1d(n[0], in0, fftout0, FFTW_MEASURE);
		pf1 = fftw_plan_dft_r2c_1d(n[1], in1, fftout1, FFTW_MEASURE);
		pf2 = fftw_plan_dft_r2c_1d(n[2], in2, fftout2, FFTW_MEASURE);
		pf3 = fftw_plan_dft_r2c_1d(n[3], in3, fftout3, FFTW_MEASURE);
		pf4 = fftw_plan_dft_r2c_1d(n[4], in4, fftout4, FFTW_MEASURE);
		pf5 = fftw_plan_dft_r2c_1d(n[5], in5, fftout5, FFTW_MEASURE);
		pb0 = fftw_plan_dft_c2r_1d(n[0], fftout0, out0, FFTW_MEASURE);
		pb1 = fftw_plan_dft_c2r_1d(n[1], fftout1, out1, FFTW_MEASURE);
		pb2 = fftw_plan_dft_c2r_1d(n[2], fftout2, out2, FFTW_MEASURE);
		pb3 = fftw_plan_dft_c2r_1d(n[3], fftout3, out3, FFTW_MEASURE);
		pb4 = fftw_plan_dft_c2r_1d(n[4], fftout4, out4, FFTW_MEASURE);
		pb5 = fftw_plan_dft_c2r_1d(n[5], fftout5, out5, FFTW_MEASURE);
		// for keeping track of incomin/outgoing samples...
		fftcountin = (int*) malloc( sizeof( int ) * numwindows);
		fftcountout = (int*) malloc( sizeof( int ) * numwindows);
		// contribution of each fft to output signal...
		fftamount = (double*) malloc( sizeof( double ) * numwindows);
		for (i = 0; i < numwindows; i++){
		  fftcountin[i] = 0.0;
		  fftcountout[i] = 0.0;
		  fftamount[i] = 0.0;
		}
		//overlap-add frame...
		halfframe0 = (double*) fftw_malloc(sizeof(double) * halfn[0]);
		halfframe1 = (double*) fftw_malloc(sizeof(double) * halfn[1]);
		halfframe2 = (double*) fftw_malloc(sizeof(double) * halfn[2]);
		halfframe3 = (double*) fftw_malloc(sizeof(double) * halfn[3]);
		halfframe4 = (double*) fftw_malloc(sizeof(double) * halfn[4]);
		halfframe5 = (double*) fftw_malloc(sizeof(double) * halfn[5]);
		for (i = 0; i < halfn[0]; ++i) halfframe0[i] = 0.0;
		for (i = 0; i < halfn[1]; ++i) halfframe1[i] = 0.0;
		for (i = 0; i < halfn[2]; ++i) halfframe2[i] = 0.0;
		for (i = 0; i < halfn[3]; ++i) halfframe3[i] = 0.0;
		for (i = 0; i < halfn[4]; ++i) halfframe4[i] = 0.0;
		for (i = 0; i < halfn[5]; ++i) halfframe5[i] = 0.0;
		//buffers for queueing input and output
		insig0 = (double*) fftw_malloc( sizeof(double) * n[0]);
		insig1 = (double*) fftw_malloc( sizeof(double) * n[1]);
		insig2 = (double*) fftw_malloc( sizeof(double) * n[2]);
		insig3 = (double*) fftw_malloc( sizeof(double) * n[3]);
		insig4 = (double*) fftw_malloc( sizeof(double) * n[4]);
		insig5 = (double*) fftw_malloc( sizeof(double) * n[5]);
		outsig0 = (double*) fftw_malloc( sizeof(double) * n[0]);
		outsig1 = (double*) fftw_malloc( sizeof(double) * n[1]);
		outsig2 = (double*) fftw_malloc( sizeof(double) * n[2]);
		outsig3 = (double*) fftw_malloc( sizeof(double) * n[3]);
		outsig4 = (double*) fftw_malloc( sizeof(double) * n[4]);
		outsig5 = (double*) fftw_malloc( sizeof(double) * n[5]);
		for (i = 0; i < n[0]; ++i) {
		  insig0[i]  = 0.0;
		  outsig0[i] = 0.0;
		}
		for (i = 0; i < n[1]; ++i) {
		  insig1[i]  = 0.0;
		  outsig1[i] = 0.0;
		}
		for (i = 0; i < n[2]; ++i) {
		  insig2[i]  = 0.0;
		  outsig2[i] = 0.0;
		}
		for (i = 0; i < n[3]; ++i) {
		  insig3[i]  = 0.0;
		  outsig3[i] = 0.0;
		}
		for (i = 0; i < n[4]; ++i) {
		  insig4[i]  = 0.0;
		  outsig4[i] = 0.0;
		}
		for (i = 0; i < n[5]; ++i) {
		  insig5[i]  = 0.0;
		  outsig5[i] = 0.0;
		}
	}
	// FFT destructor function
	virtual void FFTDestruct( void ) {
	  printf("Cleaning up...\n");
	  fftw_destroy_plan(pf0);
	  fftw_destroy_plan(pb0);
	  fftw_free(in0);
	  fftw_free(in1);
	  fftw_free(in2);
	  fftw_free(in3);
	  fftw_free(in4);
	  fftw_free(in5);
	  fftw_free(out0);
	  fftw_free(out1);
	  fftw_free(out2);
	  fftw_free(out3);
	  fftw_free(out4);
	  fftw_free(out5);
	  fftw_free(fftout0);
	  fftw_free(fftout1);
	  fftw_free(fftout2);
	  fftw_free(fftout3);
	  fftw_free(fftout4);
	  fftw_free(fftout5);
	  fftw_free(insig0);
	  fftw_free(insig1);
	  fftw_free(insig2);
	  fftw_free(insig3);
	  fftw_free(insig4);
	  fftw_free(insig5);
	  fftw_free(outsig0);
	  fftw_free(outsig1);
	  fftw_free(outsig2);
	  fftw_free(outsig3);
	  fftw_free(outsig4);
	  fftw_free(outsig5);
	  fftw_free(halfframe0);
	  fftw_free(halfframe1);
	  fftw_free(halfframe2);
	  fftw_free(halfframe3);
	  fftw_free(halfframe4);
	  fftw_free(halfframe5);
	  free(n);
	  free(nc);
	  free(halfn);
	  free(fftcountin);
	  free(fftcountout);
	  free(fftamount);
	}
	virtual void checksize( int count ) {
	    for (i = 0; i < numwindows; i++)
	    {
		//Hack: check that smallest fft buffer size is >= COUNT
		//      AND n vals must also be an integer multiple of count...
		if (2 * count > n[i]) {
		  printf("FFT window size must be at least twice the"
			 " control block size.\n");
		  printf("FFT window: %d, control block size: %d.\n",
			 n[i], count);
		  exit(1);
		}
		else if (n[i] % count) {
		  printf("FFT window size must be an integer multiple of the"
			 " control block size.\n");
		  printf("FFT window: %d, control block size: %d.\n",
			 n[i], count);
		  exit(1);
		}
	    }
	}

	~mydsp() {
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
		    &fslider2, 13.0f, 11.0f, 18.0f, 0.1f
		);
		interface->addHorizontalSlider(
		    "Level (db)",
		    &fslider0, 0.0f, -60.0f, 90.0f, 0.1f
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


		//Queue up this control block of input samples for the various
		//FFTs...
		for (i = 0; i < count; i++){
		      insig0[i+fftcountin[0]] = input0[i];
		      insig1[i+fftcountin[1]] = input0[i];
		      insig2[i+fftcountin[2]] = input0[i];
		      insig3[i+fftcountin[3]] = input0[i];
		      insig4[i+fftcountin[4]] = input0[i];
		      insig5[i+fftcountin[5]] = input0[i];
		}
		for (i = 0; i < numwindows; i++) fftcountin[i] += count;

		//check/perform the fft/quantization/iffts for the various
		//fft sizes...


		if (fftcountin[0] == n[0])
		{
		      for (i = 0; i < n[0]; i++) in0[i] = insig0[i];
		      sinewindow(in0, n[0]); //window input
		      fftw_execute(pf0);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[0]; i++){
			fftout0[i][0] = quantize(bits, fftout0[i][0], (int) n[0]/2.0);
			fftout0[i][1] = quantize(bits, fftout0[i][1], (int) n[0]/2.0);
		      }

		      fftw_execute(pb0); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[0]; i++) out0[i] /= n[0];

		      sinewindow(out0, n[0]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[0]; i++){
			outsig0[fftcountout[0]+i] = fmin(
			    1.0, fmax(-1.0, out0[i] + halfframe0[i])
			);
			insig0[i]     = insig0[i + halfn[0]];
			halfframe0[i] = out0[i + halfn[0]];
		      }

		      fftcountin[0]  -= halfn[0];
		      fftcountout[0] += halfn[0];
		}

		if (fftcountin[1] == n[1])
		{
		      for (i = 0; i < n[1]; i++) in1[i] = insig1[i];
		      sinewindow(in1, n[1]); //window input
		      fftw_execute(pf1);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[1]; i++){
			fftout1[i][0] = quantize(bits, fftout1[i][0], (int) n[1]/2.0);
			fftout1[i][1] = quantize(bits, fftout1[i][1], (int) n[1]/2.0);
		      }

		      fftw_execute(pb1); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[1]; i++) out1[i] /= n[1];

		      sinewindow(out1, n[1]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[1]; i++){
			outsig1[fftcountout[1]+i] = fmin(
			    1.0, fmax(-1.0, out1[i] + halfframe1[i])
			);
			insig1[i]     = insig1[i + halfn[1]];
			halfframe1[i] = out1[i + halfn[1]];
		      }

		      fftcountin[1]  -= halfn[1];
		      fftcountout[1] += halfn[1];
		}

		if (fftcountin[2] == n[2]) //perform the fft/quantization/ifft
		{
		      for (i = 0; i < n[2]; i++) in2[i] = insig2[i];
		      sinewindow(in2, n[2]); //window input
		      fftw_execute(pf2);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[2]; i++){
			fftout2[i][0] = quantize(bits, fftout2[i][0], (int) n[2]/2.0);
			fftout2[i][1] = quantize(bits, fftout2[i][1], (int) n[2]/2.0);
		      }

		      fftw_execute(pb2); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[2]; i++) out2[i] /= n[2];

		      sinewindow(out2, n[2]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[2]; i++){
			outsig2[fftcountout[2]+i] = fmin(
			    1.0, fmax(-1.0, out2[i] + halfframe2[i])
			);
			insig2[i]     = insig2[i + halfn[2]];
			halfframe2[i] = out2[i + halfn[2]];
		      }

		      fftcountin[2]  -= halfn[2];
		      fftcountout[2] += halfn[2];
		}

		if (fftcountin[3] == n[3]) //perform the fft/quantization/ifft
		{
		      for (i = 0; i < n[3]; i++) in3[i] = insig3[i];
		      sinewindow(in3, n[3]); //window input
		      fftw_execute(pf3);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[3]; i++){
			fftout3[i][0] = quantize(bits, fftout3[i][0], (int) n[3]/2.0);
			fftout3[i][1] = quantize(bits, fftout3[i][1], (int) n[3]/2.0);
		      }

		      fftw_execute(pb3); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[3]; i++) out3[i] /= n[3];

		      sinewindow(out3, n[3]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[3]; i++){
			outsig3[fftcountout[3]+i] = fmin(
			    1.0, fmax(-1.0, out3[i] + halfframe3[i])
			);
			insig3[i]     = insig3[i + halfn[3]];
			halfframe3[i] = out3[i + halfn[3]];
		      }

		      fftcountin[3]  -= halfn[3];
		      fftcountout[3] += halfn[3];
		}

		if (fftcountin[4] == n[4]) //perform the fft/quantization/ifft
		{
		      for (i = 0; i < n[4]; i++) in4[i] = insig4[i];
		      sinewindow(in4, n[4]); //window input
		      fftw_execute(pf4);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[4]; i++){
			fftout4[i][0] = quantize(bits, fftout4[i][0], (int) n[4]/2.0);
			fftout4[i][1] = quantize(bits, fftout4[i][1], (int) n[4]/2.0);
		      }

		      fftw_execute(pb4); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[4]; i++) out4[i] /= n[4];

		      sinewindow(out4, n[4]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[4]; i++){
			outsig4[fftcountout[4]+i] = fmin(
			    1.0, fmax(-1.0, out4[i] + halfframe4[i])
			);
			insig4[i]     = insig4[i + halfn[4]];
			halfframe4[i] = out4[i + halfn[4]];
		      }

		      fftcountin[4]  -= halfn[4];
		      fftcountout[4] += halfn[4];
		}

		if (fftcountin[5] == n[5]) //perform the fft/quantization/ifft
		{
		      for (i = 0; i < n[5]; i++) in5[i] = insig5[i];
		      sinewindow(in5, n[5]); //window input
		      fftw_execute(pf5);  //forward transform

		      // quantize fft values
		      bits = (int) fslider1; //bit depth from GUI
		      for(i=0; i<nc[5]; i++){
			fftout5[i][0] = quantize(bits, fftout5[i][0], (int) n[5]/2.0);
			fftout5[i][1] = quantize(bits, fftout5[i][1], (int) n[5]/2.0);
		      }

		      fftw_execute(pb5); // backward transform

		      // scale by n after inverse tranform
		      for(i=0; i<n[5]; i++) out5[i] /= n[5];

		      sinewindow(out5, n[5]); //window output

		      //overlap-and-add
		      for (i = 0; i < halfn[5]; i++){
			outsig5[fftcountout[5]+i] = fmin(
			    1.0, fmax(-1.0, out5[i] + halfframe5[i])
			);
			insig5[i]     = insig5[i + halfn[5]];
			halfframe5[i] = out5[i + halfn[5]];
		      }

		      fftcountin[5]  -= halfn[5];
		      fftcountout[5] += halfn[5];
		}


		//--------------
		// combine into count-segmented audio streams for each fft
		for (i = 0; i < count; ++i) output0[i] = (FAUSTFLOAT) 0.0;
		for (i = 0; i < numwindows; i++) //for crossfading
		{
		    fftamount[i] = max(1.0 - fabs(11.0+i-fslider2),0.0);
		}
		if (fftcountout[0]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig0[i] * fftamount[0]);
		      }
		      for (i = 0; i<(fftcountout[0]-count); i++){
			  outsig0[i] = outsig0[i + count];
		      }
		      fftcountout[0] -= count;
		}
		if (fftcountout[1]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig1[i] * fftamount[1]);
		      }
		      for (i = 0; i<(fftcountout[1]-count); i++){
			  outsig1[i] = outsig1[i + count];
		      }
		      fftcountout[1] -= count;
		}
		if (fftcountout[2]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig2[i] * fftamount[2]);
		      }
		      for (i = 0; i<(fftcountout[2]-count); i++){
			  outsig2[i] = outsig2[i + count];
		      }
		      fftcountout[2] -= count;
		}
		if (fftcountout[3]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig3[i] * fftamount[3]);
		      }
		      for (i = 0; i<(fftcountout[3]-count); i++){
			  outsig3[i] = outsig3[i + count];
		      }
		      fftcountout[3] -= count;
		}
		if (fftcountout[4]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig4[i] * fftamount[4]);
		      }
		      for (i = 0; i<(fftcountout[4]-count); i++){
			  outsig4[i] = outsig4[i + count];
		      }
		      fftcountout[4] -= count;
		}
		if (fftcountout[5]>0){
		      for (i = 0; i < count; ++i){
			  output0[i] += (FAUSTFLOAT)(outsig5[i] * fftamount[5]);
		      }
		      for (i = 0; i<(fftcountout[5]-count); i++){
			  outsig5[i] = outsig5[i + count];
		      }
		      fftcountout[5] -= count;
		}

		//-------------------------------------------------------------
		//Volume and mute
		for (int i=0; i<count; i++) {
			output0[i] = (FAUSTFLOAT)(fSlow1 * ((float)output0[i] * fSlow0));
		}
	}
};


