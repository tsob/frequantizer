(
o = s.options;
o.memSize = 32768; //increase server's memory size for larger FFTs
s.boot;
s.meter;
)


(
SynthDef("myCoder", { arg out=0, bits=8, fftsize=4096, mul=1.0;
  // Applies frequency quantization on input audio -- CAREFUL OF FEEDBACK!
	var in, outSig, chain;
	//in = SinOsc.ar(SinOsc.kr(SinOsc.kr(0.08, 0, 6, 6.2).squared, 0, 100, 800));
	in = SoundIn.ar(0);
	chain = FFT(LocalBuf(fftsize), in);
	chain = PV_Quantize(chain, bits);
	Out.ar(out, IFFT(chain)*mul);
}).add;
)

// Try this out with a long window length -- left channel.
// Notice the temporal artifacts.
a = Synth("myCoder", [\out, 0, \bits, 8, \fftsize, 32768]);
a.set(\mul, 15.0); //increase gain - but be careful
a.set(\bits, 2);   //very pronounced quantization effects now.
a.set(\bits, 6);
a.set(\bits, 10);
a.free;

//try with a smaller window in the right channel.
b = Synth("myCoder", [\out, 1, \bits, 2, \fftsize, 256]);
b.set(\mul, 12.0); //increase gain - but be careful
b.set(\bits, 10);
b.set(\bits, 6);
b.set(\bits, 4);
b.set(\bits, 2);
b.free;


// Mouse x position is bit depth and y position is volume.
(
SynthDef("myCoder2", { arg out=0, fftsize=4096;
	var in, outSig, chain;
	in = SoundIn.ar(0);
	chain = FFT(LocalBuf(fftsize), in);
	chain = PV_Quantize(chain, MouseX.kr(2,14) );
	Out.ar(out, IFFT(chain)*MouseY.kr(0.0,30.0) );
}).add;
)


c = Synth("myCoder2", [\out, 0, \fftsize, 32768]);
c.set(\mul, 15.0); //increase gain - but be careful
c.set(\fftsize, 4096);
c.free;

(
SynthDef("myCoder3", { arg out=0, size=3;
	var in, outSig, chain, chain1, chain2, chain3, chain4, chain5;
	in = SoundIn.ar(0);

	chain1 = FFT(LocalBuf(2048), in);
	chain2 = FFT(LocalBuf(4096), in);
	chain3 = FFT(LocalBuf(8192), in);
	chain4 = FFT(LocalBuf(16384), in);
	chain5 = FFT(LocalBuf(32768), in);

	chain1 = PV_Quantize(chain1, MouseX.kr(2,14) );
	chain2 = PV_Quantize(chain2, MouseX.kr(2,14) );
	chain3 = PV_Quantize(chain3, MouseX.kr(2,14) );
	chain4 = PV_Quantize(chain4, MouseX.kr(2,14) );
	chain5 = PV_Quantize(chain5, MouseX.kr(2,14) );

	Out.ar(
		out,
		Select.ar(size,
			[
				IFFT(chain1),
				IFFT(chain2),
				IFFT(chain3),
				IFFT(chain4),
				IFFT(chain5),
			]
		).dup
		* MouseY.kr(0.0,30.0)
	);

}).add;
)


d = Synth("myCoder3", [\out, 0, \size, 1]);
d.set(\size,1);
d.set(\size,2);
d.set(\size,3);
d.set(\size,4);
d.set(\size,5);
d.free;


(
SynthDef("myCoder4", { arg out=0, gain=30;
	var in, outSig, chain, chain1, chain2, chain3, chain4, chain5;
	in = SoundIn.ar(0);

	chain1 = FFT(LocalBuf(2048), in);
	chain2 = FFT(LocalBuf(4096), in);
	chain3 = FFT(LocalBuf(8192), in);
	chain4 = FFT(LocalBuf(16384), in);
	chain5 = FFT(LocalBuf(32768), in);

	chain1 = PV_Quantize(chain1, MouseX.kr(2,14) );
	chain2 = PV_Quantize(chain2, MouseX.kr(2,14) );
	chain3 = PV_Quantize(chain3, MouseX.kr(2,14) );
	chain4 = PV_Quantize(chain4, MouseX.kr(2,14) );
	chain5 = PV_Quantize(chain5, MouseX.kr(2,14) );

	Out.ar(
		out,
		SelectX.ar(MouseY.kr(0,5),
			[
				IFFT(chain1),
				IFFT(chain2),
				IFFT(chain3),
				IFFT(chain4),
				IFFT(chain5),
			]
		).dup
		* gain
	);

}).add;
)

// Mouse x position is bit depth (2 to 14), y position is window length
e = Synth("myCoder4", [\out, 0, \gain, 30]);
e.free;





(
SynthDef("myCoder5", { arg out=0, fftsize=4096, gain=15.0;
	var in, outSig, chain, chainA, chainB;
	in = SoundIn.ar(0);
	chainA = FFT(LocalBuf(fftsize), in);
	chainB = PV_Copy(chainA,LocalBuf(fftsize));
	chainB = PV_Quantize(chainA, MouseX.kr(2,14) );
  chain = PV_BinWipe(chainA, chainB, MouseY.kr(-1, 1));
	Out.ar(out, IFFT(chain)*gain.dup );
}).add;
)

f = Synth("myCoder5", [\out, 0, \fftsize, 8192, \gain, 20]);
f.free;


// fft bin stretching - could use some work
(
SynthDef("myCoder6", { arg out=0, fftsize=4096, gain=15.0;
	var in, outSig, chain, chainA, chainB;
	in = SoundIn.ar(0);
	chain = FFT(LocalBuf(fftsize), in);
	chain = PV_BinShift(chain, MouseY.kr(0.2, 2, \exponential) );
	chain = PV_Quantize(chain, MouseX.kr(2,14), 0, 1);
	Out.ar(out, IFFT(chain)*gain.dup );
}).add;
)

g = Synth("myCoder6", [\out, 0, \fftsize, 16384, \gain, 20]);
g.free;


(
SynthDef("myCoder7", { arg out=0, fftsize=4096, gain=15.0;
	var in, outSig, chain, chainA, chainB;
	in = SoundIn.ar(0);
	chain = FFT(LocalBuf(fftsize), in);
  chain = PV_BinShift(chain, 1, MouseY.kr(-128, 128), 1 );
	chain = PV_Quantize(chain, MouseX.kr(2,14) );
	Out.ar(out, IFFT(chain)*gain.dup );
}).add;
)

h = Synth("myCoder7", [\out, 0, \fftsize, 1024, \gain, 20]);
h.free;

(
SynthDef("myCoder8", { arg out=0, fftsize=4096, gain=15.0;
	var in, outSig, chain, chainA, chainB;
	in = SoundIn.ar(0);
	chain = FFT(LocalBuf(fftsize), in);
  chain = PV_BinShift(chain, MouseY.kr(0.2, 1.5, \exponential), MouseY.kr(-128, 128), 1 );
	chain = PV_Quantize(chain, MouseX.kr(2,14) );
	Out.ar(out, IFFT(chain)*gain.dup );
}).add;
)

i = Synth("myCoder8", [\out, 0, \fftsize, 1024, \gain, 20]);
i.free;
