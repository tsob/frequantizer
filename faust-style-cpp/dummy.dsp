declare name "Frequency Domain Quantizer";
declare version "0.1";
declare author "Tim O'Brien";

import("music.lib");

onswitch = checkbox("On");
bits = hslider(
    "Quantization bits",
    4, 2, 16, 1.0
    ); //bits for quantization
fft_power = hslider(
    "FFT block size (power of two)",
    60, 6, 18, 1.0
    ); //power of two
smooth(c) = *(1-c) : +~*(c);
level  = hslider("Level (db)",
    0, -96, 4, 0.1
    ) : db2linear : smooth(0.999);

process = _ : *(level) : *(onswitch) : *(fft_power) : *(bits);
