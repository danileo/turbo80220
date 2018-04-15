# Turbo code in IEEE 802.20
Matlab (with MEX) implementation of the source coding scheme in IEEE 802.20 Mobile Broadband Wireless Access.

## Index of MEX-files:
turbo_code.cpp: turbo modulator
turbo.cpp: turbo decoder
modulator.cpp: interleaving + QAM modulation
demodulator.cpp: deinterleaving + QAM demodulation

## Instructions
Instructions on how to use the MEX-functions are at the beginning
of the respective files.
Parameters of turbo code are at the beginning of utility.cpp
In main_turbo.m and main_bicm.m the MEX-functions are used to
simulate turbo/bicm encoding/decoding with an AWGN channel.
