# Assignment 6



In this assignment, a bits to bits OFDM system is simulated and the performance of LDPC coding is studied in the case of an AWGN channel for 16 and 64 QAM constellations. The LDPC encoding and decoding, rate matching and recovery are all implemented using functions provided in the MATLAB 5G NR toolbox. The objective is to obtain the BER vs Eb/N0 curves for different constellations at different code rates. The following specifications are used for this assignment-

- TS 38.212 
  - Section 7.2



- [Problem statement](https://drive.google.com/file/d/1xsRpaC-CvVn1NWNDqieZTxPUsDYt3H05/view?usp=sharing)  
- [Report outline](https://drive.google.com/file/d/1tUAnCLgcfeZDdLT6tbHtJJ7GxPTbk7DQ/view?usp=sharing)



## Contents 

- **```OFDM_AWGN_LDPC.m```**: a script to complete the simulation given in the problem statement above. It simulates a bits to bits communications system where random bits are generated, LDPC encoded, rate matched, modulated, transmitted through an AWGN channel, demodulated into soft bits, rate recovered, decoded and converted into bits which can be used to compute the BER. The output of this script is
  - a plot of BER vs Eb/N0 for coded and uncoded transmission over an AWGN channel for each constellation specified by the variable ```QAMorderList```. The number of bits to be transmitted and the code rate for modulation order x (expressed as a power of 2) are given by the variables ```A``` and ```R{x}```. 
- **```plotResourceGrid.m```**: function to plot a heatmap of a specified 2 dimensional array- this is useful when analysing the resource grid and resource allocation



## Useful resources

- [LDPC processing for DL-SCH and UL-SCH in MATLAB](https://in.mathworks.com/help/5g/gs/ldpc-processing-chain-for-dl-sch.html) 
