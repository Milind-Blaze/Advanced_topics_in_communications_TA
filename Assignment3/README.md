# Assignment 3



This assignment requires the simulation of a bits to bits OFDM system in the presence of AWGN and multipath fading channels. Further, the code provided in this folder can be used to analyze the performance of the system (as measured using BER) under different conditions- channels, timing offset and carrier frequency offset. Timing offset correction is also implemented. In the case of the frequency selective channel, channel estimation and equalization are performed. 



-  [Problem statement](https://drive.google.com/file/d/13ozoxXyeC3dL7V2y-aMDWoWdYGJ0KhNV/view?usp=sharing)
- [Report outline](https://drive.google.com/file/d/1XIJYuTLABmYUy9K8noZjpsXUcwShArEF/view?usp=sharing)



## Contents

- **```OFDM_5G_AWGN.m```**: script to implement Simulations 1,2 and 3. It also produces the OFDM spectrum as described in the Problem Statement. Here, random bits are generated, modulated, transmitted through an AWGN channel, demodulated and the BER is computed. Timing offset is implemented by adding zeros to the time domain signal at the receiver. CFO is implemented by multiplying with a complex exponential in the time domain. The output includes
  - a plot of BER vs Eb/N0 for the specified QAM values, specified CFO values, specified timing offset and timing offset correction (if enabled). These values can be adjusted using the ```QAMorderList```,``` CFO```, ```timing_offset``` and ```TOcorrection``` variables respectively.
  - a scatter plot of the received 16 QAM constellation (if 16 QAM is indeed transmitted) at the highest SNR value being used for the first CFO value
  - a plot of the OFDM spectrum
- **```OFDM_5G_TDL.m```**: script ot implement Simulation 4 and the second bonus question. Here, random bits are generated, modulated and transmitted through a multipath fading channel, channel estimation and equalization are performed, and the bits obtained after demodulation are used to compute the BER at different Doppler shift frequencies. The output is a 
  - plot of BER vs Eb/N0 for different QAM values and a specified Doppler shift frequency which can be adjusted using ```QAMorderList``` and ```v``` respectively. ```v``` is the relative velocity between the UE and base station.
- **```plotResourceGrid.m```**: function to plot a heatmap of a specified 2 dimensional array- this is useful when analysing the resource grid and resource allocation