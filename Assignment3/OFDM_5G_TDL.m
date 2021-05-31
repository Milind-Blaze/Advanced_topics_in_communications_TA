%{
This script is used to simulate an OFDM system with a TDL channel
according to 5G specifications. 

Data is transmitted through a TDL channel. Pilot insertion and channel
estimation are performed. ZF equalization is also performed. 

Author: Milind Kumar V
%}

clc; 
close all;
clear;

%% Regarding fd

% In contrast to the requirements of the assignment where the variable fd
% is to be provided at the top of the script, fd is instead computed in the
% "TDL channel parameters" section using the relative velocity between the
% transmitter and receiver.

%% Defining the 5G NR parameters

CPLen1 = 352;                   % length of the cyclic prefix for the first OFDM symbol
CPLen2 = 288;                   % length of the cyclic prefix for the OFDM symbols 2-14
FFTsize = 4096;                 % size of the FFT being used
SCS = 30e3;                     % Subcarrier spacing
SamplingRate = SCS*FFTsize;     % symbol/sample rate for time domain data
numSlots = 1;                   % number of slots
numSymbolSlot = 14;             % number of symbols per slot
numREperRB = 12;                % number of resource elements per resource block
QAMorderList = [4,16,64,256];   % the value of M in M-QAM, modulation used for data
QAMencoding = 'gray';           % QAM encoding type
SNR_list = -20:2:25;            % Range of SNRs simulated                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          25;             % SNR for the AWGN channel

%% TDL channel parameters

v = 0;                              % Relative speed of UE wrt basestation
c = physconst('lightspeed');        % Speed of light
fc = 3.5e9;                         % Carrier frequency
fd = fc*v/c;                        % Doppler shift
delaySpread = 30e-9;                % Delay spread of channel in s
MIMOCorrelation = 'Low';            % Correlation between UE and BS antennas
delayProfile = 'TDL-C';             % Different channel models from the spec TR 38.901, options are TDL-A/B/C/D/E
numReceiveAntennas = 1;             % number of receive antennas
numTransmitAntennas = 1;            % number of transmit antennas
numIter = 10;                       % number of times the TDL channel is initialized

%% Non-ideal effects

timing_offset = 0;
CFO = 0;

%% Defining Resource Grid (RG) parameters

numRB = 50;                             % number of resource blocks allocated 
numRE = numRB*numREperRB;               % number of REs available for transmission 
numSymbols = numSlots*numSymbolSlot;    % Total number of OFDM symbols in the slot

%% Output variables

BER_vs_SNR = zeros(length(SNR_list),1);
BER_vs_QAM = zeros(length(SNR_list), length(QAMorderList));
EbN0_vs_QAM = zeros(length(SNR_list), length(QAMorderList));
BER_vs_QAM_theory = zeros(length(SNR_list), length(QAMorderList));


%% Creating the TDL channel 

% % Creating nrCarrierConfig object
% carrierObject = nrCarrierConfig;
% carrierObject.SubcarrierSpacing = SCS/1e3;      % SCS defined in kHz
% carrierObject.CyclicPrefix = 'normal';          % 14 OFDM symbols in a slot
% carrierObject.NSizeGrid = numRB;                % Number of RBs in RG


tdl = nrTDLChannel;
tdl.DelayProfile = delayProfile;
tdl.DelaySpread = delaySpread;
tdl.MaximumDopplerShift = fd;
tdl.SampleRate = SamplingRate;
tdl.MIMOCorrelation = MIMOCorrelation;
tdl.NumTransmitAntennas = numTransmitAntennas;
tdl.NumReceiveAntennas = numReceiveAntennas;

tdlinfo = info(tdl);
pathFilters = getPathFilters(tdl);


%% Simulation

% Iterating over different Modulation orders

for M_index = 1:length(QAMorderList)
    QAMorder = QAMorderList(M_index);

    %% Generating the resource grid 

    RG = zeros(numRE, numSymbolSlot*numSlots);

    %% Generating transmit signal
    
    %% Generating pilots to transmit
    pilotGrid = zeros(numRE, numSymbols);
    pilotGrid(1:2:numRE,1:2:numSymbols) = (1/(2))*(1 + 1i);
    pilotLoc = find(pilotGrid ~= 0);
    dataLoc = find(pilotGrid == 0);
    
    % Plot the pilot grid
%     plotResourceGrid(abs(pilotGrid), "Pilot grid", "OFDM symbol", "REs");
    
    %% Generating transmit data
    
    % Computing the number of bits necessary
    numBits = log2(QAMorder)*length(dataLoc);

    % Generating the bits
    dataBits = randi([0,1], numBits, 1);

    % Generating QAM symbols 
    modulatedSymbols = qammod(dataBits, QAMorder, QAMencoding, "InputType",...
        "bit", "UnitAveragePower", true);

    % Populating the reosource grid
    RG(dataLoc) = modulatedSymbols;
    RG = RG + pilotGrid;

    % Plotting the generated resource grid
%     plotResourceGrid(abs(RG), "RG", "OFDM symbols", "REs");

    %% Generating the FFT grid

    FFTgrid = zeros(FFTsize, numSymbolSlot*numSlots);

    % Generating the suitable index mapping. This is necessary because of the 
    % IFFT implementation in MATLAB. One needs to rearrange [0,a,b,c,d,e,f,0] 
    % to [d,e,f,0,0,a,b,c]. This occurs only along the first (RE) dimension.
    orgIndexSet = 0:numRE-1;
    newIndexSet = mod(orgIndexSet - numRE/2, FFTsize) + 1;
    FFTgrid(newIndexSet,:) = RG;

    % Plotting the generated FFT grid
%     plotResourceGrid(abs(FFTgrid), "FFT grid", "OFDM symbols", "REs");

    %% Generating the time domain data

    % IFFT
    normalizingFactor = sqrt(FFTsize);
    timeDataParallel = normalizingFactor*ifft(FFTgrid, FFTsize, 1);

    % Add CP and time domain sequence
    timeDataSerial = [];

    % Add CP according to the OFDM symbol number
    for i = 1:numSymbols
        if mod(i-1,14) == 0
            CPlength = CPLen1;
        else
            CPlength = CPLen2;
        end
        timeDataSymbol = timeDataParallel(:,i);
        timeDataSerial = [timeDataSerial; timeDataSymbol(end-CPlength + 1:end,:);...
            timeDataSymbol];
    end

    % Change input shape to feed to MATLAB TDL channel
    timeDataSerial = squeeze(timeDataSerial);

    %% OFDM spectrum
    
    % Plot the spectrum of the OFDM signal
    % figure;
    % pspectrum(timeDataSerial, SamplingRate);
    
    
    %% Passing through a TDL channel
    BER_vs_iter = [];         % Re-initialized for every M-index
    for iter = 1:numIter
        
        [channelOutputTDL, pathGains, sampleTimes] = 	tdl(timeDataSerial);

        % Estimating the channel
    %     pathFilters = getPathFilters(tdl);
    %     perfectChannel = nrPerfectChannelEstimate(pathGains, pathFilters,...
    %         numRB, SCS/1e3, numSlots, 'SampleRate', SamplingRate);


        % Iterating over SNR values
        BER_vs_SNR = zeros(length(SNR_list),1); % Re-initializing for each QAMorder
        for SNR_index = 1:length(SNR_list)
            %% Adding Gaussian noise
            SNR = SNR_list(SNR_index);
            channelOutput = awgn(channelOutputTDL, SNR, "measured", 1234);

            %% Receiver

            % Introducing timing offset
            timing_offset_signal = zeros(timing_offset, 1);
            channelOutput = [timing_offset_signal;...
                channelOutput(1:end-timing_offset)];



            % Introducing carrier frequency offset
            residual_freq = exp(2*pi*1i*CFO*(1/FFTsize).*...
                (0:(length(channelOutput)-1)));
            RXinput = channelOutput.*residual_freq';

            %% CP strip and serial to parallel

            CPcumulative = 0;
            RXtimeDataParallel = zeros(FFTsize, numSymbols);
            for i = 1:numSymbols
                if mod(i-1, 14) == 0
                    CPlength = CPLen1;
                else
                    CPlength = CPLen2;
                end
                CPcumulative = CPcumulative + CPlength;
                RXtimeDataSymbol = RXinput(CPcumulative + (i-1)*FFTsize + 1: CPcumulative +...
                    i*FFTsize,:);
                RXtimeDataParallel(:,i) = RXtimeDataSymbol;
            end

            %% FFT
            % Scaling performed again because IFFT was also scaled
            RXFFTgrid = 1/normalizingFactor*fft(RXtimeDataParallel, FFTsize, 1);

            %% Correction for timing offset

            indexArray = 1:FFTsize;
            indexArray = indexArray' - 1;
            TOcorrectionVector = exp(2*pi*1i*indexArray*timing_offset/FFTsize);
            RXFFTgrid = RXFFTgrid.*TOcorrectionVector;

            %% Mapping FFT grid to RB grid
            % This is necessary to undo the mapping done before hand and to rearrange 
            % [d,e,f,0,0,a,b,c] to [0,a,b,c,d,e,f,0]. This occurs only along the 1st
            % (RE) dimension.

            RXfreqData = RXFFTgrid(newIndexSet, :);
            % plotResourceGrid(abs(RXfreqData), "RX RG", "OFDM symbols", "REs");
            %% Determining the channel 

            channelEstimate = ones(numRE, numSymbols);
            channelEstimate(pilotLoc) = RXfreqData(pilotLoc)./RG(pilotLoc);
            channelEstimate(2:2:numRE, 1:2:numSymbols) = channelEstimate(1:2:numRE,1:2:numSymbols);
            channelEstimate(:,2:2:numSymbols) = channelEstimate(:,1:2:numSymbols);

            % Plotting the estimated channel
    %         plotResourceGrid(abs(channelEstimate), "Estimated channel",...
    %             "OFDM symbol", "RE");
            %% Equalization
            RXfreqData = RXfreqData./channelEstimate;

            %% Obtaining the transmitted bits

            RXfreqSerial = reshape(RXfreqData(dataLoc),[],1);                     
            RXbits = qamdemod(RXfreqSerial, QAMorder, QAMencoding, "OutputType",...
                "bit", "UnitAveragePower", true);

            %% BER

            BER = sum(RXbits ~= dataBits)/numBits;
            BER_vs_SNR(SNR_index) = BER;
        end
        BER_vs_iter = [BER_vs_iter, BER_vs_SNR];
    end
    BER_vs_QAM(:,M_index) = mean(BER_vs_iter,2);
    EbN0_vs_QAM(:,M_index) = SNR_list ... 
        + 10*log10(FFTsize/numRE) - 10*log10(log2(QAMorder));
    BER_vs_QAM_theory(:,M_index) =  berawgn(EbN0_vs_QAM(:,M_index),...
        'qam', QAMorder);
end

legendList = {};
fig = figure;

plot_count = 1;
for M_index = 1:length(QAMorderList)
    semilogy(EbN0_vs_QAM(:,M_index),...
        BER_vs_QAM(:,M_index),...
        "o--", 'LineWidth', 2);
    legendValue = string(QAMorderList(M_index)) + " QAM empirical";
    legendList{1,plot_count} = legendValue;
    plot_count = plot_count + 1;
    hold on;
end

for M_index = 1:length(QAMorderList)
    semilogy(EbN0_vs_QAM(:,M_index),...
        BER_vs_QAM_theory(:,M_index),...
        'LineWidth', 2);
    legendValue = string(QAMorderList(M_index)) + " QAM theoretical AWGN";
    legendList{1,plot_count} = legendValue;
    plot_count = plot_count + 1;
end

grid on
ylim([1e-3,1]);
xlabel("$10\log\left(\frac{E_b}{N_0}\right)$", "Interpreter", "latex", 'FontSize', 22);
ylabel("BER", 'FontSize', 22);
titleValue = "BER curves for a TDL channel, fd: " + string(fd) + " Hz, "...
    + "iterations: " + string(numIter);

title(titleValue,'FontSize', 22);
legend(legendList);

    