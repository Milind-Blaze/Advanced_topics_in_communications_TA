%{
This script is used to simulate an OFDM system with a AWGN channel
according to 5G specifications. Further, timing offset and CFO are
implemented. Timing offset correction is also implemented.

When run, the script produces plots of BER vs Eb/N0 for a given timing
offset and CFO.

Author: Milind Kumar V
%}

clc; 
close all;
clear;

%% Non-ideal effects

timing_offset = 0;                          % Number of samples by which the signal is offset
TOcorrection = false;                       % If "true" implement TO correction, set to false otherwise
% CFO = [0,0.0005,0.001,0.0015,0.002,0.003];  % CFO = (offset frequency)/(SCS)
CFO = 0;

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

%% Defining Resource Grid (RG) parameters

numRB = 50;                             % number of resource blocks allocated 
numRE = numRB*numREperRB;               % number of REs available for transmission 
numSymbols = numSlots*numSymbolSlot;    % Total number of OFDM symbols in the slot

%% Output variables

BER_vs_SNR = zeros(length(SNR_list),1);
BER_vs_QAM = zeros(length(SNR_list), length(QAMorderList), length(CFO));
EbN0_vs_QAM = zeros(length(SNR_list), length(QAMorderList));
BER_vs_QAM_theory = zeros(length(SNR_list), length(QAMorderList));

% Iterating over different Modulation orders

for M_index = 1:length(QAMorderList)
    QAMorder = QAMorderList(M_index);

    %% Generating the resource grid 

    RG = zeros(numRE, numSymbolSlot*numSlots);

    %% Generating transmit data

    % Computing the number of bits necessary
    numBits = numRE*log2(QAMorder)*numSymbolSlot*numSlots;

    % Generating the bits
    dataBits = randi([0,1], numBits, 1);

    % Generating QAM symbols 
    modulatedSymbols = qammod(dataBits, QAMorder, QAMencoding, "InputType",...
        "bit", "UnitAveragePower", true);

    % Populating the reosource grid
    RG = reshape(modulatedSymbols, size(RG));

    % Plotting the generated resource grid
    % plotResourceGrid(abs(RG), "RG", "OFDM symbols", "REs");

    %% Generating the FFT grid

    FFTgrid = zeros(FFTsize, numSymbolSlot*numSlots);

    % Generating the suitable index mapping. This is necessary because of the 
    % IFFT implementation in MATLAB. One needs to rearrange [0,a,b,c,d,e,f,0] 
    % to [d,e,f,0,0,a,b,c]. This occurs only along the first (RE) dimension.
    orgIndexSet = 0:numRE-1;
    newIndexSet = mod(orgIndexSet - numRE/2, FFTsize) + 1;
    FFTgrid(newIndexSet,:) = RG;

    % Plotting the generated FFT grid
    % plotResourceGrid(abs(FFTgrid), "FFT grid", "OFDM symbols", "REs");

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

    % Change input shape 
    timeDataSerial = squeeze(timeDataSerial);

    % Plot the spectrum of the OFDM signal
    if M_index == 1
        figure('DefaultAxesFontSize',22);
        pspectrum(timeDataSerial, SamplingRate);
        title("OFDM spectrum",'FontSize',22)
    end
    % Iterating over different CFO values
    for CFO_index = 1:length(CFO)
        % Iterating over SNR values
        BER_vs_SNR = [];                    % Re-initializing for each QAMorder
        for SNR_index = 1:length(SNR_list)
            %% Transmitting the time domain data through the AWGN channel
            SNR = SNR_list(SNR_index);
            channelOutput = awgn(timeDataSerial, SNR, "measured", 1234);

            %% Receiver

            % Introducing timing offset
            timing_offset_signal = zeros(timing_offset, 1);
            channelOutput = [timing_offset_signal;...
                channelOutput(1:end-timing_offset)];



            % Introducing carrier frequency offset
            residual_freq = exp(2*pi*1i*CFO(CFO_index)*(1/FFTsize).*...
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

            if TOcorrection
                indexArray = 1:FFTsize;
                indexArray = indexArray' - 1;
                TOcorrectionVector = exp(2*pi*1i*indexArray*timing_offset/FFTsize);
                RXFFTgrid = RXFFTgrid.*TOcorrectionVector;
            end

            %% Mapping FFT grid to RB grid
            % This is necessary to undo the mapping done before hand and to rearrange 
            % [d,e,f,0,0,a,b,c] to [0,a,b,c,d,e,f,0]. This occurs only along the 1st
            % (RE) dimension.

            RXfreqData = RXFFTgrid(newIndexSet, :);
            % plotResourceGrid(abs(RXfreqData), "RX RG", "OFDM symbols", "REs");

            %% Obtaining the transmitted bits

            RXfreqSerial = reshape(RXfreqData, [],1);

            % Plotting the received constellation for 16 QAM for the
                % highest SNR value
                if (QAMorder == 16)
                    if (SNR_index == length(SNR_list))
                        if (CFO_index == 1)
                            fig = figure('DefaultAxesFontSize',18);
                            scatter(real(RXfreqSerial), imag(RXfreqSerial));
                            titleValue = "Received 16 QAM constellation, TO: " + ...
                                string(timing_offset) + ", SNR: " + string(SNR) + ...
                                "dB, TO correction: " + string(TOcorrection);
                            title(titleValue);
                            xlabel("Real part");
                            ylabel("Imaginary part");
                            grid on;
                        end
                    end
                end
            RXbits = qamdemod(RXfreqSerial, QAMorder, QAMencoding, "OutputType",...
                "bit", "UnitAveragePower", true);

            %% BER

            BER = sum(RXbits ~= dataBits)/numBits;
            BER_vs_SNR(SNR_index) = BER;
        end
        BER_vs_QAM(:,M_index,CFO_index) = BER_vs_SNR;
    end
    EbN0_vs_QAM(:,M_index) = SNR_list ...
        + 10*log10(FFTsize/numRE) - 10*log10(log2(QAMorder));
    BER_vs_QAM_theory(:,M_index) =  berawgn(EbN0_vs_QAM(:,M_index),...
        'qam', QAMorder);
end

%% Plotting the obtained results

legendList = {};
fig = figure;

plot_count = 1;
for M_index = 1:length(QAMorderList)
    for CFO_index = 1:length(CFO)
        semilogy(EbN0_vs_QAM(:,M_index),...
            BER_vs_QAM(:,M_index,CFO_index),...
            "o--", 'LineWidth', 2);
        legendValue = string(QAMorderList(M_index)) + " QAM empirical, " + ...
            "CFO: " + string(CFO(CFO_index));
        legendList{1,plot_count} = legendValue;
        plot_count = plot_count + 1;
        hold on;
    end
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
titleValue = "BER curves, TO: " + string(timing_offset) + ...
    " samples, TO correction: " + string(TOcorrection);

title(titleValue,'FontSize', 22);
legend(legendList);
