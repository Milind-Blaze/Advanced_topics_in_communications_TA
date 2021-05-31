%{
This script is used to simulate an OFDM system with a AWGN channel
according to 5G specifications. Further, timing offset and CFO are
implemented. Timing offset correction is also implemented. LDPC encoding
and decoding are also implemented.

When run, the script produces plots of BER vs Eb/N0 for each constellation
at the specified code rates.

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

%% Defining data to be transmitted

A = 42000;                                              % Number of bits to be transmitted
dataBits = randi([0,1], A, 1);                          % Data bits to be transmitted
QAMorderList = [16,64];                                     % the value of M in M-QAM, modulation used for data


% Code rates for different constellations from Table 5.1.3.1-2 of TS 38.214
R{4} = [120/1024, 193/1024, 308/1024, 449/1024, 602/1024];        
R{16} = [378/1024, 434/1024, 490/1024, 553/1024, 616/1024, 658/1024];
R{64} = [466/1024, 517/1024, 567/1024, 616/1024, 666/1024, 719/1024, 772/1024, 822/1024, 873/1024];
R{256} = [682.5/1024, 711/1024, 754/1024, 797/1024, 841/1024, 885/1024, 916.5/1024, 948/1024];

% Selected code rates 
R{16} = [434/1024, 616/1024];
R{64} = [466/1024, 873/1024];

% modulation correspondnig to QAMorderList
modulation{4} = 'QPSK';                         
modulation{16} = '16QAM';
modulation{64} = '64QAM';
modulation{256} = '256QAM';

% EbN0 for different modulation orders to obtain a BER of 1e-3
EbN0_list = {};
EbN0_list{4} = [-20:2:-1 -1:0.1:2 2:2:20];
EbN0_list{16} = [-20:2:1 1:0.085:5 5:2:20];
EbN0_list{64} = [-20:2:5 5:0.09:9 linspace(9,11,20) 11:2:20];
EbN0_list{256} = [-20:2:10 10:0.2:15 15:0.1:17   17:2:25];

%% Defining the 5G NR parameters

CPLen1 = 352;                                   % length of the cyclic prefix for the first OFDM symbol
CPLen2 = 288;                                   % length of the cyclic prefix for the OFDM symbols 2-14
FFTsize = 4096;                                 % size of the FFT being used
SCS = 30e3;                                     % Subcarrier spacing
SamplingRate = SCS*FFTsize;                     % symbol/sample rate for time domain data
numREperRB = 12;                                % number of resource elements per resource block


rv = 0;                         % 

%% Defining Resource Grid (RG) parameters

numRB = 50;                             % number of resource blocks allocated 
numRE = numRB*numREperRB;               % number of REs available for transmission 
nLayers = 1;                            % number of layers used for transmission

%% Output variables


BER_vs_QAM = {};
EbN0_vs_QAM = {};
BER_vs_QAM_theory = {};

% Iterating over different Modulation orders

for M_index = 1:length(QAMorderList)
    QAMorder = QAMorderList(M_index);
    BER_vs_QAM{M_index} = zeros(length(EbN0_list{QAMorder}), length(R{QAMorder}), length(CFO));
    BLER_vs_QAM{M_index} = zeros(length(EbN0_list{QAMorder}), length(R{QAMorder}), length(CFO));
    
    for R_index = 1:length(R{QAMorder})
        codeRate = R{QAMorder}(R_index);
        outlen = ceil(A/codeRate);

        %% Generating coded bits 

        cbsInfo = nrDLSCHInfo(A, codeRate);
        % Transport block CRC attachment
        tbIn = nrCRCEncode(dataBits, cbsInfo.CRC);

        % Code block segmentation and CRC attachment
        cbsIn = nrCodeBlockSegmentLDPC(tbIn, cbsInfo.BGN);

        % LDPC encoding
        enc = nrLDPCEncode(cbsIn, cbsInfo.BGN);

        % Rate matching and code block concatenation
        codedBits = nrRateMatchLDPC(enc, outlen, rv, modulation{QAMorder},...
            nLayers);

        %% Generating transmit data

        %% TODO: num(codedBits)%modOrder = 0?
        % Generating QAM symbols 
        modulatedSymbols = nrSymbolModulate(codedBits,modulation{QAMorder});
    %     size(modulatedSymbols)

        % Number of OFDM symbols needed to transmit given the frequency
        % resource grid constraints
        numSymbols = ceil(length(modulatedSymbols)/numRE);
        % Adding zero bits to fill up resource grid
        numDummySymbols = numSymbols*numRE - length(modulatedSymbols);
        dummyBits = zeros(numDummySymbols*log2(QAMorder),1);
        dummySymbols = nrSymbolModulate(dummyBits, modulation{QAMorder});
        numDummySymbols = length(dummySymbols);
        modulatedSymbols = [modulatedSymbols; dummySymbols];

        %% Generating the resource grid 

        RG = zeros(numRE, numSymbols);

        % Populating the reosource grid
        RG = reshape(modulatedSymbols, size(RG));

        % Plotting the generated resource grid
        % plotResourceGrid(abs(RG), "RG", "OFDM symbols", "REs");
        %% Generating the FFT grid

        FFTgrid = zeros(FFTsize, numSymbols);

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

        % Iterating over different CFO values
        for CFO_index = 1:length(CFO)
            % Iterating over SNR values
            BER_vs_SNR = [];                    % Re-initializing for each QAMorder
            
            for EbN0_index = 1:length(EbN0_list{QAMorder})
                %% Transmitting the time domain data through the AWGN channel
                SNR = EbN0_list{QAMorder}(EbN0_index) - 10*log10(FFTsize/numRE) + 10*log10(log2(QAMorder)) + 10*log10(R{QAMorder}(R_index));
                signal_power = rms(timeDataSerial).^2;
                noise_power = signal_power/(10^(SNR/10));
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

               
                RXfreqSerial = RXfreqSerial(1:end - numDummySymbols);
                
                RXCodedBits = nrSymbolDemodulate(RXfreqSerial, modulation{QAMorder}, noise_power);

                %% Decoding the bits
                RXCodedSoft = RXCodedBits;% double(1 - 2*RXCodedBits);
                % Rate recovery
                raterec = nrRateRecoverLDPC(RXCodedSoft, A, codeRate, rv,...
                    modulation{QAMorder}, nLayers);

                % LDPC decoding
                [decBits, actNumIter, finalParityChecks] = nrLDPCDecode(raterec,cbsInfo.BGN,25);

                % Code block desegmentation and CRC decoding
                [blk,blkErr] = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,A+cbsInfo.L);
    %             disp(['CRC error per code-block: [' num2str(blkErr) ']'])

                % Transport block CRC decoding
                [RXbits,tbErr] = nrCRCDecode(blk,cbsInfo.CRC);


                %% BER

                BER = sum(RXbits ~= dataBits)/A;
                BER_vs_SNR(EbN0_index) = BER;
            end
            BER_vs_QAM{M_index}(:,R_index,CFO_index) = BER_vs_SNR;
        end
        EbN0_vs_QAM{M_index} = EbN0_list{QAMorder};
    end
        
    BER_vs_QAM_theory{M_index} =  berawgn(EbN0_vs_QAM{M_index},...
        'qam', QAMorder);
end

%% Plotting the obtained results




for M_index = 1:length(QAMorderList)
    plot_count = 1;
    legendList = {};
    fig = figure;
    
    for CFO_index = 1:length(CFO)
        for R_index = 1:length(R{QAMorderList(M_index)})
            semilogy(EbN0_vs_QAM{M_index},...
                BER_vs_QAM{M_index}(:,R_index,CFO_index),...
                "o--", 'LineWidth', 2);
            legendValue = string(QAMorderList(M_index)) + " QAM empirical, " + ...
                "CFO: " + string(CFO(CFO_index) + ...
                ", R: " + string(R{QAMorderList(M_index)}(R_index)));
            legendList{1,plot_count} = legendValue;
            plot_count = plot_count + 1;
            hold on;
        end
    end
    
    semilogy(EbN0_vs_QAM{M_index},...
        BER_vs_QAM_theory{M_index},...
        'LineWidth', 2);
    legendValue = string(QAMorderList(M_index)) + " QAM theoretical AWGN";
    legendList{1,plot_count} = legendValue;
    plot_count = plot_count + 1;


    grid on
    ylim([1e-3,1]);
    xlabel("$10\log\left(\frac{E_b}{N_0}\right)$", "Interpreter", "latex", 'FontSize', 22);
    ylabel("BER", 'FontSize', 22);
    titleValue = "BER curves, TO: " + string(timing_offset) + ...
        " samples, TO correction: " + string(TOcorrection);

    title(titleValue,'FontSize', 22);
    legend(legendList);
    disp(legendList)
end

