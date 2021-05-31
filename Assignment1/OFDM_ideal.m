%{ 
OFDM simulation with ideal transmitter and receiver

Simulates 4 QAM, 16 QAM, 64 QAM and 256 QAM and plots the waterfall curves 
under the  assumption of an AWGN channel.

Author: Milind Kumar Vaddiraju 
%}

%% Defining variables
clc; clear; close all;

num_bits = 3*2^20           % Number of bits to be transmitted 
% more the number of bits, the smoother the curve at lower errors
N = 2048                    % number of sub carriers
cp_len = 120                % length of cyclic prefix

num_iter = 100              % number of iterations
N0_up = 1                 % upper bound on noise power
N0_lb = -2.5                % lower bound on noise power

N0s = logspace(N0_lb, N0_up, num_iter);     % list of noise powers


Ms = [4, 16, 64, 256];       % Modulation orders
X = [];                      % X coordinates for the final plots
Y = [] ;                     % Y coordinates for the final plots

timing_offset = 0;          % timing offset to the left in the received data


%% Simulation 

% Generate bits
rng("default");                     % setting the seed (sort of)
bits = randi([0,1], num_bits,1);    % Generating a column vector of bits

for mod_index = 1:length(Ms)
    
    M = Ms(mod_index);

    %% Simulation for a given modulation order
    BERs = [];           % list of bit errors rates for different noise powers
    for noise_index = 1:length(N0s)

        N0 = N0s(noise_index);
        


        %% Modulate bits 

        % specify gray coding for consistency during demodulation
        modulated_symbols = qammod(bits, M, "gray", "InputType", "bit", ...
            "UnitAveragePower", false);  

        %% IFFT 
        % There is no notion of useful subcarriers. All N subcarriers carry
        % information.

        % Converting from serial to parallel
        modulated_sym_parallel = reshape(modulated_symbols, N, []);

        % IFFT - ignoring effect of normalizing constant
        time_domain_symbols = ifft(modulated_sym_parallel, N, 1);

        %% CP addition
        
        % CP addition for all symbols at once
        transmit_signal_parallel = [time_domain_symbols(end - (cp_len -1):end, :); time_domain_symbols];

        % converting from parallel to serial
        transmit_signal = reshape(transmit_signal_parallel, [],1);

        % Add CP 
%         transmit_signal = [transmit_serial(end - (cp_len -1):end); transmit_serial];


        %% Channel
        
        E = 2/3*(M - 1);
        SNR = 10*log10(E/N0);
        received_signal = awgn(transmit_signal, SNR, "measured", 1234); % 1234 seed
        
        
        %% Timing offset 
        
        % Offset is introduced by prepending timing_offset number of zeros
        % to the signal and ignoring the last timing_offset number of
        % values
        
        offset_signal = zeros(timing_offset, 1);
        received_signal = [offset_signal; received_signal(1:end-timing_offset)];
        %% CP removal 
        
        % Converting from serial to parallel
        received_parallel = reshape(received_signal, N + cp_len, []);
        
        % CP removal for all symbols at once
        received_parallel = received_parallel(cp_len + 1: end, :);
       
        

        %% FFT
        
        % Obtaining FFT
        received_freq_domain = fft(received_parallel, N);
        
        %% Correction for timing offset
        
        indexArray = 1:N;
        indexArray = indexArray' - 1;
        correctionVector = exp(2*pi*1i*indexArray*timing_offset/N);
        received_freq_domain = received_freq_domain.*correctionVector;

        %% Demodulation 

        % Converting from parallel to serial
        received_freq_serial = reshape(received_freq_domain, [], 1);

        % Demodulating
        received_bits = qamdemod(received_freq_serial, M, "gray",...
            "OutputType", "bit"); 

        %% Comparison of transmitted and received bits

        BER = 1 - (sum(bits == received_bits)/num_bits);
        BERs = [BERs BER];
    end
    Eb = E/log2(M);

    SNRs = 10*log10(Eb*N0s.^-1);
    X = [X; SNRs];
    Y = [Y; BERs];
end

%% reshape the arrays for plotting
X = X';
Y = Y';

%% Plotting a waterfall curve
fig = figure('DefaultAxesFontSize',18)
semilogy(X, Y, "o", 'LineWidth', 2);
title("Waterfall curves")%, "FontSize", 22);
xlabel("$10\log\left(\frac{E_b}{N_0}\right)$", "Interpreter", "latex");
ylabel("BER");
ylim([1e-6, 1e1])
grid on
% legend("M_{sim} = 4", "M_{sim} = 16", "M_{sim} = 64", "M_{sim} = 256")
hold on
%% Add theoretical curves
% TODO: vectorise this
Y_theoretical = [];
for i = 1:length(Ms)
    m = Ms(i);
    x = X(:,i);
    x = x/10;
    x = 10.^x;
    y = qfunc(sqrt(3*log2(m)/(m-1)*x));
    y = 4*(1 - 1/sqrt(m))*y;
%     y = y - y.^2/4;
    y = y/log2(m);
    Y_theoretical = [Y_theoretical, y];
end
semilogy(X, Y_theoretical, "-", 'LineWidth',2)
legend("M_{sim} = 4", "M_{sim} = 16", "M_{sim} = 64", "M_{sim} = 256",...
    "M_{theory} = 4", "M_{theory} = 16", "M_{theory} = 64", "M_{theory} = 256",...
    "NumColumns", 2)


