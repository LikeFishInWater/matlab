% Demonstration of Eb/N0 Vs BER for BPSK modulation scheme
clear;
clc;
%---------Input Fields------------------------

K = 2048; % Length of Code block
Ncb = 2000;
N = K * Ncb; %Number of input bits
%EbN0dB = 0:5:30; % Eb/N0 range in dB for simulation
ModulType = 'BPSK';
ChlType = 'RAYLEIGH';
Rc=1/3; %Rc = code rate for a coded system. Since no coding is used Rc=1



% Turbo code system over AWGN channel
EbN0dB = 1.2;
Err = zeros(1,length(EbN0dB)); %Place holder for BER values for each Eb/N0
for i=1:length(EbN0dB),
    Err(i) = myTurboCodeSimu(K, Ncb, Rc, EbN0dB(i), ModulType, 'AWGN', 6);
end
BER_TC_AWGN = Err./N;
semilogy(EbN0dB,BER_TC_AWGN,'bs:','LineWidth',2);hold on;



% Turbo code system over rayleigh fading channel
EbN0dB = 30;
Err = zeros(1,length(EbN0dB)); %Place holder for BER values for each Eb/N0
for i=1:length(EbN0dB),
    Err(i) = myTurboCodeSimu(K, Ncb, Rc, EbN0dB(i), ModulType, 'RAYLEIGH', 6);
end
BER_TC_rayleigh = Err./N;
semilogy(EbN0dB,BER_TC_rayleigh,'rd:','LineWidth',2);hold on;


EbN0dB = 0:1:10;
Err = zeros(1,length(EbN0dB)); %Place holder for BER values for each Eb/N0
for i=1:length(EbN0dB),
    Err(i) = UncodedSystemSimu(N, EbN0dB(i), ModulType, 'AWGN');
end
BER_UC_AWGN = Err./N;
semilogy(EbN0dB,BER_UC_AWGN,'bo-','LineWidth',2);hold on;

EbN0dB = 0:5:30;
Err = zeros(1,length(EbN0dB)); %Place holder for BER values for each Eb/N0
for i=1:length(EbN0dB),
    Err(i) = UncodedSystemSimu(N, EbN0dB(i), ModulType, 'RAYLEIGH');
end
BER_UC_rayleigh = Err./N;
semilogy(EbN0dB,BER_UC_rayleigh,'r*-','LineWidth',2);hold on;

%EbN0=10.^(EbN0dB/10); %Eb/N0 in Linear Scale

% semilogy(EbN0dB,BER_UC_AWGN,'bo-','LineWidth',2);hold on;
% semilogy(EbN0dB,BER_UC_rayleigh,'r*-','LineWidth',2);hold on;
axis([0 30 10^-5 1.2]);
legend('TC over AWGN','TC over Rayleigh','UC over AWGN','UC over Rayleigh');
title('Eb/N0 Vs BER for BPSK over Rayleigh and AWGN Channels');
xlabel('Eb/N0(dB)');
ylabel('Bit Error Rate');
