%% 16 QAM Viterbi decoder && 2/3 Convolutional enconding && Rayleigh && Interleaver
% Author: Marcus Vinicius Bunn
% date: 03/05/2017

clc;
clear all;
N = 200e3;         
k = 2;
n = 3;
codeRate = k/n;

info = randi([0 1], N*k, 1);

% Trellis
% If the encoder diagram has k inputs and n outputs, the code generator matrix is a k-by-n matrix. 
% The element in the ith row and jth column 
% indicates how the ith input contributes to the jth output.
% 3 = 1 input + 2 shift registers adding on that input
%  1     1     0  ones representing shift-register adders
%  1     1     1
% Deseja-se n = 4 e k =7.
% 
% trellis = poly2trellis(7,[171 133]);
% tbl = 32;

trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
tbl=16;
delay = k*tbl;

% Convolutional enconding
codeword = convenc(info,trellis);

% interleaved
state = 4831;
interleaved = randintrlv(codeword,state);

% Modulacao
M = 16;
K = log2(M);
modulated = qammod(codeword,M,'InputType','bit');
%interleaved
modulatedInterleaved = qammod(interleaved,M,'InputType','bit');

% Rayleigh Fading Channel 
Ts = 1/100000;
fd = 130;
h = rayleighchan(Ts, fd);
h.ResetBeforeFiltering = 0;
h.StoreHistory=1;
txSig = filter(h,modulated);
channel_gains = h.PathGains;

%interleaved
txSigInterleaved = filter(h,modulatedInterleaved);
channel_gainsInterleaved = h.PathGains;


EbNo= -2:20;
berSoft = zeros(size(EbNo));
berHard = zeros(size(EbNo));
berSoftInterleaved = zeros(size(EbNo));
berHardInterleaved = zeros(size(EbNo));



for n = 1:length(EbNo)    
    
    %% without interleaving
    snr = EbNo(n) + 10*log10(K*codeRate);
    
    with_noise = awgn(txSig,snr,'measured');
    
    rxSig = with_noise./channel_gains;
    
    rxDataSoft = qamdemod(rxSig,M,'OutputType','llr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');
    
    dataSoft = vitdec(rxDataSoft,trellis,tbl,'cont','unquant');
    dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');

    [~,berSoft(n)] = biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHard(n)] = biterr(info(1:end-delay),dataHard(delay+1:end));
    
    %% with interleaving
        
    with_noise = awgn(txSigInterleaved,snr,'measured');
    
    rxSig = with_noise./channel_gainsInterleaved;
    
    rxDataSoft = qamdemod(rxSig,M,'OutputType','llr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');

    softDeinter = randdeintrlv(rxDataSoft,state); % Deinterleave.
    hardDeinter = randdeintrlv(rxDataHard,state); % Deinterleave.
    
    dataSoft = vitdec(softDeinter,trellis,tbl,'cont','quant','soft',3);
    dataHard = vitdec(hardDeinter,trellis,tbl,'cont','hard');
    
    [~,berSoftInterleaved(n)]= biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHardInterleaved(n)]= biterr(info(1:end-delay),dataHard(delay+1:end));
end

figure(1)
semilogy(EbNo,berSoft,EbNo,berHard,EbNo,berSoftInterleaved,EbNo,berHardInterleaved);
title('BER simulation of SDD 16 QAM 2/3 CC on Rayleigh fadding channel and interleaved')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD','SDDint','HDDint');
