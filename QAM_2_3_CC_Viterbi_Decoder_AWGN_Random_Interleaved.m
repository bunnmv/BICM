%% 16 QAM Viterbi decoder && 2/3 Convolutional enconding && Interleaver
% Author: Marcus Vinicius Bunn
% date: 03/05/2017

clc;clear all;
N = 100e3;         
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
txSig = qammod(codeword,M,'InputType','bit');

%interleaved
txSigInterleaved = qammod(interleaved,M,'InputType','bit');

EbNo= -2:20;
berSoft = zeros(size(EbNo));
berHard = zeros(size(EbNo));
berSoftInterleaved = zeros(size(EbNo));
berHardInterleaved = zeros(size(EbNo));


for n = 1:length(EbNo)    
    
    %% without interleaving
    snr = EbNo(n) + 10*log10(K*codeRate);
    
    rxSig = awgn(txSig,snr,'measured');
        
    rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');
        
    % 0 = the most confident 0 and 2^nsdec-1 = the most confident 1 
    % TO DO = func onde min(teste) = 0 e max(teste) = maximo_soft
    % meio = peso/2 como quantizar ?
    
    dataSoft = vitdec(rxDataSoft,trellis,tbl,'cont','unquant');
    dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');

    [~,berSoft(n)]= biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHard(n)]= biterr(info(1:end-delay),dataHard(delay+1:end));
    
    %% with interleaving
    
    rxSig = awgn(txSigInterleaved,snr,'measured');
        
    rxDataSoft = qamdemod(rxSig,M,'OutputType','approxllr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');
    
    softDeinter = randdeintrlv(rxDataSoft,state); % Deinterleave.
    hardDeinter = randdeintrlv(rxDataHard,state); % Deinterleave.
    
    dataSoft = vitdec(softDeinter,trellis,tbl,'cont','unquant');
    dataHard = vitdec(hardDeinter,trellis,tbl,'cont','hard');
    
    [~,berSoftInterleaved(n)]= biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHardInterleaved(n)]= biterr(info(1:end-delay),dataHard(delay+1:end));
end

figure(1)
semilogy(EbNo,berSoft,EbNo,berHard,EbNo,berSoftInterleaved,EbNo,berHardInterleaved);
title('16 QAM 2/3 CC AWGN and Interleaved')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD','SDDint','HDDint');

%     dataSoft = vitdec(rxDataSoft,trellis,tbl,'cont','unquant');
%   Quantize to prepare for soft-decision decoding.
%   qcode = quantiz(teste,[min(teste),min(teste)*peso/K,min(teste)*peso/(2*K),max(teste)/min(teste),max(teste)/(2*peso*K),max(teste)/(peso*K),max(teste)]);
    
