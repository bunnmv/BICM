%% 16 QAM Viterbi decoder && 2/3 Convolutional enconding
% Author: Marcus Vinicius Bunn
% date: 28/04/2017
clc;
clear all;
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

% Modulacao
M = 16;
K = log2(M);

% Montando vetor para modula�ao
toMod = reshape(codeword,n,[]).';
txSig = qammod(codeword,M,'InputType','bit');

EbNo= -2:10;
berSoft = zeros(size(EbNo));
berHard = zeros(size(EbNo));

for n = 1:length(EbNo)    
    
    snr = EbNo(n) + 10*log10(K*codeRate);
    
    rxSig = awgn(txSig,snr,'measured');
        
    rxDataSoft = qamdemod(rxSig,M,'OutputType','llr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');
    
    dataSoft = vitdec(rxDataSoft,trellis,tbl,'cont','unquant');
    dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');

    [~,berSoft(n)]= biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHard(n)]= biterr(info(1:end-delay),dataHard(delay+1:end));
end

ber = berawgn(EbNo,'qam',M);

figure(1)
semilogy(EbNo,berSoft,EbNo,berHard,EbNo,ber);
title('16 QAM 2/3 CC AWGN')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD','Uncoded Ber');
