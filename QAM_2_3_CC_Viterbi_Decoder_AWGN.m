%% 16 QAM Viterbi decoder && 2/3 Convolutional enconding
% Author: Marcus Vinicius Bunn
% date: 28/04/2017
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

% Modulacao
M = 16;
K = log2(M);

% Montando vetor para modulaçao
toMod = reshape(codeword,n,[]).';
txSig = qammod(codeword,M,'InputType','bit');

EbNo= -2:10;
berSoft = zeros(size(EbNo));
berHard = zeros(size(EbNo));
ndsec=3;
maximo_soft = 2^ndsec-1;
peso = maximo_soft/ndsec;

for n = 1:length(EbNo)    
    
    with_noise = awgn(txSig,EbNo(n),'measured');
    
    rxSig = with_noise./channel_gains;
    
    rxDataSoft = qamdemod(rxSig,M,'OutputType','llr'); % -1 = 1 + = 0.
    rxDataHard = qamdemod(rxSig,M,'OutputType','bit');
    
    teste = rxDataSoft*-1; % - = 0 e 1 = 1;
    
    % 0 = the most confident 0 and 2^nsdec-1 = the most confident 1 
    % TO DO = func onde min(teste) = 0 e max(teste) = maximo_soft
    % meio = peso/2 como quantizar ?
    
    qcode = quantiz(teste,[0.001,.1,.3,.5,.7,.9,.999]);

    dataSoft = vitdec(qcode,trellis,tbl,'cont','soft',ndsec);
    dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');

    [~,berSoft(n)]= biterr(info(1:end-delay),dataSoft(delay+1:end));
    [~,berHard(n)]= biterr(info(1:end-delay),dataHard(delay+1:end));
end

figure(1)
semilogy(EbNo,berSoft,EbNo,berHard);
title('BER simulation of SDD 16 QAM 2/3 CC')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD');
