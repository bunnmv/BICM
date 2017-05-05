%% 16 QAM Viterbi decoder && 4/7 Hamming && Interleaved 
% Author: Marcus Vinicius Bunn
% date: 28/04/2017

clear all; close all; clc;

n = 7;

k = 4;

matriz = [1 0 0 0 1 1 0; 0 1 0 0 1 0 1; 0 0 1 0 0 1 1; 0 0 0 1 1 1 1];
       
palavras_cod = [0 0 0 0 0 0 0; 0 0 0 1 1 1 1; 0 0 1 0 0 1 1; 0 0 1 1 1 0 0; 0 1 0 0 1 0 1; 
                0 1 0 1 0 1 0; 0 1 1 0 1 1 0;0 1 1 1 0 0 1; 1 0 0 0 1 1 0; 1 0 0 1 0 0 1; 
                1 0 1 0 1 0 1;1 0 1 1 0 1 0;1 1 0 0 0 1 1;1 1 0 1 1 0 0; 1 1 1 0 0 0 0; 1 1 1 1 1 1 1];
            
palavras_cod_polar = 2*palavras_cod-1; 

code = struct('n',n,'k',k,'matriz',matriz,'vocabulario',palavras_cod_polar);

N = 100e3;         

info = randi([0 1], N, 4);

codificado = block_enc(code, info);

% interleaved
state = 4831;
interleaved = randintrlv(codificado,state);

%% Modulacao

M = 16;

K = log2(M);

modulado = qammod(codificado,M,'InputType','bit');

%interleaved
modulatedInterleaved = qammod(interleaved,M,'InputType','bit');

%% Soft Decision
EbNo = -2:10;

for l = 1:length(EbNo)
    
    %Without Interleaver
    recebido = awgn(modulado, EbNo(l)); 
        
    demodulado_soft = qamdemod(recebido,M,'OutputType','llr'); % - = 1; + = 0.

    sdd_decodificado = block_dec_sdd(code, demodulado_soft);    
 
    [~,BER_SDD(l)]= biterr(info,sdd_decodificado);
       
    demodulado_hard = qamdemod(recebido,M,'OutputType','bit'); 
     
    hdd_decodificado = block_dec_sdd(code, -(2*demodulado_hard-1));    
     
    [~,BER_HDD(l)]= biterr(info,hdd_decodificado);
    
    % With interleaver
    recebido = awgn(modulatedInterleaved, EbNo(l)); 
        
    demodulado_soft = qamdemod(recebido,M,'OutputType','llr'); % - = 1; + = 0.

    demodulado_soft_deinter=randdeintrlv(demodulado_soft,state); % Deinterleave.
    
    sdd_decodificado = block_dec_sdd(code, demodulado_soft_deinter);  
      
    [~,BER_SDD_Interleaved(l)]= biterr(info,sdd_decodificado);
       
    demodulado_hard = qamdemod(recebido,M,'OutputType','bit'); 
    
    demodulado_hard_deinter=randdeintrlv(demodulado_hard,state); % Deinterleave.
    
    hdd_decodificado = block_dec_sdd(code, -(2*demodulado_hard_deinter-1));    
     
    [~,BER_HDD_interleaved(l)]= biterr(info,hdd_decodificado);
end

figure(1)
semilogy(EbNo,BER_SDD,EbNo,BER_HDD,EbNo,BER_SDD_Interleaved,EbNo,BER_HDD_interleaved);
title('BER simulation of SDD 16 QAM (7,4) Hamming code info')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD','SDD-Interleaved','HDD-Interleaved');

