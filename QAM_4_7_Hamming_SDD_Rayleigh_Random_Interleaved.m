%% 16 QAM - 4/7 Hamming && Rayleigh && Interleaved
% Author: Marcus Vinicius Bunn
% date: 03/05/2017

clear all; close all; clc;

n = 7;

k = 4;

matriz = [1 0 0 0 1 1 0; 0 1 0 0 1 0 1; 0 0 1 0 0 1 1; 0 0 0 1 1 1 1];
       
palavras_cod = [0 0 0 0 0 0 0; 0 0 0 1 1 1 1; 0 0 1 0 0 1 1; 0 0 1 1 1 0 0; 0 1 0 0 1 0 1; 
                0 1 0 1 0 1 0; 0 1 1 0 1 1 0;0 1 1 1 0 0 1; 1 0 0 0 1 1 0; 1 0 0 1 0 0 1; 
                1 0 1 0 1 0 1;1 0 1 1 0 1 0;1 1 0 0 0 1 1;1 1 0 1 1 0 0; 1 1 1 0 0 0 0; 1 1 1 1 1 1 1];
            
palavras_cod_polar = 2*palavras_cod-1; 

code = struct('n',n,'k',k,'matriz',matriz,'vocabulario',palavras_cod_polar);

N = 200e3;         

info = randi([0 1], N, k);

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

%% Rayleigh Fading Channel 
Ts = 1/100000;
fd = 130;
h = rayleighchan(Ts, fd);
modulado_vetor = reshape(modulado.',1,[]);
h.StoreHistory=1;
y = filter(h,modulado_vetor);
y_message = reshape(y,n,[]).';
channel_gains = reshape(h.PathGains,n,[]).';

%interleaved
modulado_vetor_interleaved = reshape(modulatedInterleaved.',1,[]);
txSigInterleaved = filter(h,modulado_vetor_interleaved);
txSigInterleaved_message = reshape(txSigInterleaved,n,[]).';
channel_gainsInterleaved = reshape(h.PathGains,n,[]).';

%% BER calculation
EbNo = -2:20;

for l = 1:length(EbNo)
    %% Without interleaver
    recebido = awgn(y_message, EbNo(l)); 
    
    equalizado = recebido./channel_gains;
    
    demodulado_soft = qamdemod(equalizado,M,'OutputType','llr'); % - = 1; + = 0.
    
    %Soft decision
    sdd_decodificado = block_dec_sdd(code, demodulado_soft);    
    
    [~,BER_SDD(l)]= biterr(info,sdd_decodificado);
       
    demodulado_hard = qamdemod(equalizado,M,'OutputType','bit'); % - = 1; + = 0.

    %hard decision
    hdd_decodificado = block_dec_sdd(code, -(2*demodulado_hard-1));    
     
    [~,BER_HDD(l)]= biterr(info,hdd_decodificado);
    
    %% With Interleaver
    recebido = awgn(txSigInterleaved_message, EbNo(l)); 
    
    equalizado = recebido./channel_gainsInterleaved;
     
    %Soft decision
    demodulado_soft = qamdemod(equalizado,M,'OutputType','llr'); % - = 1; + = 0.
    
    softDeinter = randdeintrlv(demodulado_soft,state); % Deinterleave.
       
    sdd_decodificado = block_dec_sdd(code, softDeinter);    
    
    [~,BER_SDD_Inter(l)]= biterr(info,sdd_decodificado);
       
    %hard decision
    demodulado_hard = qamdemod(equalizado,M,'OutputType','bit'); % - = 1; + = 0.
    
    hardDeinter = randdeintrlv(demodulado_hard,state); % Deinterleave.
   
    hdd_decodificado = block_dec_sdd(code, -(2*hardDeinter-1));    
     
    [~,BER_HDD_Inter(l)]= biterr(info,hdd_decodificado);
end

figure(1)
semilogy(EbNo,BER_SDD,EbNo,BER_HDD,EbNo,BER_SDD_Inter,EbNo,BER_HDD_Inter);
title('16 QAM (7,4) Hamming Rayleigh and Interleaved')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD','SDD-Inter','HDD-Inter');
