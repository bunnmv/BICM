%% 16 QAM - 4/7 Hamming && Rayleigh
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

N = 100e3;         

info = randi([0 1], N, k);

codificado = block_enc(code, info);


%% Modulacao

M = 16;

K = log2(M);

modulado = qammod(codificado,M,'InputType','bit');

%% Rayleigh Fading Channel 
Ts = 1/100000;
fd = 130;
h = rayleighchan(Ts, fd);
modulado_vetor = reshape(modulado.',1,[]);
h.StoreHistory=1;
y = filter(h,modulado_vetor);
y_message = reshape(y,n,[]).';
channel_gains = reshape(h.PathGains,n,[]).';

%% BER calculation
EbNo = -2:10;

for l = 1:length(EbNo)
    
    recebido = awgn(y_message, EbNo(l)); 
    
    equalizado = recebido./channel_gains;
    
    demodulado = qamdemod(equalizado,M,'OutputType','llr'); % - = 1; + = 0.
    
    %Soft decision
    sdd_decodificado = block_dec_sdd(code, demodulado);    
    
    [~,BER_SDD(l)]= biterr(info,sdd_decodificado);
       
    demodulado_hard = qamdemod(equalizado,M,'OutputType','bit'); % - = 1; + = 0.

    %hard decision
    hdd_decodificado = block_dec_sdd(code, -(2*demodulado_hard-1));    
     
    [~,BER_HDD(l)]= biterr(info,hdd_decodificado);
    
end

figure(1)
semilogy(EbNo,BER_SDD,EbNo,BER_HDD);
title('16 QAM (7,4) Hamming Rayleigh')
ylabel('Pb')
xlabel('Eb/No')
legend('SDD','HDD');
