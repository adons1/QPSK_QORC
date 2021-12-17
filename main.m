clear; clc; close all;

res= zeros(1,4);
N=10^6;
Tb=10^-4;
binary = randi([0 1],N,1);

tic

modemQPSK = QPSK(N, Tb, binary);
signalQPSK = modemQPSK.Modulation(false);
toc

modemQORC = QORC(N, Tb/10, Tb, binary);
signalQORC = modemQORC.Modulation(false);
toc

DensitySpectrum(modemQPSK, modemQORC, signalQORC);
toc

SNR=0:1.0:15.0;
[errsQPSK]=model(modemQPSK, signalQPSK, SNR);
[errsQORC]=model(modemQORC, signalQORC, SNR);
toc

hold on
semilogy(SNR-10*log10(2), errsQPSK, SNR-10*log10(2), errsQORC);
legend('QORC','QPSK');
hold off




