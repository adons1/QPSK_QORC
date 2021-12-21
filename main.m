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

SNR=(0:1.0:15.0);
[errsQPSK]=model(modemQPSK, signalQPSK, SNR);
SNR_for_QORC=(0:1.0:15.0)-7.0;
[errsQORC]=model(modemQORC, signalQORC, SNR_for_QORC);
toc

 figure
% hold on
semilogy(SNR-10*log10(2), errsQPSK, SNR-10*log10(2), errsQORC);
% plot(SNR-10*log10(2), log10(errsQPSK));
% plot(SNR+10*log10(2), log10(errsQORC));
legend('QPSK','QORC');
% hold off




