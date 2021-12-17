function [] = modelQORC(modem, signal, SNR)
    Errors=zeros(1,length(SNR));
    index=1;
    for snr=SNR
        [output, Errors(index)] = modem.Demodulation(noised_signal,false);
        fprintf('SNR:\t%f\nPb:\t%f\n\n', snr-10*log10(2), Errors(index));

        index=index+1;
    end
    figure;
    grid on;
    semilogy(SNR-10*log10(2),Errors)
end
