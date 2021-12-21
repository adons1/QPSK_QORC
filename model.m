function [Errors] = model(modem, signal, SNR)
    Errors=zeros(1,length(SNR));
    index=1;
    for snr=SNR
        [output, Errors(index)] = modem.Demodulation(signal, snr, false);
        fprintf('SNR:\t%f\nPb:\t%f\n\n', snr-10*log10(2), Errors(index));

        index=index+1;
    end
end

