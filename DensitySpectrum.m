function [] = DensitySpectrum(modemQPSK, modemQORC, signalQORC)
    [f_qpsk, spectrumQPSK] = modemQPSK.GetSpectrum(false);

    [f_qorc, spectrumQORC] = modemQORC.GetSpectrum(signalQORC, false);
    
    
    spectrumQORC_ds = downsample(spectrumQORC(1:end-mod(length(f_qorc),length(f_qpsk))),floor(length(f_qorc)/length(f_qpsk)));
    f_qorc_ds = downsample(f_qorc(1:end-mod(length(f_qorc),length(f_qpsk))),floor(length(f_qorc)/length(f_qpsk)));

    figure
    hold on
        xlim([-9*10^4 9*10^4]);
        plot(f_qpsk, log(spectrumQPSK));
        plot(f_qorc_ds, log(spectrumQORC_ds));
    hold off
    
    downsampled_power = downsample(spectrumQPSK(length(spectrumQPSK)/2:end), 100);
    downsampled_power1 = downsample(spectrumQORC_ds(length(spectrumQORC_ds)/2:end), 100);

    Pob_QORC = zeros(1,floor(length(downsampled_power1)/2));
    Pob_QPSK = zeros(1,floor(length(downsampled_power)/2));

    for B=1:length(downsampled_power1)
        Pob_QORC(B)=sum(downsampled_power1(B:end))./sum(downsampled_power1(1:end));
    end
    for B1=1:length(downsampled_power)
        Pob_QPSK(B1)=sum(downsampled_power(B1:end))./sum(downsampled_power(1:end));
    end
    figure;
    semilogy(Pob_QORC,'g','LineWidth',1);
    hold on
    semilogy(Pob_QPSK,'r','LineWidth',1);
    legend('QORC','QPSK')
    grid on
end

