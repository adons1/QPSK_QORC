classdef QPSK
    properties
        N
        Tb
        input
    end
    methods(Access=private)
        function output = HardDecidion(obj, symbol)
            
            resultI=0;
            resultQ=0;
            if((real(symbol)>0) && (imag(symbol)>0))
                resultI= 1;
                resultQ= 1;
            end
            if((real(symbol)>0) && (imag(symbol)<0))
                resultI= 0;
                resultQ= 1;
            end
            if((real(symbol)<0) && (imag(symbol)>0))
                resultI= 1;
                resultQ= 0;
            end
            if((real(symbol)<0) && (imag(symbol)<0))
                resultI= 0;
                resultQ= 0;
            end
            output=[resultI resultQ];
        end
        function [I_samples, Q_samples] = GetIQ(obj)
            binary=obj.input;
            binary(binary==0)=-1;
            
            binary_reshaped = reshape(binary,[2,obj.N/2]);

            I_samples = binary_reshaped(1,:);
            Q_samples = binary_reshaped(2,:);
        end
    end
    methods
        function obj = QPSK(bitLength, bitDuration, inputBits)
            obj.N=bitLength;
            obj.Tb=bitDuration;
            obj.input = inputBits;
        end
        function signal = Modulation(obj, print)
            [I_samples, Q_samples] = obj.GetIQ();

            signal = complex(Q_samples, I_samples);
            
            if(print == true)
                figure
                hold on
                xlim([-2.0 2.0]);
                ylim([-2.0 2.0]);
                scatter(Q_samples, I_samples)
                hold off
            end
        end
        function [binary, errorRate] = Demodulation(obj, signal, SNR, print)
            noised_signal=awgn(signal,SNR,'measured');
            if(print == true)
                figure
                hold on
                xlim([-2.0 2.0]);
                ylim([-2.0 2.0]);
                scatter(real(noised_signal), imag(noised_signal))
                hold off
            end
            
            offset=0;
            output=zeros(1, length(obj.N));
            for i = 1:length(noised_signal)
                decidion_result = obj.HardDecidion(noised_signal(i));
                output(offset+i:offset+i+1) = decidion_result;
                offset=offset+1;
            end
            
            binary = output;
            
            errors = xor(obj.input, output');
            errorRate=length(errors(errors==1))/obj.N;
        end
        
        function [f, spectral_density] = GetSpectrum(obj, print)
            [I_samples, Q_samples] = obj.GetIQ();
            
            upsample_N=18;
            
            bi_upsampled = upsample(I_samples, upsample_N+1, 0);
            bq_upsampled = upsample(Q_samples, upsample_N+1, 0);
            
            bi_s = zeros(1, length(bi_upsampled));
            bq_s = zeros(1, length(bq_upsampled));
                        
            for i=1:upsample_N
                bi_s= bi_s + circshift(bi_upsampled,i);
                bq_s= bq_s + circshift(bq_upsampled,i);
            end
            
            basebandQPSK=bi_s+bq_s;
            
            Fs = 1/obj.Tb;
            y = fftshift(fft(basebandQPSK));
            
            n = length(basebandQPSK);
            f = (-n/2:n/2-1)*(Fs/(obj.N/2));   
            spectral_density = abs(y).^2/n;  
            
            if(print==true)
                figure
                hold on
                grid on

                plot(f,log(spectral_density))
            end
        end
    end
end

