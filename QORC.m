classdef QORC
    properties
        N
        d_Tb
        Tb
        input
    end
    methods(Access=private)
        function output = HardDecidion(obj, Ipart, Qpart, RC)
            RC_len = length(RC);
            output = [0 0 0 0];

            buffer = zeros(2, 60);
            part1=1:RC_len+1;
            part2=RC_len/2:(RC_len+RC_len/2);
            
            buffer(1, part1) = [0 RC];
            buffer(2, part2) = [0 RC];
            f1=sum(buffer);

            buffer(1, part1) = [0 RC];
            buffer(2, part2) = -[0 RC];
            f2=sum(buffer);

            buffer(1, part1) = -[0 RC];
            buffer(2, part2) = [0 RC];
            f3=sum(buffer);

            buffer(1, part1) = -[0 RC];
            buffer(2, part2) = -[0 RC];
            f4=sum(buffer);

            correlate_F1 = sum(f1(21:RC_len).*Qpart);
            correlate_F2 = sum(f2(21:RC_len).*Qpart);
            correlate_F3 = sum(f3(21:RC_len).*Qpart);
            correlate_F4 = sum(f4(21:RC_len).*Qpart);

            corr_max= max([correlate_F1 correlate_F2 correlate_F3 correlate_F4]);

            switch corr_max
                case correlate_F1
                    L_bit = 1;
                    R_bit = 1;
                case correlate_F2
                    L_bit = 1;
                    R_bit = 0;
                case correlate_F3
                    L_bit = 0;
                    R_bit = 1;
                case correlate_F4
                    L_bit = 0;
                    R_bit = 0;
            end

            output([2 4]) = [L_bit R_bit];

            correlate_F1 = sum(f1(21:RC_len).*Ipart);
            correlate_F2 = sum(f2(21:RC_len).*Ipart);
            correlate_F3 = sum(f3(21:RC_len).*Ipart);
            correlate_F4 = sum(f4(21:RC_len).*Ipart);

            corr_max= max([correlate_F1 correlate_F2 correlate_F3 correlate_F4]);

            switch corr_max
                case correlate_F1
                    L_bit = 1;
                    R_bit = 1;
                case correlate_F2
                    L_bit = 1;
                    R_bit = 0;
                case correlate_F3
                    L_bit = 0;
                    R_bit = 1;
                case correlate_F4
                    L_bit = 0;
                    R_bit = 0;
            end

            output([1 3]) = [L_bit R_bit];
        end
        function [I_samples, Q_samples] = GetIQ(obj)
            binary=obj.input;
            binary(binary==0)=-1;
            
            binary_reshaped = reshape(binary,[2,obj.N/2]);

            I_samples = binary_reshaped(1,:);
            Q_samples = binary_reshaped(2,:);
        end
        function [RaisedCosine] = GetRaisedCosine(obj, print)
            t=0:obj.d_Tb:4*obj.Tb-obj.d_Tb;
            RaisedCosine=(sin((pi.*t)./(4.*obj.Tb))).^2;

            if (print==true)
                figure
                plot(t,RaisedCosine)
            end
        end
    end
    methods
        function obj = QORC(bitLength, d_bitDuration, bitDuration, inputBits)
            obj.N=bitLength;
            obj.d_Tb=d_bitDuration;
            obj.Tb=bitDuration;
            obj.input = inputBits;
        end
        function signal = Modulation(obj, print)
            [RC] = obj.GetRaisedCosine(print);
            [I_samples, Q_samples] = obj.GetIQ();
            
            S_I=zeros(1,length(RC)+(length(I_samples))*length(RC)/2);
            S_Q=zeros(1,length(RC)+(length(Q_samples))*length(RC)/2);
            for k=1:obj.N/2
                S_I((((k-1)*length(RC)/2)+1):(k-1)*length(RC)/2+length(RC))=S_I((((k-1)*length(RC)/2)+1):(k-1)*length(RC)/2+length(RC)) + I_samples(k).* RC;
                S_Q((((k-1)*length(RC)/2)+1):(k-1)*length(RC)/2+length(RC))=S_Q((((k-1)*length(RC)/2)+1):(k-1)*length(RC)/2+length(RC)) + Q_samples(k).* RC;
            end

            inphaseQORC=zeros(1,length(RC)+(length(I_samples))*length(RC)/2 + 20);
            quadratureQORC=zeros(1,length(RC)+(length(Q_samples))*length(RC)/2 + 20);

            inphaseQORC(1:end - 20) = S_I;
            quadratureQORC(21:end) = S_Q;
            
            signal = inphaseQORC + 1i*quadratureQORC;
            if(print == true)
                figure
                tiledlayout(2,1)
                nexttile
                hold on
                plot((1:1:length(inphaseQORC))*obj.d_Tb, inphaseQORC)
                hold off
                grid on

                nexttile
                hold on
                plot((1:1:length(quadratureQORC))*obj.d_Tb, quadratureQORC)
                hold off
                grid on
               
            end
        end
        function [binary, errorRate] = Demodulation(obj, signal, SNR, print)
            noised_signal=awgn(signal,SNR,'measured');
            
            [RC] = obj.GetRaisedCosine(false);
            RC_len = length(RC);

            output = zeros(1,length(obj.input));
            
            for i = 1:floor(length(noised_signal)/RC_len)-1
                Ipart = real(noised_signal((i-1)*RC_len + 21: i*RC_len));
                Qpart = imag(noised_signal((i-1)*RC_len + 41: i*RC_len + 20));

                [outputBits]=obj.HardDecidion(Ipart,Qpart, RC);

                output((i-1)*4+1:(i-1)*4+4) = outputBits;
            end
            
            binary = output';
            errors = xor(obj.input, binary);
            errorRate=length(errors(errors==1))/obj.N;
            
            if(print==true)
            end
        end
        function [f, spectral_density] = GetSpectrum(obj, signal, print)
            baseband=real(signal)+imag(signal);
            
            Fs = 1/obj.Tb;
            y = fftshift(fft(baseband));
            
            n = length(baseband);
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

