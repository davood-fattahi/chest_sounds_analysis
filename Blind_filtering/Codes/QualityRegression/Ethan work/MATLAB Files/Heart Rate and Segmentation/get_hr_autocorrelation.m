    %% Autocorrelation
    function [signal_autocorrelation,signal_autocorrelation_filtered]=get_hr_autocorrelation(envelope,Fs)
        y=envelope-mean(envelope);
        [c] = xcorr(y,'coeff');
        signal_autocorrelation = c(length(envelope)+1:end);
        
        % Low-pass filter the autocorrelation:
        signal_autocorrelation_filtered = butterworth_low_pass_filter(signal_autocorrelation,1,15,Fs, false);
    end