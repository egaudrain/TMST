function [AMsgram, fc, scale, step] = AMwavelet(insig, fs, kv, flags)
% AMwavelet Amplitude modulation wavelet spectrogram
%   [AMSGRAM, FC, SCALE, STEP] = AMwavelet(INSIG, FS)
%       Returns the AM wavelet spectrogram (wavelet amplitude) of signal
%       INSIG with sampling frequency FS.
%
%   Optional arguments:
%       flow: lower frequency of audio filterbank
%       fhigh: upper frequency of audio filterbank
%       mflow: lower modulation frequency
%       mfhigh: upper modulation frequency
%       modbank_Nmod: number of modulation bands
%       modbank_Qfactor: sharpness of the modulation filterbank filters
%       do_silent: true or false
%
%   Returns:
%       AMSGRAM is a N-by-M function where N is the number of modulation
%           frequencies and M is the number of audio frequencies (FC).
%       FC is a vector of center frequencies for the frequency bands.
%       SCALE the time points.
%       STEP (optional) a structure with some more details.
%
% Leo Varnet - 07/2023

    arguments
        insig (:,1) double
        fs {mustBeNumeric,mustBePositive}
        kv.flow double = NaN
        kv.fhigh double = NaN
        kv.mflow double = NaN
        kv.mfhigh double = NaN
        kv.modbank_Nmod double = NaN
        kv.modbank_Qfactor double = NaN
        flags.do_silent logical = 1
    end

    [kv, flags] = arg_varnet2017(kv, flags);

    % defines the modulation axis
    mflow  = kv.mflow;
    mfhigh = kv.mfhigh;
    % N_fsamples = kv.modbank_Nmod; 
    % f_spectra_intervals = logspace(log10(mflow), log10(mfhigh), N_fsamples+1);
    % f_spectra           = logspace(log10(sqrt(f_spectra_intervals(1)*f_spectra_intervals(2))), log10(sqrt(f_spectra_intervals(end)*f_spectra_intervals(end-1))), N_fsamples);
    % 
    % Number of steps in fractional octaves to go from mod_flow to mod_fhigh:
    % N_octave_steps = ceil(NthOct * log10(mod_fhigh/mod_flow)/log10(2));
    % mfc            = mod_flow * 2.^((0:N_octave_steps)/NthOct); % includes one extra bin (because of the ceiling)
    % cutoff_oct     = mfc*2^(-.5/NthOct); % half step down
    %%%
    
    t=(1:length(insig))/fs;
    
    %%% gammatone filtering
    % [gamma_responses,fc] = auditoryfilterbank(insig,fs,kv.flow,kv.fhigh);
    % f_bw = audfiltbw(fc);
    nERBs = ceil(diff(ERBn_number([kv.flow, kv.fhigh])));
    gammaFiltBank = gammatoneFilterBank([kv.flow, kv.fhigh], nERBs, 'SampleRate', fs);
    gamma_responses = gammaFiltBank(insig);
    fc = gammaFiltBank.getCenterFrequencies();
    f_bw = gammaFiltBank.getBandwidths();
    
    %%% AM extraction
    if ~flags.do_silent
        fprintf('E extraction\n');
    end
    E = abs(hilbert(squeeze(gamma_responses)));
    
    Nchan = length(fc);
    %%% AM spectra
    if ~flags.do_silent
        fprintf('calculating envelope wavlets\n');
    end
    fb = cwtfilterbank('SignalLength',length(t),'SamplingFrequency',fs,'FrequencyLimits',[mflow mfhigh],'Wavelet','morse');%'TimeBandwidth',120,
    
    %freqz(fb)
    for ichan=1:Nchan
        [wt,scale] = cwt(E(:,ichan),'FilterBank',fb);
        if ichan == 1
            AMsgram = abs(wt)'/Nchan;
        else
            AMsgram = AMsgram + abs(wt)'/Nchan;
        end
    end
    
    if ~flags.do_silent
        figure; h = pcolor(t,scale,abs(AMsgram(:,:)));
        set(h,'EdgeColor', 'none');
        set(gca, 'YScale', 'log');
    end
    
    if nargout>3
        step.t = t;
        step.f_bw = f_bw;
        step.gamma_responses = gamma_responses;
        step.E = E;
        step.scale = scale;
    end

end
