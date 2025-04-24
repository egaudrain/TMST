function [AMIspec, fc, mf, step] = AMIspectrum(insig, fs, kv, flags)
% AMIspectrum Amplitude modulation excitation pattern
%   [AMIspec, fc, mf, step] = AMIspectrum(INSIG, FS)
%       Returns the AMi spectrum of signal INSIG (with sampling frequency FS)
%       in excitation units (W).
%
%   Optional arguments:
%       flow: lower frequency of audio filterbank
%       fhigh: upper frequency of audio filterbank
%       mflow: lower modulation frequency
%       mfhigh: upper modulation frequency
%       modbank_Nmod: number of modulation bands
%       modbank_Qfactor: sharpness of the modulation filterbank filters
%       do_LP_150_Hz: if true, first order modulation Butterworth lowpass filter with a cut-off
%           frequency of 150 Hz. This is to remove all modulation frequencies
%           above 150 Hz. The motivation behind this filter can be found in
%           kohlrausch2000.
%       do_phase_insens_hilbert: if true, a phase insensitive Hilbert
%           transform is used.
%       do_silent: true or false
%
%   Returns:
%       AMIspec is a N-by-M function where N is the number of modulation
%           frequencies and M is the number of audio frequencies (FC).
%       fc is a vector of center frequencies for the frequency bands.
%       mf are the modulation frequencies.
%       step (optional) a structure with some more details.
%
% See Varnet et al. 2017 for more details
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
        flags.do_LP_150_Hz logical = []
        flags.do_phase_insens_hilbert logical = []
        flags.do_silent logical = 1
    end
    
    % 
    % 
    % if nargin<2
    %   error('%s: Too few input arguments.',upper(mfilename));
    % end
    % 
    % if ~isnumeric(insig) 
    %   error('%s: insig must be numeric.',upper(mfilename));
    % end
    % 
    % if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    %   error('%s: fs must be a positive scalar.',upper(mfilename));
    % end
    
    %{
    definput.import={'varnet2017'}; 
    definput.importdefaults={}; 
    
    [flags,kv]  = ltfatarghelper({'flow','fhigh'},definput,varargin);
    %}
    
    [kv, flags] = arg_varnet2017(kv, flags);
    
    % defines the modulation axis
    % mflow  = kv.mflow;
    % mfhigh = kv.mfhigh;
    % N_fsamples = kv.modbank_Nmod; 
    % f_spectra_intervals = logspace(log10(mflow), log10(mfhigh), N_fsamples+1);
    % f_spectra           = logspace(log10(sqrt(f_spectra_intervals(1)*f_spectra_intervals(2))), log10(sqrt(f_spectra_intervals(end)*f_spectra_intervals(end-1))), N_fsamples);
    
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
    
    %Nchan = length(fc);
    %%% AMi spectra
    if ~flags.do_silent
        fprintf('calculating envelope spectra\n');
    end
    
    % using the king2019_modfilterbank function
    [AMfilt, mf] = king2019_modfilterbank_updated(E,fs,kv,flags); 
    
    % % using the modfilterbank function
    % [AMfilt_temp, mf] = modfilterbank(E,fs,fc,'argimport',flags,kv);
    % AMfilt=nan(length(t),length(fc),length(mf));
    % for i=1:length(AMfilt_temp)
    %     AMfilt(:,i,1:size(AMfilt_temp{i},2)) = abs(AMfilt_temp{i});
    % end
    
    AMrms = squeeze(sqrt(mean(AMfilt.^2,1)))*sqrt(2);%squeeze(rms(AMfilt,'dim',1));
    DC = squeeze(mean(E,1));%squeeze(rms(E,'dim',1));
    AMIspec = AMrms./(DC'*ones(1,length(mf)));%(AMrms.^2*sqrt(2))./(DC'.^2*ones(1,length(mf))); % check this line
    AMIspec = AMIspec';
    
    if nargout>3
        step.t = t;
        step.f_bw = f_bw;
        step.gamma_responses = gamma_responses;
        step.E = E;
        step.mf = mf;
        step.AMrms = AMrms';
        step.DC = DC;
    end

end
