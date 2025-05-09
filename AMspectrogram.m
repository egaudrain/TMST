function [AMsgram, fc, scale, step] = AMspectrogram(insig, fs, varargin)
%AMspectrogram Amplitude modulation wavelet spectrogram
%   [AMsgram, fc, scale, step] = AMspectrogram(insig, fs, varargin)
% returns the AM wavelet spectrogram of signal insig (wavelet amplitude).
% fs: sampling frequency
% AMspec is a N-by-M function where N is the number of modulation
% frequencies (mf) and M is the number of audio frequencies (fc).
% see Varnet et al. 2017 for more details.
%
% Leo Varnet - 07/2023

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

%{
definput.import={'varnet2017'}; 
definput.importdefaults={}; 

do_silent = 1;

[flags,kv]  = ltfatarghelper({'flow','fhigh'},definput,varargin);
%}
definput = arg_varnet2017(struct());
kv = definput.keyvals;
for i = 1:length(varargin)/2
    k = varargin{2*i-1};
    v = varargin{2*i};
    kv.(k) = v;
end

do_silent = 1;

% defines the modulation axis
mflow  = kv.mflow;
mfhigh = kv.mfhigh;
N_fsamples = kv.modbank_Nmod; 
f_spectra_intervals = logspace(log10(mflow), log10(mfhigh), N_fsamples+1);
f_spectra           = logspace(log10(sqrt(f_spectra_intervals(1)*f_spectra_intervals(2))), log10(sqrt(f_spectra_intervals(end)*f_spectra_intervals(end-1))), N_fsamples);

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
if do_silent == 0
    fprintf('E extraction\n');
end
E = abs(hilbert(squeeze(gamma_responses)));

Nchan = length(fc);
%%% AM spectra
if do_silent == 0
    fprintf('calculating envelope wavlets\n');
end
fb = cwtfilterbank('SignalLength',length(t),'SamplingFrequency',fs,'FrequencyLimits',[mflow mfhigh],'Wavelet','morse');%'TimeBandwidth',120,

%freqz(fb)
for ichan=1:Nchan
    [wt,scale] = cwt(E(:,ichan),'FilterBank',fb);
    if exist('AMsgram')
        AMsgram(:,:) = AMsgram + abs(wt)'/Nchan;
    else
        AMsgram(:,:) = abs(wt)'/Nchan;
    end
end

if do_silent == 0
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
