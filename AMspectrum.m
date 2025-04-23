function [AMspec, fc, mf, step] = AMspectrum(insig, fs, varargin)
%AMspectrum Amplitude modulation spectrum
%   [AMspec, fc, mf, step] = AMspectrum(insig, fs, varargin)
% returns the AMa spectrum of signal insig (PSD measured in W/Hz unit).
% fs: sampling frequency
% AMspec is a N-by-M function where N is the number of modulation
% frequencies (mf) and M is the number of audio frequencies (fc).
% see Varnet et al. 2017 for more details
%
% Leo Varnet - 07/2023
% Etienne Gaudrain - 0

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end

% definput.import={'varnet2017'}; 
% definput.importdefaults={}; 
definput = arg_varnet2017(struct());
kv = definput.keyvals;
for i = 1:length(varargin)/2
    k = varargin{2*i-1};
    v = varargin{2*i};
    kv.(k) = v;
end

do_silent = 1;

% [flags,kv]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

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
nERBs = ceil(diff(ERBn_number([kv.flow, kv.fhigh])));
gammaFiltBank = gammatoneFilterBank([kv.flow, kv.fhigh], nERBs, 'SampleRate', fs);
gamma_responses = gammaFiltBank(insig);
fc = gammaFiltBank.getCenterFrequencies();
f_bw = gammaFiltBank.getBandwidths();

%[gamma_responses,fc] = auditoryfilterbank(insig,fs,kv.flow,kv.fhigh);
%f_bw = audfiltbw(fc);

%%% AM extraction
if do_silent == 0
    fprintf('E extraction\n');
end
E = abs(hilbert(squeeze(gamma_responses)));

Nchan = length(fc);
%%% AM spectra
if do_silent == 0
    fprintf('calculating envelope spectra\n');
end

% We are using the FFT method which is quite a bit faster
sp_meth = 2;

%sp_meth = 2; % 1: periodogram, 2: fft
if sp_meth==2
    nfft = 2^nextpow2(size(E,1));
    df = fs/nfft;

    tdf = min(f_spectra(1), f_spectra(2)-f_spectra(1));
    
    if df > tdf
        nfft = 2^nextpow2(fs / tdf);
        %df = fs/nfft;
    end
end

for ichan=1:Nchan
    %[Efft, mf] = periodogram(E(:,ichan),[],f_spectra,fs,'psd');
    if sp_meth==1
        [Efft, mf] = periodogram(E(:,ichan),[],f_spectra,fs,'psd');
    elseif sp_meth==2
        X = fft(E(:,ichan), nfft);
        X2  = abs(X).^2/fs/size(E,1);
        Xf = linspace(0, fs, nfft);
        Efft = interp1(Xf, X2, f_spectra);
        mf = f_spectra;
    end
    %[Efft, mf] = plomb(E(:,ichan),t,f_spectra);
    
    %Efft=2*Efft; % because the output is a 2-sided periodogram
    AMspec(:,ichan) = 2*Efft; % PSD : power = mean(E(:,ichan).^2) = sum(Efft.*diff(f_spectra_intervals))
end

if nargout>3
    step.t = t;
    step.f_bw = f_bw;
    step.gamma_responses = gamma_responses;
    step.E = E;
    step.mf = mf;
    step.mfb = f_spectra_intervals;
end

end
