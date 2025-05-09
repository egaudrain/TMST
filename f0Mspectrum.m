function [f0Mspec, mf, step] = f0Mspectrum(insig, fs, varargin)
%f0Mspectrum Summary of this function goes here
%   [f0Mspec, fc, mf, step] = f0Mspectrum(insig, fs, varargin)
% returns the f0M spectrum as defined in Varnet et al 2017

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
% 
% do_silent = 1;
% 
% [flags,kv]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

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

t=(1:length(insig))/fs;

%%% f0 extraction

% Parameters for the YIN algorithm (see help yin)
undersample = kv.undersample; % undersampling for speeding up f0M calculation
ap0_thres = kv.ap0_thres;
yin_thresh = kv.yin_thresh;
% Parameters for f0 extraction artifact removing (see help remove_artifacts_FM)
maxjump = kv.maxjump;
minduration = kv.minduration;
minf = kv.minf;
maxf = kv.maxf;
var_thres = kv.var_thres;

P=[]; P.hop = undersample; P.sr = fs; P.minf0 = minf; P.maxf0 = maxf; P.thresh = yin_thresh;

R = yin(insig(:), P);
f0 = 440*2.^(R.f0);
f0_withnan = f0; f0_withnan(R.ap0>ap0_thres) = NaN;
f0_withnan = remove_artifacts_FM( f0_withnan, fs/undersample, maxjump, minduration, [minf maxf], [0.4 2.5], var_thres, 0);
t_f0 = (1:length(f0_withnan))/(fs/undersample);

f0withoutnan = f0_withnan; f0withoutnan(isnan(f0_withnan))=[];
twithoutnan = t(1:undersample:end); twithoutnan=twithoutnan(1:length(f0_withnan)); twithoutnan(isnan(f0_withnan))=[];

%%% f0M spectra
if do_silent == 0
    fprintf('calculating f0M spectra\n');
end

[f0Mfft, f_periodo] = plomb(f0withoutnan,twithoutnan,f_spectra); %TODO change parameters
f0Mfft=2*f0Mfft; % because the output is a 2-sided periodogram
f0Mspec = interpmean( f_periodo, f0Mfft, f_spectra_intervals );
mf = f_spectra;
    
if nargout>=3
    step.t = t_f0;
    step.f0 = f0_withnan;
    step.mf = mf;
    step.mfb = f_spectra_intervals;
    step.flags = flags;
    step.kv = kv;
end

end
