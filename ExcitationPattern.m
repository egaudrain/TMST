function [EP,fc] = ExcitationPattern(insig, fs, varargin)
%EXCITATIONPATTERN Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end

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

% Number of steps in fractional octaves to go from mod_flow to mod_fhigh:
% N_octave_steps = ceil(NthOct * log10(mod_fhigh/mod_flow)/log10(2));
% mfc            = mod_flow * 2.^((0:N_octave_steps)/NthOct); % includes one extra bin (because of the ceiling)
% cutoff_oct     = mfc*2^(-.5/NthOct); % half step down
%%%

t=(1:length(insig))/fs;

%%% gammatone filtering

% [gamma_responses,fc] = auditoryfilterbank(insig,fs,kv.flow,kv.fhigh,'bwmul',0.1);
% f_bw = audfiltbw(fc);
nERBs = ceil(diff(ERBn_number([kv.flow, kv.fhigh])));
gammaFiltBank = gammatoneFilterBank([kv.flow, kv.fhigh], nERBs*10, 'SampleRate', fs);
gamma_responses = gammaFiltBank(insig);
fc = gammaFiltBank.getCenterFrequencies();

EP = sqrt(squeeze(mean(gamma_responses.^2,2)));

end

