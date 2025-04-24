function [kv, flags] = arg_varnet2017(kv, flags)
% function definput = arg_varnet2017(definput)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Defaults for all TMST functions:
% auditory frequency range
definput.keyvals.flow=70; % gammatone range (Hz), lowest frequency
definput.keyvals.fhigh=6700; % gammatone range (Hz), highest frequency

% modulation frequency range
definput.keyvals.mflow  = 0.5; % Hz, modbank_fmin
definput.keyvals.mfhigh = 200; % Hz, modbank_fmax
definput.keyvals.modbank_Nmod    = 200; % number of filters, for overalpped 
                               % filters choose 'modbank_Nmod'
definput.keyvals.modbank_Qfactor = 1;% 13; % Q factor for the filters. This is for king2019_modfilterbank. In Varnet et al. 2017 we used a Q factor of 13 for representation purposes
%definput.keyvals.Q_mfb = 1;% this is for modfilterbank

% auditory model
%definput.flags.modfilter_150Hz_LP = {'no_LP_150_Hz','LP_150_Hz'}; % modbank_LPfilter
%definput.flags.modfilter_150Hz_LP = 'no_LP_150_Hz';
%definput.flags.do_no_LP_150_Hz = 1;
definput.flags.do_LP_150_Hz = 0;
%definput.flags.modfilter_phase = {'no_phase_insens', 'phase_insens_hilbert'};
%definput.flags.modfilter_phase = 'no_phase_insens';
%definput.flags.do_no_phase_insens = 1;
definput.flags.do_phase_insens_hilbert = 0;

% f0 extraction
% Parameters for YIN
definput.keyvals.undersample = 20;%10;
definput.keyvals.ap0_thres = 0.8;
definput.keyvals.yin_thresh = 0.2;
% Parameters for f0 extraction artifact removing
definput.keyvals.maxjump = 10;
definput.keyvals.minduration = 0.08;
definput.keyvals.minf = 60;
definput.keyvals.maxf = 550;
definput.keyvals.var_thres = 1500;

%==========================================================================

kvfn = fieldnames(kv);
for i = 1:length(kvfn)
    k = kvfn{i};
    if isempty(kv.(k)) || isnan(kv.(k))
        kv.(k) = definput.keyvals.(k);
    end
end

flagsfn = fieldnames(flags);
for i = 1:length(flagsfn)
    f = flagsfn{i};
    if isempty(flags.(f))
        flags.(f) = definput.flags.(f);
    end
end
