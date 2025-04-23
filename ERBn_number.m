function n = ERBn_number(f, variant)

if nargin<2
    variant = 'G&M1990';
end

switch variant
    case 'G&M1990'
        n = 21.4*log10(1+4.37e-3*f);
    case 'voicebox'
        n = 11.17268*log(1+(46.06538*f)./(f+14678.49));
    case 'M&G1983'
        n = 11.17*log((f*1e-3+0.312)./(f*1e-3+14.675))+43.;
end