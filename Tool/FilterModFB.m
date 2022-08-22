%
%   Filtering by Modulation Filterbank (MFB)
%   Irino, T.
%   Created:  07 Feb 2022 from modFbank_YK_v2 in mrGEDI by YamaKatsu
%   Modified: 07 Feb 2022
%   Modified: 11 Feb 2022 Separating  calculation of filter coef. and filtering for speed up
%   Modified:  4 Aug 2022  v110  introducing 512 Hz as a default
%
%   INPUT:
%           Env:  The envelope to be filtered
%           ParamMFB.fs: sampling frequency of the envelope
%           ParamMFB.fc: center frequencies of the modulation filters 
%                                 default: [1 2 4 8 16 32 64 128 256]; --> [1 2 4 8 16 32 64 128 256 512]; 
%           ParamMFB.SwPlot:  Plot frequency response of MFB
%
%   OUTPUT:
%           OutMFB:  Temporal outputs for each of the modulation filters
%           ParamMFB: Parameter
%
%  See:
%  function x_filt = modFbank_YK_v2(Env,fsEnv,cf_mod) in mrGEDI
%  simply modified some variable names
%
%
%
function [OutMFB, ParamMFB] = FilterModFB(Env,ParamMFB)
persistent MFcoefIIR

ParamMFB.fc_default =[1, 2, 4, 8, 16, 32, 64, 128, 256, 512]; % v110 4 Aug 202

if nargin < 1, 
    % help(mfilename); 
    OutMFB = [];  
    ParamMFB.fc = ParamMFB.fc_default; % just reply the default setting
    return; 
end
if isfield(ParamMFB,'fs') ==0,  error('Specify ParamMFB.fs'); end
% if isfield(ParamMFB,'fc') ==0,  ParamMFB.fc =[1 2 4 8 16 32 64 128 256]; end
if isfield(ParamMFB,'fc') ==0,  ParamMFB.fc =ParamMFB.fc_default; end % v110 4 Aug 2022
if isfield(ParamMFB,'SwPlot') ==0, ParamMFB.SwPlot = 0; end

if isfield(MFcoefIIR,'a') == 0
    MFcoefIIR = MkCoefModFilter(ParamMFB);  % Making modulation filter
end

LenFc = length(ParamMFB.fc);
[NumEnv, LenEnv] = size(Env);
if NumEnv > 1
    error('Env should be a monoaural row vector.')
end %%%%%%%%

OutMFB = zeros(LenFc,LenEnv);
for nfc = 1:LenFc
    OutMFB(nfc,:) = filter(MFcoefIIR.b(nfc,:), MFcoefIIR.a(nfc,:), Env);
end

ParamMFB.MFcoefIIR = MFcoefIIR;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making coefficients of modulation filterbank
%  The code is the same as in modFbank_YK_v2 in mrGEDI by YamaKatsu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MFcoefIIR = MkCoefModFilter(ParamMFB)

disp('--- Making modulation filter coefficients ---')
LenFc = length(ParamMFB.fc);

IIR_b = zeros(LenFc,4);
IIR_a = zeros(LenFc,4);

for nfc = 1:LenFc

    if ParamMFB.fc(nfc) == 1  % when 1 Hz
        % Third order lowpass filter
        [b, a] = butter(3, ParamMFB.fc(nfc)/(ParamMFB.fs/2));
        b4 = b/a(1);
        a4 = a/a(1);

        IIR_b(nfc,:) = b4;
        IIR_a(nfc,:) = a4;

    else % Bandpass filter
        % Pre-warping
        w0 = 2*pi*ParamMFB.fc(nfc)/ParamMFB.fs;

        % Bilinear z-transform
        W0 = tan(w0/2);

        % Second order bandpass filter
        Q = 1;
        B0 = W0/Q;
        b = [B0; 0; -B0];
        a = [1 + B0 + W0^2; 2*W0^2 - 2; 1 - B0 + W0^2];
        b3 = b/a(1);
        a3 = a/a(1);

        IIR_b(nfc,1:3) = b3;
        IIR_a(nfc,1:3) = a3;
    end

end

MFcoefIIR.a = IIR_a;
MFcoefIIR.b = IIR_b;

if ParamMFB.SwPlot == 1   % plot & pause for confirmation
    PlotFrspMF(ParamMFB,MFcoefIIR);
    disp('Return to continue > ');
    pause;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot frequency response of the digital filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotFrspMF(ParamMFB,MFcoefIIR)

hold on
Nrsl = 1024*4; % > ParamMFB.fs;

for nfc = 1:length(ParamMFB.fc)
    [frsp, freq] = freqz(MFcoefIIR.b(nfc,:),MFcoefIIR.a(nfc,:),Nrsl,ParamMFB.fs);
    plot(freq,20*log10(abs(frsp)));
end

hold off
box on
axis([0.25 max(ParamMFB.fc)*2 -40 5]);
grid;
set(gca,'xscale','log');
set(gca,'xtick',ParamMFB.fc);
xlabel('Frequency (Hz)');
ylabel('Filter attenuation (dB)');
Str_FcMFB = num2str(ParamMFB.fc');
legend(Str_FcMFB,'location','southwest');
title('Modulation filterbank');

end

