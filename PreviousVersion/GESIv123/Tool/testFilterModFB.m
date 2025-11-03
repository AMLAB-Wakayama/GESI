%
%   test   Filtering by Modulation Filterbank (MFB)
%   Irino, T.
%   Created:  15 Oct 2022
%   Modified:  15 Oct 2022  
%
%
clear
clf

ParamMFB.fs = 2000; % sampling rate
Impulse = [0 1 zeros(1,ParamMFB.fs)];

% Envelope LPF
% ParamMFB.fcutEnv = 999; %  almost No LPF filter
% ParamMFB.fcutEnv = 512; % v110 -- v121
ParamMFB.fcutEnv = 150; % before v109 & after v122, 15 Oct 2022

[ParamMFB.bzLPF, ParamMFB.apLPF] = butter(1, ParamMFB.fcutEnv/(ParamMFB.fs/2));
ImpLP = filter(ParamMFB.bzLPF,ParamMFB.apLPF,Impulse);

[OutMFB, ParamMFB] = FilterModFB(ImpLP,ParamMFB);

[NchM, ~] = size(OutMFB);

for nch = 1:NchM
    [frsp, freq] = freqz(OutMFB(nch,:),1,pow2(ceil(log2(ParamMFB.fs))),ParamMFB.fs);
    plot(freq,20*log10(abs(frsp)))
    AmpPeakdB(nch) = max(20*log10(abs(frsp)));
    hold on
end
xlabel('Modulation Freq. (Hz)')
ylabel('Gain (dB)')
axis([0,ParamMFB.fs/2, -30,5])
grid on

plot(ParamMFB.fc,AmpPeakdB,'--')
AmpPeakdB