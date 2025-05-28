%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       test Time alighment for GESI
%       Irino, T.
%       Created:   6 Jun  2022   IT  from mrGEDI GEDI_TimeAlign
%       Modified:   6 Jun  2022   IT 
%
%
clear

[SndTest1 fs] = audioread('../wav_sample/sample_sp1.wav');
SndRef1  = audioread('../wav_sample/sample_sp_clean.wav');
rng(1234);

SndTestIn = SndTest1(:)';
if 0
    % Add additional Noise
    % 10倍のnoiseを加えても、問題なくxcorrでピークが取れる。
    % additive noiseだからか。
    RmsSndTest = rms(SndTestIn);
    SndTestIn = SndTestIn + 10*RmsSndTest*randn(size(SndTestIn));
    SndTestIn = RmsSndTest*SndTestIn/RmsSndTest; % normalize
end
SndRefIn =  SndRef1(:)';
nz = 0.1*randn(1,10);
SndTestIn = [nz -SndRefIn nz]; % 
SndTestIn = [nz SndRefIn nz]; % 

ap = audioplayer(SndTestIn,fs);
playblocking(ap);

[SndTestOut, ParamTA] = TimeAlignXcorr(SndTestIn, SndRefIn);

Error = rms(SndRefIn - SndTestOut)

figure(1);clf
plot(ParamTA.Lag,ParamTA.XcorrSnd)

figure(2);clf
bias = 0.3;
subplot(2,1,1)
plot(1:length(SndRefIn),SndRefIn,1:length(SndTestIn),SndTestIn+bias)
subplot(2,1,2)
plot(1:length(SndRefIn),SndRefIn,1:length(SndTestOut),SndTestOut+bias)




