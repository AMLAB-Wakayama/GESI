%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Time alighment of Test sound for GESI
%       Irino, T.
%       Created:     6 Jun  2022   IT  from mrGEDI GEDI_TimeAlign
%       Modified:    6 Jun  2022   IT 
%       Modified:  12 Jun  2022   IT 
%       Modified:  22 Oct  2022   IT % negative alignment
%
%
%       Alignment of signals using cross-correlation
%        --> It seems all right withing the range of SpIntel SNR 
%            in the case of additive noise
%        -->  This may be replaced by more sophisticated one
%              (eg. using GCFB representation)
%
function [SndTestOut, ParamTA] = TimeAlignXcorr(SndTestIn, SndRefIn)

disp('Time Alignmment using Xcorr.  Equalizing length(SndTest) to length(SndRef).')
if length(SndTestIn) < length(SndRefIn)  
    % To avoid zero padding in SndRef and SndTest
    error('SndTest should be longer than SndRef.')
end

LenSnd = length(SndRefIn);

[XcorrSnd, Lag] = xcorr(SndTestIn, SndRefIn);
[~,IndxMaxLag] = max(abs(XcorrSnd)); %abs is necessary
NumTimeLag = Lag(IndxMaxLag);

if NumTimeLag >= 0
    SndTest1  = [SndTestIn, zeros(1,NumTimeLag)]; % zero padding
    SndTestOut = SndTest1(NumTimeLag + (1:LenSnd));
else
    SndTest1 = [zeros(1,abs(NumTimeLag)), SndTestIn];
    SndTestOut = SndTest1(1:LenSnd);
    %error('Something strange: Negative lag. Set SndTest & SndRef properly')
end

ParamTA.XcorrSnd = XcorrSnd;
ParamTA.Lag = Lag;
ParamTA.IndxMaxLag  = IndxMaxLag;
ParamTA.NumTimeLag = NumTimeLag;

disp(['Time Lag in sample point = ' num2str(NumTimeLag)]);

end 
