%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Time alighment of GCFB output
%       Irino, T.
%       Created:   14 May 2024  IT  from TimeAlignXcorr
%       Modified:  14 May 2024   IT
%       Modified:  24 May 2024   IT  error('Specify ParamTAGC.MaxLagSec')
%
%       Input  GCoutTest,
%                 GCoutRef ParamTAGC
%                   ParamTAGC
%       Output: GCoutTestOut
%                   ParamTAGC
%
%
function [GCoutTestOut, ParamTAGC] = TimeAlignGCFB(GCoutTestIn, GCoutRefIn, ParamTAGC)

if nargin < 3
    ParamTAGC.fs = 2000;  % default value of GCFB output　= GCparamTest.DynHPAF.fs
end
if isfield(ParamTAGC,'MaxLagSec') == 0
    error('Specify ParamTAGC.MaxLagSec'); % No default here. It seems better to clarify MaxLagSec.
    % ParamTAGC.MaxLagSec = 0.010; % 10 ms seems to be good.     % 35 ms
end

ParamTAGC.MaxLag = round(ParamTAGC.MaxLagSec*ParamTAGC.fs); % it should be integer
[NumCh, LenGC] = size(GCoutTestIn);

for nch = 1:NumCh
    [XcorrGC, Lag] = xcorr(GCoutTestIn(nch,:), GCoutRefIn(nch,:), ParamTAGC.MaxLag);
    [~,IndxMaxLag] = max(abs(XcorrGC));
    NumTimeLag= Lag(IndxMaxLag);
    if NumTimeLag > 0
        GCoutTest1 = [GCoutTestIn(nch,:), zeros(1,NumTimeLag)]; % zero padding
        GCoutTestOut(nch,1:LenGC) = GCoutTest1(NumTimeLag+(1:LenGC));
    else
        GCoutTest1 = [zeros(1,abs(NumTimeLag)), GCoutTestIn(nch,:)];
        GCoutTestOut(nch,1:LenGC) = GCoutTest1(1:LenGC);
    end
    ParamTAGC.NumTimeLag(nch) = NumTimeLag;
end

return

%% Trash %%%%%

% NumTimeLag = max(NumTimeLag,0);  % possitive Time Lag only because of impulse response
% NumTimeLag = -max(-NumTimeLag,0);  % negative Time Lag only: the same as no compensation
