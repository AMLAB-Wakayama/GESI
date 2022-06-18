%
%    SSI weight function from SSI weight  --- Simplified version
%    Toshio IRINO
%    Created:    5 Feb 2022   from CalSSIweight 
%    Modified:   5 Feb 2022  
%
%   function [SSIweight, SSIparam] = F0limit2SSIweight(SSIparam),
%   INPUT:  SSIparam.Fr1 :  Filterbank channel frequency (== GCparam.Fr1==GCresp.Fr1) 
%               SSIparam.h_max = 5;
%               SSIparam.F0_limit =  F0  % spectified from adaptive F0 value
%   OUTPUT: SSIweight,  SSIparam
%
function [SSIweight, SSIparam] = F0limit2SSIweight(SSIparam)

if isfield(SSIparam, 'Fr1') == 0
    help(mfilename)
    error('SSIparam.Fr1 (== GCparam.Fr1) is essential. ');
end
if isfield(SSIparam,'F0_limit') == 0,  SSIparam.F0_limit = 150; end  % default
if isfield(SSIparam,'h_max') == 0,    SSIparam.h_max = 5; end  % default

if SSIparam.F0_limit > 0   % when F0 is possitive
    SSIparam.TI_limit = 1/SSIparam.F0_limit; % Limit of sucessive Glottal pulse
    SSIweight = min(SSIparam.Fr1*SSIparam.TI_limit, SSIparam.h_max)/SSIparam.h_max;
elseif SSIparam.F0_limit == 0  % when F0 == 0 
    SSIparam.TI_limit = Inf;
    SSIweight = ones(size(SSIparam.Fr1));
else
    error('SSIparam.F0_limit should be positive.')
end

end


%%%%%%%
%% Trash
%%%%%%
% SSIparam.ERBnum  = linspace(Freq2ERB(GCparam.FRange(1)),...
%                   Freq2ERB(GCparam.FRange(2)),GCparam.NumCh);
% SSIparam.FreqVal = ERB2Freq(SSIparam.ERBnum)';  % column vector
%DiffFr1 = rms(GCparam.Fr1 - SSIparam.FreqVal); % == 0,　Fr1を求めているだけ、、、


