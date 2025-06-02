%
%   Metric to Pcorrect
%   Irino, T.
%   Created: 30 Jan 2022 from STOI2PCorrect by KY
%   Modified: 30 Jan 2022   IT
%   Modified:  2 Feb 2022   IT % corresponding to HASPI
%   Modified: 26 May 2022  IT % Add control to Max percent correct (MapPcorrect)
%
%    INPUT: Meric:   
%           ParamSigmoid:
%             Case STOI, ESTOI, GECI objective measures
%                          ParamSigmoid: 
%                          length(ParamSigmoid) == 2 --- param for sigmoid [STOIparam.a, STOIparam.b]
%                          length(ParamSigmoid) == 3 ---  [STOIparam.a,STOIparam.b,MaxPcorrect]
%             Case HASPI or multiple metrix
%             MaxPcorrect (option) was added to the end of the parameters 23 May 2022
%                  
%    OUTPUT: Pcorrect (%)
%
% Note: 
%  STOI type
%       ParamSigmoid = [STOIparam.a, STOIparam.b, MaxPcorrect]
%       a = STOIparam.a;
%       b = STOIparam.b;
%       a = -13.1903; %Taal et al., ICASSP Proc., 2011
%       b = 6.5293;
% HASPI_v1
%       arg=bias + wgtcep*CepCorr + sum(wgtcov.*cov3');
%       Intel=1.0/(1.0 + exp(-arg)); %Logistic (logsig) function
%
function Pcorrect = Metric2Pcorrect_Sigmoid(Metric,ParamSigmoid)

LenMetric = length(Metric);
LenPS = length(ParamSigmoid);

if LenPS == LenMetric+1 % default
    MaxPcorrect = 100; 
elseif LenPS == LenMetric+2  % Setting MaxPcorrect
    MaxPcorrect = ParamSigmoid(LenPS);
    ParamSigmoid = ParamSigmoid(1:(LenPS-1));
else
    error('Lengths of Metric and ParamSigmoid are inconsistent.')
end

if LenMetric == 1  % STOI, ESTOI type
    Pcorrect = MaxPcorrect  / (1+exp(ParamSigmoid(1)*Metric+ParamSigmoid(2)));
else % HASPI type
    ValArg = ParamSigmoid(:)' * [1; Metric(:)];
    Pcorrect = MaxPcorrect/(1+exp(-ValArg));
end

end


%% %%%%%%
% Trash
%%%%%%%%
% 
% if length(ParamSigmoid)  == 2
%     if length(Metric) > 1
%         error('STOI type sigmoid: Something wrong : length(Metric) > 1'); 
%     end
%     Pcorrect = MaxPcorrect  / (1+exp(ParamSigmoid(1)*Metric+ParamSigmoid(2)));
% else
%     if length(Metric) ~= (length(ParamSigmoid)-1)
%         error('HASPI type sigmoid: Something wrong : length(Metric) is not consistent'); 
%     end
%     ValArg = ParamSigmoid(:)' * [1; Metric(:)];
%     Pcorrect = MaxPcorrect/(1+exp(-ValArg));
% end


      