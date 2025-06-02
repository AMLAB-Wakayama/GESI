%
%   Ideal observer of speech intelligibility used in GEDI (sEPSM) 
%   based on listeners' vocabulary size  see also sEPSM
%   Irino, T.
%   Created : 30 Jan 2022   IT,  from mrGEDI
%   Modified: 30 Jan 2022   IT
%
%   INPUT:  Metric :  Metric derived in the main.  SDRenv_lin when using GEDI
%               Param:  vector of [k, q, m, sigma_s]
%
%   OUTPUT:
%               Pcorrect: percent correct of speech intelligibility
%
%
%  Note in GEDI: 
%       IdealObserver: Converts the overall SDRenv to percent correct.
%
%       Usage: Pcorrect = IdealObserver(SDRenv_lin,parameters)
%       Parameters :  vector with the parameters for the ideal Observer formatted as [k q m sigma_s]
%       % d_prime = k*(SDRenv_lin).^q;
%
%       Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
%       In Signal Detection and Recognition by Human Observers,
%       edited by John A. Swets (John Wiley & Sons, New York)
%
function Pcorrect  = IdealObserver_VocabGauss(Metric,Param)

if nargin < 2
    error('You have to specify the k,q,m,sigma_s parameters for the IdealObserver')
end
k = Param(1);
q = Param(2);
m = Param(3);
sigma_s = Param(4);

% ---------- Converting from Metric (SNRenv in GEDI) to d_prime  --------------
d_prime = k*(Metric).^q;

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
Un = 1*norminv_erf(1-(1/m));
mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
Pcorrect = normcdf_erf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;  

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Un = 1*norminv(1-(1/m));   % 8 Dec 2019
% Pcorrect = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;  % 8 Dec 2019
%
% To avoide using statistics toolbox, erf function was defined
%
%    norminv_erf
%    Equivalent to norminv in Statistics Toolbox.
%    Irino, T
%    Created: 2 Dec 19
%    Modified: 2 Dec 19
%
function ni = norminv_erf(p)

ni = sqrt(2)*erfinv(2*p-1);

end

%    normcdf_erf
%    Equivalent to normcdf in Statistics Toolbox.
%    Irino, T
%    Created: 2 Dec 19
%    Modified: 2 Dec 19
%    Modified: 9 Dec 19  % debug
%
%
function nc =normcdf_erf(x,mu,sigma)

if nargin==1
    mu=0;
    sigma=1;
elseif nargin==2
    sigma=1;
end
nc = (1+erf((x-mu)./(sigma*sqrt(2))))./2;

end






















