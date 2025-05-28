%
%   test   Filtering by Modulation Filterbank (MFB)
%   Irino, T.
%   Created:  15 Oct 2022
%   Modified: 15 Oct 2022  
%   Modified: 30 Dec 2023  For Introduction of TMTF
%   Modified: 18 Nov 2024  introducing MFBparam.MaxFc,  renamed ParamMFB --> MFBpram.  
%   Modified:  7 Dec 2024  introducing MFBparam.MinFc
%
%
clear
clf

MFBparam.fs = 2000; % sampling rate
Impulse = [0 1 zeros(1,MFBparam.fs)];

% Envelope LPF
% MFBparam.fcutEnv = 999; %  almost No LPF filter
% MFBparam.fcutEnv = 512; % v110 -- v121
MFBparam.fcutEnv = 150; % before v109 & after v122, 15 Oct 2022

% 事前にenvelopeにLPFをかけて、Modulationの上限を決める。
[MFBparam.bzLPF, MFBparam.apLPF] = butter(1, MFBparam.fcutEnv/(MFBparam.fs/2));
ImpRspLPF = filter(MFBparam.bzLPF,MFBparam.apLPF,Impulse);

cpz = 0.1;
apAPF = [1, cpz ]; bzAPF = [cpz ,1];
apAPF = [1]; bzAPF = [1];
ImpRspAPF = filter(bzAPF,apAPF,Impulse);

figure(1);clf
nn = 1:50;
plot(nn,ImpRspLPF(nn),nn,ImpRspAPF(nn))


%%%%%%%%%%%
% Modulation FB処理
%%%%%%%%%%%

figure(2);clf
Env = Impulse; % ImpulseをEnvelopeとしていれる。MFBの周波数特性を直接出せる。
%EnvNoLPF = Env; % LPFをかけないEnvのまま、
EnvNoLPF =rms(Env); % LPFをかけないEnvのまま、

Env = ImpRspLPF; % LPFの特性自体をEnvelopeとして入れる。全部正なので、大丈夫。

MFBparam.MaxFc = 32 % Upper Fc
MFBparam.MinFc = 2 % Upper Fc
[OutMFB, MFBparam] = FilterModFB(Env,MFBparam);

ModDepthHIdB = -18;
ModDepthNHdB = -23;
ModDepthNH     = 10^(ModDepthNHdB/20);
ModDepthRdct   = 10^((ModDepthNHdB-ModDepthHIdB)/20);
OutMFB =ModDepthRdct*(1-ModDepthNH )*OutMFB +  ModDepthNH *EnvNoLPF; %そのまま漏れをつけても、意味なし.noise floorの導入方法は？　


[NchM, LenOut] = size(OutMFB);
NfrqRsl = pow2(ceil(log2(MFBparam.fs)));
for nch =  1:NchM
    [frsp, freq] = freqz(OutMFB(nch,:),1,NfrqRsl,MFBparam.fs);
    % plot(freq,20*log10(abs(frsp)))
    % semilogx(freq,20*log10(abs(frsp+randn(NfrqRsl,1)*1e-2)))  % 係数にnoiseは乗るが、、、いまいち。
    % semilogx(freq,max(20*log10(abs(frsp)), -25+randn(NfrqRsl,1)/3) )  % -25以下はfloorになるが、入力音圧が大きくなっても変わらないのでダメ
    semilogx(freq,20*log10(abs(frsp)))
    AmpPeakdB(nch) = max(20*log10(abs(frsp)));
    hold on
end
xlabel('Modulation Freq. (Hz)')
ylabel('Gain (dB)')
axis([0,MFBparam.fs/2, -40,5])
grid on

%%%% 
% 手計算で、ModulationのLPF gainを計算
% see also two-point method
% Takashi Morimoto, Toshio Irino, Kouta Harada, Takeshi Nakaichi, Yasuhide Okamoto, Ayako Kanno, Sho Kanzaki, and Kaoru Ogawa,
% "Two-point method for measuring the temporal modulation transfer function," Ear and Hearing,40(1),pp.55--62,Jan 2019 
%
MFBparam.fc
GainLPFdB = 10*log10( 1./(1+(MFBparam.fc/MFBparam.fcutEnv).^2));

%plot(MFBparam.fc,AmpPeakdB,'--')
semilogx(MFBparam.fc,AmpPeakdB,'--',MFBparam.fc,GainLPFdB,'-.')
[AmpPeakdB; GainLPFdB]
% LPF経由と手計算の結果はほぼ一致。


% %%%%%%
% 宮﨑／花谷 NH TMTF データ
%%%%%
NameTMTFNH = 'RsltSbjTMTF_NH_MyzkHnt23.mat';
load(NameTMTFNH)
% Fcutoffの平均の計算
geomean(geomean(RsltSbjTMTF.RsltSbj(:,:,4)))
