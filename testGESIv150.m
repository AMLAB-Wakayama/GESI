%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Test code for GESI
%       Irino, T.
%       Created :  21 Mar 2022  IT from test_mrGEDIhl
%       Modified:  21 Mar 2022  IT
%       Modified:  24 Mar 2022  IT Modified Location
%       Modified:  12 May 2022  IT Modified GESIparam.Sigmoid  --> See  GESI.m
%       Modified:  26 May 2022  IT Modification on 12 May 2022 was cancelled.
%       Modified:  29 Jul  2022  IT introduced: GESIparam.SwWeightProhibit
%       Modified:   4 Aug 2022   IT  v110  introduction of version number +  normalization of SSI weight
%       Modified:  22 Aug 2022   IT  v120  The order of input arguments was replaced
%       Modified:  31 Aug 2022   IT  v121  Introduction of time-varying SSIweight
%       Modified:  18 Oct  2022   IT  v122  adding rng()
%       Modified:  19 Oct  2022   IT  v122  using GCFBv234
%       Modified:  12 Nov 2022   IT  v123  version up. Tidy up. Renamed  from GESIv122_rPwrMulti.m (YA)
%       Modified:  18 May 2023   IT  v123  adding some comments
%       Modified:    3 Jan  2024   IT  v130  Introduction of TMTF for HI listeners.
%       Modified:  27 Aug 2024    IT  v131  For simple tests
%       Modified:  19 Sep 2024    IT  v142  introducting TMTF to CosSim output
%       Modified:  10 Oct  2024   IT  v143  clarify GCparamTest/GCparamRef
%       Modified:  24 Oct  2024   IT  v144  MFBparam.TMTFfcutoff_NH = 128;
%       Modified:  26 Oct  2024   IT  v145  debug
%       Modified:  14 Nov  2024   IT  v146  introducing weightAbvThrTest in weightGCFB : Ep above threshold should be evaluated.
%       Modified:  18 Nov  2024   IT  v147  
%       Modified:  19 Nov  2024   IT  v147  modified weightAbvThrCmpnst
%       Modified:    7 Dec  2024   IT  v148  MFBparam.MinFc = 1
%       Modified:   20 Apr  2025   IT  v148  Note on parameters
%       Modified:    3 May  2025   IT  v149  
%       Modified:   15 Jul   2025   IT  v150  %   nrPwr --> nRho, for nEta = 1:length(GESIparam.Sim.weightAbvThrCmpnst)
%
%
%  Note (v148, 20 Apr 2025):
%       GESIparam.Sim.PowerRatio == "rho" in the paper
%       GESIparam.Sim.weightAbvThrCmpnst == "eta" in the paper
%       TMTF:
%           GESIparam.MFBparam.Type = 'NH'; % TMTF is not included in the analysis
%               NH default: GESIparam.MFBparam.TMTFfcutoff = 128; %  (== fc in the paper) cutoff freq 
%               NH default: GESIparam.MFBparam.TMTFminDepth = -23; % (== Lps in the paper) modulation threshold 
%           GESIparam.MFBparam.Type = 'HI';  % in the case of degraded TMTF 
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Environment settings
% Root
DirProg = fileparts(which(mfilename)); % Directory of this m-file
DirRoot = [DirProg '/'];
%
% Essential package: dynamic compressive gammachirp filterbank
% (GCFBv233 or the later version)
% Please download and put it at the same level of this directory.
% https://github.com/AMLAB-Wakayama/gammachirp-filterbank/GCFBv234
%
%
%DirGCFB = [DirRoot '../GCFBv233/'];  % normal install
DirGCFB = [DirRoot '../../../GitHub_Public/gammachirp-filterbank/GCFBv234/'];  % local use only
%exist(DirGCFB)   % for check directory
addpath(DirGCFB)
StartupGCFB;   % startup GCFB

% Sounds
DirSnd = [DirRoot 'wav_sample/'];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GCFB & GESI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of dcGC filterbank
GCparam.Ctrl       = 'dynamic';
GCparam.NumCh  = 100;
GCparam.fs          = 48000;
GCparam.DynHPAF.StrPrc = 'frame';    % frame-base introduced by IT
GCparam.HLoss.Type = 'NH';  % NH
% Modify the Hearing loss type here  --- See usage of GCFBv234
% GCparam.HLoss.Type = 'HL2';  % 80yr preset
% GCparam.HLoss.CompressionHealth = 0.5;  % Compression health alpha
    % GCparam.HLoss.CompressionHealth = 0;  % Compression health alpha
    % GCparam.HLoss.Type = 'HL3';  % 70yr preset

% Correspondence between digital level and SPL dB
CalibToneSPLdB = 65;                          % used in the SPL calibration (1 kHz, sin wave)
CalibToneRMSDigitalLeveldB = -26;       % Digital RMS level of the calibration sound
% Parameter settings for GESI  --- Digital RMS 1 == ?? dB SPL
GESIparam.DigitalRms1SPLdB = CalibToneSPLdB - CalibToneRMSDigitalLeveldB;

% Sigmoid parameter [a, b] (temporal setting)
GESIparam.Sigmoid = [-14, 6]; % temporal sigmoid parameter values for this demo (v148)
    % GESIparam.Sigmoid = [-20, 6]; % temporal value which should be modified 
    % GESIparam.Sigmoid = [-25, 6]; % temporal value for v145

% parameter "rho"
GESIparam.Sim.PowerRatio = 0.55;  % (== "rho" in the paper.) power asymmetry valid for both NH + HI listeners
    % GESIparam.Sim.PowerRatio = 0.5;   % DO NOT USE. Valid only for NH  (== "rho" in the paper.)
    % GESIparam.Sim.PowerRatio = 0.7;  % (== "rho" in the paper.)

% parameter "eta"
GESIparam.Sim.weightAbvThrCmpnst = 0.7; % (== "eta" in the paper)  param for weightGCFB
    % GESIparam.Sim.weightAbvThrCmpnst = 0.5; % (== "eta" in the paper)  
    % GESIparam.Sim.weightAbvThrCmpnst = 1; % (== "eta" in the paper)

% Plot patterns to give useful information
GESIparam.SwPlot = 2; %  image(Result.dIntrm.GCMFBtmc*256)

% Introduction of TMTF in MFB
GESIparam.MFBparam.Type = 'NH'; % TMTF is not included 
% default NH value: GESIparam.MFBparam.TMTFfcutoff = 128; %  (== fc in the paper) cutoff freq 
%                           GESIparam.MFBparam.TMTFminDepth = -23; % (== Lps in the paper) modulation threshold 
%
if isfield(GESIparam,'MFBparam') == 0  % if NH is not specified, example TMTF of HI is used
    GESIparam.MFBparam.Type = 'HI';  % v130 and later
    % Example 
    GESIparam.MFBparam.TMTFfcutoff = 128/2;   % fc cutoff freq lower than NH (128) 
    GESIparam.MFBparam.TMTFminDepth = -10;  % Lps (modulation threshold)  higher than NH (-23) 
end

% MFB setteing: default
% GESIparam.Sim.weightMFB = ones(10,1); 
% GESIparam.Sim.weightMFB = [ones(6,1); 0; 0; 0; 0];   % 1 up to 32 Hz; 0 above it
% GESIparam.MFBparam.MaxFc = 32; % default

% Parameter settings for example speech materials
SNRList = [-6, -3, 0, 3]; %SNR between clean speech and noise

rng(12345);  % To ensure reproducibility in the simulation, although the GCFB output uses the 'randn'


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Start simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nSnd =  1:length(SNRList)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test signal (enhanced/Unprocessed speech); the speech intelligiblity is calculated
    % Name of wav-file (example: '*/GEDI_Software/wav_sample/sample_sp1')
    NameSndTest = ['sample_sp' num2str(nSnd)];
    disp(['SndTest: ' NameSndTest]);
    % Read wav-file of test speech
    [SndTest, fs] = audioread([DirSnd NameSndTest '.wav']);

    %% Reference signal (Clean speech)
    % Name of wav-file
    NameSndRef = 'sample_sp_clean';
    disp(['SndRef : ' NameSndRef]);
    % Read wav-file of clean speech
    [SndRef, fs2] = audioread([DirSnd NameSndRef '.wav']);

   if fs ~= fs2  %  IT
        error('Inconsistency of sampling rate.');
    end
    GESIparam.fs = fs;   % Samping rate of sounds.   % NG:  GESIparam.fsSnd = fs;

    % The GCout mat file is retained if the name is given.
    % These files are used for fast processing when GESI is run repeatedly.
    % GESIparam.NameSndRef  = NameSndRef;
    % GESIparam.NameSndTest = NameSndTest;

    %%%%%%%%%%%%%%%%%%%%%%
    %% Speech intelligibility prediction by GESI
    %%%%%%%%%%%%%%%%%%%%%%
    [Result, GESIparam] = GESIv150(SndRef, SndTest, GCparam, GESIparam); % v148

    Metric(nSnd)     = Result.d.GESI;
    Pcorrects(nSnd) = Result.Pcorrect.GESI; % temporal value. It should be changed by the sigmoid parameters.
    ValCosSim(nSnd) =  Result.GESI; % CosSim values

    disp('==========================================');
    disp(['Percent correct (temporally):' num2str(Pcorrects(nSnd)) '(%)']);
    if GESIparam.SwTimeAlign> 0
        disp(sprintf('TimeAlign_SndLag : %d',GESIparam.TimeAlign.NumSndLag))
    end
    disp('==========================================');

    % RsltSSI =  [Result.d.GESI, Result.Pcorrect.GESI, Result.ModIndex.Ref_SSI, Result.ModIndex.Test_SSI]
    % RsltInvSSI =  [Result.d.GESI_InvSSI, NaN, Result.ModIndex.Ref_InvSSI, Result.ModIndex.Test_InvSSRI]


    %%%%%%%%%%%%%%%%%%%
    % plot figures
    %%%%%%%%%%%%%%%%%%%
    figure(nSnd)
    subplot(1,2,1)
    image(Result.GESI.CosSim*256);
    set(gca,'YDir','normal');
    xlabel('MFB channel')
    ylabel('GCFB channel')
    title(['Metric: ' num2str(Metric(nSnd),'%5.3f') ...
        ',  Pcorrect(tmp) : ' num2str(Pcorrects(nSnd),'%4.1f')])
    drawnow

    subplot(1,2,2)
    image(GESIparam.SSIparam.weight*256*0.8);
    set(gca,'YDir','normal');
    xlabel('Frame')
    ylabel('GCFB channel')
    title(['SSIweight'])
    drawnow


end

disp(['Pcorrect : ' num2str(Pcorrects)])
disp(['Metric    : ' num2str(Metric)])

disp('==========================================');


%% Plot results
figure(nSnd+1)
plot(SNRList,Pcorrects,'o-');
xlim([-0.5+min(SNRList) max(SNRList)+0.5]);
ylim([-0.5+0 100+0.5]);
xlabel('SNR (dB)');
ylabel('Percent correct (%)')
legend('Unprocessed')
title('Results of GESI for example sounds');
grid on;


% Keep results for comparison   10 Jan 22
if 0
    NameRslt = 'Rslt_GESI';
    save([NameRslt '_Val'])
    print([NameRslt '_Fig'],'-depsc')
end
%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trash / memo
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% controlling weight matrix introduced 29 Jul 2022
% GESIparam.Sim.SwWeightProhibit = 0; % No prohibit region : conventional weight
% GESIparam.Sim.SwWeightProhibit = 1; % set prohibit region (default)  -- The result changes only slightly
%
%

% [Result, GESIparam] = GESIv120(SndRef, SndTest, GCparam, GESIparam);  % v120: SndRef, SndTest
% [Result, GESIparam] = GESIv123(SndRef, SndTest, GCparam, GESIparam);  % v123
% [Result, GESIparam] = GESIv131(SndRef, SndTest, GCparam, GESIparam);  % v131
% [Result, GESIparam] = GESIv142(SndRef, SndTest, GCparam, GESIparam); % v142
% [Result, GESIparam] = GESIv144(SndRef, SndTest, GCparam, GESIparam); % v144
%[Result, GESIparam] = GESIv145(SndRef, SndTest, GCparam, GESIparam); % v145
% [Result, GESIparam] = GESIv146(SndRef, SndTest, GCparam, GESIparam); % v146
% [Result, GESIparam] = GESIv147(SndRef, SndTest, GCparam, GESIparam); % v147


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparation of sound
% Time alignment of sounds are performed in GESI at least v130 and later
%%%%%%%%%%%%

%     SwTimeAlign_SndLag = 0; % Time alignment by apriori information
% %SwTimeAlign_SndLag = 1; % Using TimeAlignXcorr in GESI
% if SwTimeAlign_SndLag == 0 % preparation here
%     disp('-- Extraction of a speech segment in SndTest from apriori information.')
%     TimeSndBefore   = 0.35;
%     % SndTest = SndTest(fs*TimeSndBefore+(0:length(SndRef)-1));  not very different
%     SndTest = SndTest(fs*TimeSndBefore+(1:length(SndRef)));
%
%     % Taper window
%     GESIparam.DurTaperWindow = 0.02; % 20ms taper window
%     LenTaper = GESIparam.DurTaperWindow *GESIparam.fs;
%     Win = TaperWindow(length(SndRef),'han',LenTaper);
%     SndRef   = SndRef(:).* Win(:);  % column vector
%     SndTest  = SndTest(:).* Win(:);
%     GESIparam.SwTimeAlign = 1;  -- Confirmed to get the same result 20 Apr 2025
% else
%     disp('-- Extraction of a speech segment in SndTest using TimeAlignXcorr in GESI.')
% end

%%%% 26 Oct 24 check
% GESIparam.TimeAlign.GCFBmaxLagSec = 0;
% GESIparam.Sim.PowerRatio = 0.5;
% SndTest = SndRef; % the case of the same signal
