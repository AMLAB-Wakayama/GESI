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
%       Modified:  22 Oct  2022   IT  deug  GESIparam.fs <-- GESIparam.fsSnd
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
% DirGCFB = [DirRoot '../../../GitHub_Public/gammachirp-filterbank/GCFBv234/'];  % local use only
DirGCFB = [DirRoot '../gammachirp-filterbank/GCFBv234/']; 
%exist(DirGCFB)   % for check directory
addpath(DirGCFB)
StartupGCFB;   % startup GCFB

% Sounds
DirSnd = [DirRoot 'wav_sample/'];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEDI and materials
% Parameters of dcGC filterbank
GCparam.Ctrl       = 'dynamic';
GCparam.NumCh  = 100;
GCparam.fs          = 48000;
GCparam.DynHPAF.StrPrc = 'frame';    % frame-base introduced by IT
GCparam.HLoss.Type = 'NH';  % NH
% GCparam.HLoss.Type = 'HL2';  % 80yr
% GCparam.HLoss.CompressionHealth = 0.5;  % Compression health alpha


CalibToneSPLdB = 65;
CalibToneRMSDigitalLeveldB = -26;
DigitalRms1SPLdB = CalibToneSPLdB - CalibToneRMSDigitalLeveldB;

%% Parameter settings for GESI
GESIparam.DigitalRms1SPLdB = DigitalRms1SPLdB;
GESIparam.Sigmoid = [-20, 6]; % temporal value which should be modified --- 26 May 22
GESIparam.Sim.PowerRatio = 0.6;  % power asymmetry valid for both NH + HI listeners
        %GESIparam.Sim.PowerRatio = 0.5;  % do not use: valid only for NH listeners
        %GESIparam.Sim.PowerRatio = 0.55;  % intermediate : adjustment for individual listeners
        %GESIparam.Sim.PowerRatio = 0.57;  %
%GESIparam.Sim.PowerRatio = 0.5;  % do not use: valid only for NH listeners

% controlling weight matrix introduced 29 Jul 2022
% GESIparam.Sim.SwWeightProhibit = 0; % No prohibit region : conventional weight
% GESIparam.Sim.SwWeightProhibit = 1; % set prohibit region (default)  -- The result changes only slightly

GESIparam.SwPlot = 2; %  image(Result.dIntrm.GCMFBtmc*256)

% Parameter settings for materials
SNRList = [-6, -3, 0, 3]; %SNR between clean speech and noise

rng(12345);  % simulationでの再現性確保。GCFB出力でrandnを使っているので、同じにするため。

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start simulation

for nSnd = 1:length(SNRList)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    GESIparam.fsSnd = fs;

    % GCout mat file will be kept when the name is specified.
    % These files will be used when GESI is executed again for fast processing.
    % GESIparam.NameSndRef  = NameSndRef;
    % GESIparam.NameSndTest = NameSndTest;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preparation of sound
    %%%%%%%%
    SwTimeAlign = 0; % Time alignment by apriori information
    SwTimeAlign = 1; % test for TimeAlignXcorr in GESI
    if SwTimeAlign == 0 % preparation here
        disp('-- Extraction of a speech segment in SndTest from apriori information.')
        TimeSndBefore   = 0.35;
        % SndTest = SndTest(fs*TimeSndBefore+(0:length(SndRef)-1));  not very different
        SndTest = SndTest(fs*TimeSndBefore+(1:length(SndRef)));

        % Taper window
        GESIparam.DurTaperWindow = 0.02; % 20ms taper window
        LenTaper = GESIparam.DurTaperWindow *GESIparam.fsSnd;
        Win = TaperWindow(length(SndRef),'han',LenTaper);
        SndRef   = SndRef(:).* Win(:);  % column vector
        SndTest  = SndTest(:).* Win(:);
    else
        disp('-- Extraction of a speech segment in SndTest using TimeAlignXcorr in GESI.')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Speech intelligibility prediction by GESI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [Result, GESIparam] = GESIv120(SndRef, SndTest, GCparam, GESIparam);  % v120: SndRef, SndTest
    % [Result, GESIparam] = GESIv121(SndRef, SndTest, GCparam, GESIparam);  % v121
    [Result, GESIparam] = GESIv122(SndRef, SndTest, GCparam, GESIparam);  % v122
    Metric(nSnd)    = Result.d.GESI;
    Pcorrects(nSnd) = Result.Pcorrect.GESI; % temporal value. It should be changed by the sigmoid parameters.

    disp('==========================================');
    disp(['Percent correct (temporally):' num2str(Pcorrects(nSnd)) '(%)']);
    if GESIparam.SwTimeAlign> 0
        disp(sprintf('TimeAlign : %d',GESIparam.TimeAlign.NumTimeLag))
    end
    disp('==========================================');

    %%%%%%%%%%%%%%%%%%%
    % plot figures
    %%%%%%%%%%%%%%%%%%%    
    figure(nSnd)
    subplot(1,2,1)
    image(Result.dIntrm.GESI*256);
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


%%%%%%
%% Trash / memo
%%%%%%%%
% TaperWindowが、SndTestにかかるかどうかでも、%は異なってくる。
% 

