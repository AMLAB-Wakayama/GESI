%
%       Gammachirp Envelope Similarity Index (GESI)
%       Objective mesure index for speech intelligibility including hearing loss
%       Irino, T., Yamamoto, A.
%       Created : 30 Jan 2022   IT, New 1st version
%       Modified: 30 Jan 2022   IT  (v100)
%       Modified:   2 Feb 2022   IT
%       Modified:   7 Feb 2022   IT introducing ModFB
%       Modified:   8 Feb 2022   IT % keeping GCout in GESIparam.DirGCout
%       Modified:   9 Feb 2022   IT % introducing SSIweight(F0mean)
%       Modified: 19 Feb 2022   IT % Tidy up
%       Modified:   7 Mar 2022   IT using new Eqlz2MeddisHCLevel with GESIparam.DigitalRms1SPLdB, GCFBv232
%       Modified:   9 Mar 2022   IT GECI-> GESI precise naming:  cosine-similarity --  not Correlation
%       Modified:  20 Mar 2022  IT introduction of GCFBv233 --- Interspeech2022 version
%       Modified:  12 May 2022  IT  sum(weightMFB) == 1
%       Modified:  26 May 2022  IT  Quit "sum(weightMFB) == 1"  because we used mean value at the last stage. It becomes similar value with STOI.
%       Modified:   7 Jun  2022  IT  (v108)  Introduced TimeAlignXcorr　+  Taperwindow to SndRef/SndTest
%       Modified:  29 Jul  2022  IT  (v109)  introduced GESIparam.SwWeightProhibit
%       Modified:   4 Aug 2022   IT  v110  introduction of version number +  normalization of SSI weight
%       Modified:  22 Aug 2022   IT  v120  The order of input arguments (SndTest, SndRef, ... )
%                                          was replaced to (SndRef, SndTest, ... ) as the same as in stoi,estoi,& haspi
%       Modified:  31 Aug 2022   IT  v121  Introduction of time-varying SSIweight
%       Modified:  15 Oct  2022   IT  v122  return back to MFBparam.fcutEnv = 150; Error < 1 %
%       Modified:  19 Oct  2022   IT  v122  using GCFBv234
%       Modified:  22 Oct  2022   IT  deug  GESIparam.fs <--  GESIparam.fsSnd
%       Modified:  23 Oct  2022  YA, Calculation of multiple rho (rPwr)  -- OK IT 12 Nov 22
%       Modified:  12 Nov 2022   IT  v123  version up. Tidy up. Renamed  from GESIv122_rPwrMulti.m (YA)
%       Modified:  18 May 2023   IT  v123  adding some comments
%       Modified:   5 Jan  2024   IT  v130  Introduction of TMTF for HI listeners. Removed Japanese comments in the function.
%       Modified:  24 Apr  2024   IT  v130  Added some comments
%       Modified:   5  May  2024   IT  v131 Introducing SndInfo -  TimeAlign
%       Modified:  14 May  2024   IT  v140  TimeAlign at GCFB output from the idea of the SWMT
%       Modified:  24 May  2024   IT  v141  Keep tFrame
%       Modified:  19 Sep  2024   IT  v142  TMTF:  GESIparam.Sim.weightMFB = MFBparam.TMTFweight_NH;
%       Modified:  21 Oct  2024   IT  v143  clarify GCparamTest/GCparamRef
%       Modified:  24 Oct  2024   IT  v144  MFBparam.TMTFfcutoff_NH =128;      morimoto on slack Jan 2024
%       Modified:  26 Oct  2024   IT  v145  debugged line 562, default GCFBmaxLagSec = 0.030;
%       Modified:  14 Nov  2024   IT  v146  introducing weightAbvThrTest in weightGCFB : Ep above threshold should be evaluated.
%       Modified:  18 Nov  2024   IT  v147  Maximum MFB fc default 32 Hz,  Simplified CosSim, introducing weightAbvThrCmpnst
%       Modified:  19 Nov  2024   IT  v147  modified weightAbvThrCmpnst (default 0.5)
%       Modified:  23 Nov  2024   IT  v147  default weightAbvThrCmpnst==1
%       Modified:   7  Dec  2024   IT  v148  introduced MFBparam.MinFc   default 1 Hz
%       Modified:   3 May  2025   IT  v149  introduced LenAbvThr to avoid division by zero
%       Modified:  15 Jul   2025   IT  v150  %   nrPwr --> nRho, for nEta = 1:length(GESIparam.Sim.weightAbvThrCmpnst)
%       Modified:    3 Nov 2025   IT  v150  %   Reference added
%
%
%   Inputs:
%       SndRef:   input signal of speech reference  (always analyzed by NH)
%                           Clean speech or other sounds when comparing NH and HI listeners.
%       SndTest:  input signal of enhanced/unprocessed noisy speech  (NH or HI)
%       GCparam:  parameters of dynamic compressive gammachirp filterbank (dcGC-FB)
%                        GCparam.HLoss: Hearing loss settings
%       GESIparam:   parameter of GESI
%                       GESIparam.fs: Sampling frequency of sounds
%                       GESIparam.DigitalRms1SPLdB: SPL of rms(s(t)) == 1
%                            % Not use -->  GESIparam.SPL: sound pressure level of "SndRef" in dB
%                       GESIparam.NameSndTest: Name of SndTest
%                       GESIparam.NameSndRef: Name of SndRef
%                       GESIparam.MFB: Modulation filterbank parameters
%                       GESIparam.Sim: Similarity parameters
%                       GESIparam.Sigmoid:  (100/(1+exp(ax+b) as in  STOI/ESTOI)
%
%   Outputs:
%       Result.
%           d:  Metric  = mean(mean(dIntrm))
%           dIntrm:  Intermediate Metric
%           Pcorrect: percent correct of speech intelligibility calculated  from d
%           ModIndex: Modulation Index
%           ModIndexIntrm: Intermediate Modulation Index
%       GESIparam:
%           Parameters for GESI
%       SndInfo:
%           SndRef, SndTest (After Alignment), SndTestOrig (Before Alignment), TimeAlign
%
%  ------
%   Note:  GCFBv233 or the later version is required.  -- see test program  for usage
%
%   Note:  Result.d may fluctuate very slightly in every simulation
%             because of the random noise floor at the output of the GCFB
%
%   Note: 22 Aug 2022.  Important change  -- Please be careful when using the previous version
%     The order of the input arguments has been changed to the same as in STOI,ESTOI, & HASPI.
%     v110 and earlier:  [Result, GESIparam] = GESI(SndTest, SndRef, GCparam, GESIparam)
%     v120 and later:     [Result, GESIparam] = GESI(SndRef, SndTest, GCparam, GESIparam)
%
%   Note: 5 Jan 2024
%      Introduction of TMTF for NH & HI
%      Implimentation of TMTF was changed. Not using butterworth LPF (apLPF,bzLPF)
%
%
%   Reference:
%       Ayako Yamamoto, Fuki Miyazaki, and Toshio Irino , "Predicting speech intelligibility 
%       in older adults for speech enhancement using the Gammachirp Envelope Similarity 
%       Index, GESI," Speech Communication, 175, 103318, Nov. 2025. 
%       [DOI: 10.1016/j.specom.2025.103318] . [arXiv preprint, [arXiv: DOI 10.48550/arXiv.2504.14437] )
%
%
function [Result, GESIparam,SndInfo] = GESIv150(SndRef, SndTest, GCparam, GESIparam)

[DirProg, NameProg] = fileparts(which(mfilename)); % Directory of this program
addpath([DirProg '/Tool/']);     % Path to Tool
addpath([DirProg '/Tool/world-0.2.4_matlab/']);     % Path to WORLD
Result.Name = NameProg;
disp(['##### Start: ' Result.Name ' #####'])


% making directory and data for keeping GCout
if isfield(GESIparam,'DirGCout') == 0
    GESIparam.DirGCout = [getenv('HOME') '/Data/GESI/GCout/']; %default directory
    if exist(GESIparam.DirGCout) == 0
        mkdir(GESIparam.DirGCout);
    end
end

if isfield(GESIparam,'SwPlot') == 0
    GESIparam.SwPlot = 0;
end

StrSPLcond = ['_Rms1SPL' int2str(GESIparam.DigitalRms1SPLdB) 'dB'];
StrHLossCond = [GCparam.HLoss.Type '_'];
if strcmp(GCparam.HLoss.Type,'NH') == 0
    if isfield(GCparam.HLoss,'CompressionHealth') == 0
        error('GCparam.HLoss.CompressionHealth should be specified.')
    end
    StrHLossCond = [GCparam.HLoss.Type '_' int2str(GCparam.HLoss.CompressionHealth*100) '_'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters of GCFB & GESI
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters of dcGC filterbank
if isfield(GCparam,'fs')  == 0,         GCparam.fs = 48000; end % GCFB should always work at 48kHz
% Recommended to use fs = 48000 Hz independent of fs of Snd
if isfield(GCparam,'NumCh')  == 0, GCparam.NumCh = 100; end
if isfield(GCparam,'FRange')  == 0, GCparam.FRange = [100, 8000]; end % covering audiogram
if isfield(GCparam,'OutMidCrct')  == 0,  GCparam.OutMidCrct = 'FreeField'; end % default FreeField  not ELC
if isfield(GCparam, 'HLoss') == 0,   GCparam.HLoss.Type = 'NH';  end  % NH or HL
GCparam.Ctrl = 'dynamic';  % mandatory
GCparam.DynHPAF.StrPrc = 'frame-base'; % mandatory
GCparam.StrFloor = 'NoiseFloor'; % rms=1 randn
GCparam.FloorLevel = 1; % it should be 1 (rms) when using 'NoiseFloor';


%%% GESI
if isfield(GESIparam,'fs') == 0
    warning(['GESIparam.fs was not specified.']);
    disp(['Do you use fs = ' num2str(GCparam.fs) 'Hz ?  [Return to yes] > ']);
    pause;
    GESIparam.fs = GCparam.fs;
end

% GCoutRef & GCoutTest are saved for speed up.
if isfield(GESIparam,'NameSndRef') == 1
    GESIparam.NameGCoutRef = ['GCref_' StrHLossCond GESIparam.NameSndRef StrSPLcond ];
    SwSave = 1;
else
    GESIparam.NameGCoutRef = '';
    SwSave = 0;
end
if isfield(GESIparam,'NameSndTest') == 1
    GESIparam.NameGCoutTest = ['GCtest_' StrHLossCond GESIparam.NameSndTest StrSPLcond ];
    SwSave = 1;
else
    GESIparam.NameGCoutTest = '';
    SwSave = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sound check
%%%%%%%%%%%%%%%%%%%%%%%%%%
[mm1, nn1] = size(SndTest);
[mm2, nn2] = size(SndRef);
if mm1 > 1 && nn1 > 1 || mm2 > 1 && nn2 > 1
    error('SndTest and SndRef should be monaural.')
end
SndRef  = SndRef(:)';  % row vector for GCFB
SndTest = SndTest(:)';


%%%%%%%%%%%%%%%%
%% Time alignment  of input sound   7 Jun 22
%%%%%%%%%%%%%%%%%
if  length(SndTest) == length(SndRef) && isfield(GESIparam,'SwTimeAlign') == 0
    disp('GESI: No time alignment:  length(SndRef) == length(SndTest)')
    GESIparam.SwTimeAlign = 0; % No alignment
    SndTestOrig = [];
    % SndTest&SndRef are assumed to be prepared properly (e.g. using TaperWindow)
else
    if isfield(GESIparam,'SwTimeAlign') == 0 % controlable from main
        GESIparam.SwTimeAlign = 1;  % default: using Xcorr
    end

    if GESIparam.SwTimeAlign == 1
        disp('GESI: Time alignment of SndRef and SndTest');
        SndTestOrig = SndTest; % Keep original, 5 May 2024;
        [SndTest, ParamTA] = TimeAlignXcorr(SndTestOrig, SndRef); % Time alignment + length eqaulization
        SndInfo.TimeAlign = ParamTA;
        GESIparam.TimeAlign.NumSndLag =  ParamTA.NumTimeLag;
    else
        error('--- Not prepared yet: Another TimeAlign algorithm. ---')
    end
    % Taper here for stable estimation
    if isfield(GESIparam,'DurTaperWindow') ==0
        SndInfo.DurTaperWindow = 0.02; % 20ms taper window
        GESIparam.DurTaperWindow = SndInfo.DurTaperWindow; % for backward compatibility % 5 May 2024
    else
        SndInfo.DurTaperWindow = GESIparam.DurTaperWindow;
    end
    LenTaper = SndInfo.DurTaperWindow *GESIparam.fs;
    Win = TaperWindow(length(SndRef),'han',LenTaper);
    SndRef   = SndRef .* Win;  % row vector
    SndTest  = SndTest .* Win;
end
SndInfo.Ref = SndRef; % SndInfo added: 5 May 2024
SndInfo.Test = SndTest;
SndInfo.TestOrig = SndTestOrig;


%%%%%%%%%%%%%%%%
%% Time alignment at the GCFB output   24 May 24
%%%%%%%%%%%%%%%%%
if isfield(GESIparam,'TimeAlign') == 0
    error('You need time align for sound: GESIparam.SwTimeAlign should be 1')
end
if isfield(GESIparam.TimeAlign,'GCFBmaxLagSec') == 0 % controlable from main
    % Default for TimeAlignGCFB (see around line 360)
    %GESIparam.TimeAlign.GCFBmaxLagSec = 0.010; % up to v144
    GESIparam.TimeAlign.GCFBmaxLagSec = 0.030; % from v145
    % You can set GESIparam.TimeAlign.GCFBmaxLagSec = 0 to disable TimeAlignGCFB
end


%%%%%%%%%%%%%%%%%%%%%%
%% GESIparam.MFBparam
%%%%%%%%%%%%%%%%%%%%%%

if isfield(GESIparam,'MFBparam')== 0
    MFBparam.Type = 'NH'; % default : Compatibility with v123 and earlier
else
    MFBparam = GESIparam.MFBparam;
    if isfield(MFBparam,'Type')== 0
        error('Specify GESIparam.MFBparam.Type : ''NH'' or ''HI'''); % Specify type
    end
end

% default value for average NH listeners
% MFBparam.TMTFfcutoff_NH = 150; % (Hz) TMTF cutoff freq. of NH  (default value) (Previous name: MFBparam.fcutEnv)
% MFBparam.TMTFminDepth_NH = -23; % (dB)  TMTF minimum detection depth in dB of NH  (default value)
%
% slack wg_temporalresponse,  2023/12/30, 2024/1/6：
% No clear evidence of 150Hz.  TMTFminDepth is also unclear
% Morimoto (2018) Two-Point Method for Measuring the Temporal Modulation Transfer Function, Ear & Hearing, 40(1), pp. 55–62
% Results Obtained With the Conventional Method
% The average values of Lps_C and fcutoff_C were −22.6 dB and 127.3 Hz, respectively,
% when using the conventional method for　NH subjects. （26 NH participants
% These values are approximately the same as the values reported by Shen and Richards (2013),
% that is, −23.9 dB and 140.5 Hz for four NH subjects. (4 NH participants)
%
% We concluded to use -23 dB and 128 Hz for the NH TMTF parameters.
% Modified 24 Oct 2024　GESIv144
MFBparam.TMTFfcutoff_NH = 128; % (Hz) TMTF cutoff freq. of NH  (default value)
MFBparam.TMTFminDepth_NH = -23; % (dB)  TMTF minimum detection depth in dB of NH  (default value)

if strcmp(MFBparam.Type,'NH') == 1  % When NH, default setting
    MFBparam.TMTFfcutoff = MFBparam.TMTFfcutoff_NH;
    MFBparam.TMTFminDepth = MFBparam.TMTFminDepth_NH;
end

% In the case of 'HI', you need set these parameters.
if isfield(MFBparam,'TMTFfcutoff') == 0
    error('Specify GESIparam.MFBparam.TMTFfcutoff -- for HI');
end
if isfield(MFBparam,'TMTFminDepth') == 0
    error('Specify GESIparam.MFBparam.TMTFminDepth -- for HI');
end

% default MaxFc/MinFc for MFB filter
if isfield(MFBparam,'MaxFc') == 0    % 18 Nov 2024
    MFBparam.MaxFc = 32;  % 32 Hz
end
if isfield(MFBparam,'MinFc') == 0    % 7 Dec 2024  v148
    MFBparam.MinFc = 1;  % 1 Hz
end


%%%%%%%%%%%%%%%%%%%%%%
%% GESIparam.Sim
%%%%%%%%%%%%%%%%%%%%%%
% FilterModFBは、1Hz ~ MFBparam.MaxFcまですべて作る　-- 複雑になってしまうので。
[~, MFBparam] = FilterModFB([],MFBparam); % load the default settings of MFB -- always use this

LenMFB = length(MFBparam.fc);

% from v142
% If no control from the main routine in v142 (19 Sep 2024), then
%  GESIparam.Sim.weightMFB = MFBparam.TMTFweight_NH
% See line 408
%
% Before v141
% if isfield(GESIparam.Sim,'weightMFB') == 0  % To control from the main routine % up to v141
%     LenMFB = length(MFBparam.fc); % length of MFB fc see FilterModFB.m
%     GESIparam.Sim.weightMFB = [ones(LenMFB,1)];
%     % GESIparam.Sim.weightMFB = [0; 0; 0; ones(4,1); 0; 0; 0]; % --- [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
% else
%     LenMFB = length(GESIparam.Sim.weightMFB); % When the parameter is entered
% end

if isfield(GESIparam,'Sim') == 0 ||  isfield(GESIparam.Sim,'PowerRatio') ==0  % To control from the main
        % GESIparam.Sim.PowerRatio = 0.6; % 6:4
        GESIparam.Sim.PowerRatio = 0.55; % rho :   default 18 Nov 24
        disp('GESIparam.Sim.PowerRatio (rho) is set to 0.55 (default) -- OK? Return to continue > ')
        pause
        % Normalizing power of Ref / Test ーーー　rPwr == pho in Paper
        % rPwr = 1; % Ref only
        % rPwr = 0.75; % 3:1  Ref:Test
        % rPwr = 0.6; % 6:4  Ref:Test
        %   Note for Interspeech 2022
        %   There was huge difference in -20 dB condtion between lab and remote.
        %   It seems related to the reported number of Tone Pip test.
        %　Adjusting this value may  improve the prediction
end
if isfield(GESIparam.Sim,'weightAbvThrCmpnst') ==0  % Eta:  To control from the main
    GESIparam.Sim.weightAbvThrCmpnst  = 0.70; % default 15 Jul 25
    disp('GESIparam.Sim.weightAbvThrCmpnst is set to 0.70 (default) -- OK? Return to continue > ')
    pause
end


if MFBparam.MaxFc > 128
    disp(['MFBparam.MaxFc (' int2str(MFBparam.MaxFc) ' Hz) > 128 Hz' ]);
    disp('--- This is an unusual contion. No compensation nor NaN is introduced. ')
    input('If OK, return to continue > ')
    % if isfield(GESIparam.Sim,'SwWeightProhibit') == 0  % To control from the main routine
    %     GESIparam.Sim.SwWeightProhibit = 1;  % default  switch for prohibit region　 introduced 29 Jul 2022
    %     % If you want to use the conventional method, set the next line when calling this function
    %     % GESIparam.Sim.SwWeightProhibit = 0;
    % end
    % if isfield(GESIparam.Sim,'RangeWeightProhibit') == 0  % To control from the main routine
    %     GESIparam.Sim.RangeWeightProhibit  = 1; % The range of prohibit region: default 1
    % end
end


%%%%%%%%%%%%%%%%%%
%% sound  sampling rate conversion & normalization
% GCFBparam.fs may be always 48000 Hz.
% GESIparam.fs : sampling frequency of ref/test sounds
%%%%%%%%%%%%%%%%%%
if GCparam.fs ~= GESIparam.fs    % Resampling to 48 kHz
    SndTest = resample(SndTest(:)',GCparam.fs,GESIparam.fs);
    SndRef  = resample(SndRef(:)',GCparam.fs,GESIparam.fs);
end

% Calibrate input level of SndRef by using Meddis hair cell level for GCFB
% [SndRef, MdsAmpdB]   = Eqlz2MeddisHCLevel(SndRef, GESIparam.SPL);
[SndRef, MdsAmpdB]   = Eqlz2MeddisHCLevel(SndRef, [], GESIparam.DigitalRms1SPLdB); % 7 Mar 21
% SndTest should be convert with MdsAmpdB for the same ratio.
SndTest =  10^(MdsAmpdB(2)/20)*SndTest; % adjusted by the same ratio


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%    Analysis by GCFB (dynamic compressive gammachirp filterbank)
% If there are files: GESIparam.NameSndRef,  GESIparam.NameSndTest, save GCout for speed up
%       GCFB takes about 3sec for 0.9sec sound.  RTF = 3; (M1 Mac, R2022b, 14 May 2024)
%       MFB takes about 0.8 sce for 0.9sec sound. RTF < 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Test signal analized by GCFB with either NH  or HI listener
%    For test siginal, you should specify the hearling loss
DirNameTest = [GESIparam.DirGCout, GESIparam.NameGCoutTest '.mat'];
if exist(DirNameTest) == 0
    % GCFB analysis
    GCparamTest = GCparam;
    [GCoutTest, ~, GCparamTest] = GCFBv234(SndTest,GCparamTest);  % you need GCparam
    [NumCh,LenFrame] = size(GCoutTest);
    GCoutTest = EqlzGCFB2Rms1at0dB(GCoutTest, GCparam.StrFloor);  % 0dB: Abs. Thresh. level
    if SwSave == 1
        disp(['save: ' DirNameTest ])
        save(DirNameTest,'GCoutTest','GCparamTest','SndTest','SndTestOrig');
    end
else
    disp(['load: ' DirNameTest ])
    load(DirNameTest);
    [NumCh,LenFrame] = size(GCoutTest);% to avoid error 24 Oct 2024
end


MeanEpTest = mean(GCoutTest,2); % average Excitation Pattern(Ep) for test
weightAbvThrTest0 = (MeanEpTest >= GCparam.FloorLevel); % Ep above threhsold (Floor level) binary 0/1
LenAbvThr = length(find(weightAbvThrTest0)); % Number of channel above threshold
if LenAbvThr > 0 %  to avoid division by zero   3 May 2025
    for nEta = 1:length(GESIparam.Sim.weightAbvThrCmpnst)  % 15 Jul 2025
        EfficiencyVal(nEta) = (length(MeanEpTest)/LenAbvThr).^( GESIparam.Sim.weightAbvThrCmpnst(nEta) + eps );
    end                                                                                             % eps: to avoid  (anyval)^0 == 1;
else %     % GESIparam.Sim.weightAbvThrCmpnst == "eta" in Paper,    
    EfficiencyVal = zeros(size(GESIparam.Sim.weightAbvThrCmpnst));
end
weightAbvThrTest = weightAbvThrTest0(:)*EfficiencyVal(:)';



%%  Ref signal always analized by GCFB with NH
DirNameRef = [GESIparam.DirGCout, GESIparam.NameGCoutRef '.mat'];
if exist(DirNameRef) == 0
    % GCFB analysis
    GCparamRef = GCparam; % clarify Ref
    GCparamRef.HLoss = ''; % clear HLoss
    GCparamRef.HLoss.Type = 'NH';  % Always analyzed by NH. Overwrite the condition
    [GCoutRef, ~, GCparamRef] = GCFBv234(SndRef,GCparamRef);
    [NumCh, LenFrame] = size(GCoutRef);
    GCoutRef = EqlzGCFB2Rms1at0dB(GCoutRef, GCparam.StrFloor);

    if GCparamRef.DynHPAF.fs ~= GCparamTest.DynHPAF.fs %Check
        error('Something wrong: GCparamRef.DynHPAF.fs ~= GCparamTest.DynHPAF.fs ')
    end
    tFrame = (0:LenFrame-1)/GCparamRef.DynHPAF.fs; % Frame time in sec
    GCparamRef.tFrame = tFrame; % keep 24 May 2024

    if SwSave == 1
        disp(['save: ' DirNameRef ])
        save(DirNameRef,'GCoutRef','GCparamRef','SndRef');
    end
else
    disp(['load: ' DirNameRef ])
    load(DirNameRef);
    tFrame = GCparamRef.tFrame; % recover
end

% --- MeanEpRef = mean(GCoutRef,2); % Not necessary because of NH for Ref
% --- weightAbvThrRef = (MeanEpRef >= GCparam.FloorLevel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% TimeAlignment at the output of GCFB    14 May 2024, mod 24 May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  GESIparam.TimeAlign.GCFBmaxLagSec > 0
    disp(['GESI: Time alignment at the GCFB output of SndRef and SndTest: MaxLagSec = ' ...
        num2str(GESIparam.TimeAlign.GCFBmaxLagSec) ' (sec)']);
    GCoutTestOrig = GCoutTest; % Keep original
    ParamTAGC.fs  = GCparamTest.DynHPAF.fs;
    ParamTAGC.MaxLagSec = GESIparam.TimeAlign.GCFBmaxLagSec;
    [GCoutTest, ParamTAGC] = TimeAlignGCFB(GCoutTestOrig, GCoutRef, ParamTAGC);
    GESIparam.TimeAlign.GCFB = ParamTAGC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%  ModFB analysis
%    For test siginal, you should specify the hearling loss
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Test signal
MFBparam.fs = GCparamTest.DynHPAF.fs;
MFBparam.TMTFampSpec = sqrt(1./(1+(MFBparam.fc/MFBparam.TMTFfcutoff).^2)); % v130, 5 Jan 2024
MFBparam.TMTFampReduct = 10^((MFBparam.TMTFminDepth_NH - MFBparam.TMTFminDepth)/20); % Decline in depth
MFBparam.TMTFweight = [1;  MFBparam.TMTFampReduct*MFBparam.TMTFampSpec(2:end)'];
TMTFweightMtrx = MFBparam.TMTFweight*ones(1,LenFrame);

GCModEnvTest = zeros(NumCh,LenMFB,LenFrame);
for nch = 1:GCparam.NumCh
    % frame-base processing, i.e. all positive value
    % Env = filter(MFBparam.bzLPF,MFBparam.apLPF,GCoutTest(nch,:));  %  LPF  v123 and earlier

    [ModEnv, MFBparam] = FilterModFB(GCoutTest(nch,:),MFBparam);
    GCModEnvTest(nch,:,:)  = TMTFweightMtrx.*ModEnv;  % v130, 5 Jan 2024

    if  length(find(nch == [1 50 100])) > 0
        disp(['> Modulation Filterbank Analysis: GCFB ch = #' int2str(nch)]);
    end
end

%% Ref signal
% ModFB analysis by NH condition
MFBparam.TMTFampSpec_NH = sqrt(1./(1+(MFBparam.fc/MFBparam.TMTFfcutoff_NH).^2)); %  NH v130, 5 Jan 2024
MFBparam.TMTFweight_NH = [1; MFBparam.TMTFampSpec_NH(2:end)'];  %  NH　 v130, 5 Jan 2024
TMTFweightMtrx_NH = MFBparam.TMTFweight_NH*ones(1,LenFrame);
%v142 -- NG
% if isfield(GESIparam.Sim,'weightMFB') == 0  % if no control from the main routine % v142
%     % GESIparam.Sim.weightMFB = MFBparam.TMTFweight_NH;  % set to Normal hearing TMTF v142, 19 Sep 24
%
% end

GCModEnvRef = zeros(NumCh,LenMFB,LenFrame);
for nch = 1:GCparam.NumCh
    % frame-base processing, i.e. all positive value
    % Env_NH = filter(MFBparam.bzLPF_NH,MFBparam.apLPF_NH,GCoutRef(nch,:)); % LPF  v123 and earlier

    [ModEnv_NH, ~] = FilterModFB(GCoutRef(nch,:),MFBparam);
    GCModEnvRef(nch,:,:)  = TMTFweightMtrx_NH.*ModEnv_NH; % v130, 5 Jan 2024

    if  length(find(nch == [1 50 100])) > 0
        disp(['> Modulation Filterbank Analysis: GCFB ch = #' int2str(nch)]);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Calculation of SSIweight
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean Fo analysis of Ref & SSIweight using WORLD
HarvestRef = Harvest(SndRef,GCparam.fs); % every 5ms
F0Frame = interp1(HarvestRef.temporal_positions,HarvestRef.f0,tFrame,'linear','extrap');
F0Frame = max(F0Frame,0); % to avoid negative F0Frame by using interp1 22 Oct 22
F0MeanRef =  geomean(HarvestRef.f0(HarvestRef.f0>0)); %geomean of F0
disp(['Fo Mean of Ref sound: ' num2str(F0MeanRef,'%5.1f') ' Hz']);
if length(find(isnan(F0Frame))) > 0
    error('Error in F0Frame.')
end

%SSIparam.SwSSIweight = 1; % v110 fixed at F0mean (Interspeech 2022 version)
SSIparam.SwSSIweight = 2; % > v121 time-varying (Latest)
SSIparam.Fr1        = GCparamRef.Fr1;

if SSIparam.SwSSIweight == 1 % fixed value == the same as the process in v1
    warning('Fixed value of SSIparam.SwSSIweight, i.e., obsolete ---  Really OK?')
    SSIparam.F0_limit = F0MeanRef;
    [SSIweight, SSIparam] = F0limit2SSIweight(SSIparam);
    SSIparam.SwSSIweight_Norm = 1; % == v110
    if SSIparam.SwSSIweight_Norm == 1 %   == v110     % normalization
        %  Introdued normalization of SSIweight  v110    4 Aug 2022
        SSIparam.F0_limit_Norm = 125; % F0 for normalization
        SSIparamNorm = SSIparam; % copy
        SSIparamNorm.F0_limit = SSIparam.F0_limit_Norm ; % freq for normalization
        [SSIweightNorm] = F0limit2SSIweight(SSIparamNorm);
        SSIparam.weight_AmpNorm = sum(SSIweightNorm)/sum(SSIweight);
        SSIweightVal = SSIparam.weight_AmpNorm * SSIweight;
    else  % just for check
        error('Check code -- Not in use --- Obsolete')
        SSIweightVal = SSIweight/mean(SSIweight);
        % normalized for every channel -- simple.23 Aug 2022 --
        % It was not adjusted for 125 Hz. ---- Obsolete
    end
    SSIweightMtrx = SSIweightVal(:)*ones(1,LenFrame); % Matrix

elseif SSIparam.SwSSIweight == 2 % time-varying v120 and later
    for nFrame = 1:LenFrame
        SSIparam.F0_limit = F0Frame(nFrame);
        [SSIweight, SSIparam] = F0limit2SSIweight(SSIparam);
        SSIweightMtrx(:,nFrame) = SSIweight/mean(SSIweight);  % normalize by mean value. some of them > 1
    end
end
SSIparam.weight = SSIweightMtrx;

GESIparam.GCparam.Ref = GCparamRef;
GESIparam.GCparam.Test = GCparamTest;
GESIparam.MFBparam  = MFBparam;  % common in Test and Ref
GESIparam.SSIparam = SSIparam;   % only in Ref

if GESIparam.SwPlot == 1
    gcfKeep = gcf;
    figure(gcfKeep.Number+1); clf

    amp = 6;  % it seems reasonable. max(GCFB output) is nearly 30
    subplot(2,1,1);  image(amp*GCoutTest); set(gca,'YDir','normal');
    subplot(2,1,2);  image(amp*GCoutRef);  set(gca,'YDir','normal');
    drawnow
    figure(gcfKeep.Number+2); clf
    plot(SSIparam.weight)
    figure(gcfKeep.Number); % Return back to original figure number
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Cosine similarity analysis ---  Not Pearson's correlation
%           Extended cosine similarity was much better than other methods.
%           We can introduce the SSI weight with the channel-by-channel calculation.
%           You can change the order of mutiplication because it is a linear algebra.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NumCh, LenMFB, ~] = size(GCModEnvRef); % you need when reloading from mat file

if isfield(GESIparam.Sim,'weightMFB') == 0  % if no control from the main routine % v147
    GESIparam.Sim.weightMFB = ones(1,LenMFB);
end

% -------------------------------
% Modified for multiple nRho for large simulation
% -------------------------------
for nEta = 1:length(GESIparam.Sim.weightAbvThrCmpnst) % 15 Jul 25
    for nRho = 1:length(GESIparam.Sim.PowerRatio)

        for nch = 1:NumCh  % == GCparam.NumCh
            %weightGCFB = SSIparam.weight(nch,:); % time-varying frame SSI weight
            %weightGCFBinv = max(1 - SSIparam.weight(nch,:),0);  % Inverse SSI  --- not precise but good approximation
            weightGCFB    = weightAbvThrTest(nch,nEta)*SSIparam.weight(nch,:); % 14 Nov 24,  using only EP above threshold
            weightGCFBinv = weightAbvThrTest(nch,nEta)*max(1 - SSIparam.weight(nch,:),0);  % 14 Nov 24, sing only EP above threshold

            weightMFB  = GESIparam.Sim.weightMFB/mean(GESIparam.Sim.weightMFB);  % normalized to be mean == 1
            for nMFB = 1:LenMFB
                ModEnvRef  = squeeze(GCModEnvRef(nch,nMFB,:))';   % row vector
                ModEnvTest = squeeze(GCModEnvTest(nch,nMFB,:))';  % row vector

                PwrRef  = sum(ModEnvRef.^2);      % sum power
                PwrTest = sum(ModEnvTest.^2);     % sum power
                % rPwr = GESIparam.Sim.PowerRatio;  % == "rho" in paper.  Power ratio  0<=rPwr<=1
                rPwr = GESIparam.Sim.PowerRatio(nRho); % added by YA, v122 on 23 Oct 22

                % CosSim = sum(ModEnvRef.*ModEnvTest)/(PwrRef^rPwr*PwrTest^(1-rPwr)); % original
                % introduced weighted sum for time-varying SSIweight   31 Aug 22
                CosSim = sum(weightGCFB.*ModEnvRef.*ModEnvTest)/(PwrRef^rPwr*PwrTest^(1-rPwr));
                % added invese SSIweight 4 Jan 24   positive
                CosSim_InvSSI = sum(weightGCFBinv.*ModEnvRef.*ModEnvTest)/(PwrRef^rPwr*PwrTest^(1-rPwr));

                CosSim_Mtrx(nch,nMFB)             = weightMFB(nMFB)*CosSim;
                CosSim_RawMtrx(nch,nMFB)       = CosSim;
                CosSim_InvSSIMtrx(nch,nMFB)   = weightMFB(nMFB)*CosSim_InvSSI; % v130, 4 Jan 2024
                CosSim_InvSSIRawMtrx(nch,nMFB) = CosSim_InvSSI; % v130, 4 Jan 2024

            end % for nMFB = 1:LenMFB
        end     % for nch = 1:NumCh

        % Results of GESI
        Result.GESI.CosSim(:,:,nRho,nEta)               = CosSim_Mtrx;
        Result.GESI.CosSim_Raw(:,:,nRho,nEta)       = CosSim_RawMtrx;
        Result.GESI.CosSim_InvSSI(:,:,nRho,nEta)    = CosSim_InvSSIMtrx;
        Result.GESI.CosSim_InvSSIRaw(:,:,nRho,nEta) = CosSim_InvSSIRawMtrx;
        Result.GESI.CosSim_MeanGCFB(:,nRho,nEta)    = mean(CosSim_Mtrx,1); % mean value across GCFB
        Result.GESI.CosSim_MeanMFB(:,nRho,nEta)     = mean(CosSim_Mtrx,2); % mean value across MFB
        Result.GESI.dVal(nRho,nEta)                 = mean(CosSim_Mtrx,[1 2]); % mean value across GCFB

        Result.GESI.Pcorrect(nRho,nEta) = Metric2Pcorrect_Sigmoid(Result.GESI.dVal(nRho,nEta),GESIparam.Sigmoid);
        %　Preliminary output of Pcorrect -- It does not represent acculate SI

    end  % nRho = 1:length(GESIparam.Sim.PowerRatio)
end  %  nEta = 1:length(GESIparam.Sim.weightAbvThrCmpnst)


Result.d.GESI = Result.GESI.dVal; % for backward compativility
Result.Pcorrect.GESI = Result.GESI.Pcorrect; % for  backward compativility


% -------------------------------
% plot Result.dIntrm.CosSimMtrx
if GESIparam.SwPlot == 2 & length(GESIparam.Sim.PowerRatio) == 1 &  length(GESIparam.Sim.weightAbvThrCmpnst) == 1
    image(Result.GESI.CosSim*256)  %rPwr = 0.55
    set(gca,'YDir','normal');
    xlabel('MFB channel')
    ylabel('GCFB channel')
    drawnow
end



end % end of function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trash & Memo.
% Comments for the records and future studies. Please ignore them.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Removed

% Result.dIntrm.GESI(:,:,nRho) = CosSimMtrx;  % Original intermediate  representation for d
% Result.d.GESI(:,nRho)        = mean(CosSimMtrx(~isnan(CosSimMtrx))); % average non-NaN values  29 Jul 22
% It was removed for simplicity because default MFBparam.MaxFc == 32 Hz, i.e.,   MFBparam.fc <= 32 Hz
%if GESIparam.Sim.SwWeightProhibit == 1 % with prohibit region
%    if GCparamRef.Fr1(nch) < MFBparam.fc(nMFB)*GESIparam.Sim.RangeWeightProhibit
%        % The calculation for "Ch Freq < Modulation Freq" makes no sense. It is prohibitted.
%        weightMFB(nMFB) = NaN; % set NaN for such region
%    end
%end

%             if nRho == 1
%                 %%% Calculation of Raw SSIweighted Modulation Index,  v130, 4 Jan 2024
%                 ModIndexRefSSIMtrx(nch,nMFB)       = weightMFB(nMFB)* sum(weightGCFB.*ModEnvRef);
%                 ModIndexRefInvSSIMtrx(nch,nMFB)  = weightMFB(nMFB)* sum(weightGCFBinv.*ModEnvRef);
%                 ModIndexTestSSIMtrx(nch,nMFB)      = weightMFB(nMFB)* sum(weightGCFB.*ModEnvTest);
%                 ModIndexTestInvSSIMtrx(nch,nMFB) = weightMFB(nMFB)* sum(weightGCFBinv.*ModEnvTest);
%             end
% % summary of Modulation Index v130, 4 Jan 2024
% Result.ModIndex.Ref_SSI       = mean(ModIndexRefSSIMtrx(~isnan(ModIndexRefSSIMtrx)));
% Result.ModIndex.Ref_InvSSI  = mean(ModIndexRefInvSSIMtrx(~isnan(ModIndexRefInvSSIMtrx)));
% Result.ModIndex.Test_SSI      = mean(ModIndexTestSSIMtrx(~isnan(ModIndexTestSSIMtrx)));
% Result.ModIndex.Test_InvSSI  = mean(ModIndexTestInvSSIMtrx(~isnan(ModIndexTestInvSSIMtrx)));
% Result.ModIndex.Ratio_SSI      = Result.ModIndex.Test_SSI/Result.ModIndex.Ref_SSI;
% Result.ModIndex.Ratio_InvSSI = Result.ModIndex.Test_InvSSI/Result.ModIndex.Ref_InvSSI;
% GESIv144 : There is no reason why PowerRatio appears here. --- Bug
% CosSimMtrx(nch,nMFB) = GESIparam.Sim.PowerRatio(nMFB)*CosSim;
% GESIv145 : debugged


% % Intermediate Modulation Index, v130, 4 Jan 2024
% Result.ModIndexIntrm.Ref_SSI       = ModIndexRefSSIMtrx;
% Result.ModIndexIntrm.Ref_InvSSI  = ModIndexRefInvSSIMtrx;
% Result.ModIndexIntrm.Test_SSI     = ModIndexTestSSIMtrx;
% Result.ModIndexIntrm.Test_InvSSI = ModIndexTestInvSSIMtrx;
%
%
%
% Note: 5 Jan 2024 v123 --> v130
% Original LPF/TMTF implimentation, v123 and earlier
% ==========================================
% Pcorrect : 39.5188      57.0723      81.9954      93.9013
% Metric    : 0.27872     0.31424      0.3758     0.43671
% ==========================================
%
% New TMTF implimentation, v130
%==========================================
%Pcorrect : 38.798      56.4263      80.9875      93.3187
%Metric    : 0.27721     0.31292     0.37246     0.43184
%==========================================
%
% ratio in Metric
%   1.0054    1.0042    1.0090    1.0113
%  They are less than 1%.  --- OK small difference
%


%    GESIparam.IdealObserver =  [k q m sigma_s](not used in the current version)

% comments when writing Interspeech 2022
% Normalizing power of Ref / Test ーーー　rPwr == pho in Paper
% rPwr = 1; % Ref only
% rPwr = 0.75; % 3:1  Ref:Test
% rPwr = 0.6; % 6:4  Ref:Test
%     田丸実験で、 -20dB条件の SRTが人ごとに大きく異なる。これは、あきらかに tone pip数と関連しているであろう。
%　  個人ごとの違いを出すためにはこれを調整か？


%    HarvestRef = Harvest(SndRef,GCparam.fs); % every 5ms
% 次は、v121 f0をちゃんと反映させる。
% interp1(0:180,HarvestRef.f0,0:0.1:180)
% を使って、SSIweightを時点ごとで正確に反映
% similarityの計算の中に重みとしていれる。modulationの計算には入れない方が良いので。
% 子音も有声音も同じSSIweightで重み付けはやはり変。


%    v110 -- v121:  MFBparam.fcutEnv =  max(MFBparam.fc);  % v110  4 Jun 2022  extended fcut of the MFB  for using 512 Hz ModFilter ：　256にする理由もない。
% MFBparam.fcutEnv = 150; % before v109 & after v122, 15 Oct 2022
%  v122 TMTFを考えにいれておいた方が良い. muiti-resolution sESPM, 山克GEDIでは導入済み。
%   誤差は1% 程度でほとんど影響がないことを確認。rms(Metric-Metric1)/rms(Metric) = 0.0114
%   いままでの研究の流れから、これはいれておく。v122

% GCoutRef = EqlzGCFB2Rms1at0dB(GCoutRef, GCparam.StrFloor)+ alpha*randn(GCparamRef.NumCh,length(GCoutRef)); % Add noise, 28 Apr, Yamamoto


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 将来の変更に対する余地を残している。
% 最終結果は、以下のよう出力を決めうちして、GESIという名前にしたい。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Result.Note       = 'Result.d.GESI  == Result.d.GCMFBtmc';
%Result.d.GESI    = Result.d.GCMFBtmc;
%Result.Pcorrect.GESI = Result.Pcorrect.GCMFBtmc;

% Previous implementation
% dGCMFBtmc(nch,nMFB) = weightGCFB(nch)*weightMFB(nMFB)*CosSim; %



%  History of variants.
%  Basically they were not very successful. -- Comment out now
%%%%%%%%%%%%%%%%%%%%%%
% Correlation analysis of ESTOI like in spectral domain: frame by frame
% GESI_CorrFrame
% correlation anaysis in MFB spectral dimension
% GESI_CorrMFBspec
%  Memo and trash
% GESI_Memo


% if isfield(GESIparam.Sim, 'ThrshLowPower') == 0
%     GESIparam.Sim.ThrshLowPower = 1.5;
% end

% Absolute thresholdより内部の雑音の方が影響 YA
% GCoutTest = EqlzGCFB2Rms1at0dB(GCoutTest, GCparam.StrFloor) + alpha*randn(GCparamTest.NumCh,length(GCoutTest)); % Add noise, 28 Apr, Yamamoto.

% ModEnvRef  = GCModEnvRef(nch,nMFB,:);
% ModEnvTest = GCModEnvTest(nch,nMFB,:);
%  これにrPwrを入れた場合、80/70yr 反転問題は解消したものの、Unprocが一番下、LowLevelが一番上。これはダメ。
%  -->   重み関数の導入　SSIweight
%  9 Feb 22
%     SSI weightこそ、ここで使える。
%     Articulation indexの重み付け係数^0.3とほぼ同等となるので、説明できる....
%     0.3乗は、compression分。でもAIは、dB上の重み関数。
%     See  ./Tool/testArticulationIndexHL.m
%  19 Feb 22
%   ... と思ったのはこじつけかも。
%  20 Feb 22
%   計算の手順上、ModEnvRef/ModEnvTestのところに入れる。

%
% crTmp= corrcoef(ModEnvRef,ModEnvTest);  Pearson Correlation:  =="rPwr =0.5"
% dGCMFBtmpPC(nch,nMFB) = crTmp(1,2);

% rPwr = 0.6 SSIweight:   Unpro == LowL > 70yr > 80yr  ほぼOK
% rPwr = 0.6 : pwr== 1.2, Wt8, Wt9  Unpro == LowL > 70yr > 80yr  ほぼOK
%  rPwr = 0.75 : pwr== 1.5, Wt8 Unpro >> LowL > 70yr > 80yr  LowL小さい。80yrもっとあげたい
%  rPwr = 0.75 : pwr== 1.5, Wt6 Unpro >> 70yr ==  > 80yr   LowLは下がる。高い周波数重要
% pwr== 0.75, Wt7 Unpro = LowL >> 70yr > 80yr レベルがだいぶあってきた
% pwr== 1: pwr== 2, Wt7 Unpro >> LowL >> 70yr > 80yr 傾向的にはOKだが、やはり、レベルが違いすぎ
% pwr== 1: pwr== 2, Wt5 Wt6 も同じ傾向　Unpro >>> 70yr > > 80yr
% Result.dIntrm.GCMFBtmp   = dGCMFBtmp;  % Pearson Correlation Original intermediate d
% Result.d.GCMFBtmp           = mean(mean(Result.dIntrm.GCMFBtmp));
% Result.Pcorrect.GCMFBtmp = Metric2Pcorrect_Sigmoid(Result.d.GCMFBtmp,GESIparam.Sigmoid);
% if GESIparam.SwPlot == 1
%     subplot(2,1,1);  imagesc(dGCMFBtmp);
%     subplot(2,1,2);  imagesc(dGCMFBtmprPwr);
% end
%

%       以下の組み合わせはダメ
%       ModEnvRef  = SSIparam.weight(nch)*GCModEnvRef(nch,nMFB,:);
%         ModEnvTest = GCModEnvTest(nch,nMFB,:);

% こちらは、normalizeした方がband数によっての変動を吸収できる。
% weightMFB   = GESIparam.Sim.weightMFB/sum(GESIparam.Sim.weightMFB); % sum = 1;  % 12 May 22
%
% weightGCFBは、Foによって変わってしまうため、以下は使わない方が良い。
% weightGCFB = SSIparam.weight/sum(SSIparam.weight); % sum == 1    % 12 May 22



% ModEnvRef  = squeeze(GCModEnvRef(nch,nMFB,:)) + 0.5*randn(LenEnv,1); % リファレンスにもnoiseをつけてみる, 2 May, YA
% ModEnvTest = squeeze(GCModEnvTest(nch,nMFB,:)) + 0.5*randn(LenEnv,1); % Add noise, 27 Apr, YA

% Note: 15 Oct 2022   v121--> v122
% v121 Metric
% ==========================================
% Pcorrect : 38.7108      56.6441      80.9039      93.2404
% Metric    : 0.27703     0.31337     0.37219     0.43121
% ==========================================
%
% v122 Metric1
% ==========================================
% Pcorrect : 39.5992      56.8733      81.7984      94.1002
% Metric    : 0.27889     0.31383     0.37514     0.43847
% ==========================================

% [MFBparam.bzLPF, MFBparam.apLPF] = butter(1, MFBparam.TMTFfcutoff/(MFBparam.fs/2));  % v123 and earlier
% [MFBparam.bzLPF_NH, MFBparam.apLPF_NH] = butter(1,MFBparam.TMTFfcutoff_NH/(MFBparam.fs/2)); %  LPF  v123 and earlier



% v148 original, but it is difficult to interpret and also includes the bug of division by zero.   3 May 2025
% weightAbvThrTest1 = ( weightAbvThrTest0/mean(weightAbvThrTest0) ).^( GESIparam.Sim.weightAbvThrCmpnst + eps );
% if mean(abs(weightAbvThrTest-weightAbvThrTest1)) ~= 0
%    [weightAbvThrTest(:)'; weightAbvThrTest1(:)']
%    warning('Check weightAbvThrTest')
% else
%   disp('===============  OK ============================')
% end

% weightAbvThrTest = weightAbvThrTest/mean(weightAbvThrTest); %
%  ---  normalizedto be mean === 1  -- NG too much compensation



% introducing weightAbvThrTest (14 Nov 2024) weightAbvThrCmpnst (v147 18 Nov 2024)
% default値はなしにする。15 Jul 2025
% if isfield(GESIparam.Sim,'weightAbvThrCmpnst') == 0 % new param v147 18 Nov 2024
%    %  GESIparam.Sim.weightAbvThrCmpnst = 0.5; % It means the power is the same. 18 Nov 2024
% GESIparam.Sim.weightAbvThrCmpnst = 1; % Control by rho 22 Nov 2024
% end