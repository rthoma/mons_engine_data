%% Baseline test for FREout
% Trip Load at 15 min
% Line trip at 30 min
% Trip Gen at 45 min

close all; clear; clc;
addpath('../../../'); %Location of MONS software

%% Load MiniWecc simulation
%Path = '../../../../';
% Path = 'C:\Users\dtrudnowski\OneDrive - Montana Tech\Documents\2024_25\Consulting\Sandia_RealTimeAssessment\MiniWECC\TrippingCases\';
%File = 'd_minniWECC_V3C_C3_6_C_PowerPlantC_20m_freq_LoadNoise1_PMU.mat';
File = 'd_minniWECC_V3C_C3_6_C_PowerPlantC_60m_3freq_LoadNoise1_Case_3freq_PMU.mat';
% load([Path File], 'pmuData');
load(File, 'pmuData');
fs = 60; % sample rate for all PMUs

% Synchrophasors to keep:
SynMap = [...
    %Type   fbus    tbus
    1       33      0;      % Colstrip HV Voltage
    2       32      33;     % Gen 14
    2       126     33;];   % Gen 15

SessionInfo.FREreportRate = 60;
SessionInfo.ODreportRate = 1;

%% MONS Input Setup
% SynData will be InData for the engine
nS = size(SynMap, 1);
x.Data = [];
x.Flag = [];
x.Name = [];
SynData(1:nS) = x;
for k = 1:nS
    f = false;
    for kk = 1 : pmuData.n_pmu
        if max(abs(SynMap(k,:) - pmuData.pmu_con(kk,1:3))) < eps
            f = true;
            break
        end
    end
    if ~f; error('Synchro not found'); end
    if SynMap(k,1) == 1 % V vs I pmu
        s = 500/sqrt(3);
        SynData(k).Name = ['Bus ' num2str(SynMap(k,2)) ' V'];
    elseif SynMap(k,1) == 2
        s = (100e6)/(sqrt(3)*500e3);
        SynData(k).Name = ['Bus ' num2str(SynMap(k,2)) ' to ' num2str(SynMap(k,3)) ' I'];
    else
        error('PMU type error')
    end
    SynData(k).Data = s*pmuData.pmu{kk}.data;
    SynData(k).Flag = false(length(pmuData.pmu{kk}.data),1);
end
clear x k s

% Build Configuration for Synchros in the data file
x.Name = [];
x.Type = [];
x.Format = 1;
x.Units = [];
x.SampleRate = fs;
x.Bounds = [];
ConfigSyn(1:pmuData.n_pmu) = x;
for k = 1 : pmuData.n_pmu
    if pmuData.pmu_con(k,1) == 1
        ConfigSyn(k).Type = 'V';
        ConfigSyn(k).Name = ['Bus ' num2str(pmuData.pmu_con(k,2)) ' V'];
        ConfigSyn(k).Units = 'kVLN';
        ConfigSyn(k).Bounds = [99999999; 0];
    elseif pmuData.pmu_con(k,1) == 2
        ConfigSyn(k).Type = 'I';
        ConfigSyn(k).Name = ['Bus ' num2str(pmuData.pmu_con(k,2)) ' to ' num2str(pmuData.pmu_con(k,3)) ' I'];
        ConfigSyn(k).Units = 'A';
        ConfigSyn(k).Bounds = [99999999; 0];
    else
        error('PMU type error')
    end
end

% Configure derived signals
ConfigDerived.AngleJumpThreshold = 0;

ConfigDerived.MW(1).DerivedName = 'Bus 32 to 33 MW';
ConfigDerived.MW(1).SigName = {'Bus 33 V', 'Bus 32 to 33 I'};
ConfigDerived.MW(1).ScaleFac = 1;

ConfigDerived.MVAR = [];

ConfigDerived.MAG = [];

ConfigDerived.AngleDerivative(1).DerivedName = 'Bus 33 f';
ConfigDerived.AngleDerivative(1).SigName = {'Bus 33 V'};
ConfigDerived.AngleDerivative(1).ScaleFac = 1;

% ParamAMB
ParamAMB = [];

% ParamFRE
k = 1;
ParamFRE(k).Name = 'FRE 1';
ParamFRE(k).MWname = 'Bus 32 to 33 MW';
ParamFRE(k).ADname = 'Bus 33 f';
ParamFRE(k).MAshortWin = 5/60; % sec
ParamFRE(k).MAlongWin = 200/60; % sec
ParamFRE(k).MAthreshold = 0.03; % Hz
ParamFRE(k).BufferWin = 5.5 * 60; % sec
ParamFRE(k).SSwin = 140; % sec
ParamFRE(k).SSsubWin = 45; % sec
ParamFRE(k).SSoffsetPre = 1/3; % seconds
ParamFRE(k).SSoffsetPost = 10; % sec
ParamFRE(k).MWthreshold = 20; % MW
ParamFRE(k).PrctInvld = 10; % percent

clear x
%% Run MONS Matlab
DispFlag = true;
[AMBout,FREout,ODout,DerivedOut] = funRunMONSv1Engine_Matlab( ...
    SynData, ...
    SessionInfo, ...
    ConfigSyn, ...
    ConfigDerived, ...
    ParamAMB, ...
    ParamFRE, ...
    DispFlag);

%% External Solver
x.R.Val = [];
x.R.Time = [];
x.R.GoodData = [];
FREext(1:length(ParamFRE)) = x;
clear x

k = 1;

kD = find(strcmp(ConfigDerived.MW(k).DerivedName,{DerivedOut.Name}));
p_x = DerivedOut(kD).Data;
kD = find(strcmp(ConfigDerived.AngleDerivative(k).DerivedName,{DerivedOut.Name}));
f_x = DerivedOut(kD).Data;

% MA filter signals
% Short MA fitler
win_s = fs * ParamFRE(k).MAshortWin;
th_s = ones(1,win_s)/win_s;
f_x_ma_s = funMA(f_x, th_s);

% Long MA filter
win_l = fs * ParamFRE(k).MAlongWin;
th_l = ones(1,win_l)/win_l;
f_x_ma_l = funMA(f_x, th_l);

% Edge trigger
MA_threshold = ParamFRE(k).MAthreshold * 10000; % frequency scaling
MA_flag = zeros(size(p_x));
MA_flag(abs(f_x_ma_l-f_x_ma_s)>MA_threshold) = 1;
n_block = fs * ParamFRE(k).BufferWin/2; 

i = 1;
indx_search = n_block; % skip first n data points
while_flag = true;
while while_flag
    tmp = find(diff(MA_flag(indx_search:end,:)) ~= 0, 1, "first");
    if isempty(tmp)
        break;
    end
    indx_event(i) = tmp + indx_search - 1;
    indx_search = indx_event(i) + n_block;
    i = i + 1;
end
event_times = (indx_event-1)/fs;

clear i indx_search tmp

% Calculate parameter for event each event
win_ss = fs * ParamFRE(k).SSsubWin;
N = fs * ParamFRE(k).SSwin;
off_pre = fs * ParamFRE(k).SSoffsetPre; % seconds
off_post = fs * ParamFRE(k).SSoffsetPost; % sec

Nevents = length(indx_event);


for i = 1:Nevents
    n = indx_event(i);

    pre_e = n - off_pre;
    pre_s = pre_e - N;
    [f_x_ss_pre, ind_f_x_ss_pre]=funSSScan(f_x,pre_s,pre_e,win_ss);
    [p_x_ss_pre, ind_p_x_ss_pre]=funSSScan(p_x,pre_s,pre_e,win_ss);
    
    post_s = n + off_post;
    post_e = post_s + N;
    [f_x_ss_post, ind_f_x_ss_post]=funSSScan(f_x,post_s,post_e,win_ss);
    [p_x_ss_post, ind_p_x_ss_post]=funSSScan(p_x,post_s,post_e,win_ss);

    if (abs(p_x_ss_post - p_x_ss_pre) > 10)
        good_data = true; 
    else 
        good_data = false; 
    end

    droop = (p_x_ss_post - p_x_ss_pre)/((f_x_ss_post - f_x_ss_pre)/10000);

    FREext(k).R.Val(i) = droop;
    FREext(k).R.Time(i) = event_times(i);
    FREext(k).R.GoodData(i) = good_data;
end

clear k

%% Compare Results

disp('Results');

format long

k = 1;

FREind1 = (FREext.R.GoodData > 0);
FREind2 = (FREout(k).R.Freshness > 0);

disp('External Calc. Droop');
disp(FREext(k).R.Val(FREind1)');
disp('MONs Calc. Droop');
disp(FREout(k).R.Val(FREind2));

disp('External Calc. Event Times');
disp(FREext(k).R.Time(FREind1)');
disp('MONs Calc. Event Times');
disp(FREout(k).R.Time(FREind2));

clear k
format default
%% Plot DerivedOUT against External calc.
tStart = 0.1*60; 
n = length(ConfigDerived.MW);
figure
for k=1:n
    subplot(n,1,k)
    kD = find(strcmp(ConfigDerived.MW(k).DerivedName,{DerivedOut.Name}));
    T = (1/DerivedOut(kD).SampleRate) * [0:length(DerivedOut(kD).Data)-1]';
    n = find(T>=tStart);
    plot(T(n)./60,DerivedOut(kD).Data(n),'k')
    ylabel('MW')
    title(ConfigDerived.MW(k).DerivedName)
end
xlabel('Time (min.)')

n = length(ConfigDerived.AngleDerivative);
figure
for k=1:n
    subplot(n,1,k)
    kD = find(strcmp(ConfigDerived.AngleDerivative(k).DerivedName,{DerivedOut.Name}));
    T = (1/DerivedOut(kD).SampleRate) * [0:length(DerivedOut(kD).Data)-1]';
    n = find(T>=tStart);
    plot(T(n)./60,DerivedOut(kD).Data(n)./10000,'k')
    hold on;
    plot(T(n)./60,FREout(k).MA(n,1)./10000);
    plot(T(n)./60,f_x_ma_s(n)./10000,'--'); 
    plot(T(n)./60,FREout(k).MA(n,2)./10000);
    plot(T(n)./60,f_x_ma_l(n)./10000,'--'); 
    hold off;
    ylabel('Hz')
    title(ConfigDerived.AngleDerivative(k).DerivedName)
    legend('Data','MONs MA short','Ext. MA short','MONs MA long','Ext. MA long');
    legend('location','best');
end
xlabel('Time (min.)')


%% Functions
% filter signal by some window

function filtered_sig = funMA(signal, th)
    % signal
    % th is filter coefficients
    filtered_sig = NaN(size(signal));
    Nwin = max(size(th));
    ns = Nwin;
    ne = max(size(signal));
    for n = ns:ne
        phi = signal(n-(Nwin-1):n);
        filtered_sig(n) = sum(phi.'*th.');
    end
end

% scan for steady state value
function [ss_val, ind] = funSSScan(signal,start,stop,win)
    min_slope = inf;
    ind = NaN;
    ones_chunk = ones(win+1,1);
    x_chunk = [1:(win+1)].';

    for k = start : (stop - win)
        y_chunk = signal(k : (k+win));
        th = [x_chunk, ones_chunk]\y_chunk;
        slope = abs(th(1));
        if slope < min_slope
            min_slope = slope;
            ind = k;
        end
    end

    ss_val = mean(signal(ind : (ind+win)));
end