%Test for AMBout
%Gen 38 is tripped at 30 min point in simulation.
clear; close all; clc
addpath('../../../'); %Location of MONS software

%% Load MiniWECC simulation
% Path = 'C:\Users\dtrudnowski\OneDrive - Montana Tech\Documents\2024_25\Consulting\Sandia_RealTimeAssessment\MiniWECC\';
File = 'd_minniWECC_V3C_C3_6_C_PowerPlantC_Gen38trip_60m_LoadNoise1_PMU.mat';
% load([Path File],'pmuData');
load(File, 'pmuData');
fs = 60; %sample rate for all PMUs

%Synchros to keep:
SynMap = [...
    %Type  fbus tbus
    1      33   0;     %Colsrip HV voltage
    2      32   33;    %Gen 14
    2     126   33];   %Gen 38

SessionInfo.AMBreportRate = 60;
SessionInfo.ODreportRate = 1;

%% Solve using MONS
% Build Synchro data compatilble with MONS (assume all are 500 kV buses)
nS = size(SynMap,1);
x.Data = [];
x.Flag = [];
x.Name = [];
SynData(1:nS) = x;
for k=1:nS
    f = false;
    for kk=1:pmuData.n_pmu
        if max(abs(SynMap(k,:) - pmuData.pmu_con(kk,1:3))) < eps
            f = true;
            break
        end
    end
    if ~f; error('Synchro not found'); end
    if SynMap(k,1)==1
        s = 500/sqrt(3);
        SynData(k).Name = ['Bus ' num2str(SynMap(k,2)) ' V'];
    elseif SynMap(k,1)==2
        s = (100e6)/(sqrt(3)*500e3);
        SynData(k).Name = ['Bus ' num2str(SynMap(k,2)) ' to ' num2str(SynMap(k,3)) ' I'];
    else
        error(' ')
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
for k=1:pmuData.n_pmu
    if pmuData.pmu_con(k,1)==1
        ConfigSyn(k).Type = 'V';
        ConfigSyn(k).Name = ['Bus ' num2str(pmuData.pmu_con(k,2)) ' V'];
        ConfigSyn(k).Units = 'kVLN';
        ConfigSyn(k).Bounds = [99999999; 0];
    elseif pmuData.pmu_con(k,1)==2
        ConfigSyn(k).Type = 'I';
        ConfigSyn(k).Name = ['Bus ' num2str(pmuData.pmu_con(k,2)) ' to ' num2str(pmuData.pmu_con(k,3)) ' I'];
        ConfigSyn(k).Units = 'A';
        ConfigSyn(k).Bounds = [99999999; 0];
    else
        error(' ')
    end
end

% Configure derived signals
ConfigDerived.AngleJumpThreshold = 0;
ConfigDerived.MW(1).DerivedName = 'Bus 32 to 33 MW';
ConfigDerived.MW(1).SigName = {'Bus 33 V', 'Bus 32 to 33 I'};
ConfigDerived.MW(1).ScaleFac = 1;

ConfigDerived.MVAR(1).DerivedName = 'Bus 32 to 33 MVAR';
ConfigDerived.MVAR(1).SigName = {'Bus 33 V', 'Bus 32 to 33 I'};
ConfigDerived.MVAR(1).ScaleFac = 1;

ConfigDerived.MAG(1).DerivedName = 'Bus 33 VMAG';
ConfigDerived.MAG(1).SigName = {'Bus 33 V'};
ConfigDerived.MAG(1).ScaleFac = 1;

ConfigDerived.AngleDerivative(1).DerivedName = 'Bus 33 f';
ConfigDerived.AngleDerivative(1).SigName = {'Bus 33 V'};
ConfigDerived.AngleDerivative(1).ScaleFac = 1;

clear pmuData

% ParamAMB
k = 1;
ParamAMB(k).Name = 'AMB 1';
ParamAMB(k).MWname = 'Bus 32 to 33 MW';
ParamAMB(k).ADname = 'Bus 33 f';
ParamAMB(k).MVARname = 'Bus 32 to 33 MVAR';
ParamAMB(k).VMAGname = 'Bus 33 VMAG';
ParamAMB(k).TfftSynD = 50;
ParamAMB(k).TfftGov = 400;
ParamAMB(k).Twin = 15*60;
ParamAMB(k).RMSeMax.vmag = [0.8 0.6 0.1 0.1]; %kV
ParamAMB(k).RMSeMax.f = [400 100 2 1]; %10*milliHz
ParamAMB(k).RMSeMax.mw =   [6 4 0.3 0.1]; %MW
ParamAMB(k).RMSeMax.mvar = [6 4 0.3 0.1]; %MVAR
ParamAMB(k).PrctInvld = 10;

% ParamFRE
ParamFRE = [];

% Run MONS Matlab
DispFlag = true;
[AMBout,FREout,ODout,DerivedOut] = funRunMONSv1Engine_Matlab(SynData,SessionInfo,ConfigSyn,ConfigDerived,ParamAMB,ParamFRE,DispFlag);

%% Plot RMSe
for k=1:length(ODout)
    figure
    for kk=1:4
        subplot(4,1,kk)
        plot(ODout(k).Time./60,ODout(k).RMSEnergy(:,kk),'k','LineWidth',1.5)
        grid
        if kk==1; title(ODout(k).Name); end
        ylabel(['Band ' num2str(kk)])
    end
    xlabel('Time (min.)')
end
clear k kk

%% Plot Bus 33 freq and RMSe
tR = 60*30 + [-10 50];
kD = find(strcmp('Bus 33 f',{DerivedOut.Name}));
t = (1/DerivedOut(kD).SampleRate) * [0:length(DerivedOut(kD).Data)-1]';
n = find(t>=tR(1) & t<tR(2));
kO = find(strcmp('OD for Bus 33 f',{ODout.Name}));
m = find(ODout(kO).Time>=tR(1) & ODout(kO).Time<tR(2));
figure
subplot(5,1,1)
plot(t(n),60+DerivedOut(kD).Data(n)./10000,'k','LineWidth',1.5)
ylim([59.8 60.2])
grid
ylabel('Freq. (Hz)')
title('Bus 33 f')
for kk=1:4
    subplot(5,1,kk+1)
    plot(ODout(kO).Time(m),0.1*ODout(kO).RMSEnergy(m,kk),'k','LineWidth',1.5)
    x = axis;
    ylim([0 x(4)])
    grid
    ylabel(['Band ' num2str(kk) ', mHz'])
end
xlabel('Time (sec.)')
set(gcf,"Position",[680 0.5*458 560 1.7*420])
%clear tR kD t n kO m

%% Plot derived signals
tStart = 10*60; 
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
n = length(ConfigDerived.MVAR);
figure
for k=1:n
    subplot(n,1,k)
    kD = find(strcmp(ConfigDerived.MVAR(k).DerivedName,{DerivedOut.Name}));
    T = (1/DerivedOut(kD).SampleRate) * [0:length(DerivedOut(kD).Data)-1]';
    n = find(T>=tStart);
    plot(T(n)./60,DerivedOut(kD).Data(n),'k')
    ylabel('MVAR')
    title(ConfigDerived.MVAR(k).DerivedName)
end
xlabel('Time (min.)')
n = length(ConfigDerived.MAG);
figure
for k=1:n
    subplot(n,1,k)
    kD = find(strcmp(ConfigDerived.MAG(k).DerivedName,{DerivedOut.Name}));
    T = (1/DerivedOut(kD).SampleRate) * [0:length(DerivedOut(kD).Data)-1]';
    n = find(T>=tStart);
    plot(T(n)./60,DerivedOut(kD).Data(n),'k')
    ylabel('kV')
    title(ConfigDerived.MAG(k).DerivedName)
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
    ylabel('Hz')
    title(ConfigDerived.AngleDerivative(k).DerivedName)
end
xlabel('Time (min.)')

%% Plot Kd and Ks
for k=1:length(AMBout)
    figure
    subplot(211)
    plot(AMBout(k).Time./60,AMBout(k).Ksd.Ksbar,'k*')
    ylabel('Ksbar')
    subplot(212)
    plot(AMBout(k).Time./60,AMBout(k).Ksd.Kdbar,'k*')
    ylabel('Kdbar')
    xlabel('Time (min.)')
end
