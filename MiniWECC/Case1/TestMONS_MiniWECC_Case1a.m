%Baseline test for AMBout
clear; close all; clc
addpath('../../../'); %Location of MONS software

%% Load MiniWECC simulation
% Path = 'C:\Users\dtrudnowski\OneDrive - Montana Tech\Documents\2024_25\Consulting\Sandia_RealTimeAssessment\MiniWECC\';
File = 'd_minniWECC_V3C_C3_6_C_PowerPlantC_20m_LoadNoise1_PMU.mat';
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
ParamAMB(k).RMSeMax.vmag = [NaN 0.6 0.1 0.1]; %kV
ParamAMB(k).RMSeMax.f = [NaN 100 2 1]; %10*milliHz
ParamAMB(k).RMSeMax.mw = [NaN 4 0.3 0.1]; %MW
ParamAMB(k).RMSeMax.mvar = [NaN 4 0.3 0.1]; %MVAR
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

%% Solve outside of MONS for comparison
for k=1:1%length(ConfigDerived.MW)
    nv = find(strcmp(ConfigDerived.MW(k).SigName{1},{SynData.Name}));
    ni = find(strcmp(ConfigDerived.MW(k).SigName{2},{SynData.Name}));
    S = 3*SynData(nv).Data.*conj(SynData(ni).Data)./1000; %MW
    V = abs(SynData(nv).Data);
    x = (fs/(2*pi))* angle(SynData(nv).Data(2:end).*conj(SynData(nv).Data(1:end-1)));
    F = [x(1);x]; %HZ
    T = (1/fs)*[0:length(S)-1]';
    clear nv ni

    %plot derived signals
    figure
    subplot(411)
    kD = find(strcmp(ConfigDerived.MW(k).DerivedName,{DerivedOut.Name}));
    plot(T./60,real(S),'k','linewidth',2)
    hold on
    plot(T./60,DerivedOut(kD).Data,'r')
    ylabel('MW')
    title(ConfigDerived.MW(k).DerivedName)
    subplot(412)
    kD = find(strcmp(ConfigDerived.MVAR(k).DerivedName,{DerivedOut.Name}));
    plot(T./60,imag(S),'k','linewidth',2)
    hold on
    plot(T./60,DerivedOut(kD).Data,'r')
    ylabel('MVAR')
    subplot(413)
    kD = find(strcmp(ConfigDerived.MAG(k).DerivedName,{DerivedOut.Name}));
    plot(T./60,V,'k','linewidth',2)
    hold on
    plot(T./60,DerivedOut(kD).Data,'r')
    ylabel('V')
    subplot(414)
    kD = find(strcmp(ConfigDerived.AngleDerivative(k).DerivedName,{DerivedOut.Name}));
    plot(T./60,F,'k','linewidth',2)
    hold on
    plot(T(3:end)./60,DerivedOut(kD).Data(3:end)./10000,'r')
    ylabel('Hz')
    xlabel('Time (min.)')

    %Conduct analysis
    n = length(T);
    m = round(SessionInfo.AMBreportRate*fs);
    Ns = round(ParamAMB(k).Twin*fs);
    Nd = [1:m:n];
    Nd = Nd(1:round((n-Ns)/m));
    NdL = length(Nd);
    Time = T(Nd);
    clear n m P
    for kt=1:NdL
        disp(['Iteration ' num2str(k) ', ' num2str(kt) ' of ' num2str(NdL)])
        kw = [Nd(kt):Nd(kt)+Ns];
        t = T(kw) - T(kw(1));
        y = [F(kw) real(S(kw))];

        %Kdamp and Ksyn
        [Syy,Sxy,Cxy,f] = funSpectrum(t,y,[],fs,ParamAMB(k).TfftSynD,'H',1,[],[]);
        if kt==1
            z = zeros(length(f),NdL);
            Ksd.Sww = z;
            Ksd.Cw_p = z;
            Ksd.Gw_p = z;
            Ksd.Ks= z;
            Ksd.Kd = z;
            %fn_index = cell(NdL,1);
            Ksd.bar_index = zeros(NdL,2);
            Ksd.Ksbar = zeros(NdL,1);
            Ksd.Kdbar = zeros(NdL,1);
            Ksd.f = f;
            clear z
        end
        Ksd.Sww(:,kt) = Syy(:,1);
        Ksd.Cw_p(:,kt) = Cxy(:,2);
        Ksd.Gw_p(:,kt) = conj(Sxy(:,2) ./ (10.^(Syy(:,1)./10)));
        Ks = -imag(Ksd.Gw_p(:,kt)) ./ (2*pi*f);
        n = find(f>1 & f<3 & Ks<0);
        if isempty(n); error(' '); end
        [~,i] = min(abs(f-0.7*f(n(1))));
        m = find(f>=0.1 & f<=f(i) & Ksd.Cw_p(:,kt)>0.5);
        Ksd.Ks(:,kt) = Ks;
        Ksd.Ksbar(kt) = mean(Ks(m));
        Ksd.Kd(:,kt) = -real(Ksd.Gw_p(:,kt));
        Ksd.Kdbar(kt) = mean(Ksd.Kd(m,kt));
        Ksd.bar_index(kt,:) = [m(1) m(end)];
        clear m n f Syy Sxy Cxy f

        %Kgov
        [Syy,Sxy,Cxy,f] = funSpectrum(t,y,[],fs,ParamAMB(k).TfftGov,'H',1,[],[]);
        m = find(f<0.1);
        if kt==1
            Kg.f = f(m);
            z = zeros(length(f(m)),NdL);
            Kg.Sww = z;
            Kg.Cw_p = z;
            Kg.Kg = z;
            clear z
        end
        Kg.Sww(:,kt) = Syy(m,1);
        g = conj(Sxy(:,2) ./ (10.^(Syy(:,1)./10)));
        Kg.Cw_p(:,kt) = Cxy(m,2);
        Kg.Kg(:,kt) = -real(g(m));
        clear Syy Sxy Cxy f m g t y
    end

    figure
    subplot(211)
    plot((ParamAMB(k).Twin+T(Nd))./60,Ksd.Ksbar,'ko',AMBout(k).Time./60,AMBout(k).Ksd.Ksbar,'r*')
    ylabel('Ksbar')
    subplot(212)
    plot((ParamAMB(k).Twin+T(Nd))./60,Ksd.Kdbar,'ko',AMBout(k).Time./60,AMBout(k).Ksd.Kdbar,'r*')
    ylabel('Kdbar')
    xlabel('Time (min.)')

    nP = 16;
    figure
    subplot(211)
    semilogx(Ksd.f(2:end),Ksd.Ks(2:end,nP-14),'k',AMBout(k).Ksd.f(2:end,nP),AMBout(k).Ksd.Ks(2:end,nP),'r')
    ylabel('Ks')
    subplot(212)
    semilogx(Ksd.f(2:end),Ksd.Kd(2:end,nP-14),'k',AMBout(k).Ksd.f(2:end,nP),AMBout(k).Ksd.Kd(2:end,nP),'r')
    ylabel('Kd')
    xlabel('freq (Hz)')

end
