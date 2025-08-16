clear;clc;
%选择文件夹
% path='E:\Lai\2022_Mn_14365\D2_site01_TRIAL07';
% cd(path)
load Combined.mat
%%

%Set the parameters of the MTM.
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 300];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.

%speparam=params;
%speparam.fpass=[0 120];

%时间窗参数length steplength
movingwin=[0.3 0.03];

% 获取当前工作空间中的所有变量名
vars = who;
indxL=startsWith(vars, 'LFP');
indxS=startsWith(vars, 'SPKC');
nameL=vars(find(indxL==1));
nameS=vars(find(indxS==1));

%channelnum=length(nameL);

%%
%total time-independent spectrum
filepath=[path '\results\spectrum']; mkdir(filepath);
for i=1:length(nameL) 
    cuee=['Total spectrum Processing .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[Ca,phia,~,Sf,Ssa,f]=coherencycpt(',lfp,',',spk1,',','speparam);']);
    eval(['[Cb,phib,~,Sf,Ssb,f]=coherencycpt(',lfp,',',spk2,',','speparam);']);
    spectrum(2*i-1).channel=lfp;spectrum(2*i-1).neuron='a';spectrum(2*i-1).C=Ca;spectrum(2*i-1).Phi=phia;spectrum(2*i-1).Sf=Sf;spectrum(2*i-1).Ss=Ssa;spectrum(2*i-1).f=f;
    spectrum(2*i).channel=lfp;spectrum(2*i).neuron='b';spectrum(2*i).C=Cb;spectrum(2*i).Phi=phib;spectrum(2*i).Sf=Sf;spectrum(2*i).Ss=Ssb;spectrum(2*i).f=f;
    figure('Visible', 'off');
    subplot(2,4,1);
    plot(f, Ssa, 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Power (Hz)');title ('Spectrum Neuron_a');set(gca, 'FontSize', 14)
    subplot(2,4,2);
    plot(f, 10*log10(Sf), 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Power (Hz)');title(lfp);set(gca, 'FontSize', 14)
    subplot(2,4,3);
    plot(f, Ca, 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Coherence');title('Coherence');set(gca, 'FontSize', 14);
    subplot(2,4,4);
    rose(phia);
    subplot(2,4,5);
    plot(f, Ssb, 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Power (Hz)');title ('Spectrum Neuron_b');set(gca, 'FontSize', 14)
    subplot(2,4,6);
    plot(f, 10*log10(Sf), 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Power (Hz)');title(lfp);set(gca, 'FontSize', 14)
    subplot(2,4,7);
    plot(f, Cb, 'k', 'LineWidth', 2);xlim([0 120]);xlabel('Frequency');ylabel('Coherence');title('Coherence');set(gca, 'FontSize', 14);
    subplot(2,4,8);
    rose(phib);
    
    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Total spectrum done! (14%/100%)  ')

%% time-dependent spectrum
%speparam.trialave=0;
filepath=fullfile(pwd,'results'); mkdir(filepath);
cd(filepath);
filepath=fullfile(pwd,'spectrum'); mkdir(filepath);
cd(filepath);
filepath=fullfile(pwd,'mtp'); mkdir(filepath);
cd(filepath)
totalC=[];
totalSf=
for i=1:length(nameL)
    cuee=['Time spectrum .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk=nameS{i};
    eval(['[C,phi,~,Sf,Ss,t,f]=cohgramcpt(',lfp,',',spk,',','movingwin,params);']);
    figure;
    imagesc(t-1,f,Ss');colorbar;
    xlabel('Time(s)');ylabel('Frequency(Hz)');title ([lfp,' Time-frequency']);set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savePath, [lfp,'.png']));
    figure;
    imagesc(t-1,f,Sf');colorbar;
    xlabel('Time(s)');ylabel('Frequency(Hz)');title ([spk,' Spike Time-frequency']);set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savePath,[spk,'.png']));
    figure;
    imagesc(t-1,f,C');colorbar;
    xlabel('Time(s)');ylabel('Frequency(Hz)');title ([spk,'_',lfp,' Coherence']);set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savePath, [spk,'_',lfp,' Cohe.png']));
    
    imagesc(t-1,f,phia');colorbar;
    xlabel('time(s)');ylabel('Frequency(Hz)');title ('Phase of Coherence');set(gca, 'FontSize', 14)
    subplot(2,4,5);
    imagesc(t-1,f,Ssb');colorbar;
    xlabel('time(s)');ylabel('Frequency(Hz)');title ('Time-frequency Neuron b');set(gca, 'FontSize', 14)
    subplot(2,4,6);
    imagesc(t-1,f,Sf');colorbar;
    xlabel('time(s)');ylabel('Frequency(Hz)');title (lfp);set(gca, 'FontSize', 14)
    subplot(2,4,7);
    imagesc(t-1,f,Cb');colorbar;
    xlabel('time(s)');ylabel('Frequency(Hz)');title ('Coherence');set(gca, 'FontSize', 14)
    subplot(2,4,8);
    imagesc(t-1,f,phib');colorbar;
    xlabel('time(s)');ylabel('Frequency(Hz)');title ('Phase of Coherence');set(gca, 'FontSize', 14)
    
    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Time spectrum done! (28%/100%)  ')
%%
%STA & FTA analysis
%{
speparam.Fs=Fs;
speparam.trailave=0;
t=0:1/40000:3;
for i=1:length(nameL)
    lfp=nameL{i};
    for j=1:length(nameS)
        spk1=nameS{j};
        spk2=nameS{j+1};
        eval(['spk1=spk_ts2bin(',spk1,',','speparam);']);
        eval(['Analyze_spklfp_spectrum(spk1,',lfp,',speparam);']);
        %subplot(2,2,1);plot()
        eval(['spk2=spk_ts2bin(',spk2,',','speparam);']);
        eval(['Analyze_Spike_LFP(spk2,',lfp,',speparam);']);
        input('Press Enter to continue...');
    end
end
%}

%% chronux sta- time independent
smp=0:1/40000:3;
T=[0 3];
D=[-0.1 0.1];
plt='n';
filepath=[path '\results\indepsta']; mkdir(filepath);
for i=1:length(nameL)
    cuee=['STA .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[s1,t1] = sta(',spk1,',',lfp,',','smp,plt,0,T,D);']);
    eval(['[s2,t2] = sta(',spk2,',',lfp,',','smp,plt,0,T,D);']);
    indepsta(2*i-1).channel=lfp;indepsta(2*i-1).neuron='a';indepsta(2*i-1).sta=s1;indepsta(2*i-1).t=t1;
    indepsta(2*i).channel=lfp;indepsta(2*i).neuron='b';indepsta(2*i).sta=s2;indepsta(2*i).t=t2;
    figure('Visible', 'off');
    subplot(1,2,1);plot(t1,mean(s1)-mean(s1(:)) ,'k', 'LineWidth', 2);
    %xticks(linspace(1, length(t), 3)); xticklabels(linspace(-0.1,0.1,3));
    xlabel('Time (s)');ylabel('Voltage (mV)');title([lfp 'STA Neuron a']);
    subplot(1,2,2);plot(t2,mean(s2)-mean(s2(:)) ,'k', 'LineWidth', 2);
    %xticks(linspace(1, length(t), 3)); xticklabels(linspace(-0.1,0.1,3));
    xlabel('Time (s)');ylabel('Voltage (mV)');title([lfp 'STA Neuron b']);
    
    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('STA done! (42%/100%)  ')

%% chronux sta- time dependent
smp=0:1/40000:3;
Tc=[0 3];
Tinc=0.1;
Tw=0.5;
D=[-0.1 0.1];
plt='n';
filepath=[path '\results\timesta']; mkdir(filepath);
for i=1:length(nameL)
    cuee=['Time STA .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);    

    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[S1,tau1,tc1] = staogram(',spk1,',',lfp,',smp,plt,Tc,Tinc,Tw,0,D);']);
    eval(['[S2,tau2,tc2] = staogram(',spk2,',',lfp,',smp,plt,Tc,Tinc,Tw,0,D);']);
    timesta(2*i-1).channel=lfp;timesta(2*i-1).neuron='a';timesta(2*i-1).sta=S1;timesta(2*i-1).tc=tc1;timesta(2*i-1).tau=tau1;
    timesta(2*i).channel=lfp;timesta(2*i).neuron='b';timesta(2*i).sta=S2;timesta(2*i).tc=tc2;timesta(2*i).tau=tau2;
    figure('Visible', 'off');
    subplot(1,2,1);imagesc(tc1-1,tau1,squeeze(S1)');  set(gca,'ydir','normal');xlabel('time (s)');ylabel('lag(s)');title([lfp 'STA Neuron a']);colorbar;
    subplot(1,2,2);imagesc(tc2-1,tau2,squeeze(S2)');  set(gca,'ydir','normal');xlabel('time (s)');ylabel('lag(s)');title([lfp 'STA Neuron b']);colorbar;

    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Time spectrum done! (56%/100%) ')
%% chronux fta- time independent
t=0:1/40000:3;
fs=40000;
T=[0 3];
wb=[0 150];
plt='n';
filepath=[path '\results\indepfta']; mkdir(filepath);
for i=1:length(nameL)
    cuee=['FTA .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[fta1,phi1] = fta(',spk1,',',lfp,',','fs,wb,t);']);
    eval(['[fta2,phi2] = fta(',spk2,',',lfp,',','fs,wb,t);']);
    indepfta(2*i-1).channel=lfp;indepfta(2*i-1).neuron='a';indepfta(2*i-1).fta=fta1;indepfta(2*i-1).phi=phi1;
    indepfta(2*i).channel=lfp;indepfta(2*i).neuron='b';indepfta(2*i).fta=fta2;indepfta(2*i).phi=phi2;
    figure('Visible', 'off');
    subplot(1,2,1);plot(phi1, mean(fta1, 1), 'k', 'LineWidth', 2);
    xlabel('Phase');ylabel('FTA');title([lfp 'FTA Neuron a']);
    xlim([-pi pi]);ylim([0 0.05]);set(gca, 'FontSize', 14)
    subplot(1,2,2);plot(phi2, mean(fta2, 1), 'k', 'LineWidth', 2);
    xlabel('Phase');ylabel('FTA');title([lfp 'FTA Neuron b']);
    xlim([-pi pi]);ylim([0 0.05]);set(gca, 'FontSize', 14)
    
    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Time spectrum done! (70%/100%)  ')
%% chronux fta- time dependent
t=0:1/40000:3;
Tc=[0 3];
Tinc=0.05;
Tw=0.5;
wb=[0 120];
fs=40000;
plt='n';
filepath=[path '\results\timefta']; mkdir(filepath);
for i=1:length(nameL)
    cuee=['Time FTA .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[F1,phi1,tc1] = ftaogram(',spk1,',',lfp,',t,Tc,Tinc,Tw,wb,fs);']);
    eval(['[F2,phi2,tc2] = ftaogram(',spk2,',',lfp,',t,Tc,Tinc,Tw,wb,fs);']);
    timefta(2*i-1).channel=lfp;timefta(2*i-1).neuron='a';timefta(2*i-1).fta=F1;timefta(2*i-1).phi=phi1;timefta(2*i-1).tc=tc1;
    timefta(2*i).channel=lfp;timefta(2*i).neuron='b';timefta(2*i).fta=F2;timefta(2*i).phi=phi2;timefta(2*i).tc=tc2;
    figure('Visible', 'off');
    subplot(1,2,1);imagesc(tc1-1,phi1,squeeze(F1)');  set(gca,'ydir','normal');xlabel('time (s)');ylabel('phase');title([lfp 'FTA Neuron a']);colorbar;colormap(hot);
    subplot(1,2,2);imagesc(tc2-1,phi2,squeeze(F2)');  set(gca,'ydir','normal');xlabel('time (s)');ylabel('phase');title([lfp 'FTA Neuron b']);colorbar;colormap(hot);

    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Time spectrum done! (84%/100%) \n ')
%% GC same channel
filepath=[path '\results\GC']; mkdir(filepath);
for i=1:length(nameL)
    cuee=['Granger Causality .......... ' num2str(i) '/' num2str(channelnum)];
    disp(cuee);
    
    lfp=nameL{i};
    spk1=nameS{2*i-1};
    spk2=nameS{2*i};
    eval(['[~,~,~,~,tt1,f1,GCs2p_1,GCp2s_1]=spk_lfp_GC(',lfp,',',spk1,',params,movingwin);']);
    eval(['[~,~,~,~,tt2,f2,GCs2p_2,GCp2s_2]=spk_lfp_GC(',lfp,',',spk2,',params,movingwin);']);
    GC(2*i-1).channel=lfp;GC(2*i-1).neuron='a';GC(2*i-1).t=tt1;GC(2*i-1).f=f1;GC(2*i-1).s2p=GCs2p_1;GC(2*i-1).p2s=GCp2s_1;
    GC(2*i).channel=lfp;GC(2*i).neuron='b';GC(2*i).t=tt2;GC(2*i).f=f2;GC(2*i).s2p=GCs2p_2;GC(2*i).p2s=GCp2s_2;
    figure('Visible', 'off');
    subplot(1,2,1);
    plot(tt1-1, sum(GCs2p_1,2));title([lfp ' GCs2p Neuron a']);xlabel('time (s)');ylabel('GC');hold on;
    plot(tt1-1, sum(GCp2s_1,2));title([lfp ' GCp2s Neuron a']);xlabel('time (s)');ylabel('GC');hold off;
    subplot(1,2,2);
    plot(tt2-1, sum(GCs2p_2,2));title([lfp ' GCs2p Neuron b']);xlabel('time (s)');ylabel('GC');hold on;
    plot(tt2-1, sum(GCp2s_2,2));title([lfp ' GCp2s Neuron b']);xlabel('time (s)');ylabel('GC');hold off;
    
    fig = gcf; fig.WindowState = 'maximized';
    filename = [filepath '\channel' lfp(4:5) '.png'];saveas(fig, filename);
end
disp('Time spectrum done! (98%/100%) \n ')
%% Save results
save('result.mat', 'spectrum','tspectrum','indepsta','timesta','indepfta','timefta','GC');
close all;
disp('Data of results saved! (100%/100%)  ')
%% GC cross channels
%{
for i=1:length(nameL)
    tmpname=nameL{i};
    eval(['lfp=' tmpname,';']);
    for j=1:length(nameS)
        tmpname=nameS{j};
        eval(['spk=' tmpname,';']);
        [S1,S2,C,Phi,tt,f,GCs2p,GCp2s]=spk_lfp_GC(lfp,spike,params,movingwin);
        figure;
        %subplot(2,2,1);plot()
    end
end



% 在命令窗口显示筛选出的变量名
%disp(selectedVars);
%}

%% filter test

%{

[ys,f]=mtspectrumc(LFP01,params);
plot(f,10*log10(ys));xlim([0 80])

    Wn = [49 51]/20000;
    [B,A]=butter(2,Wn);
    tm=LFP01';
    for k = 1:size(tm,1)
        LFPf(k,:) = filter(B, A, tm(k, :));
    end
%     ord = 10000;
%     b = fir1(ord, Wn);
%     tm=LFP01';
%     for k = 1:size(tm,1)
%         LFPf(k,:) = filtfilt(b, 1, tm(k, :));
%     end

[hs,f]=mtspectrumc(LFPf',params);
plot(f,10*log10(hs));xlim([0 80])

%}

