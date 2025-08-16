clc;clear;close all;
cd G:\2022_Mn_14365\20221117_Expt_Week08\D4_site02\combined\
load Combined.mat

savePath = 'G:\2022_Mn_14365\20221117_Expt_Week08\D4_site02\results\channel_lfp';
mkdir(savePath);
%%
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 30];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.

movingwin=[0.15 0.01];
%%
% aquire channel
vars = who;
numbers_channel = {};
channels_name={};
for i = 1:length(vars)
    matches = regexp(vars{i}, '\d+', 'match');
    if ~isempty(matches)
        numbers_channel = [numbers_channel, matches];
        channels_name=[channels_name,['channel',matches{1}]];
    end
end
numbers_channel = unique(numbers_channel);
num_channel=[];
for i=1:numel(numbers_channel)
    num_channel=[num_channel,str2num(numbers_channel{i})];
end
num_channel=sort(num_channel);
channels_name=unique(channels_name);

%%
probeDic=[3,7,11,15,14,10,6,2,1,5,9,13,16,12,8,4];
probe_num=probeDic(num_channel);
[probe,idx]=sort(probe_num);
probe_name=[];
for pp=1:numel(probe)
    probe_name=[probe_name,{['Probe',sprintf('%02d', probe(pp));]}];
end
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s

%% 计算channel base要，但是比较channel不运行这段
% aquire 0.2s-1s 实际是1.2s-2s 
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s
timeindx=recordtime>=0.2 & recordtime <=1 | abs(recordtime-0.2)<eps | abs(recordtime-1)<eps;
needtime=recordtime(timeindx);
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
for k=1:numel(vars_spk)
    eval(['lfp=',vars_lfp{k},';']);
    eval(['spk=',vars_spk{k},';']);
    for kk=1:numel(spk)
        tmpts=spk(kk).times;
        tmpindx=tmpts>=1.2 & tmpts<=2 | abs(tmpts-1.2)<eps | abs(tmpts-2)<eps;
        tmpts=tmpts(tmpindx);
        % extratime=1.2;%onset变为0
        % spk(kk).times=tmpts-extratime;   
        spk(kk).times=tmpts;   
    end
    lfp=lfp(timeindx,:);
    eval([vars_lfp{k},'=lfp;']);
    eval([vars_spk{k},'=spk;']);
end
disp('now ts is 0.2-1s!(actually 1.2-2s) ')
%disp('now ts is 0.2-1s! and onset is 0 point!!')

%%
% coherence - neuron base&channel base
results_cohe=struct();
for i=1:numel(numbers_channel)
    disp(['channel ',numbers_channel{i},' is processing...'])
    vars_lfp = who('-regexp', ['^LFP',numbers_channel{i}]);
    vars_spk = who('-regexp', ['^SPKC',numbers_channel{i}]);
    totallfp=[];
    totalspk=[];
    num_neuron=numel(vars_spk);
    figure;
    set(gcf,'unit','normalized','position',[0,0,1,0.62]);
    hold on;
    for j=1:num_neuron
        eval(['lfp=',vars_lfp{j},';']);
        eval(['spk=',vars_spk{j},';']);
        disp(['neuron ',vars_spk{j},' is processing...'])
        % !!!change all timestamp into positive num!
        [dN,t]=binspikes(spk,Fs,[1.2 2]);
        lfp=lfp_norm(lfp);
        totallfp=[totallfp,lfp];
        totalspk=[totalspk,dN];
        [C,phi,~,Sf,Ss,f]=coherencycpb(lfp,dN,params);
        eval(['cohe',numbers_channel{i},'.',vars_spk{j},'=C;']);
        subplot(1,num_neuron+1,j);
        plot(f, C, 'k', 'LineWidth', 2);xlim([0 30]);xlabel('Frequency');ylabel('Coherence');set(gca, 'FontSize', 14);
        xticks([0 4 8 12 16 21 30])
        title(['neuron ',extractAfter(vars_spk{j},'SPKC')]);
    end
    [C,phi,~,Sf,Ss,f]=coherencycpt(totallfp,totalspk,params);
    eval(['cohe',numbers_channel{i},'.f=f;']);
    eval(['cohe',numbers_channel{i},'.combineC=C;']);
    eval(['cohe',numbers_channel{i},'.lfp=totallfp;']);
    eval(['cohe',numbers_channel{i},'.Sf=Sf;']);
    subplot(1,num_neuron+1,num_neuron+1);
    plot(f, C, 'k', 'LineWidth', 2);xlim([0 30]);xlabel('Frequency');ylabel('Coherence');set(gca, 'FontSize', 14);
    xticks([0 4 8 12 16 21 30])
    title('Combined');
    hold off
    sgtitle(['Channel ',numbers_channel{i},' Coherence']);
    saveas(gcf, fullfile(savePath, ['Channel',numbers_channel{i},'.png']));
    disp(['channel ',numbers_channel{i},' done'])
    eval(['results_cohe.Channel',numbers_channel{i},'=cohe',numbers_channel{i},';'])
end

save(fullfile(savePath,'results_cohe.mat'),'-struct','results_cohe');
disp('channel base-neuron base finished!')
%%
% %% 单独绘图上面的coherence结果
% close all
% clear cohe
% vars_cohe=who('-regexp','^cohe');
% save('cohe_channel_new.mat', vars_cohe{:});
% for n=1: numel(vars_cohe)
%     eval(['cohe=',vars_cohe{n},';']);
%     vars_ncoh=fieldnames(cohe);
%     vars_ncoh=vars_ncoh(startsWith(vars_ncoh,'SPKC'));
%     figure;
%     set(gcf,'unit','normalized','position',[0,0,1,0.62]);
%     for nn=1:numel(vars_ncoh)
%         subplot(1,numel(vars_ncoh)+1,nn);
%         f=cohe.f;
%         eval(['c=cohe.',vars_ncoh{nn},';']);
%         plot(f, mean(c,1),'k','LineWidth', 2);xlim([0 30]);xlabel('Frequency');ylabel('Coherence');set(gca, 'FontSize', 14);
%         xticks([0,4,8,12,21,30]);
%         title(['neuron ',extractAfter(vars_ncoh{nn},'SPKC')]);
%     end
%     subplot(1,numel(vars_ncoh)+1,numel(vars_ncoh)+1);
%     plot(cohe.f, mean(cohe.combineC,1),'k','LineWidth', 2);
%     xlim([0 30]);xlabel('Frequency/Hz');ylabel('Coherence');set(gca, 'FontSize', 14);
%     title('Combined');
%     xticks([0,4,8,12,21,30]);
%     sgtitle(['Channel ',extractAfter(vars_cohe{n},'cohe'),' Coherence']);
%     savePath = '/Users/lai/Desktop/AnneWen/neurolfp/Result/Combine/nocut/channel_timecoh';
%     fileName = ['Channel',extractAfter(vars_cohe{n},'cohe'),'.png'];
%     % 使用saveas保存图形到指定路径
%     saveas(gcf, fullfile(savePath, fileName));
% end
% clear cohe

%% choherence compare
clear cohe
vars_cohe=who('-regexp','^cohe');
figure;
hold on;
for p=1: numel(idx)
    eval(['cohe=',vars_cohe{idx(p)},';']);
    plot(cohe.f, mean(cohe.combineC,1),'LineWidth', 2);
    
end
legend(probe_name);
xticks([0,4,8,12,21,30]);
xlim([0 30]);
xlabel('Frequency');ylabel('Coherence');set(gca, 'FontSize', 14);

figure;
hold on;
for p=1: numel(idx)
    eval(['cohe=',vars_cohe{idx(p)},';']);
    plot(cohe.t, mean(cohe.combineC,2),'LineWidth', 2);
    
end
legend(probe_name);
%xticks([0,4,8,12,21,30]);
%xlim([0 30]);
xlabel('Time/s');ylabel('Coherence');set(gca, 'FontSize', 14);

%%
% totallfp compare
totallfp_probe=struct();
for i=1:numel(idx)
    disp(['channel ',numbers_channel{idx(i)},' is processing...'])
    vars_lfp = who('-regexp', ['^LFP',numbers_channel{idx(i)}]);
    totallfp=[];
    num_lfp=numel(vars_lfp);
    for j=1:num_lfp
        eval(['lfp=',vars_lfp{j},';']);
        disp(['lfp ',vars_lfp{j},' is processing...'])
        
        % !!!change all timestamp into positive num!
        totallfp=[totallfp,lfp];
    end
   
    probe_num= num2str(probe(i), '%02d');
    eval(['Probe',probe_num,'.LFP=totallfp;'])

    totallfp=lfp_norm(totallfp);

    spectro=[];
    for j=1:size(totallfp,2)
        tmplfp=totallfp(:,j);
        [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[0 30]);
        spectro(:,:,j)=scalogram;
    end
    spectro=mean(spectro,3);
    totalBand.S=spectro;
    totalBand.f=frequency;
    totalBand.t=t;
    eval(['Probe',probe_num,'.totalBand=totalBand;'])
    
    spectro=[];
    for j=1:size(totallfp,2)
        tmplfp=totallfp(:,j);
        [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[0 4]);
        spectro(:,:,j)=scalogram;
    end
    spectro=mean(spectro,3);
    deltaBand.S=spectro;
    deltaBand.f=frequency;
    deltaBand.t=t;
    eval(['Probe',probe_num,'.deltaBand=deltaBand;'])

    spectro=[];
    for j=1:size(totallfp,2)
        tmplfp=totallfp(:,j);
        [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[4 8]);
        spectro(:,:,j)=scalogram;
    end
    spectro=mean(spectro,3);
    thetaBand.S=spectro;
    thetaBand.f=frequency;
    thetaBand.t=t;
    eval(['Probe',probe_num,'.thetaBand=thetaBand;'])
    
    spectro=[];
    for j=1:size(totallfp,2)
        tmplfp=totallfp(:,j);
        [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[8 12]);
        spectro(:,:,j)=scalogram;
    end
    spectro=mean(spectro,3);
    alphaBand.S=spectro;
    alphaBand.f=frequency;
    alphaBand.t=t;
    eval(['Probe',probe_num,'.alphaBand=alphaBand;'])

    spectro=[];
    for j=1:size(totallfp,2)
        tmplfp=totallfp(:,j);
        [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[12 30]);
        spectro(:,:,j)=scalogram;
    end
    spectro=mean(spectro,3);
    betaBand.S=spectro;
    betaBand.f=frequency;
    betaBand.t=t;
    eval(['Probe',probe_num,'.betaBand=betaBand;'])

    eval(['totallfp_probe.Probe',probe_num,'=Probe',probe_num,';'])
end
save(fullfile(savePath,'ProbeLFP.mat'),'-struct',"totallfp_probe");


vars_probelfp = who('-regexp', '^Probe');
colors=parula(numel(vars_probelfp));
figure;
set(gcf,'unit','normalized','position',[0,0,0.5,0.5]);
hold on;
for p=1:numel(vars_probelfp)
    eval(['lfp=',vars_probelfp{p},'.totalBand;']);
    plot(lfp.t-1,(10*log10(sum(lfp.S,1))),'LineWidth', 2,'Color',colors(p,:));
end
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off;
xlabel('Time/s');ylabel('dB');
title('Total Band -Total Channels');
lg=legend(probe_name','location','BestOutside');
lg.Direction = 'reverse';
%lg.Position=[0.4205  0.87  0  0]; 
saveas(gcf, fullfile(savePath, ['totalband','.png']));

figure;
set(gcf,'unit','normalized','position',[0,0,0.6,1]);
t = tiledlayout('flow','TileSpacing','compact');
nexttile
hold on;
for p=1:numel(vars_probelfp)
    eval(['lfp=',vars_probelfp{p},'.deltaBand;']);
    plot(lfp.t-1,(10*log10(sum(lfp.S,1))),'LineWidth', 2,'Color',colors(p,:));
end
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off;
xlabel('Time/s');ylabel('dB');
title('Delta Band');


nexttile
hold on;
for p=1:numel(vars_probelfp)
    eval(['lfp=',vars_probelfp{p},'.thetaBand;']);
    plot(lfp.t-1,(10*log10(sum(lfp.S,1))),'LineWidth', 2,'Color',colors(p,:));
end
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off;
xlabel('Time/s');ylabel('dB');
title('Theta Band');


nexttile
hold on;
for p=1:numel(vars_probelfp)
    eval(['lfp=',vars_probelfp{p},'.alphaBand;']);
    plot(lfp.t-1,(10*log10(sum(lfp.S,1))),'LineWidth', 2,'Color',colors(p,:));
end
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off;
xlabel('Time/s');ylabel('dB');
title('Alpha Band');


nexttile
hold on;
for p=1:numel(vars_probelfp)
    eval(['lfp=',vars_probelfp{p},'.betaBand;']);
    plot(lfp.t-1,(10*log10(sum(lfp.S,1))),'LineWidth', 2,'Color',colors(p,:));
end
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off;
xlabel('Time/s');ylabel('dB');
title('Beta Band');

lgd=legend(probe_name);
lgd.Layout.Tile='east';
lgd.Direction = 'reverse';

% % 创建额外的坐标区域
% ax = axes('Units', 'normalized', 'Position', [1, 0, 0.05, 1], 'Visible', 'on');
% 
% % 绘制图例
% legend(ax, probe_name, 'Location', 'bestoutside');

saveas(gcf, fullfile(savePath, ['varyBand','.png']));

%%
num_colors = 9;
%colors = distinguishable_colors(num_colors);
colors=parula(9);

%totallfpcoh=totallfpcohh;

savePath = '/Users/lai/Desktop/AnneWen/neurolfp/Result/Combine/nocut/channel_lfp';
%0-30Hz
figure;
hold on
for i=1:size(totallfpcoh,2)
    plot(f,(totallfpcoh(:,i)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlim([0,30]);
xlabel('Frequency/Hz');ylabel('Power');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'total.png'));

%delta
figure;
hold on
for i=1:size(totallfpcoh,2)
    plot(f,10*log10(totallfpcoh(:,i)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlim([0,4]);
xlabel('Frequency/Hz');ylabel('Power');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'delta.png'));

%theta
figure;
hold on
for i=1:size(totallfpcoh,2)
    plot(f,10*log10(totallfpcoh(:,i)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlim([4,8]);
xlabel('Frequency/Hz');ylabel('Power');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'theta.png'));

%alpha
figure;
hold on
for i=1:size(totallfpcoh,2)
    plot(f,10*log10(totallfpcoh(:,i)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlim([8,12]);
xlabel('Frequency/Hz');ylabel('Power');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'alpha.png'));

%beta
figure;
hold on
for i=1:size(totallfpcoh,2)
    plot(f,10*log10(totallfpcoh(:,i)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlim([12,30]);
xlabel('Frequency/Hz');ylabel('Power');
lg=legend(probe_name);
l.Direction = 'reverse';
saveas(gcf, fullfile(savePath, 'beta.png'));

%% plot waveform and psth
close all;
vars_wf=who('-regexp','^Waveform');
vars_wt=who('-regexp','^Wavetime');
vars_lfp=who('-regexp','^LFP');
vars_spk=who('-regexp','^SPKC');
for i=1:numel(vars_wf)
    eval(['wf=',vars_wf{i},';']);
    eval(['wt=',vars_wt{i},';']);
    eval(['lfp=',vars_lfp{i},';']);
    eval(['spk=',vars_spk{i},';']);
    [spkc,spkct]=binspikes(spk,50,[0,2.5]);
    figure;
    set(gcf,'unit','normalized','position',[0,0,0.2,1]);
    subplot(4,1,1);
    plot(wt,wf,'red','LineWidth',2);
    xlabel('Time/s');ylabel('Volt');
    title('Waveform');
    
    subplot(4,1,2);
    plot(recordtime,mean(lfp,2),'k','LineWidth',1);
    hold on;
    %plot(recordtime,spkc,'yellow','LineWidth',1.5);
    line([0.2,0.2],ylim,'LineWidth',2,'Color','red');
    line([1,1],ylim,'LineWidth',2,'Color','Blue');
    hold off;
    xlabel('Time/s');ylabel('Volt');
    xlim([-1,1.5]);
    title('LFP');

    subplot(4,1,3);
    bar(spkct-1,sum(spkc,2),'histc');
    hold on;
    line([0.2,0.2],ylim,'LineWidth',2,'Color','red');
    line([1,1],ylim,'LineWidth',2,'Color','Blue');
    hold off;
    xlim([-1,1.5]);xlabel('Time/s');ylabel('Average Count');
    title('Spike Count (bin==0.02s)')

    [N,B]=isi(spk,[],0,100);
    subplot(4,1,4);
    if ~isempty(N)
    bar(B,N,'histc');
    end
    xlabel('Time');ylabel('Count');
    title('ISI');

    sgtitle(vars_spk{i});
    saveas(gcf, fullfile(savePath, [vars_spk{i},'.png']));
end

%% time-coherence
movingwin=[0.15 0.01];

close all;
vars_lfp=who('-regexp','^LFP');
vars_spk=who('-regexp','^SPKC');
for i=1:numel(vars_spk)
    eval(['lfp=',vars_lfp{i},';']);
    eval(['spk=',vars_spk{i},';']);
    [C,phi,~,Sf,Ss,t,f]=cohgramcpt(lfp,spk,movingwin,params);
    
end

%% time-LFP
movingwin=[0.15 0.01];

close all;
vars_lfp=who('-regexp','^LFP');
vars_spk=who('-regexp','^SPKC');
timelfp=[];
for i=1:numel(vars_spk)
    eval(['lfp=',vars_lfp{i},';']);
    [S,t,f]=mtspecgramc(lfp,movingwin,params);
    timelfp=[timelfp,S];
end
timelfpt=t;
timelfpf=f;

%%
close all;
num_colors = 9;
colors=parula(9);

savePath = '/Users/lai/Desktop/AnneWen/neurolfp/Result/Combine/nocut/channel_timelfp';

%delta
figure;
hold on
for i=1:length(vars_cohe)
    eval(['tmptimelfp=',vars_cohe{i},'.Sf;']);
    eval(['timelfpf=',vars_cohe{i},'.f;']);
    eval(['timelfpt=',vars_cohe{i},'.t;']);
    indx=timelfpf<=4;
    tmpS=tmptimelfp(:,indx);
    plot(timelfpt,zscore(mean(10*log10(tmpS),2)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlabel('time/s');ylabel('z-score');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'delta.png'));

%theta
figure;
hold on
for i=1:length(vars_cohe)
    eval(['tmptimelfp=',vars_cohe{i},'.Sf;']);
    eval(['timelfpf=',vars_cohe{i},'.f;']);
    eval(['timelfpt=',vars_cohe{i},'.t;']);
    indx=timelfpf>4 & timelfpf<=8;
    tmpS=tmptimelfp(:,indx);
    plot(timelfpt,zscore(mean(10*log10(tmpS),2)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlabel('time/s');ylabel('z-score');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'theta.png'));

%alpha
figure;
hold on
for i=1:length(vars_cohe)
    eval(['tmptimelfp=',vars_cohe{i},'.Sf;']);
    eval(['timelfpf=',vars_cohe{i},'.f;']);
    eval(['timelfpt=',vars_cohe{i},'.t;']);
    indx=timelfpf>8 &timelfpf<=12;
    tmpS=tmptimelfp(:,indx);
    plot(timelfpt,zscore(mean(10*log10(tmpS),2)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlabel('time/s');ylabel('z-score');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'alpha.png'));

%beta
figure;
hold on
for i=1:length(vars_cohe)
    eval(['tmptimelfp=',vars_cohe{i},'.Sf;']);
    eval(['timelfpf=',vars_cohe{i},'.f;']);
    eval(['timelfpt=',vars_cohe{i},'.t;']);
    indx=timelfpf>12;
    tmpS=tmptimelfp(:,indx);
    plot(timelfpt,zscore(mean(10*log10(tmpS),2)),'LineWidth', 2,'Color',colors(i,:));
end
hold off;
xlabel('time/s');ylabel('z-score');
legend(probe_name);
saveas(gcf, fullfile(savePath, 'beta.png'));
