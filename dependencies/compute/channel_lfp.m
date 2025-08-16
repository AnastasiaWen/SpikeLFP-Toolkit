function channel_lfp(sitepath)
%不同深度-不同频带LFP功率谱计算绘图
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load Combined.mat

savePath = fullfile(sitepath,'results\channel_lfp');
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


%%
% totallfp compare
totallfp_probe=struct();
powermax=[];powermin=[];
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
    powermax=max(spectro(:));powermin=min(spectro(:));
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

%%
vars_probelfp = who('-regexp', '^Probe');
%colors=parula(numel(vars_probelfp));
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

