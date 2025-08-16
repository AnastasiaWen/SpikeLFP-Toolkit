function channel_cohe(sitepath)
%计算channel base的coherence
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load p_data_cut_inc.mat

savePath = fullfile(sitepath,'results\channel_coh');
if isfolder(savePath)
    rmdir(savePath,'s');
%     savePath = fullfile(sitepath,'results\channel_coh\nocut');
%     mkdir(savePath);
% else 
%     savePath = fullfile(sitepath,'results\channel_coh\nocut');
%     mkdir(savePath);
end
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
    [C,phi,~,Sf,Ss,f]=coherencycpb(totallfp,totalspk,params);
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

%% 单独绘图上面的coherence结果
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
% clear cohe
% vars_cohe=who('-regexp','^cohe');
% figure;
% hold on;
% for p=1: numel(idx)
%     eval(['cohe=',vars_cohe{idx(p)},';']);
%     plot(cohe.f, mean(cohe.combineC,1),'LineWidth', 2);
% 
% end
% legend(probe_name);
% xticks([0,4,8,12,21,30]);
% xlim([0 30]);
% xlabel('Frequency');ylabel('Coherence');set(gca, 'FontSize', 14);
% 
% figure;
% hold on;
% for p=1: numel(idx)
%     eval(['cohe=',vars_cohe{idx(p)},';']);
%     plot(cohe.t, mean(cohe.combineC,2),'LineWidth', 2);
% 
% end
% legend(probe_name);
% %xticks([0,4,8,12,21,30]);
% %xlim([0 30]);
% xlabel('Time/s');ylabel('Coherence');set(gca, 'FontSize', 14);


