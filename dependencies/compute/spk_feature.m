function spk_feature(sitepath)
%spk单个neuro的feature图
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
% if exist("Combined_cut_new.mat")
load p_data_cut_inc.mat
% else 
%     load Combined_cut.mat
% end

savePath = fullfile(sitepath,'results\spk_feature');
if isfolder(savePath)
    rmdir(savePath,'s');
end
mkdir(savePath);
%%
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 200];
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


%% plot waveform and psth
close all;
vars_wf=who('-regexp','^Waveform');
vars_wt=who('-regexp','^Wavetime');
vars_lfp=who('-regexp','^LFP');
vars_spk=who('-regexp','^SPKC');

for i=1:numel(vars_spk)
    eval(['wf=',vars_wf{i},';']);
    eval(['wt=',vars_wt{i},';']);
    eval(['spk=',vars_spk{i},';']);

    tmpname=vars_spk{i};tmpname=tmpname(end-2:end);
    orderspk=eval(['Order',tmpname]);
    tmplfp=eval(['LFP',tmpname(1:2)]);
    lfporder=eval(['Order',tmpname(1:2)]);
    [lfp,spk] = match_trials(tmplfp, spk, lfporder, orderspk);
    %% firing rate 
        [spkc,spkct]=binspikes(spk,1/0.03,[0.151,2.341]);
    fr=spkc./0.03;
    fr=mean(fr,2);
    eval(['frRe.',vars_spk{i},'=fr;'])

    % %%
    % 
    % 
    % 
    % [spkc,spkct]=binspikes(spk,50,[0,2.5]);
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
save(fullfile(savePath,'fr.mat'),'-struct',"frRe");
