function LFP_bad_channel(trialpath)
close all;
load(fullfile(trialpath,'data.mat'));


vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPK');

hw=figure;
set(gcf,'unit','normalized','position',[0,0,0.4,1]);
ha=tight_subplot(8,2,0.04);

% figure;
hold on
cm=[];
for i=1:numel(vars_lfp)
    eval(['lfp=',vars_lfp{i},';']);
    % plot(lfp(:,1));
    cm=[cm,lfp(:,1)];

    axes(ha(i));
    plot(-1:1/1000:1.4999,lfp(:,1));
    title(vars_lfp{i});
    ylim([-0.2 0.2]);
end


% 假设你的LFP数据存储在变量 lfp_data 中，每列代表一个通道
% 比如 lfp_data 是一个 2500 x 16 的矩阵 (2500个时间点, 16个通道)
lfp_data = cm; % 这里是一个示例数据，替换为你的LFP数据
% 计算每个通道的均值或中位数
channel_max = max(lfp_data,[],1); % 计算每列的中位数
channel_min = min(lfp_data, [],1); % 计算每列的中位数

% 设定一个阈值，比如你可以设置为所有通道中位数的 10th 百分位
threshold1 = mean(channel_max)-1.5*std(channel_max);
threshold2 = mean(channel_min)+1.5*std(channel_min);
% 找出那些中位数小于阈值的通道 (坏通道)
bad_channels = find(channel_min >threshold2 & channel_max < threshold1);

%for 22数据
bad_channels=[];


vars_lfp(bad_channels)=[];
if  ~isempty(bad_channels)
for  hh=1:numel(bad_channels)
% 删除坏掉的通道
    
    axes(ha(bad_channels(hh)));
    title('delete','Color','red');
end
end

saveas(gcf,fullfile(trialpath,'orig_sig_channelcheck.png'))
savename=[vars_spk;vars_lfp];
save(fullfile(trialpath,'data_delch.mat'),savename{:});

