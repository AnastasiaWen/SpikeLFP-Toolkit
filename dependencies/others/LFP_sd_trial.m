function LFP_sd_trial(trialpath)

load(fullfile(trialpath,'data_delch.mat'));
OrigAng=load(fullfile(trialpath,'StimAngle.mat'));
OrigSpe=load(fullfile(trialpath,'StimSpeed.mat'));
savePath=fullfile(trialpath,'data_sd.mat');
%%
detail_change=0;
Fs=1000;

vars_lfp = who('-regexp', '^LFP');


for i=1:numel(vars_lfp)
    tmpname=vars_lfp{i};
    tmpname=tmpname(end-1:end);
    sname=['StimAngle',tmpname];
    eval([sname,'=reshape(OrigAng.angleRec,[],1);']);
    sname=['StimSpeed',tmpname];
    eval([sname,'=reshape(OrigSpe.speedRec,[],1);']);
end

d = designfilt('bandpassfir', 'FilterOrder', 300, ...
   'CutoffFrequency1', 1, 'CutoffFrequency2', 300,...
   'SampleRate', Fs,'Window','hamming');
for i=1:size(LFP05,2)
    for j=1:numel(vars_lfp)
    LFP_data=eval(vars_lfp{j});
    LFP_data=LFP_data(:,i);
    LFP_data=filtfilt(d,LFP_data);
    eval([vars_lfp{j},'(:,i)=LFP_data;'])
    end
end
for i=1:size(LFP02,2)
    for j=1:numel(vars_lfp)
    LFP_data=eval(vars_lfp{j});
    LFP_data=LFP_data(:,i);
    LFP_data=filtfilt(d,LFP_data);
    eval([vars_lfp{j},'(:,i)=LFP_data;'])
    end
end



vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');

outliers=[];
for j=1:numel(vars_lfp)
    outliers=[];
    removed_trials_index=[];
% 假设 LFP_data 是 2500x120 的 LFP 数据矩阵
% 2500 行代表时间点，120 列代表不同的 trial
LFP_data=eval(vars_lfp{j});
for i = 1:size(LFP_data, 2)
    % 检查每个 trial 是否超出阈值范围
    if any(LFP_data(:, i) < -0.4 | LFP_data(:, i) > 0.4 )
        outliers(i) = true;
    end
end

% 保存被删除的 trial 的索引
removed_trials_index = find(outliers);

tmpname=vars_lfp{j};tmpname=tmpname(end-1:end);tmpname=['^SPKC',tmpname];
vars_spk = who('-regexp', tmpname);
nn=zeros(2500,1)*nan;
for hh=1:numel(removed_trials_index)
idxx=removed_trials_index(hh);
for hj=1:numel(vars_spk)
eval([vars_spk{hj},'(idxx).times=[];']);
end
eval([vars_lfp{j},'(:,idxx)=nn;']);
eval([vars_ang{j},'(idxx)=NaN;']);
eval([vars_spd{j},'(idxx)=NaN;']);
end
% for hj=1:numel(vars_spk)
% eval([vars_spk{hj},'(removed_trials_index)=[];']);
% end
% eval([vars_lfp{j},'(:,removed_trials_index)=[];']);
% eval([vars_ang{j},'(removed_trials_index)=[];']);
% eval([vars_spd{j},'(removed_trials_index)=[];']);

% % 删除超过阈值范围的 trial
% filtered_LFP_data = LFP_data(:, ~outliers);
% 
% % 输出结果
% fprintf('删除了 %d 个 trial。\n', sum(outliers));
% disp('被删除的 trial 的索引：');
% disp(removed_trials_index);
end


outliers=[];
for j=1:numel(vars_lfp)
    removed_trials_index=[];
    outliers=[];
% 假设 LFP_data 是 2500x120 的 LFP 数据矩阵
% 2500 行代表时间点，120 列代表不同的 trial
LFP_data=eval(vars_lfp{j});
% 预分配变量
mean_values = mean(LFP_data, 2);             % 计算每个 trial 的均值 (1x120)
std_values = std(LFP_data, 0, 2);            % 计算每个 trial 的标准差 (1x120)

% 定义阈值范围
lower_bound = mean_values - 5 * std_values;  % 下限
upper_bound = mean_values + 5 * std_values;  % 上限

% 找出超出范围的 trial
outliers = false(1, size(LFP_data, 2));      % 创建逻辑数组标记异常trial

for i = 1:size(LFP_data, 2)
    % 检查每个 trial 是否超出阈值范围
    if any(LFP_data(:, i) < lower_bound(i) | LFP_data(:, i) > upper_bound(i))
        outliers(i) = true;
    end
end


% 保存被删除的 trial 的索引
removed_trials_index = find(outliers);


tmpname=vars_lfp{j};tmpname=tmpname(end-1:end);tmpname=['^SPKC',tmpname];
vars_spk = who('-regexp', tmpname);
nn=zeros(2500,1)*nan;
for hh=1:numel(removed_trials_index)
idxx=removed_trials_index(hh);
for hj=1:numel(vars_spk)
eval([vars_spk{hj},'(idxx).times=[];']);
end
eval([vars_lfp{j},'(:,idxx)=nn;']);
eval([vars_ang{j},'(idxx)=NaN;']);
eval([vars_spd{j},'(idxx)=NaN;']);
end
% for hj=1:numel(vars_spk)
% eval([vars_spk{hj},'(removed_trials_index)=[];']);
% end
% eval([vars_lfp{j},'(:,removed_trials_index)=[];']);
% eval([vars_ang{j},'(removed_trials_index)=[];']);
% eval([vars_spd{j},'(removed_trials_index)=[];']);
% 
% 

% % 删除超过阈值范围的 trial
% filtered_LFP_data = LFP_data(:, ~outliers);
% 
% % 输出结果
% fprintf('删除了 %d 个 trial。\n', sum(outliers));
% disp('被删除的 trial 的索引：');
% disp(removed_trials_index);
end
clearvars 'LFP_data'
vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');

savevar=[vars_lfp;vars_spk;vars_ang;vars_spd;vars_wf;vars_wft];
save(savePath, savevar{:});
end
