function LFP_sd_delete(sitepath)

combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load p_data_cut.mat
load p_data_detail.mat
savePath=fullfile(combinepath,'p_data_cut_sd.mat');
%%
detail_change=0;
vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');

outliers=[];
for j=1:numel(vars_lfp)
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

eval([vars_spk{j},'(removed_trials_index)=[];']);
eval([vars_lfp{j},'(:,removed_trials_index)=[];']);
eval([vars_ang{j},'(removed_trials_index)=[];']);
eval([vars_spd{j},'(removed_trials_index)=[];']);

tmpdata=eval(vars_lfp{j});
if isempty(tmpdata)
    clearvars(vars_spk{j},vars_lfp{j},vars_ang{j},vars_spd{j},vars_wf{j},vars_wft{j});
        increaseName(strcmp(increaseName, vars_spk{j})) = [];
        decreaseName(strcmp(decreaseName, vars_spk{j})) = [];
        detail_change=1;
end
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

eval([vars_spk{j},'(removed_trials_index)=[];']);
eval([vars_lfp{j},'(:,removed_trials_index)=[];']);
eval([vars_ang{j},'(removed_trials_index)=[];']);
eval([vars_spd{j},'(removed_trials_index)=[];']);

tmpdata=eval(vars_lfp{j});
if isempty(tmpdata)
    clearvars(vars_spk{j},vars_lfp{j},vars_ang{j},vars_spd{j},vars_wf{j},vars_wft{j});
        increaseName(strcmp(increaseName, vars_spk{j})) = [];
        decreaseName(strcmp(decreaseName, vars_spk{j})) = [];
        detail_change=1;
end
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
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');

savevar=[vars_lfp;vars_spk;vars_ang;vars_spd;vars_wf;vars_wft];
save(savePath, savevar{:});
if detail_change==1
    save(fullfile(combinepath,'p_data_detail.mat'), 'increaseName','decreaseName');
end

end
