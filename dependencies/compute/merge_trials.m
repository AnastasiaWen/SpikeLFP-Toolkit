function [merged_Angle,merged_Speed,merged_SPKC, merged_Order]  = merge_trials(SPKC, Order, Angle, Speed)
    % SPKC: 输入的 spike 数据结构体数组，长度为 trial 数，每个结构体中有一个 'times' 字段保存 spike timestamps
    % Order: 对应的 order 数组，记录每个 trial 的信息，2xN 的数组
    % 返回 merged_SPKC: 合并后的 SPKC 结构体，所有重复的 trial 将被合并 spike times

    % 获取所有 trial 的 unique order（通过第二行，即每个 trial 对应的具体标识）
    [~, unique_indices, groups] = unique(Order', 'rows');
    
    % 初始化新的 SPKC 结构体数组
    merged_SPKC = struct('times', []);
    merged_SPKC = repmat(merged_SPKC, length(unique_indices), 1);
    merged_Order = zeros(size(Order, 1), length(unique_indices));  % 初始化新的 Order 数组
    merged_Angle = zeros(length(unique_indices), 1);
    merged_Speed = zeros(length(unique_indices), 1);

    % 遍历每个唯一的 trial
    for i = 1:length(unique_indices)
        % 找到所有相同 trial 的索引
        trial_indices = find(groups == i);
        
        % 合并这些 trial 的 spike times
        merged_times = [];
        for j = 1:length(trial_indices)
            merged_times = [merged_times;SPKC(trial_indices(j)).times];
        end
        
        % 将 spike times 排序并去重
        merged_times = unique(merged_times);
        
        % 保存到新的 SPKC 结构体中
        merged_SPKC(i).times = merged_times;

        merged_Order(:, i) = Order(:, trial_indices(1));
        merged_Angle(i) = Angle(trial_indices(1));
        merged_Speed(i) = Speed(trial_indices(1));
    end
end
