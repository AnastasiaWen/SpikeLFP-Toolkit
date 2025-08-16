function [tmplfp,tmpspk] = match_trials(LFPdata, SPKCdata, OrderLFP, OrderSPKC)

    tmplfp = [];  % 初始化空结构体，用于存储匹配的 SPKC 数据
    validTrials = [];   % 初始化空数组，用于记录有效的 LFP trials

    % 遍历每个 LFP trial
    for i = 1:size(OrderSPKC, 2)  % 遍历 OrderLFP 的每一列（每个 trial）
        current_trial = OrderSPKC(:, i);  % 当前 LFP trial 的 Order 信息
        
        % 在 OrderSPKC 中查找相应的 trial
        matchIdx = find(ismember(OrderLFP', current_trial', 'rows'));

        if ~isempty(matchIdx)
            % 如果在 SPKC 中找到对应的 trial，保留 LFP trial 并将 SPKC 数据加入
            validTrials = [validTrials, i];  % 记录有效的 LFP trial
            % 将 SPKC 数据的当前匹配元素加入到 tmpspk 中
            tmplfp=[tmplfp,LFPdata(matchIdx)]; 
        end
    end

    % 删除 LFP 中没有匹配 SPKC 的 trial
    tmpspk = SPKCdata(validTrials);
end
