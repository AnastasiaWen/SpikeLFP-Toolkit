function [tt,f,p_value_fisher, real_GCs2p,real_GCp2s] = permutation_test_granger(LFP, SPKC, num_permutations,ctl)
    % LFP: 输入的 LFP 数据
    % SPKC: 输入的 Spike 数据（spike times 或其他格式）
    % num_permutations: 置换次数
    % 返回 p_value_fisher: Fisher 合并后的 p 值
    % 返回 GC_real_s2p: 实际的 Spike 到 LFP 的 Granger 因果性值 (时间 x 1)
    % 返回 GC_real_p2s: 实际的 LFP 到 Spike 的 Granger 因果性值 (时间 x 1)
    % 返回 GC_perm_s2p: 置换生成的 Spike 到 LFP 的 Granger 因果性分布 (时间 x 置换次数)
    % 返回 GC_perm_p2s: 置换生成的 LFP 到 Spike 的 Granger 因果性分布 (时间 x 置换次数)

    % 设置参数
    Fs = 1000;                  % 采样频率
    TW = 3;						% Time-bandwidth product 设为 3
    ntapers = 2*TW-1;			% 使用 #tapers
    fpass = [0 55];             % 频率范围
    pad = 2;                    % Padding 参数
    
    params.Fs = Fs;				
    params.tapers = [TW, ntapers];
    params.fpass = fpass;      
    params.pad = pad;             
    params.trialave = 1;		
    
    movingwin = [0.3 0.03];     % 时间窗口参数 (0.3秒窗口，步长0.03秒)

    % 1. 计算实际的 Granger 因果性值
    [~,~,~,~,tt,f,real_GCs2p,real_GCp2s] = spk_lfp_GC(LFP, SPKC, params, movingwin);
    if ctl==1
        p_value_fisher=[1,1];
        return
    end
    
    % 对频率进行平均
    GC_real_s2p = mean(real_GCs2p, 2);  % 频率平均后的 Spike 到 LFP 的 Granger 因果性
    GC_real_p2s = mean(real_GCp2s, 2);  % 频率平均后的 LFP 到 Spike 的 Granger 因果性

    % 2. 初始化存储置换的 GC 值
    num_time_points = size(GC_real_s2p, 1);  % 时间点数
    GC_perm_s2p = zeros(num_time_points, num_permutations);  % 置换后的 GC (spike to LFP)
    GC_perm_p2s = zeros(num_time_points, num_permutations);  % 置换后的 GC (LFP to spike)

    % 3. 进行置换检验
    for i = 1:num_permutations
        % 随机打乱 spike 数据
        permuted_SPKC = SPKC(randperm(size(SPKC, 2)) );

        % 重新计算 Granger 因果性
        [~,~,~,~,~,~,perm_GCs2p,perm_GCp2s] = spk_lfp_GC(LFP, permuted_SPKC, params, movingwin);
        
        % 对频率进行平均并存储
        GC_perm_s2p(:, i) = mean(perm_GCs2p, 2);  % 置换后的 Spike 到 LFP 的 Granger 因果性
        GC_perm_p2s(:, i) = mean(perm_GCp2s, 2);  % 置换后的 LFP 到 Spike 的 Granger 因果性
        disp([num2str(i),'/1000']);
    end

    % 4. 计算 p 值：对于每个时间点计算显著性 p 值
    p_values_s2p = zeros(num_time_points, 1);
    p_values_p2s = zeros(num_time_points, 1);
    for t = 1:num_time_points
        % Spike 到 LFP 的 p 值
        p_values_s2p(t) = mean(GC_perm_s2p(t, :) >= GC_real_s2p(t));
        % LFP 到 Spike 的 p 值
        p_values_p2s(t) = mean(GC_perm_p2s(t, :) >= GC_real_p2s(t));
    end

    % 5. 使用 Fisher 合并检验来计算每个方向的总体 p 值
    p_value_fisher_s2p = fisher_combined_p(p_values_s2p);
    p_value_fisher_p2s = fisher_combined_p(p_values_p2s);

    % 返回 Fisher 合并后的 p 值
    p_value_fisher = [p_value_fisher_s2p, p_value_fisher_p2s];

    % 输出结果
    disp('Fisher combined p-values computed:');
    disp(['Spike to LFP p-value: ', num2str(p_value_fisher_s2p)]);
    disp(['LFP to Spike p-value: ', num2str(p_value_fisher_p2s)]);
end


