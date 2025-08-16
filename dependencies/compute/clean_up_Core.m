function clean_up_Core(sitepath)
%清理coherence的结果（保留increase的
if exist(fullfile(sitepath,'results/newGC/totalGC.mat'))
        load(fullfile(sitepath,'results/newGC/data_detail.mat'))
repath=fullfile(sitepath,'results/channel_cohe_ch/');
cd(repath)

    matData = load('totalCoh.mat');
    
    % 获取totalGC结构体中的字段名
    fields_s2p = fieldnames(matData.totalC);
    fields_p2s = fieldnames(matData.totalphi);

    % 获取所有 SPKC 通道
    spkc_channels = {}; % 存储所有 SPKC 通道的基本编号（如 01、02 等）
    for i = 1:length(increaseName)
        spkc_channels{end+1} = increaseName{i}(5:6); % 提取 SPKC 名字中的通道号 (如 '01')
    end
    spkc_channels = unique(spkc_channels); % 保证通道号唯一

    
    % 遍历当前文件夹中的文件
    fileList = dir('*.png'); % 获取所有.jpg图片
    for i = 1:length(fileList)
        fileName = fileList(i).name;
        
        % 提取文件名中的SPKC部分
        [lfp_part, spkc_part] = strtok(fileName, '_');
        spkc_part = strtok(spkc_part, '.'); % 去掉.jpg后缀
        spkc_part=spkc_part(2:end);
        
        % 检查该SPKC神经元是否在increaseName中
        if ~ismember(spkc_part, increaseName)
            % 删除不在increaseName中的图片
            delete(fileName);
            
            % 删除totalGCs2p和totalGCp2s中的相关字段
            if isfield(matData.totalC, [lfp_part '_' spkc_part])
                matData.totalC = rmfield(matData.totalC, [lfp_part '_' spkc_part]);
            end
            if isfield(matData.totalphi, [lfp_part '_' spkc_part])
                matData.totalphi = rmfield(matData.totalphi, [lfp_part '_' spkc_part]);
            end
        end

         % 检查对应的 LFP 是否有没有对应的 SPKC 通道
        lfp_channel = lfp_part(4:5); % 提取 LFP 名字中的通道号 (如 '03')
        if ~ismember(lfp_channel, spkc_channels)
            % 如果没有对应的 SPKC 通道，则删除该 LFP 图片和 MAT 文件中的相关数据
            delete(fileName);
            
            if isfield(matData.totalC, [lfp_part '_' spkc_part])
                matData.totalC = rmfield(matData.totalC, [lfp_part '_' spkc_part]);
            end
            if isfield(matData.totalphi, [lfp_part '_' spkc_part])
                matData.totalphi = rmfield(matData.totalphi, [lfp_part '_' spkc_part]);
            end
        end

    end
    
    % 保存修改后的totalGC.mat文件
    save('totalCoh2.mat', '-struct', 'matData');
else
    return
end
end
